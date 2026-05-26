%% LAB 1: ACOUSTIC SOURCE LOCALIZATION (FINAL CODE)
% Group Members: Shahariar Ryehan, Marti Arjona, Pranav Prakasan
clear all; close all; clc;

% =========================================================================
% PART 1: DATA LOADING & PARAMETERS
% =========================================================================
% *** CHANGE THIS PATH TO YOUR DATA FOLDER ***
FileDir = 'D:\france semester\3_Acoustic Imaging\project3\shahariar\1_ReceivedData\';
BaseName = 'PCB'; 
AcqNumber = '1'; 
Nbmic = 14;

% Parameters
c = 343;          % Speed of sound (m/s)
d = 0.05;         % Pitch (5 cm)
mic_pos = (-(Nbmic-1)/2 : (Nbmic-1)/2) * d; % Mic positions (centered at 0)
fs = 144000;      % Default Sampling Frequency (will be updated from file)

% Load Data
y_all = []; 
for mic = 1:Nbmic
    FileName = fullfile(FileDir, [BaseName, num2str(mic), 'Acq', AcqNumber, '.mat']);
    if isfile(FileName)
        data = load(FileName);
        if mic == 1
            t = data.A.tr;
            if isfield(data.A, 'fsample')
                fs = data.A.fsample;
            else
                fs = 1/mean(diff(t)); 
            end
        end
        % Ensure column vector
        y_all(:, mic) = data.A.yr(:); 
    else
        error(['File not found: ', FileName]);
    end
end
% Normalize Data (Important for plots)
y_all = y_all / max(max(abs(y_all)));

fprintf('Data Loaded: %d samples x %d mics. Fs = %.0f Hz\n', size(y_all), fs);

% =========================================================================
% PART 2: SIGNAL ANALYSIS (Colorful Waterfall & Spectrum)
% =========================================================================

% 1. Frequency Spectrum (Received Signal)
Nfft = 2^12; 
f_vec = linspace(0, fs, Nfft);
Y_spectrum = fft(y_all, Nfft);

figure('Name', 'Fig 1: Frequency Spectrum', 'Color', 'w');
plot(f_vec, abs(Y_spectrum(:,1)), 'b', 'LineWidth', 1.5); 
xlim([0 10000]); xlabel('Frequency (Hz)'); ylabel('Magnitude');
title('Frequency Content of Received Signal (Mic 1)'); 
grid on;

% 2. Waterfall Plot (Colorful - Like Friend's Report)
figure('Name', 'Fig 2: Time Domain Waterfall', 'Color', 'w');
offset = 0.8; 
myColors = jet(Nbmic); % Generate 14 distinct colors

hold on;
for i = 1:Nbmic
    plot(t*1000, y_all(:,i) + (i-1)*offset, 'LineWidth', 1.2, 'Color', myColors(i,:)); 
end
xlabel('Time (ms)'); ylabel('Microphone Index');
title('Time Domain Signals across Array (Waterfall)');
axis tight; grid on; box on;

% =========================================================================
% PART 3: PLANE WAVE BEAMFORMING (Bartlett)
% =========================================================================
% Frequency Band for Processing
f_min = 3000; f_max = 7000;
[~, idx_min] = min(abs(f_vec - f_min));
[~, idx_max] = min(abs(f_vec - f_max));
freqs = f_vec(idx_min:idx_max);

theta_scan = -90:0.5:90; 
P_plane = zeros(length(theta_scan), 1);

% Incoherent Averaging Loop
for fi = 1:length(freqs)
    k = 2*pi*freqs(fi)/c;
    X_f = Y_spectrum(idx_min+fi-1, :).'; 
    
    % Steering Vector Matrix for all angles
    % V_plane size: (Nbmic x Angles)
    V_plane = exp(-1j * k * mic_pos' * sind(theta_scan)); 
    
    % Bartlett Power: |v'*x|^2
    % (V_plane' * X_f) gives a vector of responses for all angles
    P_plane = P_plane + abs(V_plane' * X_f).^2;
end

P_plane_dB = 10*log10(P_plane / max(P_plane));

figure('Name', 'Fig 3: Plane Wave Bartlett', 'Color', 'w');
plot(theta_scan, P_plane_dB, 'LineWidth', 2);
xlabel('Angle (degrees)'); ylabel('Normalized Power (dB)');
title('Plane Wave Beamforming (Bartlett)');
ylim([-20 0]); xline(0, '--r'); grid on;

% =========================================================================
% PART 4: POINT SOURCE BEAMFORMING (2D Maps)
% =========================================================================
% Grid Setup (Scanning Area)
x_grid = -1:0.02:1;     % X from -1m to 1m
y_grid = 0.5:0.02:2.5;  % Y from 0.5m to 2.5m
[X_map, Y_map] = meshgrid(x_grid, y_grid);
P_point = zeros(size(X_map));

fprintf('Computing Point Source Map (This may take a few seconds)...\n');

for r = 1:size(X_map, 1)
    for c_idx = 1:size(X_map, 2)
        target_pos = [X_map(r,c_idx); Y_map(r,c_idx)];
        
        % Calculate Distances
        dists = sqrt((mic_pos - target_pos(1)).^2 + (0 - target_pos(2)).^2);
        
        temp_pow = 0;
        for fi = 1:length(freqs)
            k = 2*pi*freqs(fi)/c;
            X_f = Y_spectrum(idx_min+fi-1, :).';
            
            % Steering Vector (Phase alignment only)
            v = exp(-1j * k * dists).';
            
            % Incoherent Sum
            temp_pow = temp_pow + abs(v' * X_f)^2;
        end
        P_point(r, c_idx) = temp_pow;
    end
end

% Plotting 2D Map (Like Friend's "Incoherent beamformer")
figure('Name', 'Fig 4: Point Source Map', 'Color', 'w');
imagesc(x_grid, y_grid, 10*log10(P_point/max(P_point(:))));
axis xy; colorbar; colormap jet; % Using JET colormap as requested
xlabel('Lateral Position X (m)'); ylabel('Range Y (m)');
title('Point Source Localization (Incoherent Bartlett)');
clim([-12 0]); % Adjust dynamic range for better contrast
hold on;
plot(mic_pos, zeros(size(mic_pos)), 'wv', 'MarkerFaceColor','k'); % Show Array

% Calculate Max Position
[~, max_idx] = max(P_point(:));
[row_max, col_max] = ind2sub(size(P_point), max_idx);
est_x = X_map(row_max, col_max);
est_y = Y_map(row_max, col_max);
plot(est_x, est_y, 'w+', 'MarkerSize', 15, 'LineWidth', 2);
text(est_x+0.1, est_y, sprintf('Est: [%.2f, %.2f]m', est_x, est_y), 'Color','w', 'FontWeight','bold');

fprintf('Estimated Location: X=%.2f m, Y=%.2f m\n', est_x, est_y);

%% PART 5: TIME GATING & EIGENVECTOR (FIXED)
% =========================================================================

% 1. Time Gating (Removing Reflections)
% FIX: t(:) ensures it is a column vector to match y_all dimensions
t_start = 0.003; % 3 ms
t_end   = 0.007; % 7 ms

% Create Window (Force column vector)
win = (t(:) >= t_start) & (t(:) <= t_end);

% Apply Window (Element-wise multiplication)
y_gated = y_all .* win; 

% Re-Run Plane Wave on Gated Data
Y_gated_spec = fft(y_gated, Nfft);
P_gated = zeros(length(theta_scan), 1);

for fi = 1:length(freqs)
    k = 2*pi*freqs(fi)/c;
    X_f = Y_gated_spec(idx_min+fi-1, :).'; 
    
    % Steering Vector
    V_plane = exp(-1j * k * mic_pos' * sind(theta_scan)); 
    
    % Bartlett Power
    P_gated = P_gated + abs(V_plane' * X_f).^2;
end
P_gated_dB = 10*log10(P_gated / max(P_gated));

figure('Name', 'Fig 5: Time Gating Effect', 'Color', 'w');
plot(theta_scan, P_plane_dB, 'b', 'LineWidth', 1.5); hold on;
plot(theta_scan, P_gated_dB, 'r--', 'LineWidth', 2);
legend('Original (With Reflections)', 'Gated (Direct Path Only)');
xlabel('Angle (deg)'); ylabel('Power (dB)');
title('Effect of Time Gating on Beamforming'); grid on;


% 2. Eigenvector Beamforming (Noise Filtering)
P_eigen = zeros(length(theta_scan), 1);

for fi = 1:length(freqs)
    k = 2*pi*freqs(fi)/c;
    X_f = Y_spectrum(idx_min+fi-1, :).'; 
    
    % CSM and Eigendecomposition
    R = X_f * X_f'; 
    [V, D] = eig(R);
    [d_vals, idx_sort] = sort(diag(D), 'descend');
    V = V(:, idx_sort);
    
    % Filter Noise (Project onto Signal Subspace - Top 1 Eigenvector)
    Vs = V(:, 1); 
    R_clean = Vs * Vs' * d_vals(1); 
    
    % Steering Vector
    V_plane = exp(-1j * k * mic_pos' * sind(theta_scan));
    
    % Calculate Power
    for ti=1:length(theta_scan)
         v = V_plane(:,ti);
         P_eigen(ti) = P_eigen(ti) + abs(v' * R_clean * v);
    end
end

figure('Name', 'Fig 6: Eigenvector Beamforming', 'Color', 'w');
plot(theta_scan, 10*log10(P_eigen/max(P_eigen)), 'k', 'LineWidth', 2);
xlabel('Angle (deg)'); ylabel('Power (dB)');
title('Eigenvector Beamforming (Noise Subspace Filtered)');
grid on; ylim([-30 0]);

fprintf('Processing Complete.\n');