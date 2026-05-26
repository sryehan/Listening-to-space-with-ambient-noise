%% S. Rakotonarivo & Student
% FINAL REPORT GENERATOR
% This code generates the acoustic images for the Lab Report.
% It performs Delay-and-Sum Beamforming and Differential Imaging.

clearvars; close all;

%------------------------------
%% 1. SETTINGS
%------------------------------
% --- UPDATE PATHS ---
Generaldir = 'D:\france semester\3_Acoustic Imaging\lab\Shahariar\DATA\';
FileDir = strcat(Generaldir, 'Data_MAT\');

FileCyl = 'CC50_Cylinder.mat';
FileEmp = 'CC50_Empty.mat'; % Jodi eta na thake, code automatic handle korbe

% Imaging Grid Settings (Resolution 1cm)
step = 0.01; 
fov = 0.4;   % -0.4m to +0.4m
x_vec = -fov:step:fov;
y_vec = -fov:step:fov;
[GridX, GridY] = meshgrid(x_vec, y_vec);

%------------------------------
%% 2. LOAD DATA
%------------------------------
disp('Loading Cylinder Data...');
if ~isfile(strcat(FileDir, FileCyl))
    error(['Cylinder File not found! Check path: ' FileDir]);
end
load(strcat(FileDir, FileCyl)); % Loads I
I_Cyl = I; 

% Load Empty Data (If exists)
hasEmpty = false;
if isfile(strcat(FileDir, FileEmp))
    disp('Loading Empty Data (For Comparison)...');
    load(strcat(FileDir, FileEmp));
    I_Emp = I;
    hasEmpty = true;
else
    warning('Empty Data not found. Differential Image will be skipped.');
end

%------------------------------
%% 3. SETUP GEOMETRY (Mic Positions)
%------------------------------
% Inner Ring (17-32): R=0.26m | Outer Ring (1-16): R=0.27m
Nbmic = I_Cyl.Nbmic;
theta = linspace(0, 2*pi, Nbmic+1); theta(end) = [];
R_in = 0.26; R_out = 0.27;

X_L = R_out * cos(theta); Y_L = R_out * sin(theta); % Outer
X_R = R_in * cos(theta);  Y_R = R_in * sin(theta);  % Inner
XRec = [X_L, X_R]; 
YRec = [Y_L, Y_R];

c0 = 343;
fs = I_Cyl.fs;
t_Gt = I_Cyl.t_Gt;

%------------------------------
%% 4. BEAMFORMING FUNCTION (Fast Loop)
%------------------------------
disp('Running Beamforming Algorithm...');

compute_image = @(Gt) run_das(Gt, fs, c0, XRec, YRec, t_Gt, x_vec, y_vec);

% Compute Images
tic;
Img_Cyl = compute_image(I_Cyl.Gt);
t_end = toc;
fprintf('Cylinder Image computed in %.2f seconds.\n', t_end);

if hasEmpty
    Img_Emp = compute_image(I_Emp.Gt);
end

%------------------------------
%% 5. PLOTTING FOR REPORT
%------------------------------

% --- FIGURE 1: GREEN'S FUNCTION (Physics) ---
figure('Name', 'Fig 1: Signal Physics', 'Color', 'w', 'Position', [100, 100, 600, 400]);
plot(t_Gt*1000, squeeze(I_Cyl.Gt(1, 9, :)), 'LineWidth', 1.5);
title('Reconstructed Green''s Function (Mic 1 vs Mic 9)');
xlabel('Time (ms)'); ylabel('Amplitude (a.u.)');
grid on; xlim([0 5]);
subtitle('Peak represents direct flight time between sensors');
% Save for Report
% saveas(gcf, 'Fig1_GreensFunction.png');

% --- FIGURE 2: THE ACOUSTIC IMAGE (Main Result) ---
figure('Name', 'Fig 2: Cylinder Detection', 'Color', 'w', 'Position', [150, 150, 700, 600]);
imagesc(x_vec, y_vec, Img_Cyl);
set(gca, 'YDir', 'normal');
colormap('jet'); colorbar; axis equal;
title('Passive Acoustic Image (With Cylinder)');
xlabel('X (m)'); ylabel('Y (m)');
hold on;
% Draw Cylinder Circle (White)
viscircles([0, 0], 0.045, 'Color', 'w', 'LineStyle', '--');
text(0.15, 0.35, 'Bright spot = Object', 'Color', 'white', 'BackgroundColor', 'black');

% --- FIGURE 3: DIFFERENTIAL IMAGE (Advanced Analysis) ---
if hasEmpty
    figure('Name', 'Fig 3: Differential Imaging', 'Color', 'w', 'Position', [200, 200, 700, 600]);
    
    % Subtraction (Removes background noise)
    Img_Diff = abs(Img_Cyl - Img_Emp);
    
    % Normalize for better look
    Img_Diff = Img_Diff / max(Img_Diff(:));
    
    imagesc(x_vec, y_vec, Img_Diff);
    set(gca, 'YDir', 'normal');
    colormap('jet'); colorbar; axis equal;
    title('Differential Image (Cylinder - Empty)');
    subtitle('Background noise removed using subtraction');
    xlabel('X (m)'); ylabel('Y (m)');
    
    % Automatic Detection
    [~, idx] = max(Img_Diff(:));
    [r, c] = ind2sub(size(Img_Diff), idx);
    x_peak = x_vec(c); y_peak = y_vec(r);
    
    hold on;
    viscircles([0, 0], 0.045, 'Color', 'w', 'LineStyle', '-');
    plot(x_peak, y_peak, 'r+', 'MarkerSize', 20, 'LineWidth', 2);
    
    legend('True Position', 'Detected Center');
    
    fprintf('\n----------------------------\n');
    fprintf('DETECTED COORDINATES: X=%.2f m, Y=%.2f m\n', x_peak, y_peak);
    fprintf('----------------------------\n');
end

%------------------------------
%% HELPER FUNCTION (Local)
%------------------------------
function Image = run_das(Gt, fs, c0, XRec, YRec, t_Gt, x_vec, y_vec)
    Image = zeros(length(y_vec), length(x_vec));
    % Pre-calculate distances for speed
    for ix = 1:length(x_vec)
        for iy = 1:length(y_vec)
            d_all = sqrt((x_vec(ix) - XRec).^2 + (y_vec(iy) - YRec).^2);
            val = 0;
            for i = 1:32
                for j = 1:32
                    tau = (d_all(i) + d_all(j)) / c0;
                    idx = round(tau * fs) + 1;
                    if idx > 0 && idx < length(t_Gt)
                        % Summing envelope (abs)
                        val = val + abs(Gt(i, j, idx));
                    end
                end
            end
            Image(iy, ix) = val;
        end
    end
end