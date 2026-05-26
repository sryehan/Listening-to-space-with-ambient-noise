%% Exercise 1: Forward Modeling
clear all; close all; clc;

% 1. Setup Domain (Same as your main script)
fmax = 10;
nx = ceil(10*fmax);
dom = domain([0 1 0 1], [nx nx]);

% 2. Setup Velocity Model (Perturbed Model for Ex 1.1)
c_background = 1 * ones(nx^2, 1);       % Uniform background
c_perturbed  = c_background * 1.1;      % +10% perturbation
m_background = 1./c_background.^2;      % Convert to squared slowness
m_perturbed  = 1./c_perturbed.^2;

% 3. Define a single source in the middle
source_loc = [0.5; 0.5];
b = generate_sources(dom, source_loc);

% 4. Solve Forward Problem at specific Frequency
freq = 5; % Change this for Ex 1.2 (e.g., 3Hz, 5Hz, 10Hz)

% Solve for Background
fprintf('Solving for Background...\n');
H_back = helmholtz_2d(m_background, freq, dom);
A_back = invertA(H_back, 1);
u_back = A_back.apply(b);

% Solve for Perturbed
fprintf('Solving for Perturbed...\n');
H_pert = helmholtz_2d(m_perturbed, freq, dom);
A_pert = invertA(H_pert, 1);
u_pert = A_pert.apply(b);

% 5. Visualize Wavefields (Real part)
figure;
subplot(1,3,1);
dom.imagesc(real(u_back)); title(['Background Wavefield (f=' num2str(freq) ')']);
axis square; colorbar;

subplot(1,3,2);
dom.imagesc(real(u_pert)); title('Perturbed Wavefield');
axis square; colorbar;

subplot(1,3,3);
dom.imagesc(real(u_pert - u_back)); title('Difference (Scattered Field)');
axis square; colorbar;
%% 6. Visualize 1D Cross-section (Slice) to see Phase Shift
% Convert the vector u to a 2D matrix
U_back_2D = reshape(real(u_back), nx, nx);
U_pert_2D = reshape(real(u_pert), nx, nx);

% Select the middle row index
mid_row = ceil(nx/2);

% Extract the slice
slice_back = U_back_2D(:, mid_row); % Slicing through the middle
slice_pert = U_pert_2D(:, mid_row);

% Plot comparison
figure('Name', '1D Cross-Section Comparison');
plot(slice_back, 'b', 'LineWidth', 1.5); hold on;
plot(slice_pert, 'r--', 'LineWidth', 1.5);
legend('Background (c=1)', 'Perturbed (c=1.1)');
xlabel('Grid Points (x)');
ylabel('Amplitude (Real part of p)');
title(['Wavefield Slice at f = ' num2str(freq) ' Hz']);
grid on;


%% 7. Time Domain Reconstruction (Using the Professor's Hint)
% Hint: f(t) = A * cos(w*t + phi)

% 1. Select a receiver position (e.g., a point near the middle of the domain or far from the source)
% We are taking a point in the middle row (where we took the slice)
rec_idx = mid_row + (mid_row - 1)*nx; % This is an index near the center of the domain

% 2. Extract the Amplitude (A) and Phase (phi)
% For the background model
A_back = abs(u_back(rec_idx));
phi_back = angle(u_back(rec_idx));

% For the perturbed model
A_pert = abs(u_pert(rec_idx));
phi_pert = angle(u_pert(rec_idx));

% 3. Create a time vector (to see a few cycles)
T_period = 1/freq;
t = linspace(0, 3*T_period, 100); % Show 3 periods
omega = 2 * pi * freq;

% 4. Create the time-domain signal according to the hint
p_t_back = A_back * cos(omega * t + phi_back);
p_t_pert = A_pert * cos(omega * t + phi_pert);

% 5. Plot the signals
figure('Name', 'Time Domain Signal Reconstruction');
plot(t, p_t_back, 'b', 'LineWidth', 2); hold on;
plot(t, p_t_pert, 'r--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Pressure Amplitude p(t)');
title(['Time Domain Signal at Receiver (f = ' num2str(freq) ' Hz)']);
legend('Background (c=1.0)', 'Perturbed (c=1.1)');
grid on;