%% Exercise 3.2: Multi-Scale Strategy Comparison
clear all; close all; clc;

% --- 1. SETUP ---
fmax = 10;
nx = ceil(10 * fmax);
dom = domain([0 1 0 1], [nx nx]);

% Source/Receiver Setup
source_info.type = 'ring'; source_info.radius = 0.35; source_info.center = [0.5, 0.5];
sources = sources_and_receivers(8, source_info); % 8 sources
receiver_info.type = 'ring'; receiver_info.radius = 0.35; receiver_info.center = [0.5, 0.5];
receivers = sources_and_receivers(16, receiver_info); % 16 receivers

% Define Frequency Bands
freqs_low = [1, 2, 3];        % Low frequencies (f1 to f2)
freqs_high = [8, 9, 10];      % High frequencies (f3 to f4)

% Models
ctr = 0.05;
phan = phantom2(dom, 0.35, ctr, 0.5);
c_true = 1 * dom.mat2vec(phan);     % True Model
c0 = ones(nx^2, 1);                 % Smooth Initial Model (Homogeneous)

window_info.type = 'disk'; window_info.radius = 0.3; window_info.center = [0.5, 0.5];
sigma = 0; maxit = 15;

% --- 2. RUN INVERSIONS ---

% Case 1: Low Frequencies Only (Start from c0)
fprintf('Running Case 1: Low Frequencies Only...\n');
[m_low, ~] = adjoint_state_2d(dom, freqs_low, sources, receivers, window_info, c_true, c0, sigma, maxit);

% Case 2: High Frequencies Only (Start from c0) -> This should fail!
fprintf('Running Case 2: High Frequencies Only...\n');
[m_high_only, ~] = adjoint_state_2d(dom, freqs_high, sources, receivers, window_info, c_true, c0, sigma, maxit);

% Case 3: Progressive Strategy (Start High freq using result of Low freq)
fprintf('Running Case 3: Progressive Strategy...\n');
% Convert m_low back to velocity to use as initial model for high freq
c_start_progressive = sqrt(1 ./ m_low); 
[m_progressive, ~] = adjoint_state_2d(dom, freqs_high, sources, receivers, window_info, c_true, c_start_progressive, sigma, maxit);

% --- 3. PLOT COMPARISON ---
figure('Name', 'Exercise 3.2: Multi-Scale Comparison', 'Color', 'w', 'Position', [100 100 1200 400]);

% Plot 1: Low Only
subplot(1, 3, 1);
dom.imagesc(sqrt(1./m_low));
title('1. Low Freq Only [1-3 Hz]');
axis square; colorbar;
xlabel('Result: Low Resolution, Stable');

% Plot 2: High Only
subplot(1, 3, 2);
dom.imagesc(sqrt(1./m_high_only));
title('2. High Freq Only [8-10 Hz]');
axis square; colorbar;
xlabel('Result: Cycle Skipping (Failed)');

% Plot 3: Progressive
subplot(1, 3, 3);
dom.imagesc(sqrt(1./m_progressive));
title('3. Progressive Strategy');
axis square; colorbar;
xlabel('Result: High Resolution, Accurate');

colormap(jet);