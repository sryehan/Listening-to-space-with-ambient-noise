%% Exercise 2: Sensitivity Kernels
clear all; close all; clc;

% Setup
fmax = 10;
nx = ceil(10*fmax);
dom = domain([0 1 0 1], [nx nx]);
c0 = ones(nx^2, 1);
m0 = 1./c0.^2;
freq = 5; % Chosen frequency

% Define 1 Source and 1 Receiver
src_loc = [0.2; 0.5]; % Left side
rec_loc = [0.8; 0.5]; % Right side

% Generate Source Vector (b) and Receiver Vector (bq)
b = generate_sources(dom, src_loc);
bq = generate_sources(dom, rec_loc); % Receiver acts as a source for adjoint

% Solve Forward Wavefield (u)
fprintf('Solving Forward...\n');
H = helmholtz_2d(m0, freq, dom);
A = invertA(H, 1);
u = A.apply(b);

% Solve Adjoint Wavefield (q) - Backpropagate from receiver
fprintf('Solving Adjoint...\n');
q = A.applyt(bq); % applyt is transpose (adjoint)

% Calculate Sensitivity Kernel: K = w^2 * u * q
omega = 2 * pi * freq;
K = omega^2 * real(u .* conj(q)); % Eq (7) in manual

% Plot
figure;
dom.imagesc(K);
hold on; plot(src_loc(1), src_loc(2), 'go'); plot(rec_loc(1), rec_loc(2), 'rx');
title('Sensitivity Kernel (Banana-Doughnut)');
axis square; colorbar;