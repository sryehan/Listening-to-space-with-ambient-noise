%% Exercise 3.1: Full Inversion & Cost Function Plot
clear; close all; clc;

% --- 1. SETUP PARAMETERS ---
fmax = 10;                      % Max frequency
nx = ceil(10 * fmax);           % Grid points
dom = domain([0 1 0 1], [nx nx]);

% Source & Receiver Setup
source_info.type = 'ring';
source_info.radius = 0.35;
source_info.center = [0.5, 0.5];
sources = sources_and_receivers(6, source_info); % 6 Sources

receiver_info.type = 'ring';
receiver_info.radius = 0.35;
receiver_info.center = [0.5, 0.5];
receivers = sources_and_receivers(5, receiver_info); % 5 Receivers

% Frequencies for Inversion
freqs = linspace(1, fmax, 10);

% --- 2. CREATE MODELS ---
% True Model (Shepp-Logan Phantom)
ctr = 0.05; % Contrast
phan = phantom2(dom, 0.35, ctr, 0.5);
c_true = 1 * dom.mat2vec(phan);

% Initial Model (Smooth/Homogeneous)
c0 = ones(nx^2, 1); 

% Window & Noise Info
window_info.type = 'disk';
window_info.center = [0.5, 0.5];
window_info.radius = 0.3;
sigma = 0;      % No noise (Perfect data)
maxit = 20;     % Max iterations

% --- 3. RUN INVERSION (Generates 'out') ---
fprintf('Running FWI... Please wait.\n');
% এই লাইনটিই 'out' ভেরিয়েবল তৈরি করে
[m_inv, out] = adjoint_state_2d(dom, freqs, sources, receivers, ...
                                window_info, c_true, c0, sigma, maxit);

% --- 4. PLOT RESULTS & COST FUNCTION ---
figure('Name', 'Exercise 3.1 Results', 'Color', 'w');

% Plot a) Inverted Model
subplot(1, 2, 1);
dom.imagesc(sqrt(1./m_inv));
title('a) Inverted Velocity Model');
axis square; colorbar;
colormap(jet);

% Plot b) Cost Function Evolution
subplot(1, 2, 2);
if isfield(out, 'J') && ~isempty(out.J)
    % Filter out zeros from initialization
    J_history = out.J(out.J > 0); 
    iterations = 1:length(J_history);
    
    semilogy(iterations, J_history, '-or', 'LineWidth', 2, 'MarkerFaceColor', 'r');
    xlabel('Iterations');
    ylabel('Cost Function J(m)');
    title('b) Cost Function Evolution');
    grid on; axis square;
    
    % Add initial vs final cost text
    legend(['Final Cost: ' num2str(J_history(end), '%.2e')]);
else
    text(0.5, 0.5, 'No Cost Data Available', 'HorizontalAlignment', 'center');
end