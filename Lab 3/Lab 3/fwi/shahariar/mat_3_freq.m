%% Exercise 1.2: Frequency Sensitivity Analysis
% Goal: Compute synthetic data for three frequencies f1 > f2 > f3
clear; clc; close all;

% --- 1. Setup Parameters ---
% Define frequencies (High, Medium, Low)
f1 = 15; % High Hz
f2 = 8;  % Medium Hz
f3 = 3;  % Low Hz
freqs = [f1, f2, f3];

% Domain setup (1km x 1km)
nx = 150; % Grid points
dom = domain([0 1 0 1], [nx nx]);

% Constant velocity model (c = 1 km/s)
c0 = 1.0 * ones(nx^2, 1);
m0 = 1 ./ c0.^2;

% Source in the center
source_info.type = 'point';
source_info.loc = [0.5; 0.5];
b = generate_sources(dom, [0.5; 0.5]); 

% --- 2. Compute and Plot ---
figure('Name', 'Ex 1.2 Frequency Analysis', 'Color', 'w');
colormap(jet);

for i = 1:length(freqs)
    f = freqs(i);
    
    % Solve Helmholtz Equation: (nabla^2 + w^2*m) u = -s
    % 'invertA' computes the inverse/solver for the Helmholtz operator
    u = invertA(helmholtz_2d(m0, f, dom), 1).apply(b);
    
    % Plotting
    subplot(1, 3, i);
    dom.imagesc(real(u)); % Visualize the real part of the wavefield
    axis square;
    
    % Calculate theoretical wavelength lambda = c / f
    lambda = 1.0 / f;
    
    title(sprintf('f_{%d} = %d Hz\n\\lambda = %.3f km', i, f, lambda));
    xlabel('x (km)'); 
    if i==1, ylabel('z (km)'); end
end

sgtitle('Exercise 1.2: Wavefield vs Frequency');