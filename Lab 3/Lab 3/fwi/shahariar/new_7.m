%% Initial Condition Plot (Reproducing Figure 1)
clear all; close all; clc;

%% 1. Problem Setup and Parameters
fmax = 10;                      % Maximum frequency
nx = ceil(10 * fmax);           % Number of grid points (100x100 for fmax=10)
dom = domain([0 1 0 1], [nx nx]); % Square domain 1x1km

% Parameters
ns = 6;                         % Number of sources (Image shows 6 green circles)
nr = 6;                         % Number of receivers
ctr = 0.05;                     % Contrast
wpml = 0.1;                     % Width of PML
smoothness = 0.5;               % Smoothness of phantom

%% 2. Configuration: Sources, Receivers, and Model
% Sources (Green circles)
source_info.type = 'ring';
source_info.center = [0.5, 0.5];
source_info.radius = 0.35;
sources = sources_and_receivers(ns, source_info);

% Receivers (Red crosses)
receiver_info.type = 'ring';
receiver_info.center = [0.5, 0.5];
receiver_info.radius = 0.35;
receivers = sources_and_receivers(nr, receiver_info);

% True Velocity Model (Shepp-Logan Phantom)
phan = phantom2(dom, 0.35, ctr, smoothness);
c_true = 1 * dom.mat2vec(phan); % Vectorized true velocity

% Initial Background Model (Uniform Blue)
c0 = ones(nx^2, 1);

% PML (Perfectly Matched Layer) - Yellow Border
pml_info.type = 'pml';
pml_info.width = wpml;
[~, PML] = dom.window(pml_info);

% Window (Region of Interest) - Yellow Circle
window_info.type = 'disk';
window_info.center = [0.5, 0.5];
window_info.radius = 0.3;
[~, W] = dom.window(window_info);

%% 3. Generate Figure 1 (Four Subplots)
figure('Name', 'Figure 1: Initial Condition', 'Color', 'w');

% Subplot 1: Problem Setup
subplot(2, 2, 1);
dom.imagesc(c_true);                % Plot true velocity
hold on;
dom.plot(sources(1,:), sources(2,:), 'go', 'LineWidth', 1.5); % Sources
dom.plot(receivers(1,:), receivers(2,:), 'rx', 'LineWidth', 1.5); % Receivers
hold off;
axis square;
title('Problem setup');
colorbar;
c_vec = caxis; % Save color limits for consistency

% Subplot 2: Initial Background Velocity
subplot(2, 2, 2);
dom.imagesc(c0);                    % Plot uniform initial model
axis square;
title('Initial background velocity');
caxis(c_vec);                       % Apply same color scale

% Subplot 3: With PML
subplot(2, 2, 3);
% Visualize PML by adding it to the velocity matrix
imagesc(dom.vec2mat(c_true) + PML); 
hold on;
dom.plot(sources(1,:), sources(2,:), 'go');
dom.plot(receivers(1,:), receivers(2,:), 'rx');
hold off;
axis square;
axis xy; % Ensure correct orientation
title('With PML');
caxis([c_vec(1) c_vec(2)]); % Adjust color limits

% Subplot 4: With Window
subplot(2, 2, 4);
% Visualize Window (W) by adding it
% Note: Using flipud might be needed depending on your domain definition
imagesc(dom.vec2mat(c_true) + flipud(W)); 
hold on;
dom.plot(sources(1,:), sources(2,:), 'go');
dom.plot(receivers(1,:), receivers(2,:), 'rx');
hold off;
axis square;
axis xy;
title('With window');
caxis([c_vec(1) c_vec(2)+1]); % Adjust color limit to show the yellow mask