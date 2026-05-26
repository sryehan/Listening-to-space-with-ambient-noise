%% Exercise 1.2 / 2.1: Frequency Sensitivity Kernels (Figure 5 & 6)
% CORRECTED VERSION
clear all; close all; clc;

% --- 1. SETUP PARAMETERS ---
fmax = 10;
nx = ceil(10 * fmax);
dom = domain([0 1 0 1], [nx nx]);

% Source & Receiver Setup (Ring Geometry)
source_info.type = 'ring';
source_info.center = [0.5, 0.5];
source_info.radius = 0.35;
sources = sources_and_receivers(6, source_info); % 6 Sources

receiver_info.type = 'ring';
receiver_info.center = [0.5, 0.5];
receiver_info.radius = 0.35;
receivers = sources_and_receivers(16, receiver_info); % 16 Receivers

% Frequencies to test
test_freqs = [1, 3, 6, 10]; 

% --- 2. MODELS ---
% True Model (Phantom)
ctr = 0.05;
phan = phantom2(dom, 0.35, ctr, 0.5);
c_true = 1 * dom.mat2vec(phan);
m_true = 1./c_true.^2;

% Initial Model (Homogeneous)
c0 = ones(nx^2, 1);
m0 = 1./c0.^2;

% Window (To mask the result like the image)
window_info.type = 'disk';
window_info.center = [0.5, 0.5];
window_info.radius = 0.3;
[win_inds, W] = dom.window(window_info);

% --- 3. LOOP OVER FREQUENCIES ---
figure('Name', 'Figure 5 & 6: Frequency Kernels', 'Color', 'w');
colormap(parula); 

for i = 1:length(test_freqs)
    f = test_freqs(i);
    fprintf('Computing Kernel for f = %d Hz...\n', f);
    
    % A. Generate Observed Data (d_obs) from True Model
    H_true = helmholtz_2d(m_true, f, dom);
    A_true = invertA(H_true, 1);
    b_src = generate_sources(dom, sources);
    u_true = A_true.apply(b_src);
    
    % B. Generate Modeled Data (d_calc) from Initial Model
    H_init = helmholtz_2d(m0, f, dom);
    A_init = invertA(H_init, 1);
    u_init = A_init.apply(b_src);
    
    % C. Compute Residual & Adjoint Source
    r_ind = dom.loc2ind(receivers);
    
    % --- FIX IS HERE: Added (:, :) to ensure matrix dimensions match ---
    % u_init is (Nx6), r_ind is (16x1). u_init(r_ind, :) gives (16x6).
    resid = u_init(r_ind, :) - u_true(r_ind, :); 
    
    % Adjoint Source vector
    bq = zeros(dom.N, size(sources,2)); % Size: N x 6
    bq(r_ind, :) = resid;               % Assign 16x6 residual to receiver locations
    
    % D. Adjoint Solve (Backpropagation)
    q = A_init.applyt(bq);
    
    % E. Compute Gradient/Kernel: K = Re[w^2 * u * q*]
    omega = 2 * pi * f;
    % Sum over all 6 sources (columns) to get one single image
    Kernel = sum(real(omega^2 * u_init .* conj(q)), 2);
    
    % F. PLOTTING
    subplot(2, 2, i);
    
    % Apply Window Mask for visualization
    K_plot = zeros(size(Kernel));
    K_plot(win_inds) = Kernel(win_inds);
    
    dom.imagesc(K_plot);
    title(['Kernel ' num2str(f) ' Hz']);
    axis square; axis off;
    
    % Adjust contrast
    clim = max(abs(K_plot(:)));
    if clim > 0
        caxis([-clim clim] * 0.8);
    end
end