% Exercise 1.1: Forward Modeling
% Assumes 'dom', 'sources', 'receivers', 'c_true' are already defined (see Ex 3.1 setup)

% 1. Modify Velocity Model (+10%)
c_perturbed = c_true * 1.1; 
m_perturbed = 1 ./ c_perturbed.^2;

% 2. Solve Forward Problem for a specific frequency (e.g., 5 Hz)
freq = 5; 
% Generate source vector
b = generate_sources(dom, sources); 
% Solve Helmholtz: (nabla^2 + w^2*m) u = -s
u = invertA(helmholtz_2d(m_perturbed, freq, dom), 1).apply(b(:,1)); % Solve for 1st source

% 3. Extract Data at Receivers
r_ind = dom.loc2ind(receivers);
d_obs = u(r_ind);

% 4. Time Domain Reconstruction (Single Frequency Approximation)
% f(t) = |d| * cos(w*t + phase)
A = abs(d_obs);
phi = angle(d_obs);
w = 2 * pi * freq;
t = linspace(0, 1/freq * 5, 200); % Plot 5 periods

% Select 1st receiver to plot
rec_idx = 1;
signal = A(rec_idx) * cos(w * t + phi(rec_idx));

figure;
plot(t, signal);
xlabel('Time (s)'); ylabel('Amplitude');
title(['Reconstructed Signal at Receiver ' num2str(rec_idx) ' (c + 10%)']);
grid on;