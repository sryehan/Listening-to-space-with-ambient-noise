%% Exercise 3.1 (b): Plot Evolution of Cost Function
% এই কোডটি রান করার আগে নিশ্চিত করুন যে আপনার workspace-এ 'out' ভেরিয়েবলটি আছে
% [m, out] = adjoint_state_2d(...) কমান্ডটি আগে রান করা থাকতে হবে।

figure('Name', 'Exercise 3.1b: Cost Function Evolution', 'Color', 'w');

% 1. Extract non-zero cost values (LBFGS আউটপুটে অনেক সময় বাড়তি শূন্য থাকে)
J_history = out.J(out.J > 0); 
iters = 1:length(J_history);

% 2. Plot on Semi-Log Scale (Standard for convergence plots)
semilogy(iters, J_history, '-o', ...
    'LineWidth', 2, ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', 'b', ...
    'Color', 'b');

% 3. Labels and Formatting
xlabel('Iteration Number', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Cost Function J(m)', 'FontSize', 12, 'FontWeight', 'bold');
title('Convergence History: Cost Function Evolution', 'FontSize', 14);
grid on;
axis square;

% 4. Add Normalize Cost annotation (Optional)
% Shows how much the error decreased relative to the start
J_reduction = (J_history(1) - J_history(end)) / J_history(1) * 100;
text_str = sprintf('Total Reduction: %.2f%%', J_reduction);
text(iters(end)*0.6, J_history(1)*0.5, text_str, ...
     'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k');

fprintf('Initial Cost: %.4e\n', J_history(1));
fprintf('Final Cost:   %.4e\n', J_history(end));