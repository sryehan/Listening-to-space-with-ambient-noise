function [m,out]=adjoint_state_2d(dom,all_freqs,sources,receivers,window_info,c_true,c_0,sigma,maxit)
    % [Standard Header comments...]
    N = dom.N;
    [win_inds,~,W] = dom.window(window_info);
    nf = length(all_freqs);
    ns = size(sources,2);
    nw = length(win_inds);
    r_ind = dom.loc2ind(receivers);
    I = speye(N);
    E_rec = I(:,r_ind);
    m_true = 1./c_true.^2;
    m0 = m_true;
    m0(win_inds) = 1./c_0(win_inds).^2;
    b = generate_sources(dom,sources);
    opt.Niter = maxit;

    % Generate Data
    d = generate_seismic_data(dom,sources,receivers,all_freqs,m_true,sigma);

    % Run LBFGS
    [dm,out] = lbfgs(@adjoint_state_gradient,zeros(nw,1),opt);
    m = m0+W*dm;

    function [J,DJ] = adjoint_state_gradient(dm)
        fprintf('Solving Helmholtz, frequency: ')
        DJ = zeros(nw,1);
        J = 0;
        for ii = 1:nf
            fprintf('%d, ',ii)
            % Compute m for current iteration
            m_curr = m0 + W*dm;
            A = invertA(helmholtz_2d(m_curr,all_freqs(ii),dom),1);
            for jj = 1:ns
                u = A.apply(b(:,jj));               % Forward
                res = E_rec' * (u(r_ind)-d(:,ii,jj)); % Residual at receivers
                bq = E_rec * res;                   % Pad residual to full domain
                q = A.applyt(bq);                   % Adjoint

                % Gradient Update (Sum over freqs and sources)
                omega = 2*pi*all_freqs(ii);
                % Gradient is Real part of u * conj(q) restricted to window
                grad_update = omega^2 * real(u(win_inds).*conj(q(win_inds)));
                DJ = DJ + grad_update;

                % Objective Update
                J = J + 0.5 * norm(u(r_ind)-d(:,ii,jj))^2;
            end
        end
        fprintf('\n')
    end
end