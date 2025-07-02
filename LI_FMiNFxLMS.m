classdef LI_FMiNFxLMS < handle
    %#ok<*PROPLC>
    properties
        N, Nx, mu, mus, tau, lambda_w, alpha, beta, gamma, eta, lambda, q, M, R, log_vars, qr_off, anc_on
        wv, sh, lfsr, xv, xv1, uv, rv, time, qe, qd, qrx, qdx, gu, qf, qrT, log_values
    end
    methods
        function o = LI_FMiNFxLMS(param)
            % Read Parameters
            N = param{1}; Nx = param{2}; mu = param{3}; mus = param{4}; tau = param{5}; lambda_w = param{6}; alpha = param{7}; beta = param{8}; gamma = param{9}; eta = param{10}; lambda = param{11}; q = param{12}; M = param{13}; R = param{14}; log_vars = param{15}; qr_off = param{16}; anc_on = param{17};

            wv = zeros(N, 1);
            sh = zeros(N, 1);

            lfsr_b = ceil(log2(4*N));
            lfsr = rand(1, lfsr_b) > 0.5;

            xv = zeros(N,1);
            xv1 = zeros(N,1);
            uv = zeros(N,1);
            rv = zeros(N,1);

            time = 0;
            gu = 1;
            qe = 0;
            qd = 0;
            qf = 0;
            qrT = 0;
            log_values = nan;

            qrx = nan;
            qdx = nan;

            % Write properties
            o.N = N; o.Nx = Nx; o.mu = mu; o.mus = mus; o.tau = tau; o.lambda_w = lambda_w; o.alpha = alpha; o.beta = beta; o.gamma = gamma; o.eta = eta; o.lambda = lambda; o.q = q; o.M = M; o.R = R; o.log_vars = log_vars; o.qr_off = qr_off; o.anc_on = anc_on; o.wv = wv; o.sh = sh; o.lfsr = lfsr; o.xv = xv; o.xv1 = xv1; o.uv = uv; o.rv = rv; o.time = time; o.qe = qe; o.qd = qd; o.qrx = qrx; o.qdx = qdx; o.gu = gu; o.qf = qf; o.qrT = qrT; o.log_values = log_values;
        end
        function u = step(o, e)
            % Read properties
            N = o.N; Nx = o.Nx; mu = o.mu; mus = o.mus; tau = o.tau; lambda_w = o.lambda_w; alpha = o.alpha; beta = o.beta; gamma = o.gamma; eta = o.eta; lambda = o.lambda; q = o.q; M = o.M; R = o.R; log_vars = o.log_vars; qr_off = o.qr_off; anc_on = o.anc_on; wv = o.wv; sh = o.sh; lfsr = o.lfsr; xv = o.xv; xv1 = o.xv1; uv = o.uv; rv = o.rv; time = o.time; qe = o.qe; qd = o.qd; qrx = o.qrx; qdx = o.qdx; gu = o.gu; qf = o.qf; qrT = o.qrT; log_values = o.log_values;

            time = time + 1;

            uv = [0; uv(1:end-1)];
            lfsr = [lfsr(2:end), xor(lfsr(1), lfsr(2))];
            %norm_sh_sq = max(norm(sh)^2, 1/R^2); % worst case norm(sh)^2 is 1/R^2
            norm_sh_sq = norm(sh)^2;
            qr = alpha*qf/norm_sh_sq;

            % reference calculation
            dh = e - sh'*uv;
            x = dh; % x = d;
            xv = [x ; xv(1:end-1)];
            x1 = sh'*xv;

            qe = lambda*qe + (1-lambda)*e^2;
            qd = lambda*qd + (1-lambda)*dh^2;

            if time >= anc_on
                qrx = gamma*qrx + (1-gamma)*min(2*qrx, qr);    % power limit
                qdx = gamma*qdx + (1-gamma)*min(2*qdx, qd);    % power limit

                % calculate controler
                xv1 = [x1; xv1(1:end-1)];
                mu1 = mu/(xv1'*xv1 + q);
                eh = x + wv'*xv1;  % eh != e = x + sh u
                wv = lambda_w*wv - mu1*xv1*eh;

                u0 = wv'*xv;

                if qe > beta*qdx + norm_sh_sq*qrx
                    gu = 0.01;
                    qrT = 0;
                else
                    gu = min(gu*eta, 1);
                end

                r = sqrt(qrx)*(2*lfsr(1)-1);
                u1 = gu * u0;
            else
                qrx = qr_off;   % power limit
                qdx = qd;       % power limit
                u1 = 0;
                r = sqrt(qr_off)*(2*lfsr(1)-1);
            end

            u = max(-M, min(M, u1 + r));

            rv = [r; rv(1:end-1)];
            rvx = rv(2:end);

            f = e - sh'*rv;
            qf = lambda*qf + (1-lambda)*f^2;

            qrT = qrT + r^2;
            qs = tau^2/R^2;
            mu10 = max(1/qrT, 2*qs/(N*qf));
            mu1 = min(mu10, mus/(rvx'*rvx + q));

            if time >= Nx
                sh(2:end) = sh(2:end) + mu1*rvx*f;
            end

            uv(1) = u;

            log_values = eval(log_vars);
            % Write properties
            o.N = N; o.Nx = Nx; o.mu = mu; o.mus = mus; o.tau = tau; o.lambda_w = lambda_w; o.alpha = alpha; o.beta = beta; o.gamma = gamma; o.eta = eta; o.lambda = lambda; o.q = q; o.M = M; o.R = R; o.log_vars = log_vars; o.qr_off = qr_off; o.anc_on = anc_on; o.wv = wv; o.sh = sh; o.lfsr = lfsr; o.xv = xv; o.xv1 = xv1; o.uv = uv; o.rv = rv; o.time = time; o.qe = qe; o.qd = qd; o.qrx = qrx; o.qdx = qdx; o.gu = gu; o.qf = qf; o.qrT = qrT; o.log_values = log_values;
        end
    end
end

