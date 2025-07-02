% CFANVC
classdef CFANVC < handle
    %#ok<*PROPLC>
    properties
        N, L, M, P, R, Lu, wd, wu, deltax, delta, alpha, eta, Q, qr_off, anc_on
        u, ev, uv, wuv, wev, bc, ac, bc0, ac0, udv, qr, uM, time
    end
    methods
        function o = CFANVC(param)
            % Read Parameters
            N = param{1}; L = param{2}; M = param{3}; P = param{4}; R = param{5}; Lu = param{6}; wd = param{7}; wu = param{8}; deltax = param{9}; delta = param{10}; alpha = param{11}; eta = param{12}; Q = param{13}; qr_off = param{14}; anc_on = param{15};
            u = 0;                   % anti-noise signal
            ev = zeros(M+N,1);       % residual noise buffer
            wuv = wu*ones(L+N,1);
            wev = ones(L,1);         % cost weigth of each residual noise sample
            bc = zeros(1,N);
            ac = zeros(1,N);
            bc0 = zeros(1,N);
            ac0 = zeros(1,N);
            udv = zeros(L+N,1);
            uv = zeros(max(M+N+1, 2*L), 1);        % anti-noise buffer
            qr = 0;
            uM = 0;
            time = 0;
            % Write properties
            o.N = N; o.L = L; o.M = M; o.P = P; o.R = R; o.Lu = Lu; o.wd = wd; o.wu = wu; o.deltax = deltax; o.delta = delta; o.alpha = alpha; o.eta = eta; o.Q = Q; o.qr_off = qr_off; o.anc_on = anc_on; o.u = u; o.ev = ev; o.uv = uv; o.wuv = wuv; o.wev = wev; o.bc = bc; o.ac = ac; o.bc0 = bc0; o.ac0 = ac0; o.udv = udv; o.qr = qr; o.uM = uM; o.time = time;
       end
       function u = step(o, e)
            % Read properties
            N = o.N; L = o.L; M = o.M; P = o.P; R = o.R; Lu = o.Lu; wd = o.wd; wu = o.wu; deltax = o.deltax; delta = o.delta; alpha = o.alpha; eta = o.eta; Q = o.Q; qr_off = o.qr_off; anc_on = o.anc_on; u = o.u; ev = o.ev; uv = o.uv; wuv = o.wuv; wev = o.wev; bc = o.bc; ac = o.ac; bc0 = o.bc0; ac0 = o.ac0; udv = o.udv; qr = o.qr; uM = o.uM; time = o.time;

            time = time + 1;
            ev = [e; ev(1:end-1)];
            uv = [u; uv(1:end-1)]; % simulation and algorithm

            if time >= anc_on
                if mod(time, Q) == 0

                    % Calculate the parameters
                    Em = ev((1:M)'+(1:N));
                    Um = uv((0:M-1)'+(1:N+1));
                    H = [-Em, Um];
                    % eb = H x
                    Rxx = H'*H;
                    Deltax = deltax*eye(2*N+1);
                    x = (Rxx+Deltax)\H'*ev(1:M);
                    ab = x(1:N); bb=x(N+1:end);

                    qvM = max(mean(reshape(abs(ev(1:M)-H*x).^2, P, [])));

                    Delta = delta*eye(2*N+1);
                    Sxx = qvM*inv(Rxx + Delta);

                    % good values for ah and bh
                    %a0 = conv(a,pz); a0 = a0(2:end); a0 = [a0; zeros(N-length(a0),1)];
                    %b0 = conv(b,pz); b0 = [b0; zeros(N+1-length(b0),1)];

                    Saa = Sxx(1:N, 1:N);
                    Saax = [zeros(N+1,1), [zeros(1,N); Saa]];
                    Sbb = Sxx(N+1:end, N+1:end);
                    %Sab = Sxx(1:N, N+1:end);
                    Sba = Sxx(N+1:end, 1:N);
                    Sbax = [zeros(N+1,1), Sba];
                    Qbb = bb*bb'+ Sbb;
                    Qba = bb*[1;ab]' + Sbax;
                    Qaa = [1;ab]*[1;ab]' + Saax;
                    %Qab = Qba';

                    Imx = diag([1./wev; zeros(N,1)]);
                    Raax = conv2(Imx, flip(flip(Qaa),2), 'valid');
                    RaaI = inv(Raax);
                    Raa = conv2(RaaI, Qaa);
                    Rbb = conv2(RaaI, Qbb) + diag(wuv);
                    Rba = conv2(RaaI, Qba);

                    wv = [zeros(L,1); ones(N,1)];
                    Wm = diag(wv);
                    Rmd = decomposition(Rbb + wd*Wm);
                    Bc = Rmd\Rba;
                    Ac = Rmd\Wm*wd;
                    bc0 = bc; ac0 = ac;
                    bc = Bc(L, L+1:end);
                    ac = - Ac(L, L+1:end);

                    edv = zeros(L+N, 1);
                    edv(L+1:end) = ev(1:N);
                    u0v = [zeros(L,1); uv(1:N)];
                    %udv0 = udv;
                    udv = Bc*edv + Ac*u0v;
                    xi = edv'*Raa*edv + udv'*Rbb*udv - 2*udv'*Rba*edv;
                    %xi0 = edv'*Raa*edv + u0v'*Rbb*u0v - 2*u0v'*Rba*edv;
                    qr = alpha*xi/trace(Rbb);
                    uM = eta*max(abs(uv(1:R)));
                end
                % u = udv0(L-mod(k,Q)-Q);
                u = bc0*ev(1:N) - ac0*uv(1:N);
                u = u + sqrt(qr)*randn;
                u = min(Lu, max(-Lu,u));
                u = min(uM, max(-uM,u));

                %log_qv(k,n_sim) = qvM;
                %log_xi(k,n_sim) = xi;
                %log_xi0(k,n_sim) = xi0;
            else
                u = sqrt(qr_off)*randn;
            end
            % Write properties
            o.N = N; o.L = L; o.M = M; o.P = P; o.R = R; o.Lu = Lu; o.wd = wd; o.wu = wu; o.deltax = deltax; o.delta = delta; o.alpha = alpha; o.eta = eta; o.Q = Q; o.qr_off = qr_off; o.anc_on = anc_on; o.u = u; o.ev = ev; o.uv = uv; o.wuv = wuv; o.wev = wev; o.bc = bc; o.ac = ac; o.bc0 = bc0; o.ac0 = ac0; o.udv = udv; o.qr = qr; o.uM = uM; o.time = time;
       end
    end
end
