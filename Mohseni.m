% Mohseni
classdef Mohseni < handle
    %#ok<*PROPLC>
    properties
        N, L, tau_n, tau_d, eta, lambda_cc, Eid, Eccy, Eccu, Ct, Lu, qr_off, anc_on
        ah, bh, xh, Theta_id, Pid, Theta_cc, Pcc, Gammah, uh, lambda_id, zidv, uv, yv, time
    end
    methods
        function o = Mohseni(param)
            % Read Parameters
            N = param{1}; L = param{2}; tau_n = param{3}; tau_d = param{4}; eta = param{5}; lambda_cc = param{6}; Eid = param{7}; Eccy = param{8}; Eccu = param{9}; Ct = param{10}; Lu = param{11}; qr_off = param{12}; anc_on = param{13};

            ah = zeros(N+1,1); ah(1) = 1;
            bh = zeros(N+1,1); bh(1) = 1;
            xh = zeros(N,1);
            Theta_id = [ah(2:end);bh];
            Pid = eye(2*N+1);
            Theta_cc = zeros(L,1);
            Pcc = eye(L);
            Gammah = zeros(L,N);
            uh = zeros(L,1);
            lambda_id = 1;
            zidv = zeros(tau_d,1);
            yv = zeros(L,1);
            uv = zeros(L,1);
            time = 0;

            % Write properties
            o.N = N; o.L = L; o.tau_n = tau_n; o.tau_d = tau_d; o.eta = eta; o.lambda_cc = lambda_cc; o.Eid = Eid; o.Eccy = Eccy; o.Eccu = Eccu; o.Ct = Ct; o.Lu = Lu; o.qr_off = qr_off; o.anc_on = anc_on; o.ah = ah; o.bh = bh; o.xh = xh; o.Theta_id = Theta_id; o.Pid = Pid; o.Theta_cc = Theta_cc; o.Pcc = Pcc; o.Gammah = Gammah; o.uh = uh; o.lambda_id = lambda_id; o.zidv = zidv; o.uv = uv; o.yv = yv; o.time = time;
       end
       function u = step(o, y)
            % Read properties
            N = o.N; L = o.L; tau_n = o.tau_n; tau_d = o.tau_d; eta = o.eta; lambda_cc = o.lambda_cc; Eid = o.Eid; Eccy = o.Eccy; Eccu = o.Eccu; Ct = o.Ct; Lu = o.Lu; qr_off = o.qr_off; anc_on = o.anc_on; ah = o.ah; bh = o.bh; xh = o.xh; Theta_id = o.Theta_id; Pid = o.Pid; Theta_cc = o.Theta_cc; Pcc = o.Pcc; Gammah = o.Gammah; uh = o.uh; lambda_id = o.lambda_id; zidv = o.zidv; uv = o.uv; yv = o.yv; time = o.time;

            time = time + 1;

            yv = [y; yv(1:end-1)];
            % calculate the parameters
            Yid = Eid*y;
            phi = [-yv(2:N+1); uv(1:N+1)]';
            Phi_id = Eid*phi;
            Lid = Pid/lambda_id;
            Pid = Lid - Lid*Phi_id'*(1+Phi_id*Lid*Phi_id')^-1*Phi_id*Lid;
            zid = Yid - Phi_id*Theta_id;
            Theta_id = Theta_id + Pid*Phi_id'*zid;
            ah(2:end) = Theta_id(1:N);
            bh = Theta_id(N+1:end);

            % forgetting factor
            zidv = [zid; zidv(1:end-1)];
            zidvx = zidv(1:tau_n);
            xi = sqrt(mean(zidvx.^2)/mean(zidv.^2));
            f = max(0,xi-1);
            lambda_id = 1/(1+eta*f);

            % calculate the state
            xh(1) = y - bh(1)*uv(1);
            for i=2:N
                xh(i) = - ah(i+1:end)'*yv(2:N-i+2) + bh(i+1:end)'*uv(2:N-i+2);
            end

            if time > anc_on
                % control
                Ah = zeros(N);
                Ah(N+1:N+1:end) = 1;
                Ah(:,1) = -ah(2:end);
                Bh = bh(2:end) - ah(2:end)*bh(1);
                Ch = zeros(1,N); Ch(1) = 1;
                Dh = bh(1);

                Cx = Ch;
                for i=1:L
                    Gammah(i,:) = Cx;
                    Cx = Cx*Ah;
                end

                Hh = Gammah*Bh;
                Th = toeplitz([Dh;Hh(1:end-1)], [Dh;zeros(L-1,1)]);

                xh = Ah*xh + Bh*uv(1);
                Ycc = - Eccy*Gammah*xh;
                Phi_cc = Eccy*Th + Eccu;

                Lcc = Pcc/lambda_cc;
                Pcc = Lcc - Lcc*Phi_cc'*((eye(L) + Phi_cc*Lcc*Phi_cc')\Phi_cc)*Lcc;
                Theta_cc = Theta_cc + Pcc*Phi_cc'*(Ycc - Phi_cc*Theta_cc);

                u = Theta_cc(1);
                u = min(Lu, max(-Lu,u));
            else
                u = sqrt(qr_off)*randn;
            end

            uv = [u; uv(1:end-1)];

            % Write properties
            o.N = N; o.L = L; o.tau_n = tau_n; o.tau_d = tau_d; o.eta = eta; o.lambda_cc = lambda_cc; o.Eid = Eid; o.Eccy = Eccy; o.Eccu = Eccu; o.Ct = Ct; o.Lu = Lu; o.qr_off = qr_off; o.anc_on = anc_on; o.ah = ah; o.bh = bh; o.xh = xh; o.Theta_id = Theta_id; o.Pid = Pid; o.Theta_cc = Theta_cc; o.Pcc = Pcc; o.Gammah = Gammah; o.uh = uh; o.lambda_id = lambda_id; o.zidv = zidv; o.uv = uv; o.yv = yv; o.time = time;
       end
    end
end
