% CFANVC
classdef Carini < handle
    %#ok<*PROPLC>
    properties
        Nw, D, Lu, R, lambda_h, lambda, mu_scale, Nsh, mu_s_min, qr_off, anc_on
        G, w, sh1, mh, xv, xfv, yv, vv, pe, pf, pg, ps, ns, nw, mu, time
    end
    methods
        function o = Carini(param)
            % Read Parameters
            Nw = param{1}; D = param{2}; Lu = param{3}; R = param{4}; lambda_h = param{5}; lambda = param{6}; mu_scale = param{7}; Nsh = param{8}; mu_s_min = param{9}; qr_off = param{10}; anc_on = param{11};

            G = 1;
            w = zeros(Nw, 1);
            sh1 = zeros(Nsh+D, 1); sh1(1:D+1) = 1;
            mh = zeros(Nw,1);
            xv = zeros(Nw, 1);
            xfv = zeros(Nw, 1);
            yv = zeros(Nw, 1);
            vv = zeros(Nsh+D, 1);
            pe=1;
            pf=1;
            pg=1;
            ps=1;
            ns=0;
            nw=0;
            mu=nan;
            time = 0;

            % Write properties
            o.Nw = Nw; o.D = D; o.Lu = Lu; o.R = R; o.lambda_h = lambda_h; o.lambda = lambda; o.mu_scale = mu_scale; o.Nsh = Nsh; o.mu_s_min = mu_s_min; o.qr_off = qr_off; o.anc_on = anc_on; o.G = G; o.w = w; o.sh1 = sh1; o.mh = mh; o.xv = xv; o.xfv = xfv; o.yv = yv; o.vv = vv; o.pe = pe; o.pf = pf; o.pg = pg; o.ps = ps; o.ns = ns; o.nw = nw; o.mu = mu; o.time = time;
       end
       function y = step(o, e)
            % Read properties
            Nw = o.Nw; D = o.D; Lu = o.Lu; R = o.R; lambda_h = o.lambda_h; lambda = o.lambda; mu_scale = o.mu_scale; Nsh = o.Nsh; mu_s_min = o.mu_s_min; qr_off = o.qr_off; anc_on = o.anc_on; G = o.G; w = o.w; sh1 = o.sh1; mh = o.mh; xv = o.xv; xfv = o.xfv; yv = o.yv; vv = o.vv; pe = o.pe; pf = o.pf; pg = o.pg; ps = o.ps; ns = o.ns; nw = o.nw; mu = o.mu; time = o.time;

            time = time + 1;

            pe = lambda*pe + (1-lambda)*e^2;
            sh = sh1(D+1:end);
            ps = lambda*ps + (1-lambda)*(sh'*sh);

            yv = [0; yv(1:end-1)];
            x = e - sh'*yv;
            % x = d(n);
            xv = [x; xv(1:end-1)];

            v = G*randn;
            vv = [v; vv(1:end-1)];
            f = e - vv'*sh1;
            pf = lambda*pf + (1-lambda)*f^2;

            dh = e - yv'*sh;
            xf = sh'*xv(1:Nsh);
            xfv = [xf; xfv(1:end-1)];
            g = dh - w'*xfv;
            pg = lambda*pg + (1-lambda)*g^2;

            mh = lambda_h*mh + (1-lambda_h)*g*xfv/(xfv'*xfv);
            nw = lambda*nw + (1-lambda)*g*mh'*xfv;

            if time >= anc_on
                mu = nw/pg;
                w = w + mu_scale*mu/(xfv'*xfv)*xfv*g; % variante

                G = sqrt(pe/((R+1)*ps));
                y1 = w'*xv(1:Nw);
                y = min(Lu, max(-Lu, -y1 + vv(D+1)));
            else
                G = sqrt(qr_off);
                y = vv(D+1);
            end
            yv(1) = y;

            sh0 = sh1(1:D);
            ns = lambda*ns + (1-lambda)*((sh0'*sh0)*(vv'*vv))/D;
            mu_s = max(mu_s_min, ns/pf);
            %mu_s = mu_s_min; % variant, no delayed coefficient technique

            sh1 = sh1 + mu_s/(vv'*vv+1e-9)*vv*f;

            % Write properties
            o.Nw = Nw; o.D = D; o.Lu = Lu; o.R = R; o.lambda_h = lambda_h; o.lambda = lambda; o.mu_scale = mu_scale; o.Nsh = Nsh; o.mu_s_min = mu_s_min; o.qr_off = qr_off; o.anc_on = anc_on; o.G = G; o.w = w; o.sh1 = sh1; o.mh = mh; o.xv = xv; o.xfv = xfv; o.yv = yv; o.vv = vv; o.pe = pe; o.pf = pf; o.pg = pg; o.ps = ps; o.ns = ns; o.nw = nw; o.mu = mu; o.time = time;
       end
    end
end
