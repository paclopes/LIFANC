% simulate_all
%
% Source Code for the work:
% Limiting Instabilities in Feedback Active Noise Control
% Paulo A. C. Lopes
% to be published.
% 
% Paulo Alexandre Crisóstomo Lopes, 2 July 2025
% Instituto de Engenharia de Sistemas e Computadores - Investigação e Desenvolvimento
% Instituto Superior Técnico
% Universidade de Lisboa

seed = 218;
fprintf(1, 'seed: %d\n', seed);
rng(seed);
seed = seed + 1;

N_algorithms = 4;
N_simulations = 3;

plots = {'e', 'd', 'u'};
logs = unique(flatten(plots));
log_vars = join_strings(logs);
        
fs = 2000;
Ln = 20*fs;
log_values = zeros(N_algorithms, N_simulations, Ln, length(logs));
plot_handle = 1;
plot_names = {'LI-FMiNFxLMS', 'Carini', 'Mohseni', 'CFANVC'};

for alg_index = 1:N_algorithms
    for simulation = 1:N_simulations
        
        % simulations parameters and signals
        Ms = 3;         % number of sins
        As = rand(1,Ms);
        ff = (100*rand + 50)*[1,2,3]';
        phase = 2*pi*rand(Ms,1);
        qv = 0.1;       % background noise power
        R = 3;          % 1/R < abs(S(f)) < R
        M = 10;         % output stauration
        Ns = 128;       % actual secondary path length
        dM = 8;
        dC = 3;
        anc_on = Ln/10;
        qr_off = 0.1;
        
        n = 1:Ln;
        d_signal = As*sin(2*pi*ff/fs*n + phase);
        d_signal = d_signal/std(d_signal) + sqrt(qv)*randn(1,Ln);
        
        [sv1, Sa, Sd] = calc_random_path(Ns, dM, R);
        sv = sv1;
        uvx = zeros(Ns, 1);
        
        switch alg_index
            case 1
                %system('python.exe .\class_helper.py LI_FMiNFxLMS.m');
            
                % algorithm parameters
                N = 32;             % filter length
                Nx = 2*N;           % wait for transients
                mu  = 0.01;         % normalized step size
                mus = 0.1;          % secondary path normalized step size
                tau = 0.3;          % proportional to max abs(ds w)
                PS = 0.1;           % min_j p_j/mean p_i, p_i power of sin i
                % mu/(sum_i pi) pi > mu/N PS >> mu/N PS/10
                lambda_w = 1 - mu/N*PS/10;
                alpha = 1/4;        % to calculate the residual noise power
                beta = 1.2;         % to calculate the instability threshold
                gamma = 1 - 0.1/fs;     % to limit the change rate
                eta = 1 + 1/fs;     % increase in gu
                lambda = 1 - 1/100; % to calculate average power
                q = 1e-3;
        
                parameters = {N, Nx, mu, mus, tau, lambda_w, alpha, beta, gamma, eta, lambda, q, M, R, 'e', qr_off, anc_on};
                algorithm = LI_FMiNFxLMS(parameters);
            case 2
                %system('python.exe .\class_helper.py Carini.m');
        
                % algorithm parameters
                Nw = 32;            % control filter length
                D = 16;             % delay of delayed coefficient
                Lu = 10;            % saturation of the anti-noise signal
                R = 3;
                lambda_h = 0.9;
                lambda = 0.99;      % forgeting factor to calculate the power
                
                mu_scale = 1;       % factor that multiplies by the calculated step size.
                Nsh = 32;           % Nh: Secondary path estimate length
                mu_s_min = 0.01;    % due to the delay coefficients thechnique this will
                % become the actual step size after some time
        
                parameters = {Nw, D, Lu, R, lambda_h, lambda, mu_scale, Nsh, mu_s_min, qr_off, anc_on};
                algorithm = Carini(parameters);
            case 3
                %system('python.exe .\class_helper.py Mohseni.m');
        
                % algorithm parameters
                N = 16;             % model order (size-1)
                L = 64;             % horizont size
                tau_n = 100;
                tau_d = 200;
                eta = 0.1;
                lambda_cc = 0.1; %  This low value results in a simple matrix inversion
                Eid = 1;
                Eccy = 100*eye(L);
                Eccu = 1*eye(L);
                Ct = 1;             % This is hardwired in the code.
                Lu = 10;
        
                parameters = {N, L, tau_n, tau_d, eta, lambda_cc, Eid, Eccy, Eccu, Ct, Lu, qr_off, anc_on};
                algorithm = Mohseni(parameters);
            case 4
                % system('python.exe .\class_helper.py CFANVC.m');
        
                % simulations parameters and signals
                N = 31;                     % model order (size-1)
                L = 64;                     % the equations are 2*L, from i = -L+1 to L
                M = 32*(N+1);               % algorithm memory
                P = M;                      % size of qv estimation blocks
                R = N;                      % number of samples used to obtain uM
                Lu = 10;                    % saturation of the antinoise signal
                wd = 10000;                 % cost weigth of past control signal changes
                wu = 0;                     % cost weigth of control effort
                deltax = 1e-3;              % prior on x variance
                delta = 1e-9;               % actual value added Rxx to calc Sxx (<deltax)
                alpha = 0;                  % controls auxiliary noise power
                eta = 1.01;                 % max grouth rate of u (per Q samples)
                Q = 8;                      % controller update period
        
                parameters = {N, L, M, P, R, Lu, wd, wu, deltax, delta, alpha, eta, Q, qr_off, anc_on};
                algorithm = CFANVC(parameters);
        end
       
        % simulations
        tic;
        for n = 1:Ln
        
            if n == Ln/2
                sv = -sv;
            end
        
            d = d_signal(n);
            uvx = [0; uvx(1:end-1)]; % sv(1) is taken as zero
            e = d + sv'*uvx;
        
            u = algorithm.step(e);
            uvx(1) = u;
        
            %log_values_n(sim_logs) = [d, e, u];
            log_values_n = eval(log_vars);
            log_values(alg_index, simulation, n, :) = log_values_n;
        end
        toc;
    end

    cool_fig(plot_handle);
    N_means = 100;
    t = (1:N_means:Ln)/fs;
    values = reshape(log_values(alg_index, 1, :, 1), N_means , Ln/N_means); % d
    plot(t, mean(values.^2));
    hold on;
    values = reshape(log_values(alg_index, :, :, 2), N_simulations, N_means, Ln/N_means); % e
    plot(t, squeeze(mean(values.^2,2)));
    legend('d', 'sim 1', 'sim 2', 'sim 3');
    grid on;
    xlabel('time (s)')
    ylabel('Noise power')
    set(gcf, 'Name', [num2str(plot_handle), ' - ', plot_names{plot_handle}]);
    plot_handle = plot_handle + 1;
    hold off;
    drawnow;
end

figs = get(0, 'Children');
for j = 1:length(figs)
    set(figs(j),'NumberTitle','off');
    set(get(figs(j), 'CurrentAxes'), 'Ylim', [0,6]);
end

% save2pdf('FMiFxLMS', 1, 600);
% save2pdf('Carini', 2, 600);
% save2pdf('Mohseni', 3, 600);
% save2pdf('CFANC', 4, 600);
