function [p, A, d] = calc_random_path(N, dM, R)
P = 3; % pulses
r = inf;
while r > R
    d = 2 + rand(1,P)*dM;
    A = randn(P,1);
    n = (0:N-1)';
    p = sinc(n-d)*A;
    p(1) = 0;
    h = freqz(p, 1, 1024);
    r = max(max(abs(h)), 1/min(abs(h)));
end