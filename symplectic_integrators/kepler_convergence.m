% KEPLER_CONVERGENCE  Observed order of convergence of four integrators
%                     on the Kepler problem (eccentricity e = 0.7).
%
% The initial point is the pericentre (1 - e, 0). The exact position at the final
% time t = 3 is obtained from Kepler's equation M = E - e*sin(E), solved with
% Newton's method. The observed order is the slope of the global error against
% the step size in a log-log plot.
%
% The final time is not a whole period: at t = 2*pi the orbit returns to its
% starting point and the symmetry of the orbit hides the first-order error of
% symplectic Euler (the observed order would be 2 instead of 1).
%
% Methods: explicit Euler (order 1), symplectic Euler (order 1),
%          implicit midpoint (order 2), Stormer-Verlet (order 2).

clear; close all
addpath(fileparts(mfilename('fullpath')));

e = 0.7;
t0 = 0;
tstop = 3;
y0 = [1 - e; 0; 0; sqrt((1 + e) / (1 - e))];     % [q; p]
q_exact = kepler_position(e, tstop);             % exact position at t = tstop

r3 = @(q) sqrt(q(1)^2 + q(2)^2)^3;
f  = @(t, y) [y(3); y(4); -y(1) / r3(y(1:2)); -y(2) / r3(y(1:2))];
fq = @(p, q) [p(1); p(2)];
fp = @(p, q) [-q(1) / r3(q); -q(2) / r3(q)];

flag = 0;                                        % fp does not depend on p
ns = [250 500 1000 2000 4000];                   % numbers of steps
names = {'Explicit Euler', 'Symplectic Euler', 'Implicit midpoint', 'Stormer-Verlet'};
expected = [1 1 2 2];

h = (tstop - t0) ./ ns;
err = zeros(4, numel(ns));
for k = 1:numel(ns)
    n = ns(k);
    [~, y1] = Euler_Explicit(f, t0, tstop, y0, n);
    [~, y2] = Symplectic_euler(t0, tstop, y0, n, flag, fp, fq);
    [~, y3] = Midpoint_Implicit(f, t0, tstop, y0, n);
    [~, y4] = Stormer_Verlet(t0, tstop, y0, n, flag, fp, fq);
    err(1, k) = norm(y1(1:2, end) - q_exact);
    err(2, k) = norm(y2(1:2, end) - q_exact);
    err(3, k) = norm(y3(1:2, end) - q_exact);
    err(4, k) = norm(y4(1:2, end) - q_exact);
end

fprintf('Observed order of convergence (e = %.1f)\n', e);
for j = 1:4
    c = polyfit(log(h), log(err(j, :)), 1);
    fprintf('  %-18s %.2f (expected %d)\n', names{j}, c(1), expected(j));
    assert(abs(c(1) - expected(j)) < 0.25, ...
           'Observed order of %s is not close to %d', names{j}, expected(j));
end

figure
loglog(h, err', 'o-', 'LineWidth', 1.5)
hold on
loglog(h, err(2, end) * (h / h(end)), 'k--', h, err(4, end) * (h / h(end)).^2, 'k:')
grid on
xlabel('step size h'); ylabel('global error at t = 3')
legend([names, {'h', 'h^2'}], 'Location', 'southeast')
title('Order of convergence on the Kepler problem')
