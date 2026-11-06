% KEPLER_MULTISTEP_CONVERGENCE  Observed order of convergence of the Adams
%                               methods on the Kepler problem (e = 0.7).
%
% Same problem, final time, step counts and reference as kepler_convergence.m
% in symplectic_integrators: the global error is the position error at t = 3
% against the exact solution of Kepler's equation (kepler_position.m).
% The starting values of the multistep methods come from ode45 with tight
% tolerances (see the Adams*.m files).
%
% Methods: AB1-AB4 (orders 1-4), AM1-AM4 (orders 2-5), PECE (order 5).

clear; close all
here = fileparts(mfilename('fullpath'));
addpath(here, fullfile(here, '..', 'symplectic_integrators'));

e = 0.7;
t0 = 0;
tstop = 3;
y0 = [1 - e; 0; 0; sqrt((1 + e) / (1 - e))];     % [q; p]
q_exact = kepler_position(e, tstop);             % exact position at t = tstop

r3 = @(q) sqrt(q(1)^2 + q(2)^2)^3;
f  = @(t, y) [y(3); y(4); -y(1) / r3(y(1:2)); -y(2) / r3(y(1:2))];

ns = [250 500 1000 2000 4000];                   % numbers of steps
names = {'AB1','AB2','AB3','AB4','AM1','AM2','AM3','AM4','PECE'};
expected = [1 2 3 4 2 3 4 5 5];
kmax = 10000;                                    % fixed-point iterations (Adams-Moulton)
dummy = nan(4, 1);   % the Adams*.m functions need an exact vector for their own
                     % error output; here the error is computed below instead

h = (tstop - t0) ./ ns;
err = zeros(numel(names), numel(ns));
for k = 1:numel(ns)
    n = ns(k);
    [~, y] = AdamsB1(f, t0, tstop, y0, n, dummy);         err(1, k) = norm(y(1:2, end) - q_exact);
    [~, y] = AdamsB2(f, t0, tstop, y0, n, dummy);         err(2, k) = norm(y(1:2, end) - q_exact);
    [~, y] = AdamsB3(f, t0, tstop, y0, n, dummy);         err(3, k) = norm(y(1:2, end) - q_exact);
    [~, y] = AdamsB4(f, t0, tstop, y0, n, dummy);         err(4, k) = norm(y(1:2, end) - q_exact);
    [~, y] = AdamsMoulton1(f, t0, tstop, y0, n, kmax, dummy); err(5, k) = norm(y(1:2, end) - q_exact);
    [~, y] = AdamsMoulton2(f, t0, tstop, y0, n, kmax, dummy); err(6, k) = norm(y(1:2, end) - q_exact);
    [~, y] = AdamsMoulton3(f, t0, tstop, y0, n, kmax, dummy); err(7, k) = norm(y(1:2, end) - q_exact);
    [~, y] = AdamsMoulton4(f, t0, tstop, y0, n, kmax, dummy); err(8, k) = norm(y(1:2, end) - q_exact);
    [~, y] = AdamsPECE(f, t0, tstop, y0, n, dummy);       err(9, k) = norm(y(1:2, end) - q_exact);
end

fprintf('Observed order of convergence (e = %.1f), error at t = %g\n', e, tstop);
fprintf('%-6s', 'n'); fprintf('%12d', ns); fprintf('\n');
for j = 1:numel(names)
    c = polyfit(log(h), log(err(j, :)), 1);
    fprintf('%-6s', names{j}); fprintf('%12.3e', err(j, :));
    fprintf('   observed %.2f (theory %d)\n', c(1), expected(j));
    assert(abs(c(1) - expected(j)) < 0.4, ...
           'Observed order of %s is not close to %d', names{j}, expected(j));
end

figure
cols = jet(numel(names));
for j = 1:numel(names)
    loglog(h, err(j, :), 'o-', 'Color', cols(j, :), 'LineWidth', 1.2); hold on
end
grid on
xlabel('step size h'); ylabel('global error at t = 3')
legend(names, 'Location', 'southeast')
title('Adams methods on the Kepler problem')
