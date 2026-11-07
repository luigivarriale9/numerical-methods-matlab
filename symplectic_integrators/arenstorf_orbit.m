% ARENSTORF_ORBIT  Comparison of ODE45, implicit midpoint, explicit Euler,
% Stormer-Verlet and symplectic Euler on the Arenstorf periodic orbit of the
% planar circular restricted three-body problem (rotating frame).
%
% Problem and initial data as in
%   E. Hairer, S. P. Norsett, G. Wanner, "Solving Ordinary Differential
%   Equations I" (Springer, 1987), after R. F. Arenstorf, Am. J. Math. 85
%   (1963), 27-35; see also https://www.johndcook.com/blog/2020/02/08/arenstorf-orbit/
% The high-precision values of y0 and T are verified numerically: the orbit
% closes after one period T.
%
% The energy error of each method is measured against the ODE45 solution
% (tolerances 1e-8). The exact solution is periodic with period T.
% Requires: Midpoint_Implicit.m, Euler_Explicit.m, Stormer_Verlet.m,
% Symplectic_euler.m.

close all
clear
clc

% ---------- Data ----------
mu1 = 0.012277471;                % mass of the secondary
mu2 = 1 - mu1;                    % mass of the primary
T   = 17.0652165601579625588917206249;          % period
y0  = [0.994, 0.0, 0.0, -2.00158510637908252240537862224];   % [x y vx vy]
t0 = 0;  tstop = T;

% Number of steps. 300000 is the value used in the thesis; the implicit midpoint
% rule solves a nonlinear system at every step, so use e.g. 30000 for a quick test.
n = 300000;
flagSV = 0;                       % Stormer-Verlet (same result as flag = 1, faster)
flagSE = 1;                       % symplectic Euler

% ---------- Equations of motion ----------
d1 = @(x,y) sqrt(((x + mu1).^2 + y.^2).^3);     % |distance to primary|^3
d2 = @(x,y) sqrt(((x - mu2).^2 + y.^2).^3);     % |distance to secondary|^3

f  = @(t,y) [y(3); y(4);
    y(1) + 2*y(4) - mu2*(y(1) + mu1)/d1(y(1),y(2)) - mu1*(y(1) - mu2)/d2(y(1),y(2));
    y(2) - 2*y(3) - mu2*y(2)/d1(y(1),y(2))          - mu1*y(2)/d2(y(1),y(2))];

fq = @(p,q) [p(1); p(2)];
fp = @(p,q) [q(1) + 2*p(2) - mu2*(q(1) + mu1)/d1(q(1),q(2)) - mu1*(q(1) - mu2)/d2(q(1),q(2));
             q(2) - 2*p(1) - mu2*q(2)/d1(q(1),q(2))          - mu1*q(2)/d2(q(1),q(2))];

% Energy (Hamiltonian in the rotating frame), vectorised over columns
H = @(Y) 0.5*(Y(3,:).^2 + Y(4,:).^2) - 0.5*(Y(1,:).^2 + Y(2,:).^2) ...
         - mu1./sqrt((Y(1,:) - mu2).^2 + Y(2,:).^2) ...
         - mu2./sqrt((Y(1,:) + mu1).^2 + Y(2,:).^2);

% ---------- Integration ----------
h  = (tstop - t0)/n;
tt = linspace(t0, tstop, n+1);

disp('ODE45')
opts = odeset('RelTol',1e-8,'AbsTol',1e-8,'Stats','on');
[~, y1] = ode45(f, tt, y0, opts);       y1 = y1';

disp('Implicit midpoint')
[~, y2] = Midpoint_Implicit(f, t0, tstop, y0, n);

disp('Explicit Euler')
[~, y3] = Euler_Explicit(f, t0, tstop, y0, n);

disp('Stormer-Verlet')
[~, y4] = Stormer_Verlet(t0, tstop, y0, n, flagSV, fp, fq);

disp('Symplectic Euler')
[~, y5] = Symplectic_euler(t0, tstop, y0, n, flagSE, fp, fq);

% ---------- Energy errors ----------
Y  = {y1, y2, y3, y4, y5};
names = {'ODE45', 'Implicit midpoint', 'Explicit Euler', 'Stormer-Verlet', 'Symplectic Euler'};
E  = cellfun(H, Y, 'UniformOutput', false);
fprintf('\nMaximum energy error with respect to ODE45 (h = %.3e):\n', h);
for k = 2:5
    fprintf('  %-18s %.3e\n', names{k}, max(abs(E{k} - E{1})));
end
fprintf('Distance from y0 after one period:\n');
for k = 1:5
    fprintf('  %-18s %.3e\n', names{k}, norm(Y{k}(:,end) - y0'));
end

% ---------- Figures ----------
outDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~isfolder(outDir), mkdir(outDir); end
slugs = {'ode45', 'implicit_midpoint', 'explicit_euler', 'stormer_verlet', 'symplectic_euler'};

for k = 1:5
    plot_orbit(k, Y{k}, names{k});
    exportgraphics(figure(k), fullfile(outDir, ['arenstorf_' slugs{k} '.png']), 'Resolution', 200);
end

figure(6); clf; hold on
plot(tt, E{1} - E{1}(1), 'r-',  'DisplayName', names{1})
plot(tt, E{2} - E{1}(1), 'm-',  'DisplayName', names{2})
plot(tt, E{4} - E{1}(1), 'b-',  'DisplayName', names{4})
plot(tt, E{5} - E{1}(1), 'k-',  'DisplayName', names{5})
grid on; xlabel('t'); ylabel('H(t) - H(0)'); legend('Location','best')
title('Energy error')
exportgraphics(figure(6), fullfile(outDir, 'arenstorf_energy.png'), 'Resolution', 200);

figure(7); clf; hold on
plot(tt, E{1} - E{1}(1), 'r-', 'DisplayName', names{1})
plot(tt, E{3} - E{1}(1), 'b-', 'DisplayName', names{3})
grid on; xlabel('t'); ylabel('H(t) - H(0)'); legend('Location','best')
title('Energy error: explicit Euler')
exportgraphics(figure(7), fullfile(outDir, 'arenstorf_energy_explicit_euler.png'), 'Resolution', 200);
fprintf('Figures saved in %s\n', outDir);

% =====================================================================
function plot_orbit(n, Y, ttl)
figure(n); clf; hold on
plot(Y(1,:), Y(2,:), 'LineWidth', 2)
plot(Y(1,1),   Y(2,1),   'g.', 'MarkerSize', 25)    % start
plot(Y(1,end), Y(2,end), 'r.', 'MarkerSize', 25)    % end
grid on; axis equal
xlabel('q_1'); ylabel('q_2'); title(ttl)
end