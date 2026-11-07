close all
clear
clc

% ---------- Problem data ----------
% Test case (mu, y0) taken from Rubinsztejn, Matlab-Astrodynamics-Library
mu = 1.215e-2;
y0 = [1-mu, .0455, 0, -.5322, .2, 0];
t0 = 0;  tstop = 20;  tspan = [t0 tstop];
opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
flag = 1;
n  = 15394;      % steps: implicit midpoint, Stormer-Verlet
n1 = 120000;     % steps: symplectic Euler
n2 = 700000;     % steps: explicit Euler

% ---------- System (variables q, v) ----------
f  = @(t,y) [y(4:6); accel(y(1:3), y(4:6), mu)];
fq = @(p,q) p(:);
fp = @(p,q) accel(q(:), p(:), mu);

% ---------- Jacobi constant of the initial condition ----------
d1 = norm([y0(1)+mu,   y0(2), y0(3)]);
d2 = norm([y0(1)-1+mu, y0(2), y0(3)]);
jc = y0(1)^2 + y0(2)^2 + 2*(1-mu)/d1 + 2*mu/d2 - sum(y0(4:6).^2);

% ---------- Choice of the energy level (Jacobi constant) ----------
jc0  = jc;
Cmax = jc0 + sum(y0(4:6).^2);          % maximum value at the initial point (v = 0)
ans_C = input(sprintf('Jacobi constant C [Enter = %.4f, max %.4f]: ', jc0, Cmax));
if ~isempty(ans_C)
    if ans_C > Cmax
        error('C too high: no real velocity exists at the initial point (C <= %.4f).', Cmax);
    end
    v2 = y0(1)^2 + y0(2)^2 + 2*(1-mu)/d1 + 2*mu/d2 - ans_C;
    y0(4:6) = y0(4:6) * sqrt(v2) / norm(y0(4:6));   % same direction, rescaled magnitude
    jc = ans_C;
end
fprintf('Jacobi constant used: C = %.6f\n', jc);

% ---------- Choice of the methods ----------
disp('1 ode45 default | 2 Midpoint | 3 Explicit Euler | 4 ode45 tight | 5 Stormer-Verlet | 6 Symplectic Euler')
chosen = input('Methods to run, e.g. [2 5 6] (Enter = all): ');
if isempty(chosen), chosen = 1:6; end
names = {'ode45 default','Implicit midpoint','Explicit Euler','ode45 tight tolerances','Stormer-Verlet','Symplectic Euler'};
SOL = cell(1,6);  TT = cell(1,6);

if ismember(1,chosen)
    [TT{1},SOL{1}] = ode45(f, tspan, y0);          SOL{1} = SOL{1}';
end
if ismember(2,chosen)
    [TT{2},SOL{2},~] = Midpoint_Implicit(f, t0, tstop, y0, n);
end
if ismember(3,chosen)
    [TT{3},SOL{3},~] = Euler_Explicit(f, t0, tstop, y0, n2);
end
if ismember(4,chosen)
    [TT{4},SOL{4}] = ode45(f, tspan, y0, opts);    SOL{4} = SOL{4}';
end
if ismember(5,chosen)
    [TT{5},SOL{5},~] = Stormer_Verlet(t0, tstop, y0, n, flag, fp, fq);
end
if ismember(6,chosen)
    [TT{6},SOL{6},~] = Symplectic_euler(t0, tstop, y0, n1, flag, fp, fq);
end

% ---------- Forbidden region (C < jc) ----------
x  = -1.5:0.005:1.5;
y1 = -1.5:0.005:1.5;
[X,Y] = meshgrid(x, y1);
D1 = sqrt((X+mu).^2   + Y.^2);
D2 = sqrt((X-1+mu).^2 + Y.^2);
Cgrid = X.^2 + Y.^2 + 2*(1-mu)./D1 + 2*mu./D2;

% ---------- Lagrange points ----------
dOm = @(x) x - (1-mu)*(x+mu)./abs(x+mu).^3 - mu*(x-1+mu)./abs(x-1+mu).^3;
xL1 = fzero(dOm, [-mu+0.05, 1-mu-0.05]);
xL2 = fzero(dOm, [1-mu+0.05, 2]);
xL3 = fzero(dOm, [-2, -mu-0.05]);
xL4 = 0.5 - mu;
xLPos = [xL1, xL2, xL3, xL4, xL4]';
yLPos = [0, 0, 0, sqrt(3)/2, -sqrt(3)/2]';

% ---------- Trajectory plots ----------
for k = chosen
    plot_orbit(k, SOL{k}, names{k}, x, y1, Cgrid, jc, xLPos, yLPos, mu)
end

% ---------- Change of the Jacobi constant ----------
figure(7); clf; hold on
for k = chosen
    plot(TT{k}, jacobi_drift(SOL{k}, mu), 'DisplayName', names{k})
end
title('Change of the Jacobi constant')
xlabel('Time (non-dimensional)'); ylabel('\Delta C')
grid on; legend('Location','best')

% ---------- Save figures ----------
outDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figures');
if ~isfolder(outDir), mkdir(outDir); end
slugs = {'ode45_default','implicit_midpoint','explicit_euler','ode45_tight','stormer_verlet','symplectic_euler'};
for k = chosen
    exportgraphics(figure(k), fullfile(outDir, ['cr3bp_' slugs{k} '.png']), 'Resolution', 200);
end
exportgraphics(figure(7), fullfile(outDir, 'cr3bp_jacobi_drift.png'), 'Resolution', 200);
fprintf('Figures saved in %s\n', outDir);

% =====================================================================
function a = accel(q, v, mu)
% Acceleration in the rotating frame: q = position (3), v = velocity (3)
d1 = norm([q(1)+mu,   q(2), q(3)]);
d2 = norm([q(1)-1+mu, q(2), q(3)]);
a = [ 2*v(2) + q(1) - (1-mu)*(q(1)+mu)/d1^3 - mu*(q(1)-1+mu)/d2^3;
     -2*v(1) + q(2) - (1-mu)*q(2)/d1^3      - mu*q(2)/d2^3;
                    - (1-mu)*q(3)/d1^3      - mu*q(3)/d2^3 ];
end

function plot_orbit(n, Y, ttl, x, y1, Cgrid, jc, xL, yL, mu)
figure(n); clf; hold on
contourf(x, y1, Cgrid, [0 jc], 'LineColor','k', 'DisplayName','Forbidden region')
colormap([0.85 0.85 0.85])
plot(Y(1,:), Y(2,:), 'b', 'DisplayName','Trajectory')
plot(xL, yL, 'kx', 'MarkerSize',8, 'DisplayName','Lagrange points')
plot(-mu,  0, 'o', 'MarkerFaceColor',[0.8 0.2 0.2], 'DisplayName','Primary')
plot(1-mu, 0, 'o', 'MarkerFaceColor',[0.2 0.4 0.8], 'DisplayName','Secondary')
axis equal; axis([-1.5 1.5 -1.5 1.5]); grid on
title(ttl); xlabel('q_1'); ylabel('q_2'); legend('Location','best')
end