function [t, y, h] = Midpoint_Implicit(f, t0, tstop, y0, n)
%MIDPOINT_IMPLICIT  Implicit midpoint method for y' = f(t, y).
%
%   [t, y, h] = Midpoint_Implicit(f, t0, tstop, y0, n)
%
%   At each step the implicit equation for the midpoint value yh is solved
%   with fsolve, then y(i+1) = 2*yh - y(i).
%
%   Input:   f      function handle f(t, y), returns a column vector
%            t0, tstop  initial and final time
%            y0     initial state
%            n      number of steps
%   Output:  t      time grid (n+1 points)
%            y      solution, one column per time step (m x (n+1))
%            h      step size
%   Order of convergence: 2.

m = length(y0);
h = (tstop - t0) / n;
t = zeros(n + 1, 1);
y = zeros(m, n + 1);
t(1) = t0;
y(:, 1) = y0(:);
options = optimset('Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);

for i = 1:n
    to = t(i);
    yo = y(:, i);
    th = to + 0.5 * h;
    yh = yo + 0.5 * h * f(to, yo);               % initial guess
    yh = fsolve(@(yh) midpoint_residual(f, to, yo, th, yh), yh, options);
    t(i + 1) = to + h;
    y(:, i + 1) = 2.0 * yh - yo;
end
end

function r = midpoint_residual(f, to, yo, th, yh)
% Residual of the implicit equation yh = yo + (th - to) * f(th, yh)
r = yh - yo - (th - to) * f(th, yh);
end
