function [t, y, h] = Euler_Explicit(f, t0, tstop, y0, n)
%EULER_EXPLICIT  Explicit Euler method for y' = f(t, y).
%
%   [t, y, h] = Euler_Explicit(f, t0, tstop, y0, n)
%
%   Input:   f      function handle f(t, y), returns a column vector
%            t0     initial time
%            tstop  final time
%            y0     initial state
%            n      number of steps
%   Output:  t      time grid (column vector, n+1 points)
%            y      solution, one column per time step (m x (n+1))
%            h      step size
%   Order of convergence: 1.

m = length(y0);
h = (tstop - t0) / n;
t = zeros(n + 1, 1);
y = zeros(m, n + 1);
t(1) = t0;
y(:, 1) = y0(:);

for i = 1:n
    t(i + 1) = t(i) + h;
    y(:, i + 1) = y(:, i) + h * f(t(i), y(:, i));
end
end
