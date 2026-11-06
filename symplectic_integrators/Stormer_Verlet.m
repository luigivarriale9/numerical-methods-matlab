function [t, y, h] = Stormer_Verlet(t0, tstop, y0, n, flag, fp, fq)
%STORMER_VERLET  Stormer-Verlet method for a Hamiltonian system
%                q' = fq(p, q),  p' = fp(p, q).
%
%   [t, y, h] = Stormer_Verlet(t0, tstop, y0, n, flag, fp, fq)
%
%   Input:   t0, tstop  initial and final time
%            y0         initial state [q0; p0]
%            n          number of steps
%            flag       0: explicit version (fq must not depend on q)
%                       1: version with an implicit update of q (fsolve)
%            fp, fq     function handles fp(p, q) and fq(p, q), column vectors
%   Output:  t          time grid (n+1 points)
%            y          solution [q; p], one column per time step
%            h          step size
%   Order of convergence: 2.
%
%   Note: the flag = 1 branch uses the opposite sign in the half steps (as in the
%   original thesis code), so it integrates the time-reversed flow: for the Kepler
%   problem started at the pericentre the trajectory is mirrored (y -> -y).
%   Conservation properties (energy, angular momentum, Jacobi constant) are not
%   affected by this.

m = length(y0);
h = (tstop - t0) / n;
q = zeros(m/2, n + 1);
p = zeros(m/2, n + 1);
q(:, 1) = y0(1:m/2);
p(:, 1) = y0(m/2 + 1:m);
t = zeros(n + 1, 1);
t(1) = t0;
options = optimset('Display', 'off', 'TolFun', 1e-12, 'TolX', 1e-12);

if flag == 0
    for i = 1:n
        t(i + 1) = t(i) + h;
        g = p(:, i) + (h/2) * fp(p(:, i), q(:, i));
        q(:, i + 1) = q(:, i) + h * fq(g, q(:, i));
        p(:, i + 1) = g + (h/2) * fp(g, q(:, i + 1));
    end
elseif flag == 1
    for i = 1:n
        t(i + 1) = t(i) + h;
        g = p(:, i) - (h/2) * fp(p(:, i), q(:, i));
        q(:, i + 1) = fsolve(@(x) -x + q(:, i) - (h/2) * (fq(g, q(:, i)) + fq(g, x)), ...
                             q(:, i), options);
        p(:, i + 1) = g - (h/2) * fp(g, q(:, i + 1));
    end
end
y = [q; p];
end
