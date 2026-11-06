function [t, y, h] = Symplectic_euler(t0, tstop, y0, n, flag, fp, fq)
%SYMPLECTIC_EULER  Symplectic Euler method for a Hamiltonian system
%                  q' = fq(p, q),  p' = fp(p, q).
%
%   [t, y, h] = Symplectic_euler(t0, tstop, y0, n, flag, fp, fq)
%
%   Input:   t0, tstop  initial and final time
%            y0         initial state [q0; p0] (first half positions, second half momenta)
%            n          number of steps
%            flag       0: p is updated explicitly
%                       1: p is updated implicitly (fsolve), needed when fp depends on p
%            fp, fq     function handles fp(p, q) and fq(p, q), column vectors
%   Output:  t          time grid (n+1 points)
%            y          solution [q; p], one column per time step
%            h          step size
%   Order of convergence: 1.

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
        p(:, i + 1) = p(:, i) + h * fp(p(:, i), q(:, i));
        q(:, i + 1) = q(:, i) + h * fq(p(:, i + 1), q(:, i));
    end
elseif flag == 1
    for i = 1:n
        t(i + 1) = t(i) + h;
        p(:, i + 1) = fsolve(@(x) -x + p(:, i) + h * fp(x, q(:, i)), p(:, i), options);
        q(:, i + 1) = q(:, i) + h * fq(p(:, i + 1), q(:, i));
    end
end
y = [q; p];
end
