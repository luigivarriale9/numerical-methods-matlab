function q = kepler_position(e, t)
%KEPLER_POSITION  Exact position on a Kepler orbit (semi-major axis 1, mu = 1).
%
%   q = kepler_position(e, t) returns the position [x; y] at time t after the
%   pericentre, for eccentricity e. Kepler's equation M = E - e*sin(E) is
%   solved with Newton's method (mean motion = 1, so M = t).

M = t;
E = M;
for k = 1:50
    dE = (E - e * sin(E) - M) / (1 - e * cos(E));
    E = E - dE;
    if abs(dE) < 1e-15, break, end
end
q = [cos(E) - e; sqrt(1 - e^2) * sin(E)];
end
