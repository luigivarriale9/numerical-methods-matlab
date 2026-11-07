function dC = jacobi_drift(Y, mu)
% Change of the Jacobi constant with respect to its initial value.
% Y: 6 x N matrix (x, y, z, vx, vy, vz per column)
x = Y(1,:); y = Y(2,:); z = Y(3,:);
v2 = sum(Y(4:6,:).^2, 1);
d1 = sqrt((x+mu).^2   + y.^2 + z.^2);   % distance from the primary
d2 = sqrt((x-1+mu).^2 + y.^2 + z.^2);   % distance from the secondary
C  = x.^2 + y.^2 + 2*(1-mu)./d1 + 2*mu./d2 - v2;
dC = C - C(1);
end