function [t,y,err,global_err,h] = AdamsB4(f,t0,t_end,y0,n,y_exact)
% Purpose: implement the Adams-Bashforth method with k=4 (order 4).
%
% Input:
%   f      : right-hand side of the ODE, f(t,y),
%   t0     : left end of the integration interval,
%   t_end  : right end of the integration interval,
%   y0     : initial condition,
%   n      : number of steps,
%   y_exact: exact solution. Scalar problem: function handle y_exact(t).
%            System of ODEs: column vector with the exact solution at t_end.
% Output:
%   t         : time grid (column vector, n+1 points),
%   y         : approximate solution, one column per time (size l-by-(n+1)),
%   err       : error at the final time, norm(y_exact - y(:,end), 2),
%   global_err: error indicator, equal to err (kept for compatibility),
%   h         : step size.
%
% Example on a test problem:
%   f=@(t,y) -y;
%   y_exact=@(t) exp(-t);
%   [t,y,err,global_err,h] = AdamsB4(f,0,1,1,5,y_exact);
% Result: err = 1.0345e-04, h = 2.0000e-01

h=(t_end-t0)/n;          % step size
l=length(y0);
t=t0+(0:n)'*h;           % time grid
y=zeros(l,n+1);          % solution
y(:,1)=y0(:);            % initial value
% starting values computed with ode45 (tight tolerances: they must not limit the order)
opts=odeset('RelTol',1e-12,'AbsTol',1e-14);
tspan=[t(1),t(2),t(3),t(4)];
[~,y_es]=ode45(f,tspan,y0(:),opts);
y(:,2:4)=y_es(end-2:end,:).';

for i=1:n-3   % Adams-Bashforth
    y(:,i+4)=y(:,i+3)+((h/24)*(55*f(t(i+3),y(:,i+3))-59*f(t(i+2),y(:,i+2))+37*f(t(i+1),y(:,i+1))-9*f(t(i),y(:,i))));
end

% error at the final time
if l==1
    err=abs(y_exact(t(end))-y(end));
else
    err=norm(y_exact-y(:,end),2);
end
global_err=err;
end
