function [t,y,err,global_err,h] = AdamsMoulton3(f,t0,t_end,y0,n,kmax,y_exact)
% Purpose: implement the Adams-Moulton method with k=3 (order 4).
% The implicit equation of each step is solved by fixed-point iteration,
% which converges if h*L*|gamma_0| < 1 (L = Lipschitz constant of f).
%
% Input:
%   f      : right-hand side of the ODE, f(t,y),
%   t0     : left end of the integration interval,
%   t_end  : right end of the integration interval,
%   y0     : initial condition,
%   n      : number of steps,
%   kmax   : maximum number of fixed-point iterations per step,
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
%   kmax=10000;
%   [t,y,err,global_err,h] = AdamsMoulton3(f,0,1,1,5,kmax,y_exact);
% Result: err = 1.0997e-05, h = 2.0000e-01

h=(t_end-t0)/n;          % step size
l=length(y0);
t=t0+(0:n)'*h;           % time grid
y=zeros(l,n+1);          % solution
y(:,1)=y0(:);            % initial value
tol=1e-12;               % tolerance of the fixed-point iteration
% starting values computed with ode45 (tight tolerances: they must not limit the order)
opts=odeset('RelTol',1e-12,'AbsTol',1e-14);
tspan=[t(1),t(2),t(3)];
[~,y_es]=ode45(f,tspan,y0(:),opts);
y(:,2:3)=y_es(end-1:end,:).';

not_converged=0;
for i=1:n-2   % Adams-Moulton
    g=y(:,i+2)+(h/24)*(19*f(t(i+2),y(:,i+2))-5*f(t(i+1),y(:,i+1))+f(t(i),y(:,i)));   % known part of the step
    y_old=y(:,i+2);  % initial guess for the fixed point
    res=tol+1;
    it=0;
    while (res>tol) && (it<kmax)   % fixed point: y(i+3) is built starting from y(i+2)
        p=g+((9/24)*h)*f(t(i+3),y_old);
        res=norm(p-y_old,2);
        y_old=p;
        it=it+1;
    end
    if res>tol
        not_converged=not_converged+1;
    end
    y(:,i+3)=p;
end
if not_converged>0
    warning('AdamsMoulton3: fixed-point iteration did not converge in %d steps (reduce h or increase kmax).',not_converged);
end

% error at the final time
if l==1
    err=abs(y_exact(t(end))-y(end));
else
    err=norm(y_exact-y(:,end),2);
end
global_err=err;
end
