function [S_l, S_m, S_r, xx, rho_xx, u_xx, p_xx] = Euler_exact_Riemann(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol, xx_in)
% 
%   solve the exact Riemann problem of the 1D Euler equation (for ideal polytopic gas only) and compute the primitive variables at nodes xx at time T=1. 
%   Here, S_l, S_m, S_r are the wave velocities. (If there is a rarefaction fan, the output is the min and max of its velocity.)
%   The output xx is refined so that critical points are included. 
% 
% references:
% [1] Eleuterio F. Toro (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction, 3rd eds. Springer-Verlag Berlin Heidelberg. https://doi.org/10.1007/b79761

[S_l, S_r, rho_ml, rho_mr, u_m, p_m] = Euler_exact_Riemann_kernel(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol);

S_m = u_m;

% xx for computation
xx = unique(sort([xx_in(:); S_l(:); S_m(:); S_r(:)])); % normalize and add critical points

rho_xx = nan(numel(xx), 1);
u_xx = nan(numel(xx), 1);
p_xx = nan(numel(xx), 1);

a_l = sqrt(gamma*p_l/rho_l);
a_r = sqrt(gamma*p_r/rho_r);

for k = 1: numel(xx)
    if xx(k) <= min(S_l)
        % left most
        rho_xx(k) = rho_l;
        u_xx(k) = u_l;
        p_xx(k) = p_l;
    elseif xx(k) < max(S_l)
        % left-rarefaction
        u_xx(k) = ((gamma-1.0)/(gamma+1.0))*u_l + (2.0/(gamma+1.0))*(xx(k) + a_l);
        a_xx = u_xx(k) - xx(k);
        rho_xx(k) = (((rho_l^gamma)*a_xx^2)/(gamma*p_l))^(1.0/(gamma-1.0));
        p_xx(k) = rho_xx(k)*a_xx^2/gamma;
    elseif xx(k) <= S_m
        % middle-left
        rho_xx(k) = rho_ml;
        u_xx(k) = u_m;
        p_xx(k) = p_m;
    elseif xx(k) <= min(S_r)
        % middle-right
        rho_xx(k) = rho_mr;
        u_xx(k) = u_m;
        p_xx(k) = p_m;
    elseif xx(k) < max(S_r)
        % right-rarefaction
        u_xx(k) = ((gamma-1.0)/(gamma+1.0))*u_r + (2.0/(gamma+1.0))*(xx(k) - a_r);
        a_xx = xx(k) - u_xx(k);
        rho_xx(k) = (((rho_l^gamma)*a_xx^2)/(gamma*p_l))^(1.0/(gamma-1.0));
        p_xx(k) = rho_xx(k)*a_xx^2/gamma;
    else
        % right-most
        rho_xx(k) = rho_r;
        u_xx(k) = u_r;
        p_xx(k) = p_r;
    end
end

end

function [S_l, S_r, rho_ml, rho_mr, u_m, p_m] = Euler_exact_Riemann_kernel(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol)
% 
%   solve the exact Riemann problem of the 1D Euler equation (for ideal polytopic gas only)
% 
% references:
% [1] Eleuterio F. Toro (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction, 3rd eds. Springer-Verlag Berlin Heidelberg. https://doi.org/10.1007/b79761
% [2] Virginia Notaro (2024). RiemannExact(p1,rho1,u1,p4,rho4,u4,tol) (https://www.mathworks.com/matlabcentral/fileexchange/48734-riemannexact-p1-rho1-u1-p4-rho4-u4-tol), MATLAB Central File Exchange.

% check the inputs

assert(isfloat(rho_l) && isfinite(rho_l) && rho_l > 0.0, 'Euler_exact_Riemann: invalid rho_l input!');
assert(isfloat(u_l) && isfinite(u_l), 'Euler_exact_Riemann: invalid u_l input!');
assert(isfloat(p_l) && isfinite(p_l) && p_l > 0.0, 'Euler_exact_Riemann: invalid p_l input!');

assert(isfloat(rho_r) && isfinite(rho_r) && rho_r > 0.0, 'Euler_exact_Riemann: invalid rho_r input!');
assert(isfloat(u_r) && isfinite(u_r), 'Euler_exact_Riemann: invalid u_r input!');
assert(isfloat(p_r) && isfinite(p_r) && p_r > 0.0, 'Euler_exact_Riemann: invalid p_r input!');

assert(isfloat(gamma) && isfinite(gamma) && gamma > 1.0, 'Euler_exact_Riemann: invalid gamma input!');
assert(isfloat(tol) && isfinite(tol) && tol > 0.0 && tol <= 1.0, 'Euler_exact_Riemann: invalid tol input!');

% begin computation: pressure-based

a_l = sqrt(gamma*p_l/rho_l);
a_r = sqrt(gamma*p_r/rho_r);
delta_u = u_r - u_l;
gmt = (gamma-1.0)/2.0;
assert(delta_u < (a_l + a_r)/gmt, 'Euler_exact_Riemann: vacuum middle states unsupported!');

% Two-rarefaction initialization (positivity guaranteed)
p = ((a_l + a_r - 0.5*(gamma-1.0)*(u_r - u_l))/(a_l/p_l^((gamma-1.0)/(2*gamma)) + a_r/p_r^((gamma-1.0)/(2*gamma))))^(2*gamma/(gamma-1.0));

% use an iteration to enforce that p grows in the next iteration, according to the characteristics of Newton-Raphson iterations on concave functions
while F(p, rho_l,u_l,p_l,a_l, rho_r,u_r,p_r,a_r, gamma) > 0.0
    p = 0.8*p;
end

% Newton-Raphson iterations on concave functions, guaranteed to converge as long as p keeps positive
flag = true;
while flag
    if p <= 0.0
        error('Euler_exact_Riemann_kernel: iteration error: negative pressure!');
    end
    [Fv, dF] = F(p, rho_l,u_l,p_l,a_l, rho_r,u_r,p_r,a_r, gamma);
    delta_p = Fv/dF;
    p_new = p - delta_p;
    rel_chg = abs(delta_p)/abs(0.5*(p + p_new));
    disp(rel_chg);
    if rel_chg < tol
        flag = false;
    end
    p = p_new;
end
p_m = p;
[fl, ~] = f(p, rho_l, p_l, a_l, gamma);
[fr, ~] = f(p, rho_r, p_r, a_r, gamma);
u_m = 0.5*(u_l + u_r) + 0.5*(fr - fl);

if p_m >= p_l
    % left-shock
    theta = p_m/p_l;
    rat2 = (gamma-1.0)/(2*gamma);
    rat3 = 1.0 - rat2;
    rat1 = rat2 / rat3;
    rho_ml = ((rat1 + theta)/(rat1 * theta + 1.0)) * rho_l;
    a_ml = sqrt(gamma*p_m/rho_ml);
    S_l = u_l - sqrt(rat2 + rat3*theta) * a_l;
else
    % left-rarefaction
    theta = p_m/p_l;
    rho_ml = rho_l*(theta^(1.0/gamma));
    a_ml = sqrt(gamma*p_m/rho_ml);
    S_l = [u_l - a_l, u_m - a_ml];
end

if p_m >= p_r
    % right-shock
    theta = p_m/p_r;
    rat2 = (gamma-1.0)/(2*gamma);
    rat3 = 1.0 - rat2;
    rat1 = rat2 / rat3;
    rho_mr = ((rat1 + theta)/(rat1 * theta + 1.0)) * rho_r;
    a_mr = sqrt(gamma*p_m/rho_mr);
    S_r = u_m + sqrt(rat2 + rat3/theta) * a_mr;
else
    % right-rarefaction
    theta = p_m/p_r;
    rho_mr = rho_r*(theta^(1.0/gamma));
    a_mr = sqrt(gamma*p_m/rho_mr);
    S_r = [u_m + a_mr, u_r + a_r];
end

end

function [val, drv] = f(p, rho0, p0, a0, gamma)
% compute the difference of velocity of two states that are connected by a 1- or 3- wave
% return both the value and the derivative
% should be robust
AA = (2.0/(gamma+1.0)) / rho0;
BB = ((gamma-1.0)/(gamma+1.0)) * p0;
if p >= p0
    % shock
    dp = p - p0;
    pp = p + BB;
    tt = sqrt(AA/pp);
    val = dp*tt;
    drv = tt*(1.0 - 0.5*dp/pp);
else
    % rarefaction
    pr = p/p0;
    theta = (pr)^((gamma-1.0)/(2*gamma));
    val = (2.0/(gamma-1.0)) * a0 *(theta - 1.0);
    drv = a0/(gamma*p) * theta;
end
end

function [val, drv] = F(p, rho_l,u_l,p_l,a_l, rho_r,u_r,p_r,a_r, gamma)
[fl, dfl] = f(p, rho_l, p_l, a_l, gamma);
[fr, dfr] = f(p, rho_r, p_r, a_r, gamma);
val = fl + fr + (u_r - u_l);
drv = dfl + dfr;
end