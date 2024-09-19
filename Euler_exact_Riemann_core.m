function [S_l,S_r, rho_l,rho_ml,rho_mr,rho_r, u_l,u_m,u_r, p_l,p_m,p_r, a_l,a_ml,a_mr,a_r] = Euler_exact_Riemann_core(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol)
% [S_l,S_r, rho_l,rho_ml,rho_mr,rho_r, u_l,u_m,u_r, p_l,p_m,p_r, a_l,a_ml,a_mr,a_r] 
% = Euler_exact_Riemann_core(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol) 
%   
%   solve the exact Riemann problem of the 1D Euler equation (for ideal polytopic gas only) 
%   Vacuum states are supported. Zero temperature cases are not supported. 
% 
% input: 
%   rho_l, u_l, p_l:    the density, velocity and pressure on the x<0 side
%   rho_r, u_r, p_r:    the density, velocity and pressure on the x>0 side
%   gamma:              the specific heat ratio
%   tol:                relative tolerance in the Newton-Raphson iteration
% 
% output: 
%   S_l:                the 1-rarefaction speed range (if have size 2), or the 1-shock or 1-contact speed (if have size 1) 
%   u_m:                the 2-contact speed or speed of the free boundary (in the presence of vacuum) 
%   S_r:                the 3-rarefaction speed range (if have size 2), or the 3-shock or 3-contact speed (if have size 1) 
%   rho_l, u_l, p_l, a_l:   the quantities left to the 1-wave (vacuum correction may be applied) 
%   rho_ml, u_m, p_m, a_ml: the quantities between the 1-wave and the 2-wave 
%   rho_mr, u_m, p_m, a_mr: the quantities between the 2-wave and the 3-wave 
%   rho_r, u_r, p_r, a_r:   the quantities right to the 3-wave (vacuum correction may be applied) 
% 
% references:
% [1] Eleuterio F. Toro (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction, 3rd eds. Springer-Verlag Berlin Heidelberg. 

% check the inputs

assert(isfloat(rho_l) && isfinite(rho_l) && rho_l >= 0.0, 'Euler_exact_Riemann: invalid rho_l input!');
assert(isfloat(rho_r) && isfinite(rho_r) && rho_r >= 0.0, 'Euler_exact_Riemann: invalid rho_r input!');
assert(isfloat(u_r) && isfinite(u_r), 'Euler_exact_Riemann: invalid u_r input!');
assert(isfloat(u_l) && isfinite(u_l), 'Euler_exact_Riemann: invalid u_l input!');
assert(isfloat(p_l) && isfinite(p_l) && p_l >= 0.0, 'Euler_exact_Riemann: invalid p_l input!');
assert(isfloat(p_r) && isfinite(p_r) && p_r >= 0.0, 'Euler_exact_Riemann: invalid p_r input!');
assert(isfloat(gamma) && isfinite(gamma) && gamma > 1.0, 'Euler_exact_Riemann: invalid gamma input!');
assert(isfloat(tol) && isfinite(tol) && tol > 0.0 && tol <= 1.0, 'Euler_exact_Riemann: invalid tol input!');

% begin computation: pressure-based

% pre-process vacuum cases (zero density) 
if rho_l == 0.0
    p_l = 0.0;
    a_l = 0.0; % according to the isentropic limit 
end
if rho_r == 0.0
    p_r = 0.0;
    a_r = 0.0; % according to the isentropic limit 
end

% We rule out the cases where the density is zero on any side. 

% both vacuum states
if rho_l == 0.0 && rho_r == 0.0
    rho_ml = 0.0;
    rho_mr = 0.0;
    p_m = 0.0;
    a_ml = 0.0;
    a_mr = 0.0;

    % speed of free boundary 
    u_m = 0.5*(u_l + u_r);
    if u_l < u_r
        % artificial rarefaction 
        S_l = [u_l, u_m];
        S_r = [u_m, u_r];
    else
        % artificial shock 
        S_l = u_m;
        S_r = u_m;
    end

    return;
end

% vacuum right state 
if rho_r == 0.0
    rho_ml = 0.0;
    rho_mr = 0.0;
    p_m = 0.0;
    a_ml = 0.0;
    a_mr = 0.0;

    if p_l == 0.0
        a_l = 0.0;
        u_m = u_l; % speed of free boundary 
        S_l = u_l;
    else
        a_l = sqrt(gamma*p_l/rho_l);
        u_m = u_l + (2.0/(gamma-1.0)) * a_l; % speed of free boundary 
        S_l = [u_l - a_l, u_m]; % main rarefaction
    end

    if u_m < u_r
        % artificial rarefaction 
        S_r = [u_m, u_r];
    else
        % artificial shock 
        S_r = u_m;
    end

    return;
end

% vacuum left state
if rho_l == 0.0
    rho_ml = 0.0;
    rho_mr = 0.0;
    p_m = 0.0;
    a_ml = 0.0;
    a_mr = 0.0;

    if p_r == 0.0
        a_r = 0.0;
        u_m = u_r; % speed of free boundary 
        S_r = u_r;
    else
        a_r = sqrt(gamma*p_r/rho_r);
        u_m = u_r - (2.0/(gamma-1.0)) * a_r; % speed of free boundary 
        S_r = [u_m, u_r + a_r]; % main rarefaction
    end

    if u_l < u_m
        % artificial rarefaction 
        S_l = [u_l, u_m];
    else
        % artificial shock 
        S_l = u_m;
    end

    return;
end

% From now on, the densities are positive from both sides. 

% We don't allow zero pressure when the density is positive, because it means absolute zero temperature. 

assert(p_l > 0.0, 'Euler_exact_Riemann: invalid p_l input! Absolute zero temperature detected!');
assert(p_r > 0.0, 'Euler_exact_Riemann: invalid p_r input! Absolute zero temperature detected!');

% From now on, the densities and pressures are positive on both sides. 

a_l = sqrt(gamma*p_l/rho_l);
a_r = sqrt(gamma*p_r/rho_r);

% We need to rule out the case where vacuum is generated in the middle states. 

% generation of vaccum 
if 0.5*(gamma-1.0)*(u_r - u_l) >= a_l + a_r
    % two full rarefaction fans 
    S_l = [u_l - a_l, u_l + (2.0/(gamma-1.0)) * a_l];
    S_r = [u_r - (2.0/(gamma-1.0)) * a_r, u_r + a_r];
    u_m = 0.5*(S_l(2) + S_r(1));
    S_l(2) = min(S_l(2), u_m);
    S_r(1) = max(S_r(1), u_m);

    rho_ml = 0.0;
    rho_mr = 0.0;
    p_m = 0.0;
    a_ml = 0.0;
    a_mr = 0.0;
    return;
end

% From now on, the density and pressure are positive anywhere in the solution. 

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
        error('Euler_exact_Riemann_kernel: iteration error: non-positive pressure!');
    end

    [Fv, dF] = F(p, rho_l,u_l,p_l,a_l, rho_r,u_r,p_r,a_r, gamma);
    delta_p = Fv/dF;
    p_new = p - delta_p;
    rel_chg = abs(delta_p)/abs(0.5*(p + p_new));
    p = p_new;

    disp(rel_chg);
    if rel_chg < tol
        flag = false;
    end
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