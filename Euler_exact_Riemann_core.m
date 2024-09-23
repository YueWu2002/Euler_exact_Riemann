function [S_l,S_r, rho_l,rho_ml,rho_mr,rho_r, u_l,u_m,u_r, p_l,p_m,p_r, a_l,a_ml,a_mr,a_r] = Euler_exact_Riemann_core(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol)
% [S_l,S_r, rho_l,rho_ml,rho_mr,rho_r, u_l,u_m,u_r, p_l,p_m,p_r, a_l,a_ml,a_mr,a_r] 
% = Euler_exact_Riemann_core(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol) 
%   solve the exact Riemann problem of the 1D Euler equation (for ideal polytopic gas only) 
%   
%   Vacuum states are supported. Zero temperature states are not supported. 
%   The physical quantities are: 
%       rho: density, u: velocity, p: pressure, a: sonic speed 
%   The equation of state (EOS) is: p = (gamma-1)*rho*e, where e is the internal energy per unit mass 
% 
% input: 
%   rho_l, u_l, p_l:    the density, velocity and pressure on the x<0 side 
%   rho_r, u_r, p_r:    the density, velocity and pressure on the x>0 side 
%   gamma:              the specific heat ratio 
%   tol:                relative tolerance in the Newton-Raphson iteration 
% 
% output: 
%   S_l:                the 1-rarefaction speed range (size=2), or the 1-shock or 1-contact speed (size=1) 
%   u_m:                the 2-contact speed (size=1), or the speed of the only free boundary in the presence of any vacuum (size=1), or the speed range of a vacuum fan (size=2)
%   S_r:                the 3-rarefaction speed range (size=2), or the 3-shock or 3-contact speed (size=1) 
%   rho_l, u_l, p_l, a_l:   the quantities left to the 1-wave (vacuum correction may be applied) 
%   rho_ml, u_m, p_m, a_ml: the quantities between the 1-wave and the 2-wave 
%   rho_mr, u_m, p_m, a_mr: the quantities between the 2-wave and the 3-wave 
%   rho_r, u_r, p_r, a_r:   the quantities right to the 3-wave (vacuum correction may be applied) 
% 
% guaranteed relation:
%   S_l(1) (< S_l(2)) <= u_m(1) (< u_m(2)) <= S_r(1) (< S_r(2))
%
% references:
% [1] Eleuterio F. Toro (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction, 3rd eds. Springer-Verlag Berlin Heidelberg. 
% 
% author: Yue Wu (URL: https://yuewu2002.github.io/) (Email: yue_wu3@brown.edu)

% check the inputs

assert(isfloat(rho_l) && isfinite(rho_l) && rho_l >= 0.0, 'Euler_exact_Riemann_core: invalid rho_l input!');
assert(isfloat(rho_r) && isfinite(rho_r) && rho_r >= 0.0, 'Euler_exact_Riemann_core: invalid rho_r input!');
assert(isfloat(u_r) && isfinite(u_r), 'Euler_exact_Riemann_core: invalid u_r input!');
assert(isfloat(u_l) && isfinite(u_l), 'Euler_exact_Riemann_core: invalid u_l input!');
assert(isfloat(p_l) && isfinite(p_l) && p_l >= 0.0, 'Euler_exact_Riemann_core: invalid p_l input!');
assert(isfloat(p_r) && isfinite(p_r) && p_r >= 0.0, 'Euler_exact_Riemann_core: invalid p_r input!');
assert(isfloat(gamma) && isfinite(gamma) && gamma > 1.0, 'Euler_exact_Riemann_core: invalid gamma input!');
assert(isfloat(tol) && isfinite(tol) && tol > 0.0 && tol <= 1.0, 'Euler_exact_Riemann_core: invalid tol input!');

% pre-process vacuum cases (zero density) 
if rho_l == 0.0
    % according to the isentropic limit 
    p_l = 0.0;
    a_l = 0.0;
else
    % sonic speed
    a_l = sqrt(gamma*p_l/rho_l);
end
if rho_r == 0.0
    % according to the isentropic limit 
    p_r = 0.0;
    a_r = 0.0;
else
    % sonic speed
    a_r = sqrt(gamma*p_r/rho_r);
end

% left-most speed of the left rarefaction
S_lmin = u_l - a_l;

% right-most speed of the left rarefaction with a zero-pressure right end
S_lmax = u_l + (2.0/(gamma-1.0)) * a_l;

% right-most speed of the right rarefaction
S_rmax = u_r + a_r;

% left-most speed of the right rarefaction part with a zero-pressure left end
S_rmin = u_r - (2.0/(gamma-1.0)) * a_r;

% We rule out the cases where the density is zero anywhere in the solution. 
% (either in the left-initial value or the right-initial value, 
%  or GENERATED in the NON-TRIVIAL middle states from initial states with positive densities)
if rho_l == 0.0 || rho_r == 0.0 || (rho_l > 0.0 && rho_r > 0.0 && S_lmax < S_rmin)

    % vacuum middle states
    rho_ml = 0.0;
    rho_mr = 0.0;
    p_m = 0.0;
    a_ml = 0.0;
    a_mr = 0.0;

    if (rho_l == 0.0 && rho_r == 0.0) || (rho_l > 0.0 && rho_r > 0.0 && S_lmax < S_rmin)
        % both of the initial states are vacuum, or non-vacuum initial states generate a vacuum middle state (STRICT inequality for the 2nd case)

        if S_lmax < S_rmin
            % Case 1. Burgers' rarefaction (when rho_l == 0.0 && rho_r == 0.0)
            % Case 2. speeds of the free boundaries (when rho_l > 0.0 && rho_r > 0.0 && S_lmax < S_rmin)
            u_m = [S_lmax, S_rmin]; % THE ONLY CASE WHERE THERE ARE TWO VALUES!!!
        else
            % Case 1. Burgers' shock or contact (when rho_l == 0.0 && rho_r == 0.0)
            % Case 2. speed of the free boundary (when rho_l > 0.0 && rho_r > 0.0 && S_lmax < S_rmin)
            u_m = 0.5*(S_lmax + S_rmin);
        end
    elseif rho_r == 0.0
        % a non-vacuum left state and a vacuum right initial state

        % speed of the free boundary
        u_m = S_lmax;
    else
        % a non-vacuum right state and a vacuum left initial state

        % speed of the free boundary
        u_m = S_rmin;
    end

    if S_lmin < u_m(1)
        % rarefaction
        S_l = [S_lmin, u_m(1)];
    else
        % shock or contact
        S_l = u_m(1);
    end

    if S_rmax > u_m(end)
        % rarefaction
        S_r = [u_m(end), S_rmax];
    else
        % shock or contact
        S_r = u_m(end);
    end

    disp('Euler_exact_Riemann_core: vacuum detected! Analytical solution applied.');
    return;
end

% From now on, the densities are POSITIVE ANYWHERE in the solution. 

% Notice that on the two ends of a rarefaction or contact wave, if one has zero pressure, then the other one also has zero pressure. 
if p_l == 0.0 && p_r == 0.0
    if u_l <= u_r
        % In this case, we have already ruled out the subcase where u_l < u_r where a NON-TRIVIAL region of vacuum will be generated. 
        % We use '<=' here for the sake of stability. 

        % three contact waves with p_m == 0 (u_l == u_r)
        S_l = u_l;
        S_r = u_r;
        u_m = 0.5*(S_l + S_r); % for stability

        rho_ml = rho_l;
        rho_mr = rho_r;
        p_m = 0.0;
        a_ml = a_l;
        a_mr = a_r;

        disp('Euler_exact_Riemann_core: inexact initialization! Iteration needed.');
        return;
    else
        % two shock waves with p_m > 0.0
        p = sqrt(0.5*(gamma+1.0)) * (u_l - u_r) / (1.0/sqrt(rho_l) + 1.0/sqrt(rho_r));
        should_correct = false;
    end
elseif (p_l == 0.0) || (p_r == 0.0)
    % at least one is positive, which implies p_m > 0.0
    % initial guess
    p = 0.5*max(p_l, p_r);
    should_correct = true;
else
    % p_l > 0.0 && p_r > 0.0
    % two-rarefaction initialization (validity and positivity guaranteed if p_l > 0.0 && p_r > 0.0)
    p = ((a_l + a_r - 0.5*(gamma-1.0)*(u_r - u_l))/(a_l/p_l^((gamma-1.0)/(2*gamma)) + a_r/p_r^((gamma-1.0)/(2*gamma))))^(2*gamma/(gamma-1.0)); % will be more rebost than the one below (e.x. the double rarefaction problem)
    % p = ((0.5*(gamma-1.0)*(S_lmax - S_rmin))/(a_l/p_l^((gamma-1.0)/(2*gamma)) + a_r/p_r^((gamma-1.0)/(2*gamma))))^(2*gamma/(gamma-1.0));
    if p <= p_l && p <= p_r
        should_correct = false;
    else
        should_correct = true;
    end
end

if should_correct
    disp('Euler_exact_Riemann_core: inexact initialization! Iteration needed.');

    % use an iteration to enforce that p grows in the next iteration so that it won't be non-positive in the future, 
    % because the Newton-Raphson iterations on concave functions will provide a monotone sequence 
    while F(p, rho_l,u_l,p_l,a_l, rho_r,u_r,p_r,a_r, gamma) > 0.0
        p = 0.8*p;
    end

    % Newton-Raphson iterations on concave functions, guaranteed to converge as long as p keeps positive
    while true
        if p <= 0.0
            error('Euler_exact_Riemann_core: iteration error: non-positive pressure!');
        end

        [Fv, dF] = F(p, rho_l,u_l,p_l,a_l, rho_r,u_r,p_r,a_r, gamma);
        delta_p = Fv/dF;
        p_new = p - delta_p;
        rel_chg = abs(delta_p)/abs(0.5*(p + p_new));
        p = p_new;

        disp(['Relative change: ', num2str(rel_chg)]);
        if rel_chg < tol
            break;
        end
    end
else
    disp('Euler_exact_Riemann_core: exact initialization! No iteration needed.');
end
fprintf('\n');

p_m = p;
[fl, ~] = f(p_m, rho_l, p_l, a_l, gamma);
[fr, ~] = f(p_m, rho_r, p_r, a_r, gamma);
u_m = 0.5*((u_l - fl) + (u_r + fr)); % for stability

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
    S_l = [S_lmin, u_m - a_ml];
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
    S_r = [u_m + a_mr, S_rmax];
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