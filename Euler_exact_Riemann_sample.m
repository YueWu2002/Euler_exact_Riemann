function [S_l,S_m,S_r, xa, rho_xx,u_xx,p_xx,a_xx] = Euler_exact_Riemann_sample(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol, xx_in)
% [S_l,S_m,S_r, xxangle_refined, rho_xx,u_xx,p_xx,a_xx] = Euler_exact_Riemann_sample(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol, xx_in)
%   solve the exact Riemann problem of the 1D Euler equation (for ideal polytopic gas only) 
%   and compute the physical quantities at nodes xx_in at time T=1. 
%   Here, S_l, S_m, S_r are the wave velocities. (If there is a rarefaction fan, the output is the min and max of its velocity.)
%   The output xa is the array of plotting points at T=1 and is refined so that critical points are included. 
% 
% references:
% [1] Eleuterio F. Toro (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction, 3rd eds. Springer-Verlag Berlin Heidelberg. https://doi.org/10.1007/b79761

[S_l,S_r, rho_l,rho_ml,rho_mr,rho_r, u_l,u_m,u_r, p_l,p_m,p_r, a_l,a_ml,a_mr,a_r] = Euler_exact_Riemann_core(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol);
S_m = u_m;

% xx for computation
xa = xx_in(:);
if numel(S_l) > 1
    xa = [xa(:); S_l(:)];
else
    % This step can be optimized. 
    xa = [xa(:); S_l; S_l + eps(1.0)*S_l];
end
if numel(S_m) > 1
    xa = [xa(:); S_m(:)];
else
    % This step can be optimized. 
    xa = [xa(:); S_m; S_m + eps(1.0)*S_m];
end
if numel(S_r) > 1
    xa = [xa(:); S_r(:)];
else
    % This step can be optimized. 
    xa = [xa(:); S_r; S_r + eps(1.0)*S_r];
end
xa = unique(sort(xa(:))); % add critical points


rho_xx = nan(numel(xa), 1);
u_xx = nan(numel(xa), 1);
p_xx = nan(numel(xa), 1);
a_xx = nan(numel(xa), 1);

for k = 1: numel(xa)
    if xa(k) <= S_l(1)
        % left most
        rho_xx(k) = rho_l;
        u_xx(k) = u_l;
        p_xx(k) = p_l;
        a_xx(k) = a_l;
    elseif xa(k) < S_l(end) % this implies that S_l(1) < S_l(end)
        % left-rarefaction (u is linear and S is constant)
        theta = (xa(k) - S_l(1)) / (S_l(end) - S_l(1));
        u_xx(k) = (1.0-theta) * u_l + theta * (u_m(1));
        % u_xx(k) = ((gamma-1.0)/(gamma+1.0))*u_l + (2.0/(gamma+1.0))*(xa(k) + a_l); % linear in x
        a_xx(k) = (1.0-theta) * a_l + theta*a_ml;
        % a_xx(k) = u_xx(k) - xa(k); % linear in x
        rho_xx(k) = ((1.0-theta) * rho_l^(0.5*(gamma-1.0)) + theta*rho_ml^(0.5*(gamma-1.0)))^(2.0/(gamma-1.0));
        % rho_xx(k) = (((rho_l^gamma)*a_xx(k)^2)/(gamma*p_l))^(1.0/(gamma-1.0)); % (gamma-1)/2-power-linear
        p_xx(k) = rho_xx(k)*a_xx(k)^2/gamma;
    elseif xa(k) <= S_m(1)
        % middle-left
        rho_xx(k) = rho_ml;
        u_xx(k) = u_m(1);
        p_xx(k) = p_m;
        a_xx(k) = a_ml;
    elseif xa(k) < S_m(end) % this implies that S_m(1) < S_m(end)
        % the artificial zero-density rarefaction fan
        rho_xx(k) = 0.0;
        u_xx(k) = xa(k);
        p_xx(k) = 0.0;
        a_xx(k) = 0.0;
    elseif xa(k) <= S_r(1)
        % middle-right
        rho_xx(k) = rho_mr;
        u_xx(k) = u_m(end);
        p_xx(k) = p_m;
        a_xx(k) = a_mr;
    elseif xa(k) < S_r(end) % this implies that S_r(1) < S_r(end)
        % right-rarefaction (u is linear and S is constant)
        theta = (xa(k) - S_r(1)) / (S_r(end) - S_r(1));
        u_xx(k) = (1.0-theta) * u_m(end) + theta * u_r;
        % u_xx(k) = ((gamma-1.0)/(gamma+1.0))*u_r + (2.0/(gamma+1.0))*(xa(k) - a_r);
        a_xx(k) = (1.0-theta) * a_mr + theta*a_r;
        % a_xx(k) = xa(k) - u_xx(k); % linear in x
        rho_xx(k) = ((1.0-theta) * rho_mr^(0.5*(gamma-1.0)) + theta*rho_r^(0.5*(gamma-1.0)))^(2.0/(gamma-1.0));
        % rho_xx(k) = (((rho_l^gamma)*a_xx(k)^2)/(gamma*p_l))^(1.0/(gamma-1.0)); % (gamma-1)/2-power-linear
        p_xx(k) = rho_xx(k)*a_xx(k)^2/gamma;
    else
        % right-most
        rho_xx(k) = rho_r;
        u_xx(k) = u_r;
        p_xx(k) = p_r;
        a_xx(k) = a_r;
    end
end

end