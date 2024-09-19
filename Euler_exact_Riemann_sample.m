function [S_l, S_m, S_r, xx, rho_xx, u_xx, p_xx, a_xx] = Euler_exact_Riemann_sample(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol, xx_in)
% 
%   solve the exact Riemann problem of the 1D Euler equation (for ideal polytopic gas only) and compute the primitive variables at nodes xx at time T=1. 
%   Here, S_l, S_m, S_r are the wave velocities. (If there is a rarefaction fan, the output is the min and max of its velocity.)
%   The output xx is refined so that critical points are included. 
% 
% references:
% [1] Eleuterio F. Toro (2009). Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction, 3rd eds. Springer-Verlag Berlin Heidelberg. https://doi.org/10.1007/b79761

[S_l,S_r, rho_l,rho_ml,rho_mr,rho_r, u_l,u_m,u_r, p_l,p_m,p_r, a_l,a_ml,a_mr,a_r] = Euler_exact_Riemann_core(rho_l,u_l,p_l, rho_r,u_r,p_r, gamma, tol);
S_m = u_m;

% xx for computation
xx = unique(sort([xx_in(:); S_l(:); S_m(:); S_r(:)])); % normalize and add critical points

rho_xx = nan(numel(xx), 1);
u_xx = nan(numel(xx), 1);
p_xx = nan(numel(xx), 1);
a_xx = nan(numel(xx), 1);

for k = 1: numel(xx)
    if xx(k) <= min(S_l)
        % left most
        rho_xx(k) = rho_l;
        u_xx(k) = u_l;
        p_xx(k) = p_l;
        a_xx(k) = a_l;
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
        a_xx(k) = a_ml;
    elseif xx(k) <= min(S_r)
        % middle-right
        rho_xx(k) = rho_mr;
        u_xx(k) = u_m;
        p_xx(k) = p_m;
        a_xx(k) = a_mr;
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
        a_xx(k) = a_r;
    end
end

end