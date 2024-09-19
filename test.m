name = {'Sod shock tube', 'Lax shock tube', 'LeBlanc shock tube', 'Double rarefaction', 'an arbitrary problem'};
init_l = {[1.0,0.0,1.0], [0.445,0.698,3.528], [2.0,0.0,1e+9], [7.0,-1.0,0.2], [5.99924, 19.5975, 460.894]};
init_r = {[0.125,0.0,0.1], [0.5,0.0,0.571], [1e-3,0.0,1.0], [7.0,1.0,0.2], [5.99242, -6.19633, 46.0950]};
tol = {1e-14, 1e-14, 1e-14, 1e-14, 1e-14};
gamma = 1.4;
T = {0.5, 0.3, 1e-5, 0.8, 0.07};

for k = 1: 5
    xx = -1:0.001:1; % the x/t value of evaluation points
    [S_l, S_m, S_r, xx, rho_xx, u_xx, p_xx] = Euler_exact_Riemann_sample(init_l{k}(1),init_l{k}(2),init_l{k}(3), init_r{k}(1),init_r{k}(2),init_r{k}(3), gamma, tol{k}, xx/T{k});

    figure(k);
    subplot(2,3,1);

    plot(T{k}*xx, rho_xx, 'LineWidth',1);
    hold on;
    xline(T{k}*[S_l(:); S_m(:); S_r(:)], '--');
    hold off;
    grid on;
    xlabel('x');
    ylabel('\rho');

    subplot(2,3,2);
    plot(T{k}*xx, u_xx, 'LineWidth',1);
    hold on;
    xline(T{k}*[S_l(:); S_m(:); S_r(:)], '--');
    hold off;
    grid on;
    xlabel('x');
    ylabel('u');

    subplot(2,3,3);
    plot(T{k}*xx, p_xx, 'LineWidth',1);
    hold on;
    xline(T{k}*[S_l(:); S_m(:); S_r(:)], '--');
    hold off;
    grid on;
    xlabel('x');
    ylabel('p');

    subplot(2,3,4);
    plot(T{k}*xx, u_xx + (2/(gamma-1.0))*sqrt(gamma*p_xx./rho_xx), 'LineWidth',1);
    hold on;
    xline(T{k}*[S_l(:); S_m(:); S_r(:)], '--');
    hold off;
    grid on;
    xlabel('x');
    ylabel('u + 2c/(\gamma-1)');

    subplot(2,3,5);
    plot(T{k}*xx, p_xx./rho_xx.^gamma, 'LineWidth',1);
    hold on;
    xline(T{k}*[S_l(:); S_m(:); S_r(:)], '--');
    hold off;
    grid on;
    xlabel('x');
    ylabel('p/\rho^\gamma');

    subplot(2,3,6);
    plot(T{k}*xx, u_xx - (2/(gamma-1.0))*sqrt(gamma*p_xx./rho_xx), 'LineWidth',1);
    hold on;
    xline(T{k}*[S_l(:); S_m(:); S_r(:)], '--');
    hold off;
    grid on;
    xlabel('x');
    ylabel('u - 2c/(\gamma-1)');

    sgtitle([name{k}, ', T=', num2str(T{k})]);
end