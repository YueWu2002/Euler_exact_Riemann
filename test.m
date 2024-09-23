name = {'Sod shock tube', 'Lax shock tube', 'LeBlanc shock tube', 'Double rarefaction I', 'Double shock', 'Double rarefaction II', 'Vacuum middle state', 'Vacuum left state'};
init_l = {[1.0,0.0,1.0], [0.445,0.698,3.528], [2.0,0.0,1e+9], [7.0,-1.0,0.2], [5.99924, 19.5975, 460.894], [1.0, -2.0, 0.4], [1.0, -20.0, 0.4], [0.0, 0.0, 0.0]};
init_r = {[0.125,0.0,0.1], [0.5,0.0,0.571], [1e-3,0.0,1.0], [7.0,1.0,0.2], [5.99242, -6.19633, 46.0950], [1.0, 2.0, 0.4], [1.0, 20.0, 0.4], [1.0, 0.0, 1.0]};
tol = {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14};
gamma = 1.4;
T = {0.5, 0.3, 1e-5, 0.8, 0.07, 0.3, 0.03, 0.12};

plot_type = 'tight'; % 'tight', 'compact', 'loose'

for k = 1: 8
    xx = -1:0.01:1; % the x/t value of evaluation points
    [S_l, S_m, S_r, xx, rho_xx, u_xx, p_xx, a_xx] = Euler_exact_Riemann_sample(init_l{k}(1),init_l{k}(2),init_l{k}(3), init_r{k}(1),init_r{k}(2),init_r{k}(3), gamma, tol{k}, xx/T{k});

    figure(k);
    tiledlayout(2,3, "TileSpacing",'tight', "Padding",'compact');

    nexttile
    plot(T{k}*xx, rho_xx, 'LineWidth',1);
    hold on;
    xline(T{k}*S_l(:), '--', '1-wave');
    xline(T{k}*S_m(:), '--', '2-wave');
    xline(T{k}*S_r(:), '--', '3-wave');
    hold off;
    grid on;
    xlim([-1,1]);
    xlabel('x');
    ylabel('\rho');

    nexttile
    plot(T{k}*xx, u_xx, 'LineWidth',1);
    hold on;
    xline(T{k}*S_l(:), '--', '1-wave');
    xline(T{k}*S_m(:), '--', '2-wave');
    xline(T{k}*S_r(:), '--', '3-wave');
    hold off;
    grid on;
    xlim([-1,1]);
    xlabel('x');
    ylabel('u');

    nexttile
    plot(T{k}*xx, p_xx, 'LineWidth',1);
    hold on;
    xline(T{k}*S_l(:), '--', '1-wave');
    xline(T{k}*S_m(:), '--', '2-wave');
    xline(T{k}*S_r(:), '--', '3-wave');
    hold off;
    grid on;
    xlim([-1,1]);
    xlabel('x');
    ylabel('p');

    nexttile
    plot(T{k}*xx, u_xx + (2/(gamma-1.0))*a_xx, 'LineWidth',1);
    hold on;
    xline(T{k}*S_l(:), '--', '1-wave');
    xline(T{k}*S_m(:), '--', '2-wave');
    xline(T{k}*S_r(:), '--', '3-wave');
    hold off;
    grid on;
    xlim([-1,1]);
    xlabel('x');
    ylabel('u + 2c/(\gamma-1)');

    nexttile
    plot(T{k}*xx, p_xx./rho_xx.^gamma, 'LineWidth',1);
    hold on;
    xline(T{k}*S_l(:), '--', '1-wave');
    xline(T{k}*S_m(:), '--', '2-wave');
    xline(T{k}*S_r(:), '--', '3-wave');
    hold off;
    grid on;
    ytickformat("%.2f");
    xlim([-1,1]);
    xlabel('x');
    ylabel('p/\rho^\gamma');

    nexttile
    plot(T{k}*xx, u_xx - (2/(gamma-1.0))*a_xx, 'LineWidth',1);
    hold on;
    xline(T{k}*S_l(:), '--', '1-wave');
    xline(T{k}*S_m(:), '--', '2-wave');
    xline(T{k}*S_r(:), '--', '3-wave');
    hold off;
    grid on;
    xlim([-1,1]);
    xlabel('x');
    ylabel('u - 2c/(\gamma-1)');

    sgtitle([name{k}, ', T=', num2str(T{k})]);
    
    print(gcf, '-dpng', '-vector', name{k}, '-r300');
end