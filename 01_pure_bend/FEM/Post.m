% 加载误差数据
data = load('LLL1.mat');

% 提取数据
xy = data.xy;
u_error = data.u_error;
v_error = data.v_error;
sigma_x_error = data.sigma_x_error;
sigma_y_error = data.sigma_y_error;
tau_xy_error = data.tau_xy_error;

% 绘制误差图
figure;

subplot(3,2,1);
scatter(xy(:,1), xy(:,2), 2, u_error, 'filled');
axis equal;
colorbar;
title('u Error');
xlabel('x (m)');
ylabel('y (m)');
colormap('jet');

subplot(3,2,2);
scatter(xy(:,1), xy(:,2), 2, v_error, 'filled');
axis equal;
colorbar;
title('v Error');
xlabel('x (m)');
ylabel('y (m)');
colormap('jet');

subplot(3,2,3);
scatter(xy(:,1), xy(:,2), 2, sigma_x_error, 'filled');
axis equal;
colorbar;
title('\sigma_{x} Error');
xlabel('x (m)');
ylabel('y (m)');
colormap('jet');

subplot(3,2,4);
scatter(xy(:,1), xy(:,2), 2, sigma_y_error, 'filled');
axis equal;
colorbar;
title('\sigma_{y} Error');
xlabel('x (m)');
ylabel('y (m)');
colormap('jet');

subplot(3,2,5);
scatter(xy(:,1), xy(:,2), 2, tau_xy_error, 'filled');
axis equal;
colorbar;
title('\tau_{xy} Error');
xlabel('x (m)');
ylabel('y (m)');
colormap('jet');

% 保存图像
saveas(gcf, fullfile(save_path, 'Error_MATLAB.png'));
