% -------------------------------------------------------------------------
% 150MHz 磁偶极子瞬时磁场 (XZ切面) 动态辐射与近远区分界线
% -------------------------------------------------------------------------
clear; clc; close all;

%% 1. 物理参数设定
mu_0 = 4 * pi * 1e-7;
v = 3e8; 
f = 150e6;              % 频率: 150MHz
omega = 2 * pi * f;
k = omega / v;
I_amp = 1.0;
a = 0.1;

%% 2. 空间网格 (使用高密度网格 + imagesc 极速渲染)
grid_density = 400; 
max_range = 3.0;
x_vec = linspace(-max_range, max_range, grid_density);
z_vec = linspace(-max_range, max_range, grid_density);
[X, Z] = meshgrid(x_vec, z_vec);

%% 3. 核心阵列预计算
R = sqrt(X.^2 + Z.^2);
cosTheta = Z ./ R;
sinTheta = abs(X) ./ R;
signX = sign(X); signX(signX == 0) = 1;
R(R < 0.4*a) = NaN; % 屏蔽中心奇点

% 计算固定的复数相量部分 (Br 和 Btheta)
delay = exp(-1i * k * R);
Br_comp = (mu_0*a^2/2) .* ( (1./R.^3) + (1i*omega)./(R.^2 * v) ) .* I_amp .* delay .* cosTheta;
Btheta_comp = (mu_0*a^2/4) .* ( (1./R.^3) + (1i*omega)./(R.^2 * v) + (1i*omega)^2./(R .* v^2) ) .* I_amp .* delay .* sinTheta;

%% 4. 初始化图形窗口
fig = figure('Name', 'Magnetic Field Animation', 'Color', 'w', 'Position', [100, 100, 800, 700]);
colormap jet;

% 初始 t=0 的数据
Br_init = real(Br_comp);
Btheta_init = real(Btheta_comp);
BX_init = (Br_init .* sinTheta + Btheta_init .* cosTheta) .* signX;
BZ_init = Br_init .* cosTheta - Btheta_init .* sinTheta;
B_log_init = log10(sqrt(BX_init.^2 + BZ_init.^2) + 1e-12);
finite_B_log = B_log_init(isfinite(B_log_init));
color_limits = [min(finite_B_log(:)), max(finite_B_log(:))];

% 使用底层图像引擎，抛弃所有复杂曲面
h_img = imagesc(x_vec, z_vec, B_log_init);
set(gca, 'YDir', 'normal'); % 修正坐标系方向
hold on;
main_ax = gca;
caxis(color_limits);
h_colorbar_ax = axes('Position', [0.87, 0.16, 0.035, 0.72]);
colorbar_values = linspace(color_limits(1), color_limits(2), 256)';
imagesc(h_colorbar_ax, 1, colorbar_values, colorbar_values);
set(h_colorbar_ax, 'YDir', 'normal', 'XTick', [], ...
    'YAxisLocation', 'right', 'Box', 'on', 'LineWidth', 1.0);
ylabel(h_colorbar_ax, 'log_{10}|B| (T)', 'FontWeight', 'bold');
axes(main_ax);

% 绘制中心线圈
plot([-a, a], [0, 0], 'ws', 'MarkerSize', 6, 'MarkerFaceColor', 'k');


lambda = v / f;                      
r_boundary = lambda / (2 * pi);      % 计算分界线半径
theta_circle = linspace(0, 2*pi, 100);

% 绘制白色虚线
plot(r_boundary * cos(theta_circle), r_boundary * sin(theta_circle), ...
    'w--', 'LineWidth', 2.5);

% 添加带半透明黑底的文字标注
text_str = sprintf('近区/远区分界线 \n r \\approx %.2fm', r_boundary);
text(r_boundary + 0.1, r_boundary + 0.2, text_str, ...
    'Color', 'w', 'FontSize', 12, 'FontWeight', 'bold', ...
    'BackgroundColor', [0 0 0 0.5], 'EdgeColor', 'w');
% ======================================================================

h_title = title('150MHz 磁场动态扩散 (XZ 切面) | t = 0.00 ns', 'FontSize', 15);
xlabel('X 轴 (m)', 'FontWeight', 'bold');
ylabel('Z 轴 (m)', 'FontWeight', 'bold');
axis equal; axis tight;
set(gca, 'Box', 'on', 'LineWidth', 1.2);


%% 5. 动画循环
T_period = 1 / f;
frames = 40; 
cycles = 1; 
gif_file = fullfile(fileparts(mfilename('fullpath')), 'AC_150MHz_B_XYZ_animation.gif');
gif_delay = 0.06;

for i = 1:(frames * cycles)
    if ~ishandle(fig), break; end 
    
    t = (i-1) * (T_period / frames);
    time_factor = exp(1i * omega * t);
    
    % 在循环内直接乘上时间因子提取实部
    Br_inst = real(Br_comp .* time_factor);
    Btheta_inst = real(Btheta_comp .* time_factor);
    BX_inst = (Br_inst .* sinTheta + Btheta_inst .* cosTheta) .* signX;
    BZ_inst = Br_inst .* cosTheta - Btheta_inst .* sinTheta;
    B_log_new = log10(sqrt(BX_inst.^2 + BZ_inst.^2) + 1e-12);
    
    % 极速替换画面数据
    set(h_img, 'CData', B_log_new);
    set(h_title, 'String', sprintf('150MHz 磁场动态扩散 (XZ 切面) | t = %.2f ns', t*1e9));
    
    drawnow;

    frame = getframe(fig);
    [gif_frame, gif_map] = rgb2ind(frame2im(frame), 256);
    if i == 1
        imwrite(gif_frame, gif_map, gif_file, 'gif', ...
            'LoopCount', Inf, 'DelayTime', gif_delay);
    elseif i <= frames
        imwrite(gif_frame, gif_map, gif_file, 'gif', ...
            'WriteMode', 'append', 'DelayTime', gif_delay);
    end

    pause(0.01); 
end

fprintf('PPT-ready animated GIF saved to: %s\n', gif_file);
