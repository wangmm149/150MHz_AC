% 磁偶极子时变电场 (时谐场 AC) 2D与3D空间分布可视化
clear; clc; close all;

%% 1. 物理参数设定
mu_0 = 4 * pi * 1e-7;   % 真空磁导率 (H/m)
c = 3e8; v = c;         % 波速等于光速 (m/s)
f = 150e6;              % 频率设为 150MHz (为了在图中清晰看到空间波动的驻点)
omega = 2 * pi * f;     % 角频率
k = omega / v;          % 波数
I_amp = 1.0;            % 电流幅值 (A)
a = 0.1;                % 磁偶极子(圆环)半径 (m)

%% 2. 空间网格构建
grid_density = 35;
max_range = 2.0; % 扩大观测范围，以观察电磁波的空间辐射分布
x_vec = linspace(-max_range, max_range, grid_density);
y_vec = linspace(-max_range, max_range, grid_density);
z_vec = linspace(-max_range, max_range, grid_density);
[X, Y, Z] = meshgrid(x_vec, y_vec, z_vec);

%% 3. 坐标转换与奇点处理
R = sqrt(X.^2 + Y.^2 + Z.^2);
Rho = sqrt(X.^2 + Y.^2); 
Rho(Rho == 0) = eps; % 避开 Z 轴奇点

sinTheta = Rho ./ R;
cosPhi = X ./ Rho;
sinPhi = Y ./ Rho;

% 屏蔽极近区奇点
R(R < 0.5*a) = NaN; 

%% 4. 计算时变电场相量 (依据课件 P17/18)
% 延迟相位因子: e^(-jkr)
delay_term = exp(-1i * k * R);

% 电场只有 E_phi 分量 (复数相量)
% E_phi = - (mu_0*a^2 / 4) * [ (j*omega/R^2) + (j*omega)^2/(R*v) ] * I * exp(-jkR) * sin(theta)
E_phi_complex = -(mu_0 * a^2 * I_amp / 4) .* ...
    ( (1i * omega ./ R.^2) - (omega^2 ./ (R .* v)) ) .* delay_term .* sinTheta;

% 提取 t=0 时刻的瞬时物理场 (即取实部)
E_phi = real(E_phi_complex);

%% 5. 投影至直角坐标系 (Ex, Ey, Ez)
EX = -E_phi .* sinPhi;
EY = E_phi .* cosPhi;
EZ = zeros(size(E_phi)); % 电场在 z 方向始终为 0

%% 6. 可视化 (左图3D，右图2D)
figure('Color', 'w', 'Position', [100, 100, 1200, 550]);
theta_coil = linspace(0, 2*pi, 100); % 线圈坐标

% ================= 左侧：3D 电场流线图 =================
subplot(1, 2, 1); hold on; grid on; view(3);

% 绘制线圈
plot3(a * cos(theta_coil), a * sin(theta_coil), zeros(1, 100), 'k-', 'LineWidth', 4);

% 设定 3D 种子点 (沿 X 轴正半轴放置一排点即可，因为流线是圆环)
sx = linspace(0.2, max_range-0.2, 8); 
sy = zeros(size(sx));
sz_levels = [-1.0, -0.5, 0.5, 1.0]; % 在不同的高度层画圈
[SX, SZ] = meshgrid(sx, sz_levels);
SY = zeros(size(SX));

hLines = streamline(X, Y, Z, EX, EY, EZ, SX, SY, SZ, [0.1, 5000]);
set(hLines, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.5); % 经典MATLAB红/橙

title('时变磁偶极子 3D 瞬时电场线 (t=0)', 'FontSize', 14);
xlabel('X 轴 (m)', 'FontWeight', 'bold'); 
ylabel('Y 轴 (m)', 'FontWeight', 'bold'); 
zlabel('Z 轴 (m)', 'FontWeight', 'bold');
axis equal; axis([-max_range max_range -max_range max_range -max_range max_range]);
set(gca, 'Box', 'on', 'LineWidth', 1.2);

% ================= 右侧：2D XY平面截面图 =================
subplot(1, 2, 2); hold on; grid on;

% 提取 z=0 平面 (即矩阵中间的一层)
mid_idx = round(grid_density/2);
X_2D = squeeze(X(:,:,mid_idx));
Y_2D = squeeze(Y(:,:,mid_idx));
EX_2D = squeeze(EX(:,:,mid_idx));
EY_2D = squeeze(EY(:,:,mid_idx));

% 计算电场强度幅值 (用于云图颜色)
E_mag = sqrt(EX_2D.^2 + EY_2D.^2);

% 使用 pcolor 绘制底部的电场强度云图
h = pcolor(X_2D, Y_2D, E_mag); 
set(h, 'EdgeColor', 'none'); % 去除网格线让过渡更平滑
colormap parula;
%colorbar;

% 使用 quiver 绘制电场矢量箭头 (白色的箭头在彩色底图上更清晰)
% 减少箭头数量避免拥挤
skip = 2; 
quiver(X_2D(1:skip:end, 1:skip:end), Y_2D(1:skip:end, 1:skip:end), ...
       EX_2D(1:skip:end, 1:skip:end), EY_2D(1:skip:end, 1:skip:end), ...
       1.5, 'w', 'LineWidth', 1.2);

% 绘制线圈
plot(a * cos(theta_coil), a * sin(theta_coil), 'k-', 'LineWidth', 3);

title('XY 平面 (z=0) 瞬时电场强度云图与矢量方向', 'FontSize', 14);
xlabel('X 轴 (m)', 'FontWeight', 'bold'); 
ylabel('Y 轴 (m)', 'FontWeight', 'bold');
axis equal; axis([-max_range max_range -max_range max_range]);
set(gca, 'Box', 'on', 'LineWidth', 1.2);
