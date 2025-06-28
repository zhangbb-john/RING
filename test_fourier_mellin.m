% 读取图像
img1 =  imread('cameraman.tif'); % 使用非对称图像

img2 = imrotate(img1, 30, 'crop');  % 模拟旋转30度的图像

figure(1);
subplot(1,2,1); imshow(img1);
subplot(1,2,2); imshow(img2);
% 估计旋转角度
angle = fourier_mellin_rotation_estimation(img1, img2);
fprintf('Estimated rotation angle: %.2f degrees\n', angle);

function estimated_angle = fourier_mellin_rotation_estimation(img1, img2)
    % 输入: img1 和 img2 是两幅灰度图像（img2 是 img1 的旋转版本）
    % 输出: estimated_angle 是估计的旋转角度（单位：度）

    % 转换为灰度图像（如果输入是RGB）
    if size(img1, 3) == 3
        img1 = rgb2gray(img1);
    end
    if size(img2, 3) == 3
        img2 = rgb2gray(img2);
    end

    % 步骤1: 计算图像的傅里叶幅度谱（并移到中心）
    F1 = fft2(double(img1));
    F2 = fft2(double(img2));
    F1_mag = abs(fftshift(F1));
    F2_mag = abs(fftshift(F2));

    % 步骤2: 对幅度谱进行对数极坐标变换
    [rows, cols] = size(img1);
    center = [rows/2, cols/2];  % 图像中心
    radius = min(rows, cols) / 2;  % 极坐标半径
    num_angles = 360;  % 角度分辨率（1度）
    num_radii = radius;  % 径向分辨率

    % 创建对数极坐标网格
    log_polar1 = log_polar_transform(F1_mag, center, radius, num_radii, num_angles);
    log_polar2 = log_polar_transform(F2_mag, center, radius, num_radii, num_angles);
    figure(3);
    subplot(1,2,1); imshow(log_polar1 / max(max(log_polar1)));
    subplot(1,2,2); imshow(log_polar2 / max(max(log_polar2)));
    % 步骤3: 相位相关法估计平移量（即旋转角度）
    [delta_theta, ~] = phase_correlation(log_polar1, log_polar2);

    % 反推旋转角度（注意方向）
    estimated_angle = delta_theta * (360 / num_angles);
    if estimated_angle > 180
        estimated_angle = estimated_angle - 360;
    end
end

%% 对数极坐标变换函数
function log_polar = log_polar_transform(image, center, radius, num_radii, num_angles)
    [rows, cols] = size(image);
    log_polar = zeros(num_radii, num_angles);

    % 对数极坐标采样
    for r = 1:num_radii
        for theta = 1:num_angles
            % 极坐标 -> 笛卡尔坐标
            rho = log(r + 1e-5);  % 对数径向坐标（避免log(0)）
            phi = 2 * pi * (theta - 1) / num_angles;
            x = center(2) + exp(rho) * cos(phi);
            y = center(1) + exp(rho) * sin(phi);

            % 双线性插值
            if x >= 1 && x <= cols && y >= 1 && y <= rows
                log_polar(r, theta) = interp2(image, x, y, 'linear', 0);
            end
        end
    end
end

%% 相位相关法函数
function [delta_x, delta_y] = phase_correlation(img1, img2)
    F1 = fft2(double(img1));
    F2 = fft2(double(img2));
    figure(2);
    subplot(1,2,1); imshow(F1);
    subplot(1,2,2); imshow(F2);
    cross_power_spectrum = (F1 .* conj(F2)) ./ abs(F1 .* conj(F2) + 1e-10);
    correlation = abs(ifft2(cross_power_spectrum));

    % 找到峰值位置
    [max_val, max_idx] = max(correlation(:));
    [delta_y, delta_x] = ind2sub(size(correlation), max_idx);

    % 处理周期边界
    delta_x = delta_x - 1;
    delta_y = delta_y - 1;
    if delta_x > size(img1, 2)/2
        delta_x = delta_x - size(img1, 2);
    end
    if delta_y > size(img1, 1)/2
        delta_y = delta_y - size(img1, 1);
    end
end