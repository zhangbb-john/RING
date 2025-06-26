% 生成测试图像（改用非对称图像）
close all;
im1 = imread('cameraman.tif'); % 使用非对称图像
angle = 45;
im2 = imrotate(im1, angle, 'bilinear', 'crop');
im1 = im1(ceil(size(im1, 1) * 0.2) : floor(size(im1, 1) * 0.8), ceil(size(im1, 2) * 0.2) : floor(size(im1, 2) * 0.8));
im2 = im2(ceil(size(im2, 1) * 0.2) : floor(size(im2, 1) * 0.8), ceil(size(im2, 2) * 0.2) : floor(size(im2, 2) * 0.8));

figure;
subplot(1, 2, 1); imshow(im1);
subplot(1, 2, 2); imshow(im2)
% 拉东变换参数
theta = 0:1:179; % 更小的角度步长（0.1°）

% 计算拉东变换
R1 = radon(im1, theta);
R2 = radon(im2, theta);

% 可视化拉东变换
figure;
subplot(2,2,1); imagesc(theta, 1:size(R1,1), R1); title('原图拉东变换');
subplot(2,2,2); imagesc(theta, 1:size(R2,1), R2); title('旋转图拉东变换');

corr_strengths = zeros(1, length(theta));
theta = 0:2:358;
for i = 1 : length(theta) / 2
    R2_new = R2(:, [(i+1):length(theta), 1 : i]);
    corr_strengths(i) = sum(sum(R1.* R2_new));
    if (mod(i, 10) == 0)
        figure;
        subplot(1,2,1);
        imagesc(theta, 1:size(R1,1), R1); title('原图拉东变换');
        
        subplot(1,2,2)
        imagesc(theta, 1:size(R2_new,1), R2_new);
        title(num2str(theta(i)));
    end
end
[maxv, maxi] = max(corr_strengths);
figure;
plot(theta, corr_strengths, 'r.');%todo: this does not match well with imagination
disp('succeed ?')
theta(maxi)
figure;
% 对角度轴做FFT（加窗减少频谱泄漏）
window = hann(size(R1,1)) * hann(size(R1,2))';
F1 = fft(R1 .* window, [], 2);
F2 = fft(R2 .* window, [], 2);

% 计算互功率谱（增加正则化项）
cross_power = F1 .* conj(F2) ./ (abs(F1 .* F2) + 1e-5);
corr = ifft(cross_power, [], 2);

% 可视化互相关结果
subplot(2,2,3); imagesc(theta, 1:size(corr,1), abs(corr)); 
title('角度轴互相关'); xlabel('角度(°)');

% 找到峰值位置（避免边界干扰）
valid_range = 20:160; % 忽略边界区域
values = abs(corr(:, valid_range));
sum_values = sum(values, 1);
[~, idx] = max(sum_values);

angle_idx = valid_range(idx(1));
estimated_angle = 180 - theta(angle_idx);%todo: i do not know why i need 180-this

% 显示结果
disp(['真实旋转角度: ', num2str(angle), '°']);
disp(['估计旋转角度: ', num2str(estimated_angle), '°']);

% 可视化图像对齐结果
tform = affine2d([cosd(estimated_angle) sind(estimated_angle) 0; ...
                 -sind(estimated_angle) cosd(estimated_angle) 0; 0 0 1]);
im2_registered = imwarp(im2, tform, 'OutputView', imref2d(size(im1)));
subplot(2,2,4); imshowpair(im1, im2_registered, 'blend'); 
title('对齐结果（原图 vs 校正后）');