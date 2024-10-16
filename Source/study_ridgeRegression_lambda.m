% 生成随机数据
clear
d = 100;  % 特征维度
r = 10;   % 子空间维度
n = 100;  % 样本数量

% U = randn(d, r);  % 子空间表示矩阵
% z = U(:,r);  % 生成观测向量，加入噪声
A = randn(d, r)*randn(r, d);  % 子空间表示矩阵
z = A(:,d);  % 生成观测向量，加入噪声
[U,~,~] = svd(A);
U = U(:,1:r);
alpha_true = inv(U' * U)* U' * z;
% 0.1 * randn(d, 1)
% 极大值干扰
z(10) = z(10)+100;
z(40) = z(40)+100;
z(80) = z(80)+100;
% 极小值干扰
z(30) = z(30)-100;
z(60) = z(60)-100;
z(100) = z(100)-100;

% 设置不同的 lambda 值
lambda_values = [0, 0.01, 0.1, 1, 10, 100, 200, 300, 400, 500, 600,1000,2000];

% 初始化结果存储变量
alpha_results = zeros(r, length(lambda_values));
alpha_results_o = zeros(r, length(lambda_values));
mse_results = zeros(1, length(lambda_values));

% 进行实验
for i = 1:length(lambda_values)
    lambda = lambda_values(i);
    
    % 使用 Ridge 回归求解
    % alpha = (U' * U + lambda * eye(r)) \ (U' * z);
    alpha = inv(U' * U + lambda * eye(r))* U' * z;
    alpha_results(:, i) = alpha;
    alpha_results_o(:, i) = alpha_true;

    % 计算预测值
    z_hat = U * alpha;
    % 采样的极值，研究对其他正常值的恢复精度影响
%     z_hat([10,40,80]) = z([10,40,80]);
%     z_hat([30,60,100]) = z([30,60,100]);
    z_hat([10,30,40,60,80,100]) = z([10,30,40,60,80,100]);
    
    % 计算均方误差（MSE）
    mse = mean((z - z_hat).^2);
    mse_results(i) = mse;
end

% 绘制结果
figure;
subplot(1, 3, 1);
plot(lambda_values, alpha_results_o, 'o-');
xlabel('\lambda');
ylabel('\alpha_o');
legend('Dimension 1', 'Dimension 2', 'Dimension 3', 'Dimension 4', 'Dimension 5','Dimension 6', 'Dimension 7', 'Dimension 8', 'Dimension 9', 'Dimension 10');
title('\alpha without interupt');

subplot(1, 3, 2);
plot(lambda_values, alpha_results, 'o-');
xlabel('\lambda');
ylabel('\alpha');
legend('Dimension 1', 'Dimension 2', 'Dimension 3', 'Dimension 4', 'Dimension 5','Dimension 6', 'Dimension 7', 'Dimension 8', 'Dimension 9', 'Dimension 10');
title('Effect of \lambda on \alpha');

subplot(1, 3, 3);
plot(lambda_values, mse_results, 'o-');
xlabel('\lambda');
ylabel('MSE');
title('Effect of \lambda on MSE');