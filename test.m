clc
clear all
close all

delta = 0.1;
r = 1;

n = 100; % 主体の数
t = 200; % 時間の期間
variance = 0.05; % ホワイトノイズの分散
variance_s = variance * 2;
variance_s1 = variance * 100;

my_list = 1:t;


alpha = 0.5;
% 初期値の設定
size = zeros(t, n);
pie = zeros(t, n);
k = zeros(t, n);
I = zeros(t, n);
theta0 = zeros(t, n);
a0 = zeros(t, n);
s0 = zeros(t, n);
s1 = zeros(t, n);
y = zeros(t, n);
stddev = zeros(t, n);
z0 = ones(t, n);
d_a = zeros(t, n);
rho = 1;

% モデルのシミュレーション
for i = 1:t
    for j = 1:n
        size(i, j) = 1;
        eta = normrnd(0, sqrt(variance_s));
        eta1 = normrnd(0, sqrt(variance_s1));
        epsilon = normrnd(0, sqrt(variance));

        if i > 2
            size(i, j) = size(i-1, j) + 1;
            x = real(d_a(i-1, :));
            x = x(:);
            kde = fitdist(x, 'Kernel');
            density_values = pdf(kde, x);
            density_values = density_values + 1e-10;
            w = dot(log(density_values), density_values) / sum(density_values);

            if j > 50
                if  z0(i-1, j) < -0.2 && z0(i-2, j) < -0.2
                    size(i, j) = size(i-1,j);
                    theta0(i, j) = real(log(-w) + epsilon);
                    s0(i, j) = rho * theta0(i, j) + (1-rho) * theta0(i-1, j) + eta;
                    s1(i, j) = s0(i, j);
                    stddev(i, j) = var(s0(1:i, j), 1);
                    a0(i, j) = (size(i, j) * (stddev(i, j)^-2) * mean(s0(1:i, j)) + variance * a0(i-1, j)) / ((stddev(i, j)^-2) * size(i, j) + variance);
                    z0(i, j) = exp(theta0(i, j) - (theta0(i, j) - a0(i, j))^2);
                    d_a(i, j) = theta0(i, j) - a0(i, j);


                    k(i, j) = 0;
                    pie(i, j) = 0;
                    I(i,j) = 0;
                else
                    theta0(i, j) = log(-w) + epsilon;
                    s0(i, j) = rho * theta0(i, j) + (1-rho) * theta0(i-1, j) + eta;
                    s1(i, j) = 0;

%                     for i = 1:i
%                         s0(i, j) = s0(i, j) - s1(i, j);
%                     end

                    stddev(i, j) = var(s0(1:i, j), 1);
                    a0(i, j) = (size(i, j) * (stddev(i, j)^-2) * mean(s0(1:i, j)) + variance * a0(i-1, j)) / ((stddev(i, j)^-2) * size(i, j) + variance);
                    z0(i, j) = exp(theta0(i, j) - (theta0(i, j) - a0(i, j))^2);
                    d_a(i, j) = theta0(i, j) - a0(i, j);

                    k(i, j) = (1 - delta) * k(i-1, j) + I(i-1,j);
                    I(i,j) = (alpha*z0(i,j))^(1/(1-alpha)) - k(i,j);
                    pie(i, j) = z0(i, j) * k(i, j)^alpha - r*k(i, j)-I(i,j);
                end
            
        else
            if  z0(i-1, j) < -0.2 && z0(i-2, j) < -0.2
                size(i, j) = size(i-1, j);
                theta0(i, j) = real(log(-w) + epsilon);
                s0(i, j) = rho * theta0(i, j) + (1-rho) * theta0(i-1, j) + eta1;
                s1(i, j) = s0(i, j);
                stddev(i, j) = var(s0(1:i, j), 1)/size(i,j);
                a0(i, j) = (size(i, j) * (stddev(i, j)^-2) * mean(s0(1:i, j)) + variance * a0(i-1, j)) / ((stddev(i, j)^-2) * size(i, j) + variance);
                z0(i, j) = exp(theta0(i, j) - (theta0(i, j) - a0(i, j))^2);
                d_a(i, j) = theta0(i, j) - a0(i, j);


                k(i, j) = 0;
                pie(i, j) = 0;
                I(i,j) = 0;
            else
                theta0(i, j) = real(log(-w) + epsilon);
                s0(i, j) = rho * theta0(i, j) + (1-rho) * theta0(i-1, j) + eta1;
                s1(i, j) = 0;

                %for i = 1:i
                    %s0(i, j) = s0(i, j) - s1(i, j);
                %end

                stddev(i, j) = var(s0(1:i, j), 1);
                a0(i, j) = (size(i, j) * (stddev(i, j)^-2) * mean(s0(1:i, j)) + variance * a0(i-1, j)) / ((stddev(i, j)^-2) * size(i, j) + variance);
                z0(i, j) = exp(theta0(i, j) - (theta0(i, j) - a0(i, j))^2);
                d_a(i, j) = theta0(i, j) - a0(i, j);

                k(i, j) = (1 - delta) * k(i-1, j) + I(i-1,j);
                I(i,j) = (alpha*z0(i,j))^(1/(1-alpha)) - k(i,j);
                pie(i, j) = z0(i, j) * k(i, j)^alpha - r*k(i, j)-I(i,j);
            end
            end


        else
            if j < 100
                theta0(i, j) = 0;
                s0(i, j) = real(theta0(i, j) + eta);
                stddev(i, j) = 10 * sqrt(variance_s);
                a0(i, j) = theta0(i, j);
                a0(i, j) = (size(i, j) * (stddev(i, j)^-2) * mean(s0(1:i, j)) + variance * a0(i, j)) / ((stddev(i, j)^-2) * size(i, j) + variance);
                a0(i, j) = s0(i, j);
                z0(i, j) = exp(theta0(i, j) - (theta0(i, j) - a0(i, j))^2);
                d_a(i, j) = theta0(i, j) - a0(i, j);
                k(i, j) = 0;
                pie(i, j) = 0;
                I(i,j) = 0;

            else
                theta0(i, j) = 0;
                s0(i, j) = real(theta0(i, j) + eta1);
                stddev(i, j) = 10 * sqrt(variance_s);
                a0(i, j) = theta0(i, j);
                a0(i, j) = (size(i, j) * (stddev(i, j)^-2) * mean(s0(1:i, j)) + variance * a0(i, j)) / ((stddev(i, j)^-2) * size(i, j) + variance);
                a0(i, j) = s0(i, j);
                z0(i, j) = exp(theta0(i, j) - (theta0(i, j) - a0(i, j))^2);
                d_a(i, j) = theta0(i, j) - a0(i, j);
                k(i, j) = 0;
                pie(i, j) = 0;
                I(i,j) = 0;
            end
        end
        y(i,j) = z0(i, j) * k(i, j)^alpha;
    end
end





%%

time = 1:t;

figure('Position', [0, 0, 1500, 600]);
for j = 1:n
    plot(time, theta0(:, j));
    hold on;
end
hold off;
xlabel('Time');
ylabel('theta0');
title('theta0 Evolution');
legend;

figure('Position', [0, 0, 1500, 600]);
for j = 1:n
    plot(time, s0(:, j));
    hold on;
end
hold off;
xlabel('Time');
ylabel('S0');
title('S0 Evolution');
legend;

figure('Position', [0, 0, 1500, 600]);
for j = 1:n
    plot(time, a0(:, j));
    hold on;
end
hold off;
xlabel('Time');
ylabel('a0');
title('a0 Evolution');
legend;

figure('Position', [0, 0, 1500, 600]);
for j = 1:n
    plot(time, z0(:, j));
    hold on;
end
hold on;
ylabel('z0');
title('z0 Evolution');
legend;
figure('Position', [0, 0, 1500, 600]);

for j = 1:n
    plot(time, sum(k,2));
    hold on;
end
hold on;
ylabel('z0');
title('z0 Evolution');
legend;

for j = 1:n
    plot(time, y(:,j));
    hold on;
end
ylabel('y');
title('y Evolution');
legend;

%%


time = 1:t;

plot(time, z0(:,1));

ylabel('y');
title('y Evolution');
legend;
hold on;

plot(time, exp(theta0(:,1)));

%%


time = 1:t;

plot(time, mean(z0,2));

ylabel('y');
title('y Evolution');
legend;
hold on;

plot(time, mean(exp(theta0),2));

%%
plot(time, sum(k,2));

ylabel('y');
title('y Evolution');
legend;

%%
plot(time, sum(y,2));

ylabel('y');
title('y Evolution');
legend;


%%
q = linspace(min(x), max(x), 100);  % プロットする範囲を指定
pdf = pdf(kde, x);  % 確率密度関数を計算


plot(q, pdf, 'LineWidth', 2);  % プロット
title('Probability Density Function');  % グラフのタイトル
xlabel('x');  % x軸のラベル
ylabel('PDF');  % y軸のラベル


%%

model = arima(1,0,0);

estModel = estimate(model, mean(z0,2));
ARCoefficient = estModel.AR{1};
Constant = estModel.Constant;
disp(['AR Coefficient: ', num2str(ARCoefficient)]);
disp(['Constant: ', num2str(Constant)]);

%%

model = arima(1,0,0);

estModel = estimate(model, z0(:,11));
ARCoefficient = estModel.AR{1};
Constant = estModel.Constant;
disp(['AR Coefficient: ', num2str(ARCoefficient)]);
disp(['Constant: ', num2str(Constant)]);

%%

model = arima(1,0,0);

estModel = estimate(model, real(sum(y,2)));
ARCoefficient = estModel.AR{1};
Constant = estModel.Constant;
disp(['AR Coefficient: ', num2str(ARCoefficient)]);
disp(['Constant: ', num2str(Constant)]);

%%


plot(time, mean(stddev,2));



%%


phi = 0.25598;  % AR係数
initial_condition = 0 ;  % 初期条件
constant = 0.45097;

n = 600;  % データ点の数

% ホワイトノイズの生成
stddev = 0.0069249^(1/2);
mu = 0;
epsilon = stddev.*randn(n,1) + mu;


% AR(1)モデルのデータ生成
data = zeros(n, 1);
data(1) = initial_condition;
for i = 2:n
    data(i) = constant + phi * data(i-1) + epsilon(i) ;
end


% グラフのプロット
plot(data)
title('AR(1)モデル')
xlabel('時刻')
ylabel('データ')


%%

fs = 10; % サンプリング周波数
t = 0:(1/fs):1; % 時間軸

fft_data = fft(mean(z0,2)); % フーリエ変換


N = length(mean(z0,2)); % データ点数
frequency = (0:N-1)*(fs/N); % 周波数軸の生成

A = abs(fft_data(2))/N; % 振幅（ピーク周波数の成分の振幅）
f = frequency(2); % 周波数（ピーク周波数）
phi = angle(fft_data(2)); % 位相（ピーク周波数の位相）

% 近似データの構築
approx_data = A * sin(2*pi*f*t + phi);


plot(approx_data) % プロット
title('Frequency Spectrum') % タイトル
xlabel('Frequency (Hz)') % x軸ラベル
ylabel('Magnitude') % y軸ラベル
ylim([-0.1 10])



