%% 1. Параметры системы
clear; clc; close all;

m = 0.5;
l = 0.3;
g = 9.81;
b = 0.01;
c = 0.005;
I = m*l^2;

u_max = 1.0;

%% 2. Параметры регулятора
lambda = 2.0;     % коэффициент скользящей поверхности
k = 3.0;          % усиление
rho = 0.5;        % 0<rho<1 (финитность)
eps_s = 0.02;     % сглаживание tanh

%% 3. Время
dt = 1e-3;
T_end = 20;
t = 0:dt:T_end;

%% 4. Начальные условия
x = zeros(2,length(t));
x(:,1) = [pi/3; 0];

u_phys = zeros(size(t));
F_hat = zeros(size(t));

u_hist = zeros(size(t));

sat = @(v) max(-u_max,min(u_max,v));

%% 5. Основной цикл
for k_step = 2:length(t)

    x1 = x(1,k_step-1);
    x2 = x(2,k_step-1);

    %% Оценка обобщённой неизвестной динамики
    if k_step > 2
        ddy_est = (x2 - x(2,k_step-2))/dt;
    else
        ddy_est = 0;
    end

    F_hat(k_step) = ddy_est - u_hist(k_step-1);

    %% Скользящая переменная
    s = x2 + lambda*x1;

    %% Финитный однородный регулятор
    u = -F_hat(k_step) ...
        - lambda*x2 ...
        - k*abs(s)^rho*tanh(s/eps_s);

    u_phys(k_step) = sat(u);
    u_hist(k_step) = u_phys(k_step);

    %% Реальная динамика маятника
    F_true = -(g/l)*sin(x1) - (b/I)*x2 - (c/I)*sign(x2);
    ddy = F_true + u_phys(k_step)/I;

    x(:,k_step) = x(:,k_step-1) + dt*[x2; ddy];
end

%% 6. Графики
figure;
subplot(3,1,1);
plot(t,x(1,:),'LineWidth',1.2);
grid on;
ylabel('y (rad)');
title('Угол');

subplot(3,1,2);
plot(t,x(2,:),'LineWidth',1.2);
grid on;
ylabel('\dot y');

subplot(3,1,3);
plot(t,u_phys,'LineWidth',1.2);
grid on;
ylabel('u');
xlabel('t (s)');
title('Управление');

figure;
plot(t, F_true, 'r', 'LineWidth', 1.2);
hold on;
plot(t, F_hat, '--b', 'LineWidth', 1.2);
grid on;
xlabel('t (s)');
ylabel('F');
legend('F_{true}', '\hat{F}', 'Location', 'best');
title('Обобщённая неизвестная динамика и её оценка');

