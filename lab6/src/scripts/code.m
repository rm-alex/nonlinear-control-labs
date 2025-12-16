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
T_end = 5;
t = 0:dt:T_end;

%% 4. Начальные условия
x = zeros(2,length(t));
x(:,1) = [pi/3; 0];

u_phys = zeros(size(t));
F_hat = zeros(size(t));

u_hist = zeros(size(t));

sat = @(v) max(-u_max,min(u_max,v));

%% 5. Основной цикл
F_true_hist = zeros(size(t));
s_hist = zeros(size(t));
for k_step = 2:length(t)

    x1 = x(1,k_step-1);
    x2 = x(2,k_step-1);

    % Оценка обобщённой неизвестной динамики
    % --- dirty derivative для скорости ---
    tau1 = 0.06;     % настройка сглаживания
    tau2 = 0.06;
    
    if k_step == 2
        y_prev  = x(1,1);
        yd_prev = x(2,1);
        yd_hat_prev = 0;
        ydd_hat_prev = 0;
    end
    
    % коэффициенты фильтров
    a1 = (2*tau1 - dt)/(2*tau1 + dt);
    b1 = 2/(2*tau1 + dt);
    a2 = (2*tau2 - dt)/(2*tau2 + dt);
    b2 = 2/(2*tau2 + dt);
    
    % ydot_hat
    yd_hat = a1*yd_hat_prev + b1*(x1 - y_prev);
    % yddot_hat
    ydd_hat = a2*ydd_hat_prev + b2*(yd_hat - yd_prev);
    
    % сохраняем для следующего шага
    y_prev      = x1;
    yd_prev     = yd_hat;
    yd_hat_prev = yd_hat;
    ydd_hat_prev= ydd_hat;
    
    % F_hat с LPF
    tauF  = 0.15;            % сглаживание оценки
    betaF = dt/(tauF+dt);
    if k_step == 2
        F_hat(k_step) = ydd_hat - u_hist(k_step-1);
    else
        F_raw = ydd_hat - u_hist(k_step-1);
        F_hat(k_step) = (1-betaF)*F_hat(k_step-1) + betaF*F_raw;
    end

    % Скользящая переменная
    beta_s = (1 - rho) / 1;   % т.к. r1=1, r2=1-rho → x2 ~ x1^{1-rho}
    s = x2 + lambda * sign(x1) * abs(x1)^beta_s;
    s_hist(k_step) = s;

    % Финитный однородный регулятор
    u = -F_hat(k_step) - k*abs(s)^rho*tanh(s/eps_s);

    u_phys(k_step) = sat(u);
    u_hist(k_step) = u_phys(k_step);

    % Реальная динамика маятника
    F_true_hist(k_step) = -(g/l)*sin(x1) - (b/I)*x2 - (c/I)*sign(x2);
    ddy = F_true_hist(k_step) + u_phys(k_step)/I;

    x(:,k_step) = x(:,k_step-1) + dt*[x2; ddy];
end

t_settle = t(find(abs(x(1,:)) < 1e-2 & abs(x(2,:)) < 1e-2, 1, 'first'));
fprintf('Settling time (2%%): %.3f s\n', t_settle);

%% 6. Графики
figure;
subplot(3,1,1);
plot(t,x(1,:),'LineWidth',1.2);
grid on;
ylabel('y(t)');
title('Угол');

subplot(3,1,2);
plot(t,x(2,:),'LineWidth',1.2);
grid on;
ylabel('y''(t)');
title('Скорость изменения угла');

subplot(3,1,3);
plot(t,u_phys,'LineWidth',1.2);
grid on;
ylabel('u(t)');
xlabel('t (s)');
title('Управление');

figure;
plot(t, F_true_hist, 'LineWidth', 1.2);
hold on;
plot(t, F_hat, '--', 'LineWidth', 1.2);
grid on;
xlabel('t (s)');
ylabel('F(t)');
legend('F_{true}', 'F_{hat}', 'Location', 'best');
title('Обобщённая неизвестная динамика и её оценка');

figure;
plot(t, F_true_hist-F_hat, 'LineWidth', 1.2);
grid on;
xlabel('t (s)');
ylabel('e(t)');
title('Ошибка оценки неизвестной динамики');

figure;
plot(t, s_hist, 'LineWidth', 1.2);
grid on;
xlabel('t (s)');
ylabel('s(t)');
title('Скользящая поверхность');

