%% Параметры системы
a0 = 0.5;
a1 = 1.3;
b  = 2;

%% Рабочие точки (x1, x2=0)
x1_eq = [1 2 3];
x2_eq = zeros(size(x1_eq));

%% Для хранения PID
Kp = zeros(size(x1_eq));
Ki = zeros(size(x1_eq));
Kd = zeros(size(x1_eq));

%% Линеаризация и синтез PID
for i = 1:length(x1_eq)
    
    % Рабочая точка
    x1_0 = x1_eq(i);
    x2_0 = x2_eq(i);
    u0   = b*x1_0;
    
    % Коэффициент a(x1)
    a = a0 + a1*x1_0^2;
    
    % Линеаризованная система: dx/dt = A x + B u
    A = [0 1;
        -b -a];
    B = [0; 1];
    
    C = [1 0];
    D = 0;
    sys = ss(A,B,C,D);
    
    pid_ctrl = pidtune(sys,'PID'); % PID/+F
    
    Kp(i) = pid_ctrl.Kp;
    Ki(i) = pid_ctrl.Ki;
    Kd(i) = pid_ctrl.Kd;
    
    % Вывод в консоль
    fprintf('x1 = %.2f: Kp=%.2f, Ki=%.2f, Kd=%.2f\n', x1_0, Kp(i), Ki(i), Kd(i));
end
