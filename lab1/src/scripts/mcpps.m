%% Точки равновесия, якобиан, собственные числа
syms x1 x2

% change sys (1)
dx1 = x2;
dx2 = -x1+x2*(1-x1^2+0.1*x1^4);

eqns = [dx1 == 0, dx2 == 0];
S = solve(eqns, [x1, x2], 'Real', true);

x1_eq = double(S.x1);
x2_eq = double(S.x2);

J = jacobian([dx1; dx2], [x1, x2]);

for k = 1:length(x1_eq)
    Jk = double(subs(J, [x1, x2], [x1_eq(k), x2_eq(k)]));
    lambdak = eig(Jk);
    fprintf('Точка равновесия (%.3f, %.3f):\n', x1_eq(k), x2_eq(k));
    disp(Jk)
    disp(lambdak)
end

% sys 4 only
% x1=r*cos(theta);
% x2=r*sin(theta);
% 
% dx1=-x1+2*x1^3+x2;
% dx2=-x1-x2;
%  
% r_dot     = simplify((x1*dx1 + x2*dx2)/r);
% theta_dot = simplify((x1*dx2 - x2*dx1)/r^2);
% 
% disp('r_dot = ')
% pretty(r_dot)
% 
% disp('theta_dot = ')
% pretty(theta_dot)

%% Численное построение фазового портрета
clc; clear; close all;

syms x1 x2

% change sys (2)
dx1_sym = x2;
dx2_sym = -x1+x2*(1-x1^2+0.1*x1^4);

eqns = [dx1_sym == 0, dx2_sym == 0];
S = solve(eqns, [x1, x2], 'Real', true);

x1_eq = double(S.x1);
x2_eq = double(S.x2);

% change sys (3)
[x1g, x2g] = meshgrid(-3:0.2:3, -3:0.2:3);
dx1 = x2g;
dx2 = -x1g+x2g.*(1-x1g.^2+0.1.*x1g.^4);

L = sqrt(dx1.^2 + dx2.^2);
dx1n = dx1 ./ (L + eps);
dx2n = dx2 ./ (L + eps);

figure;
h1 = quiver(x1g, x2g, dx1n, dx2n, 'r');
hold on;
grid on;
xlabel('x1');
ylabel('x2');
title('Фазовый портрет системы');
% axis equal;

% change sys (4)
f = @(t,X) [X(2);
            -X(1)+X(2)*(1-X(1)^2+0.1*X(1)^4)];

inits = [-0.1 0];
tspan = [0 16];

for i = 1:size(inits,1)
    [t,X] = ode15s(f, tspan, inits(i,:));
    h2 = plot(X(:,1), X(:,2), 'b', 'LineWidth', 0.5);
end

h3 = plot(x1_eq, x2_eq, 'ko', 'MarkerFaceColor','g', 'MarkerSize',6);

legend([h1, h2, h3], ...
       {'Поле направлений','Траектории','Точки равновесия'}, ...
       'Location','best');
axis([-3 3 -3 3]);
