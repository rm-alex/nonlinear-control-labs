%% Точки равновесия
syms x1 x2

dx1 = -3*x1-x2;
dx2 = 2*x1-x2^3;

eqns = [dx1 == 0, dx2 == 0];
S = solve(eqns, [x1, x2], 'Real', true)
