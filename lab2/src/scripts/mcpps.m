%% Точки равновесия
syms x1 x2

dx1 = -3*x1-x2;
dx2 = 2*x1-x2^3;

eqns = [dx1 == 0, dx2 == 0];
S = solve(eqns, [x1, x2], 'Real', true)

%% 3

A=[0 1;2 0];
B=[0;1];

[P, J] = jordan(A)
B_J = inv(P) * B

a = 2;

cvx_begin sdp
variable P1(2,2) symmetric
variable Y1(1,2)
P1 > 0.0001*eye(2);
P1*A' + A*P1 + 2*a*P1 + Y1'*B'+ B*Y1 <= 0;
cvx_end

K=Y1*inv(P1)

eig(A+B*K)

set_param('sim1/A', 'Gain', mat2str(A));
set_param('sim1/B', 'Gain', mat2str(B));
set_param('sim1/K', 'Gain', mat2str(K));