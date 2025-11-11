import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


a = 1.0
beta0 = 0.1
eps = 0.01
Tmax = 10.0
x0 = [1.0, 0.0]

a1 = 0.5
a2 = 1.5


def sat(x):
    return np.clip(x, -1, 1)

def control(x1, x2, s):
    u = -1/3*(a*x2+a*x1*np.sin(x1)+x1*x2)-1/3*(a*x1**2+abs(x1*x2)+beta0)*sat(s/eps) # np.sign for discrete
    return u

def system(t, x):
    x1, x2 = x
    s = x2 + a * x1

    u = control(x1, x2, s)

    dx1 = x2 + a1*x1*np.sin(x1)
    dx2 = a2*x1*x2 + 3*u
    return [dx1, dx2]

def compute_s_u(sol):
    x1, x2 = sol.y
    s = x2 + a * x1
    u = control(x1, x2, s)
    return s, u


t_eval = np.linspace(0, Tmax, 2000)
sol_cont = solve_ivp(system, [0, Tmax], x0, t_eval=t_eval)

s_cont, u_cont = compute_s_u(sol_cont)

plt.figure(figsize=(12, 9))

plt.subplot(4, 1, 1)
plt.plot(sol_cont.t, sol_cont.y[0], label='x_1(t) (continuous)')
plt.ylabel('x_1(t)')
plt.grid(True)
plt.title('System state x_1(t)')
plt.legend()

plt.subplot(4, 1, 2)
plt.plot(sol_cont.t, sol_cont.y[1], label='x_2(t) (continuous)')
plt.ylabel('x_2(t)')
plt.grid(True)
plt.title('System state x_2(t)')
plt.legend()

plt.subplot(4, 1, 3)
plt.plot(sol_cont.t, s_cont, label='s(t) (continuous)')
plt.ylabel('s(t)')
plt.grid(True)
plt.title('Slipping s(t)')
plt.legend()

plt.subplot(4, 1, 4)
plt.plot(sol_cont.t, u_cont, label='u(t) (continuous)')
plt.xlabel('t, [s]')
plt.ylabel('u(t)')
plt.grid(True)
plt.title('Control u(t)')
plt.legend()

plt.tight_layout()
plt.show()
