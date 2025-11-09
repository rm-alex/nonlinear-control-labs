import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


k1 = 2.0
k2 = 3.0

def controller(x1, x2):
    # virtual control
    phi = -np.sin(x1) - x1**2 - k1 * x1
    z = x2 - phi  # = x2 + sin(x1) + x1^2 + k1*x1

    denom = 2.0 + np.sin(x1)
    u = (-x1 - x1**2 - (np.cos(x1) + 2*x1 + k1)*(x2 + np.sin(x1) + x1**2) - k2*z) / denom
    return u

def controller2(x1, x2):
    # z=x2+k1*x1
    # return -2*x1-k1*z+k1**2*x1+k1*x1**3-k2*z
    return -2*x1-k1*x2+k1*x1**3-k2*x2-k2*k1*x1

def dynamics(t, x):
    x1, x2 = x
    u = controller(x1, x2)
    dx1 = x2 + np.sin(x1) + x1**2
    dx2 = x1**2 + (2 + np.sin(x1)) * u
    return [dx1, dx2]

def dynamics2(t, x):
    x1, x2 = x
    u=controller2(x1,x2)
    dx1=x2-x1**3
    dx2=x1+u
    return [dx1, dx2]


t_span = (0, 10)
t_eval = np.arange(0, 10.001, 0.001)
x0 = [1.0, -0.5]

# sol = solve_ivp(dynamics, t_span, x0, t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-10)
sol = solve_ivp(dynamics2, t_span, x0, t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-10)

t = sol.t
x1 = sol.y[0]
x2 = sol.y[1]

# u_arr = np.array([controller(x1[i], x2[i]) for i in range(len(t))])
u_arr = np.array([controller2(x1[i], x2[i]) for i in range(len(t))])

plt.figure(figsize=(12, 8))

plt.subplot(3, 1, 1)
plt.plot(t, x1, label=r'$x_1(t)$', color='tab:blue')
plt.grid(True)
plt.ylabel(r'$x_1$')
plt.legend()
plt.title('System state x_1(t)')

plt.subplot(3, 1, 2)
plt.plot(t, x2, label=r'$x_2(t)$', color='tab:orange')
plt.grid(True)
plt.ylabel(r'$x_2$')
plt.legend()
plt.title('System state x_2(t)')

plt.subplot(3, 1, 3)
plt.plot(t, u_arr, label=r'$u(t)$', color='tab:red')
plt.grid(True)
plt.xlabel('Time $t$, s')
plt.ylabel(r'$u$')
plt.legend()
plt.title('Control u(t)')

plt.tight_layout()
plt.show()

x1[-1], x2[-1], u_arr[-1]
