import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


g = 9.81

m_hat = 0.5
l_hat = 0.8
k_hat = 0.2

a = 2.0
beta0 = 15.0
eps = 0.05

x0 = [1.0, 0.0]
Tmax = 10.0


def h_func(t):
    return 0.5*np.sin(2*t)

def sat(x):
    return np.clip(x, -1, 1)

def control(theta, dtheta, s, t):
    Phi_hat = -(k_hat/m_hat)*dtheta - (g/l_hat)*np.sin(theta) + a*dtheta
    u_s = beta0 * sat(s / eps)
    T = m_hat * (l_hat**2) * ( -Phi_hat - u_s )
    return T

def system(t, x):
    theta, dtheta = x

    s = dtheta + a * theta

    T = control(theta, dtheta, s, t)

    m = 0.8
    l = 1.0
    k = 0.15

    h = h_func(t)

    ddtheta = (T/l + m*h*np.cos(theta) - k*l*dtheta - m*g*np.sin(theta)) / (m*l)

    return [dtheta, ddtheta]

def compute_s_u(sol):
    theta, dtheta = sol.y
    s = dtheta + a * theta
    u = np.array([control(theta[i], dtheta[i], s[i], sol.t[i]) for i in range(len(sol.t))])
    return s, u


t_eval = np.linspace(0, Tmax, 4000)
sol = solve_ivp(system, [0, Tmax], x0, t_eval=t_eval)

s_vals, u_vals = compute_s_u(sol)

plt.figure(figsize=(12, 9))

plt.subplot(4, 1, 1)
plt.plot(sol.t, sol.y[0], label='θ(t) (continuous)')
plt.grid()
plt.ylabel("θ(t)")
plt.title("State θ(t)")
plt.legend()

plt.subplot(4, 1, 2)
plt.plot(sol.t, sol.y[1], label='dθ/dt(t) (continuous)')
plt.grid()
plt.ylabel("dθ/dt")
plt.title("Angular velocity dθ/dt(t)")
plt.legend()

plt.subplot(4, 1, 3)
plt.plot(sol.t, s_vals, label='s(t) (continuous)')
plt.grid()
plt.ylabel("s(t)")
plt.title("Slipping s(t)")
plt.legend()

plt.subplot(4, 1, 4)
plt.plot(sol.t, u_vals, label='u(t) (continuous)')
plt.grid()
plt.ylabel("T(t)")
plt.xlabel("t, [s]")
plt.title("Control torque T(t)")
plt.legend()

plt.tight_layout()
plt.show()
