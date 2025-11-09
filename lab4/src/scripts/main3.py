import numpy as np
import matplotlib.pyplot as plt

k1, k2, k3, k4 = 1.0, 1.0, 1.0, 1.0

def controller3(x1, x2, x3, x4, phi3_prev, dt):
    z2 = x3 - (x2*np.sin(x1) - np.sin(x1)*np.cos(x1)
               - k1*x2 + k1*np.cos(x1) - k2*x2 + k2*np.cos(x1)
               + k2*k1*x1)
    
    phi3 = (-x1*x3 + x1*np.sin(x1) + x3*np.sin(x1)
            + x2*np.cos(x1)**2 - x2**2*np.cos(x1)
            - np.cos(x1)*np.cos(2*x1) + x2*np.cos(2*x1)
            - k1*x1 - k1*x3 - (k1/2)*np.sin(2*x1)
            + k1*x2*np.sin(x1) - k2*x1 - k2*x3
            - (k2/2)*np.sin(2*x1) + k2*x2*np.sin(x1)
            + k2*k1*np.cos(x1) - k2*k1*x2 - k3*z2)
    
    if phi3_prev is None:
        phi3_dot = 0.0
    else:
        phi3_dot = (phi3 - phi3_prev) / dt
    
    z3 = x4 - phi3
    u = 0.5 * (-x2 * x3 + phi3_dot - k4 * z3)
    
    return u, phi3

def dynamics3(t, x, dt, phi3_prev):
    x1, x2, x3, x4 = x
    u, phi3 = controller3(x1, x2, x3, x4, phi3_prev, dt)
    
    dx1 = np.cos(x1) - x2
    dx2 = x1 + x3
    dx3 = x1 * x3 + (2 - np.sin(x3)) * x4
    dx4 = x2 * x3 + 2 * u
    
    return np.array([dx1, dx2, dx3, dx4]), phi3, u

t0, tf, dt = 0, 20, 0.001
t_span = np.arange(t0, tf, dt)
x0 = np.array([0.5, 0.0, 0.0, 0.0])

x_hist = np.zeros((len(t_span), 4))
u_hist = np.zeros(len(t_span))
x = x0.copy()
phi3_prev = None

for i, t in enumerate(t_span):
    dx, phi3, u = dynamics3(t, x, dt, phi3_prev)
    x = x + dx * dt
    x_hist[i, :] = x
    u_hist[i] = u
    phi3_prev = phi3

plt.figure(figsize=(10,8))
plt.subplot(5,1,1)
plt.plot(t_span, x_hist[:,0], label=r'$x_1(t)$', color='tab:blue')
plt.ylabel(r"$x_1$")
plt.title('System state x_1(t)')
plt.grid(True)
plt.legend()
plt.subplot(5,1,2)
plt.plot(t_span, x_hist[:,1], label=r'$x_2(t)$', color='tab:orange')
plt.ylabel(r"$x_2$")
plt.title('System state x_2(t)')
plt.grid(True)
plt.legend()
plt.subplot(5,1,3)
plt.plot(t_span, x_hist[:,2], label=r'$x_3(t)$', color='tab:green')
plt.ylabel(r"$x_3$")
plt.title('System state x_3(t)')
plt.grid(True)
plt.legend()
plt.subplot(5,1,4)
plt.plot(t_span, x_hist[:,3], label=r'$x_4(t)$', color='tab:purple')
plt.ylabel(r"$x_4$")
plt.title('System state x_4(t)')
plt.grid(True)
plt.legend()
plt.subplot(5,1,5)
plt.plot(t_span, u_hist, label=r'$u(t)$', color='tab:red')
plt.ylabel(r"$u$")
plt.xlabel("Time $t$, s")
plt.title('Control u(t)')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
