import sympy as sp


t = sp.Symbol('t', real=True)
k1, k2, k3 = sp.symbols('k1 k2 k3', real=True, positive=True)

x1 = sp.Function('x1')(t)
x2 = sp.Function('x2')(t)
x3 = sp.Function('x3')(t)

z2 = (
    x3
    - (x2*sp.sin(x1) - sp.sin(x1)*sp.cos(x1)
       - k1*x2 + k1*sp.cos(x1)
       - k2*x2 + k2*sp.cos(x1)
       + k2*k1*x1)
)

phi3 = (1 / (2 - sp.sin(x3))) * (
    -x1*x3
    + x1*sp.sin(x1) + x3*sp.sin(x1)
    + x2*sp.cos(x1)**2 - x2**2*sp.cos(x1)
    - sp.cos(x1)*sp.cos(2*x1) + x2*sp.cos(2*x1)
    - k1*x1 - k1*x3 - (k1/2)*sp.sin(2*x1) + k1*x2*sp.sin(x1)
    - k2*x1 - k2*x3 - (k2/2)*sp.sin(2*x1) + k2*x2*sp.sin(x1)
    + k2*k1*sp.cos(x1) - k2*k1*x2 - k3*z2
)

phi3_dot = sp.diff(phi3, t)

phi3_dot_simplified = sp.simplify(phi3_dot)
sp.print_latex(phi3_dot_simplified)
