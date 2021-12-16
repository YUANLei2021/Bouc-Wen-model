import numpy as np
import matplotlib.pyplot as plt

# ------------------ input section ------------
omega = 1.0  # omega of the sdof
gamma = 0.5  # ratio of post-yield to pre-yield
zeta = 0.02  # damping ratio
bouc_wen_a = 1.0  # A, n, a, b are parameters which regulate the shape of the hysteresis loop.
alpha = 0.5   # dz_/dt =  A*dx/dt - alpha*abs(dx/dt)*z*abs(dz/dt)^(n-1)- beta*(dx/dt)*abs(dz/dt)^n
beta = -0.5
n = 1.0
force_amp = 1.0  # amplitude of the force
force_omega = 1.0  # omega of the force

t_total = 200.0  # total time to calculate
dt = 0.01  # delta time

newmark_gamma = 1/2   # parameters of newmark method
newmark_delta = 1/4

# initial conditions
z = 0.0
x = 0.0
v = 0.0
a = 0.0

# ------------------ calculate section ---------------------
z0 = 1/newmark_delta/dt**2
z1 = newmark_gamma/newmark_delta/dt
z2 = 1/newmark_delta/dt
z3 = 1/2/newmark_delta - 1
z4 = newmark_gamma/newmark_delta - 1
z5 = dt/2*(newmark_gamma/newmark_delta-2)
z6 = dt*(1-newmark_gamma)
z7 = newmark_gamma * dt

k = gamma*omega**2
c = 2*zeta*omega
t = 0

t_list = [t]
z_list = [z]
x_list = [x]
v_list = [v]
a_list = [a]

while t <= t_total:
    u = force_amp * np.sin(force_omega*t)
    q = u - (1-gamma)*omega**2*z
    q_ = q + (z0*x + z2*v + z3*a) + c * (z1*x + z4*v + z5*a)
    k_ = k + z0 + z1*(2*zeta*omega)
    x_ = q_/k_
    a_ = z0*(x_ - x) - z2 * v - z3 * a
    v_ = v + z6 * a + z7 * a_

    x_list.append(x_)
    v_list.append(v_)
    a_list.append(a_)

    z_tem = z
    k1 = bouc_wen_a * v - alpha*np.abs(v)*z_tem*(np.abs(z_tem)**(n-1)) - beta*v*(np.abs(z_tem)**n)
    z_tem = z + dt/2*k1
    k2 = bouc_wen_a * v - alpha*np.abs(v)*z_tem*(np.abs(z_tem)**(n-1)) - beta*v*(np.abs(z_tem)**n)
    z_tem = z + dt/2*k2
    k3 = bouc_wen_a * v - alpha * np.abs(v) * z_tem * (np.abs(z_tem) ** (n - 1)) - beta * v * (np.abs(z_tem) ** n)
    z_tem = z + dt * k3
    k4 = bouc_wen_a * v - alpha * np.abs(v) * z_tem * (np.abs(z_tem) ** (n - 1)) - beta * v * (np.abs(z_tem) ** n)
    z_ = z + dt/6*(k1 + 2*k2 + 2*k3 + k4)
    z_list.append(z_)
    x = x_
    a = a_
    v = v_
    z = z_
    t += dt
    t_list.append(t)
# ----------------- plot figure -------------------
plt.figure()
plt.plot(t_list, x_list)
plt.xlabel("t")
plt.ylabel("x")
plt.xlim(0, t_total)
plt.show()

max_point = True
i = len(x_list)
while max_point:
    x_2 = x_list[i - 1]
    x_1 = x_list[i - 2]
    x_0 = x_list[i - 3]

    if x_1 > x_0 and x_1 > x_2:
        print("Dmax: ", x_1)
        max_point = False
    i -= 1
