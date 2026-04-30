import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import linregress

def F(X, sigma, rho, beta):
    x = X[0]
    y = X[1]
    z = X[2]
    return np.array([sigma*(y - x), x*(rho - z) - y, x*y - beta*z])

def rk4_step(X, h, sigma, rho, beta):
    k_1 = F(X, sigma, rho, beta)
    k_2 = F(X + (h/2)*k_1, sigma, rho, beta)
    k_3 = F(X + (h/2)*k_2, sigma, rho, beta)
    k_4 = F(X + h*k_3, sigma, rho, beta)
    return X + h/6*(k_1 + 2*k_2 + 2*k_3 + k_4)

# Following Lorenz following Saltzman, we shall let our constants be the following (we also pick the same time step of 0.01 used by Lorenz), and integrate between 0 and 50.

sigma = 10
rho = 28
beta = 8/3
h = 0.01
t_end = 50

N = int(t_end / h)

traj = np.zeros((N+1, 3))
traj[0] = np.array([1, 1, 1])

for i in range(N):
    traj[i + 1] = rk4_step(traj[i], h, sigma, rho, beta)

valeur_test = F(np.array([1, 1, 1]), sigma=10, rho=28, beta=8/3)
print("F(1,1,1) =", valeur_test)

x = traj[:, 0]
y = traj[:, 1]
z = traj[:, 2]
t = np.linspace(0, 50, len(traj)) # time vector

# 3D Figure
    
fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection='3d')
ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], lw=0.5, color='b')

fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# (x , y) projection

axes[0].plot(traj[:, 0], traj[:, 1], lw=0.5, color='teal')
axes[0].set_title("(x, y) projection")
axes[0].set_xlabel("x")
axes[0].set_ylabel("y")

# (x, z) projection

axes[1].plot(traj[:, 0], traj[:, 2], lw=0.5, color='crimson')
axes[1].set_title("(x, z) projection")
axes[1].set_xlabel("x")
axes[1].set_ylabel("z")

# (y, z) projection

axes[2].plot(traj[:, 1], traj[:, 2], lw=0.5, color='darkorange')
axes[2].set_title("(y, z) projection")
axes[2].set_xlabel("y")
axes[2].set_ylabel("z")


# Time series 

fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

axes[0].plot(t, traj[:, 0], color='teal', lw=1)
axes[0].set_title("Temporal evolution of state variables")
axes[0].set_ylabel("x(t)")


axes[1].plot(t, traj[:, 1], color='crimson', lw=1)
axes[1].set_ylabel("y(t)")

axes[2].plot(t, traj[:, 2], color='darkorange', lw=1)
axes[2].set_ylabel("z(t)")
axes[2].set_xlabel("Time t")

# Let's dive into the Butterfly effect, where two very close trajectories diverge

t_end = 30
N = int(t_end / h)
t = np.linspace(0, t_end, N+1)
epsilon = 1e-8
x0_1 = np.array([1.0, 1.0, 1.0])
x0_2 = np.array([1.0 + epsilon, 1.0, 1.0])

traj1 = np.zeros((N+1, 3))
traj1[0] = x0_1
traj2 = np.zeros((N+1, 3))
traj2[0] = x0_2

for i in range(N):
    traj1[i+1] = rk4_step(traj1[i], h, sigma, rho, beta)
    traj2[i+1] = rk4_step(traj2[i], h, sigma, rho, beta)

# We then calculate the euclidian distance for each line, i.e. each time t

d = np.linalg.norm(traj2 - traj1, axis=1)

# Let's have a look at it :

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(t, d, color='purple', lw=1.5)
ax.set_yscale('log') # Semi-logarithmic scale on y axis
ax.set_xlabel("Time t")
ax.set_ylabel("d(t) (log scale)")
ax.set_title("Butterfly effect : the divergence of two (very) close trajectories")
ax.grid(True, which="both", ls="--", alpha=0.5)


# Now, let's estimate the Lyapunov exponent via linear regression.

t_start_reg = 14
t_end_reg = 26

idx_start = int(t_start_reg / h)
idx_end = int(t_end_reg / h)

t_reg = t[idx_start:idx_end]
log_d_reg = np.log(d[idx_start:idx_end])

slope, intercept, r_value, p_value, std_err = linregress(t_reg, log_d_reg)
lambda_1 = slope

d_fit = np.exp(intercept + slope * t_reg)
ax.plot(t_reg, d_fit, 'k--', lw=2,
        label=fr'Linear fit: $\lambda_1 \approx {lambda_1:.4f}$')
ax.legend(loc='lower right', fontsize=11)

plt.savefig('butterfly_effect.png', dpi=140, bbox_inches='tight')
plt.show()


print(f"Estimated Lyapunov exponent : {lambda_1:.4f}" )
