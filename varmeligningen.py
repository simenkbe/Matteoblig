import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def solve_heat_equation(alpha, dx, dy, dt, xmax, ymax, tmax):
    nx, ny = int(xmax/dx) + 1, int(ymax/dy) + 1
    nt = int(tmax/dt)
    u = np.zeros((ny, nx))
    u[ny//2, nx//2] = 1000  #hotspot

    #Grensebetingeler
    u[:, 0] = 0
    u[:, -1] = 0
    u[0, :] = 0
    u[-1, :] = 0

    cx = alpha * dt / dx**2
    cy = alpha * dt / dy**2

    history = []  #lagrer verdier

    for _ in range(nt):
        u_new = np.copy(u)
        for j in range(1, ny-1):
            for i in range(1, nx-1):
                u_new[j, i] = u[j, i] + cx * (u[j, i+1] - 2 * u[j, i] + u[j, i-1]) + cy * (u[j+1, i] - 2 * u[j, i] + u[j-1, i])
        u = u_new
        history.append(u.copy())  

    return history

# Parameters
alpha = 0.01  # termisk diffusivitet
dx = dy = 0.1  # romlig steglengde
dt = 0.01  # time step
xmax = ymax = 10  
tmax = 2.0  #tiden for simulering


history = solve_heat_equation(alpha, dx, dy, dt, xmax, ymax, tmax)


fig, ax = plt.subplots()
cax = ax.imshow(history[0], cmap='hot', interpolation='nearest', origin='lower', vmin=0, vmax=1000)
fig.colorbar(cax)

def update(frame):
    cax.set_data(history[frame])
    return cax,


ani = FuncAnimation(fig, update, frames=len(history), interval=50, blit=True)
plt.title('Temperature distribution over time')
plt.show()