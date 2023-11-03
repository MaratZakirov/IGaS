import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from matplotlib import cm
import pandas as pd

np.random.seed(3)
g_accel = -0.5
enable_g = True
timestep = 0.05

#total energy should be constant for any time index
def total_Energy(particles):
    # mass_1 radius_1 position_2 velocity_2
    if enable_g:
        g = np.array([0, g_accel])
        E_k = (0.5 * particles[:, 0] * np.linalg.norm(particles[:, [4, 5]], axis=1) ** 2).sum()
        E_p = (-particles[:, 0] * g[1] * particles[:, 3]).sum()
        return E_k + E_p
    else:
        return (0.5 * particles[:, 0] * np.linalg.norm(particles[:, [4, 5]], axis=1) ** 2).sum()


################################################################################################################################
# Visualization of the solution with matplotlib. It use a slider to change the time
# Compute 2d Boltzmann distribution

def update(i):
    #i = int(np.rint(time/timestep))
    i = int(np.rint(i))
    vel_mod = np.linalg.norm(history[i][:, [4, 5]], axis=1)
    vel_color = np.clip(vel_mod/15, 0, 0.9)

    # Draw Particles as circles
    for j in range(particle_number):
        # mass_1 radius_1 position_2 velocity_2
        circle[j].center = tuple(history[i][j, [2, 3]].tolist())
        circle[j].set_color(cm.hot(vel_color[j]))

    hist.clear()
    hist.hist(vel_mod, bins=30, density=True, label="Simulation Data")
    hist.set_xlabel('Speed')
    hist.set_ylabel('Frecuency Density')

    height = history[i][:, 3]
    E_k = 0.5 * history[i][:, 0] * np.linalg.norm(history[i][:, [4, 5]], axis=1) ** 2
    Grad = np.corrcoef(height, E_k)[0, 1]
    hist.set_title('Corr: %0.5f' % Grad)

    # Compute 2d Boltzmann distribution
    E = total_Energy(history[i])
    ax.set_title('Energy =' + str(E))
    Average_E = E/len(history[i])
    k = 1.38064852e-23
    T = 2*Average_E/(2*k)
    m = history[0][0, 0]
    v = np.linspace(0,14,120)
    fv = m*np.exp(-m*v**2/(2*T*k))/(2*np.pi*T*k)*2*np.pi*v
    hist.plot(v, fv, label = "Maxwell–Boltzmann distribution")
    hist.legend(loc ="upper right")

history = pd.read_csv('simulation.csv', sep='\t')
particle_number = history['pnum'].max() + 1
tfin = history['t'].max()
boxsize = (int(round(history['x'].quantile(q=0.999)) + 1), int(round(history['y'].quantile(q=0.999)) + 1))
history = history[['m', 'r', 'x', 'y', 'vx', 'vy', 'col']].values.reshape(-1, particle_number, 7)

# clculate main statistics
hstat = history[len(history)//2:]
speed_stat = np.linalg.norm(hstat[..., [4, 5]], axis=2)
height_stat = hstat[..., 3]
print('Main C:\n', np.corrcoef(height_stat.reshape(-1), speed_stat.reshape(-1)))

plt.hist(speed_stat.reshape(-1), bins=100)
plt.show()

E = total_Energy(history[0])
Average_E = E/len(history[0])
k = 1.38064852e-23
T = 2*Average_E/(2*k)
m = history[0][0, 0]
v = np.linspace(0,10,120)
fv = m*np.exp(-m*v**2/(2*T*k))/(2*np.pi*T*k)*2*np.pi*v

fig = plt.figure(figsize=(12, 6))
ax = fig.add_subplot(1,2,1)

hist = fig.add_subplot(1,2,2)
hist.plot(v, fv, label = "Maxwell–Boltzmann distribution")
hist.legend(loc ="upper right")

plt.subplots_adjust(bottom=0.2,left=0.15)

ax.axis('equal')
ax.axis([-1, 30, -1, 30])
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
ax.set_xlim([0, 200])#boxsize[0]])
ax.set_ylim([0, 200])#boxsize[1]])

# Draw Particles as circles
circle = [None]*particle_number
for i in range(particle_number):
    circle[i] = plt.Circle((history[0][i, 2], history[0][i, 3]), np.clip(history[0][i, 1], a_min=0.5, a_max=2.0),
                           ec="black", lw=1.5, zorder=20, color='r')
    ax.add_patch(circle[i])

# history
border = 80.0
grads = []
for i in range(len(history)):
    # mass_1 radius_1 position_2 velocity_2
    height = history[i][:, 3]
    E_k = 0.5 * history[i][:, 0] * np.linalg.norm(history[i][:, [4, 5]], axis=1) ** 2
    Grad = np.corrcoef(height, E_k)[0, 1]
    grads.append(Grad)

# Graph Particles speed histogram
vel_mod = np.linalg.norm(history[0][:, [4, 5]], axis=1)
hist.hist(vel_mod, bins= 30, density = True, label = "Simulation Data")
hist.set_xlabel("Speed")
hist.set_ylabel("Frecuency Density")

slider_ax = plt.axes([0.1, 0.05, 0.8, 0.05])
slider = Slider(slider_ax,    # the axes object containing the slider
                  't',   # the name of the slider parameter
                  0,    # minimal value of the parameter
                  particle_number,       # maximal value of the parameter
                  valinit=0,  # initial value of the parameter
                  color = '#5c05ff'
                 )

slider.on_changed(update)
plt.show()