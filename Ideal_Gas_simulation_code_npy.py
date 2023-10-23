import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

np.random.seed(3)

#total energy should be constant for any time index
def total_Energy(particles, enable_g=False):
    # mass_1 radius_1 position_2 velocity_2
    if enable_g:
        g = np.array([0, -1.0])
        E_k = (0.5 * particles[:, 0] * np.linalg.norm(particles[:, [4, 5]], axis=1) ** 2).sum()
        E_p = (-particles[:, 0] * g[1] * particles[:, 3]).sum()
        return E_k + E_p
    else:
        return (0.5 * particles[:, 0] * np.linalg.norm(particles[:, [4, 5]], axis=1) ** 2).sum()

def compute_refl_npy(particles, step, size):
    # mass_1 radius_1 position_2 velocity_2
    r = particles[:, 1]
    v = particles[:, 4:6]
    x = particles[:, 2:4]
    projx = step * np.abs(np.dot(v, np.array([1., 0.])))
    projy = step * np.abs(np.dot(v, np.array([0., 1.])))
    Ix = np.logical_or(np.abs(x[:, 0]) - r < projx, np.abs(size - x[:, 0]) - r < projx)
    Iy = np.logical_or(np.abs(x[:, 1]) - r < projy, np.abs(size - x[:, 1]) - r < projy)
    v[Ix, 0] *= -1.
    v[Iy, 1] *= -1.

def compute_coll_npy(particles, step, look2future=False):
    """Compute velocity after collision with another particle."""
    # mass_1 radius_1 position_2 velocity_2
    N = len(particles)
    Mass = particles[:, [0]]
    Size = particles[:, 1]
    Vel = particles[:, 4:6]
    Pos = particles[:, 2:4]

    if not look2future:
        Dist = np.linalg.norm(Pos[:, None] - Pos[None], axis=2) + 1000 * np.eye(N)
        I_fil = Dist < np.abs(Size[:, None] + Size[None])
    else:
        di = Pos[:, None] - Pos[None]
        norm = np.linalg.norm(di, axis=2)
        I_fil = norm - (Size[:, None] + Size[None]) * 1.1 < step * abs(((Vel[:, None] - Vel[None]) * di).sum(axis=2))/norm

    # filter 3-collisions
    I_3col = I_fil.sum(axis=1) >= 2
    call_again = False
    if I_3col.sum() > 0:
        I_fil[I_3col, :] = False
        I_fil[:, I_3col] = False
        call_again = True

    ij = np.stack(np.meshgrid(np.arange(N), np.arange(N)), axis=2)[I_fil]

    m1 = Mass[ij[:, 0]]
    m2 = Mass[ij[:, 1]]
    v1 = Vel[ij[:, 0]]
    v2 = Vel[ij[:, 1]]
    x1 = Pos[ij[:, 0]]
    x2 = Pos[ij[:, 1]]
    di = x2 - x1

    affected = particles[ij[:, 0]]
    rest = particles[np.logical_not(I_fil.any(axis=1))]

    E_a_aff = total_Energy(affected)
    E_a_rest = total_Energy(rest)
    E_a_full = total_Energy(particles)

    affected[:, 4:6] = v1 - 2. * m2 / (m1 + m2) * ((v1 - v2) * di).sum(axis=1, keepdims=True) / (np.linalg.norm(di, axis=1, keepdims=True) ** 2.) * di
    particles[ij[:, 0]] = affected

    E_b_aff = total_Energy(affected)
    E_b_rest = total_Energy(rest)
    E_b_full = total_Energy(particles)

    # The problem here is 3 collision
    if (E_b_aff > E_a_aff * 1.0001) or (E_a_aff > E_b_aff * 1.0001) or (E_b_full > E_a_full * 1.0001) or (E_a_full > E_b_full * 1.0001):
        print(E_a_aff, E_b_aff)
        print(E_a_rest, E_b_rest)
        print(E_a_full, E_b_full)
        assert 0

    return call_again

def compute_step_npy(particles, step, enable_g=False):
    # mass_1 radius_1 position_2 velocity_2
    if enable_g:
        g = np.array([0, -1.])
        particles[:, 2:4] += step * particles[:, 4:6] + 0.5 * g * step ** 2
        particles[:, 4:6] += step * g
    else:
        particles[:, 2:4] += step * particles[:, 4:6]

def solve_step_npy(particles, step, size):
    compute_refl_npy(particles, step, size)
    compute_coll_npy(particles, step)
    compute_step_npy(particles, step)

def init_list_random_npy(N, radius, mass, boxsize):
    # mass_1 radius_1 position_2 velocity_2
    particle_list_npy = np.zeros((N, 6), dtype=np.float64)
    particle_list_npy[:, 0] = mass
    particle_list_npy[:, 1] = radius

    v_mag = np.random.rand(N, 1).astype(np.float64) * 6
    v_ang = np.random.rand(N, 1).astype(np.float64) * 2 * np.pi
    particle_list_npy[:, 4:6] = np.concatenate([v_mag * np.cos(v_ang), v_mag * np.sin(v_ang)], axis=1)
    particle_list_npy[:, 2:4] = radius + np.random.rand(N, 2)*(boxsize-2*radius)

    # filter out all particles with collision (2 particles simulteniusly)
    Size = particle_list_npy[:, 1]
    Pos = particle_list_npy[:, 2:4]
    Dist = np.linalg.norm(Pos[:, None] - Pos[None], axis=2)
    mDist = (Dist + 1000 * np.eye(N)).min(axis=1)
    return particle_list_npy[mDist > Size]

################################################################################################################################
# Visualization of the solution with matplotlib. It use a slider to change the time

# Compute 2d Boltzmann distribution

def update(time):
    i = int(np.rint(time/timestep))

    # Draw Particles as circles
    for j in range(particle_number):
        # mass_1 radius_1 position_2 velocity_2
        circle[j].center = tuple(history[i][j, [2, 3]].tolist())
    hist.clear()
    
    # Graph Particles speed histogram
    vel_mod = np.linalg.norm(history[i][:, [4, 5]], axis=1)
    hist.hist(vel_mod, bins=30, density=True, label="Simulation Data")
    hist.set_xlabel("Speed")
    hist.set_ylabel("Frecuency Density")
    
    # Compute 2d Boltzmann distribution
    E = total_Energy(history[i])
    ax.set_title('Energy =' + str(E))
    Average_E = E/len(history[i])
    k = 1.38064852e-23
    T = 2*Average_E/(2*k)
    m = history[0][0, 0]
    v = np.linspace(0,10,120)
    fv = m*np.exp(-m*v**2/(2*T*k))/(2*np.pi*T*k)*2*np.pi*v
    hist.plot(v, fv, label = "Maxwell–Boltzmann distribution")
    hist.legend(loc ="upper right")

particle_number = 250
boxsize = 200.

# You need a larger tfin and stepnumber to get the equilibrium state. But the computation takes more time.
tfin = 100
stepnumber = 1500
timestep = tfin/stepnumber
particle_list_npy = init_list_random_npy(particle_number, radius = 2, mass = 1, boxsize = 200)
particle_number = len(particle_list_npy)
history = []

# Compute simulation (It takes some time if stepnumber and particle_number are large)
for i in range(stepnumber):
    solve_step_npy(particle_list_npy, timestep, boxsize)
    history.append(np.copy(particle_list_npy))

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
ax.set_xlim([0,boxsize])
ax.set_ylim([0,boxsize])

# Draw Particles as circles
circle = [None]*particle_number
for i in range(particle_number):
    circle[i] = plt.Circle((history[0][i, 2], history[0][i, 3]), history[0][i, 1], ec="black", lw=1.5, zorder=20)
    ax.add_patch(circle[i])

# Graph Particles speed histogram
vel_mod = np.linalg.norm(history[0][:, [4, 5]], axis=1)
hist.hist(vel_mod, bins= 30, density = True, label = "Simulation Data")
hist.set_xlabel("Speed")
hist.set_ylabel("Frecuency Density")


slider_ax = plt.axes([0.1, 0.05, 0.8, 0.05])
slider = Slider(slider_ax,    # the axes object containing the slider
                  't',   # the name of the slider parameter
                  0,    # minimal value of the parameter
                  tfin,       # maximal value of the parameter
                  valinit=0,  # initial value of the parameter
                  color = '#5c05ff'
                 )

slider.on_changed(update)
plt.show()