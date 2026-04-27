import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation
import matplotlib.cm as cm

# ==========================================================
# GEOMETRY PARAMETERS
# ==========================================================
R = 116.1
r = 44.1
n_paletas = 3

gamma_ini = 75
gamma_fin = 105

# ==========================================================
# DERIVED GEOMETRY
# ==========================================================
R_ext = math.sqrt(R**2 + 3*r**2)
R_int = math.sqrt(R**2 - r**2)

x_int = (R_ext**2 - R**2 + r**2)/(2*r)
y_int = math.sqrt(R_ext**2 - x_int**2)

theta = math.degrees(math.atan2(y_int, x_int))
alpha = math.degrees(math.atan2(y_int, x_int-r))
beta  = math.degrees(math.atan(R_int/r))
eta   = math.degrees(math.asin((r/R_ext)*math.sin(math.radians(180-beta))))

desfase = 360/n_paletas

# ==========================================================
# MESH PARAMETERS
# ==========================================================
m = 10   # rows
n = 15   # columns

# ==========================================================
# FIGURE
# ==========================================================
fig, ax = plt.subplots(figsize=(8,8))
ax.set_aspect('equal')
ax.set_xlim(-R_ext-R, R_ext+R)
ax.set_ylim(-R_ext-R, R_ext+R)

# ==========================================================
# STATIC STATOR
# ==========================================================
stator = [
    patches.Arc((0,0),2*R_ext,2*R_ext,theta1=theta,theta2=180-theta,color='black'),
    patches.Arc((0,0),2*R_ext,2*R_ext,theta1=180+theta,theta2=360-theta,color='black'),
    patches.Arc(( r,0),2*R,2*R,theta1=-alpha,theta2=alpha,color='black'),
    patches.Arc((-r,0),2*R,2*R,theta1=180-alpha,theta2=180+alpha,color='black'),
    patches.Arc(( r,0),2*R,2*R,
                theta1=math.degrees(math.acos(-r/R)),
                theta2=360-math.degrees(math.acos(-r/R)),color='black'),
    patches.Arc((-r,0),2*R,2*R,
                theta1=360-math.degrees(math.acos(r/R)),
                theta2=math.degrees(math.acos(r/R)),color='black')
]

for s in stator:
    ax.add_patch(s)

# ==========================================================
# PALETAS + ROTORS
# ==========================================================
paletas = []
rotores_int = []
rotores_ext = []

colormap = cm.get_cmap('tab10', n_paletas)

for i in range(n_paletas):
    color = colormap(i)

    p_der = patches.Arc((0,0),2*r,2*r,color=color,lw=1.3)
    p_izq = patches.Arc((0,0),2*r,2*r,color=color,lw=1.3)
    p_sup = patches.Arc((0,0),2*R,2*R,color=color,lw=1.0)
    p_inf = patches.Arc((0,0),2*R,2*R,color=color,lw=1.0)

    ax.add_patch(p_der); ax.add_patch(p_izq)
    ax.add_patch(p_sup); ax.add_patch(p_inf)
    paletas.append((p_der,p_izq,p_sup,p_inf))

    rotor_i = patches.Arc((0,0),2*R_int,2*R_int,color=color,lw=1.4)
    rotor_e = patches.Arc((0,0),2*R_ext,2*R_ext,color=color,lw=1.4)

    ax.add_patch(rotor_i); ax.add_patch(rotor_e)
    rotores_int.append(rotor_i)
    rotores_ext.append(rotor_e)

# ==========================================================
# MESH LINES
# ==========================================================
mesh_i = [ax.plot([],[],'r-',lw=0.7)[0] for _ in range(m)]
mesh_j = [ax.plot([],[],'b-',lw=0.7)[0] for _ in range(n)]

# ==========================================================
# MESH GENERATION (YOUR VERIFIED LOGIC)
# ==========================================================
def generate_mesh(gamma):

    x = np.zeros((m,n))
    y = np.zeros((m,n))

    # Bottom boundary
    for i in range(n):
        th = 90 + (2*i/(n-1)-1)*math.degrees(math.asin(r/R))
        x[m-1,i] = R*math.cos(math.radians(gamma)) + R*math.cos(math.radians(th))
        y[m-1,i] = R*math.sin(math.radians(gamma)) - R_int + R*math.sin(math.radians(th))

    x[m-1] = x[m-1][::-1]
    y[m-1] = y[m-1][::-1]

    # Top boundary
    for i in range(n):
        th = (gamma-eta) + i*(2*eta/(n-1))
        x[0,i] = R_ext*math.cos(math.radians(th))
        y[0,i] = R_ext*math.sin(math.radians(th))

    x[0] = x[0][::-1]
    y[0] = y[0][::-1]

    # Right boundary
    for j in range(1,m-1):
        th = j/(m-1)*(gamma-beta)
        x[m-1-j,n-1] = R*math.cos(math.radians(gamma)) + r*math.cos(math.radians(th))
        y[m-1-j,n-1] = R*math.sin(math.radians(gamma)) + r*math.sin(math.radians(th))

    # Left boundary
    for j in range(1,m-1):
        th = beta + gamma + (180-(beta+gamma))/(m-1)*j
        x[j,0] = R*math.cos(math.radians(gamma)) + r*math.cos(math.radians(th))
        y[j,0] = R*math.sin(math.radians(gamma)) + r*math.sin(math.radians(th))

    # Initial interior
    for i in range(1,m-1):
        s = i/(m-1)
        x[i,1:-1] = (1-s)*x[0,1:-1] + s*x[m-1,1:-1]
        y[i,1:-1] = (1-s)*y[0,1:-1] + s*y[m-1,1:-1]

    # Laplace solver
    for _ in range(800):
        for i in range(1,m-1):
            for j in range(1,n-1):
                x[i,j] = 0.25*(x[i+1,j]+x[i-1,j]+x[i,j+1]+x[i,j-1])
                y[i,j] = 0.25*(y[i+1,j]+y[i-1,j]+y[i,j+1]+y[i,j-1])

    return x,y

# ==========================================================
# UPDATE FUNCTION
# ==========================================================
def update(frame):
    gamma = gamma_ini + frame
    artists = []

    # Engine geometry
    for i,(p_der,p_izq,p_sup,p_inf) in enumerate(paletas):
        g = gamma + i*desfase
        xp = R*math.cos(math.radians(g))
        yp = R*math.sin(math.radians(g))

        p_der.center = (xp,yp)
        p_der.theta1 = 180+beta+g
        p_der.theta2 = g-beta

        p_izq.center = (xp,yp)
        p_izq.theta1 = beta+g
        p_izq.theta2 = 180+(g-beta)

        p_sup.center = (xp,yp-R_int)
        p_sup.theta1 = 90-math.degrees(math.asin(r/R))
        p_sup.theta2 = 90+math.degrees(math.asin(r/R))

        p_inf.center = (xp,yp+R_int)
        p_inf.theta1 = 270-math.degrees(math.asin(r/R))
        p_inf.theta2 = 270+math.degrees(math.asin(r/R))

        rotor_i = rotores_int[i]
        rotor_e = rotores_ext[i]

        rotor_i.theta1 = g-360/n_paletas+math.degrees(math.atan(r/R_int))
        rotor_i.theta2 = g-math.degrees(math.atan(r/R_int))

        rotor_e.theta1 = g-360/n_paletas+eta
        rotor_e.theta2 = g-eta

        artists += [p_der,p_izq,p_sup,p_inf,rotor_i,rotor_e]

    # Mesh (paleta 0 chamber)
    X,Y = generate_mesh(gamma)

    for i in range(m):
        mesh_i[i].set_data(X[i,:],Y[i,:])
    for j in range(n):
        mesh_j[j].set_data(X[:,j],Y[:,j])

    return artists + mesh_i + mesh_j

# ==========================================================
# ANIMATION
# ==========================================================
frames = gamma_fin - gamma_ini + 1
ani = FuncAnimation(fig, update, frames=frames, interval=120, blit=False)

plt.show()
