import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.linalg import norm

plt.close("all")

# parameters

# distances of orbits from the sun in [AU]
sun_radius = 0.0  # lays in the center
jupiter_radius = 5.2
mars_radius = 1.52
r1 = 3.0  # planetoid 1
r2 = 3.276  # planetoid 2
r3 = 3.7  # planetoid 3

R = r1  # variable choosing the planetoid to test

GMs = 1.327e20  # heliocentric gravitational constant [m^3/s^2]
G = 6.674e-11  # gravitational constant [m^3/kg*s^2]
m_sun = 1.989e30  # mass of sun [kg]

# info Jupiter
ecc_jup = 0.0489 #0
semimajor_jup = 5.2044 
a_jup = semimajor_jup
perihelion_jup = a_jup*(1-ecc_jup)
m_jup = 1.8982e27 # mass [kg]

vp_jup = math.sqrt(GMs) * math.sqrt((1+ecc_jup) / (a_jup*(1-ecc_jup)) * (1+(m_jup/m_sun)))

# info Mars
ecc_mars = 0.0934 #0
semimajor_mars = 1.5237 # [10^6 km]
a_mars = semimajor_mars
perihelion_mars = a_mars * (1 - ecc_mars)
m_mars = 6.4171e23 # mass [kg]

vp_mars = math.sqrt(GMs) * math.sqrt((1+ecc_mars) / (a_mars*(1-ecc_mars)) * (1+(m_mars/m_sun)))

# info planetoid
ecc_pl = 0  # circular
semimajor_pl = 2.5 # [10^6 km]
a_pl = semimajor_pl
perihelion_pl = a_pl*(1-ecc_pl)
m_pl = 1e10  # mass is arbitarary since << Msun

# initial parameters
# for planetoid
vp_pl = math.sqrt(GMs) * math.sqrt((1+ecc_pl) / (a_pl*(1-ecc_pl)) * (1+(m_pl/m_sun)))
vx = 0
vy = vp_pl
x = -perihelion_pl
y = 0
# for Jupiter
vx2 = 0
vy2 = vp_jup
x2 = -perihelion_jup
y2 = 0
# for Mars
vx3 = 0
vy3 = vp_mars
x3 = -perihelion_mars
y3 = 0

starting = np.array([
    x2,y2,  #position of jup
    x3,y3,  #position of mars
    x,y,    #position of planetoid
    
    vx2,vy2,  #velocity of jup
    vx3,vy3,  #velocity of mars
    vx,vy     #velocity of planetoid
])

def system(t, Y, G, m1, m2, m3):

    # position vectors of all 3 planets
    x1 = Y[:2]
    x2 = Y[2:4]
    x3 = Y[4:6]

    # calculcating the square of distance vecotrs
    d12 = (x1 - x2)/np.power(np.linalg.norm(x1 - x2),3)
    d13 = (x1 - x3)/np.power(np.linalg.norm(x1 - x3),3)
    d23 = (x2 - x3)/np.power(np.linalg.norm(x2 - x3),3)

    d1s = x1/np.power(np.linalg.norm(x1),3)
    d2s = x2/np.power(np.linalg.norm(x2),3)
    d3s = x3/np.power(np.linalg.norm(x3),3)

    m_sun = 1.989e30  # mass of sun [kg]

    # net forces 
    f1 = G*(-m2*d12 - m3*d13-m_sun*d1s)
    f2 = G*(m1*d12 - m3*d23-m_sun*d2s)
    f3 = G*(m1*d13 + m2*d23-m_sun*d3s)

    # return the right hand side of the equation
    return np.hstack((Y[6:],f1,f2,f3))

t = 2*np.math.pi/(-vy/x)
n_steps = 1000      # the number of steps
n_periods = 100     # number of full circles
ts = np.linspace(0,n_periods*t,n_steps)

y = solve_ivp(lambda t,Y : system(t, Y, G, m_jup, m_mars, m_pl), (0,t*n_periods), starting, t_eval=ts, max_step=t/n_steps)

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(y.y[0],y.y[1],'b')
ax.plot(y.y[2],y.y[3],'g')
ax.plot(y.y[4],y.y[5],'r')
ax.legend(['jupiter','mars','planetoid'])
ax.scatter(0,0,s=200,c='y')
ax.set_aspect("equal")
ax.set(title='Planetoid: r = {} AU after {} cycle(s)'.format(R, n_periods), xlabel='Distance (AU)', ylabel='Distance (AU)')
plt.show()