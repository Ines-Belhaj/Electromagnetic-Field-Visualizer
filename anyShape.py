import numpy as np
import scipy as sp
from scipy.integrate import quad
import matplotlib.pyplot as plt
import sympy as smp
import plotly.graph_objects as go
from IPython.display import HTML

L = 0
# define the variable where to calculate the electric field on the charged shape
t = smp.symbols('t', positive=True)
# define our variables
x, y, z = smp.symbols('x y z')
# define the vector position where you are 
r = smp.Matrix([x, y, z])
# define where the charged shape is located
print('Enter the coordinates of the location where is the charged shape')
xpos = input("enter the coeficient of x position\n")
ypos = input("enter the coeficient of y position\n")
zpos = input("enter the coeficient of z position\n")

cos = input("does the x position use the cos or sin ?\n")
if cos =="cos":
    xpos1 = smp.cos(float(xpos)*t)
elif cos =="sin":
    xpos1 = smp.sin(float(xpos)*t)
else:
    xpos1 = float(xpos)*t

cos = input("does the y position use the cos or sin or none?\n")
if cos =="cos":
    ypos1 = smp.cos(float(ypos)*t)
elif cos =="sin":
    ypos1 = smp.sin(float(ypos)*t)
else:
    ypos1 = float(ypos)*t

cos = input("does the z position use the cos or sin or none ?\n")
if cos =="cos":
    zpos1 = smp.cos(float(zpos)*t)
elif cos =="sin":
    zpos1 = smp.sin(float(zpos)*t)
else:
    zpos1 = float(zpos)*t

# initialize the vector of the charged shape location
r_p = smp.Matrix([xpos1, ypos1, zpos1])
# calculate the distance and initilaize it in a vector
sep = r - r_p

# get the total charge
q = input("Enter the total charge of the shape")
Q = float(q)

# diferntiate the r' over t 
dr_pdt = smp.diff(r_p, t).norm().simplify()

# calculate the charge density lamda for Q = 1
if (xpos1 == float(xpos)*t) and (ypos1 == float(ypos)*t) and (zpos1 == float(zpos)*t):
    L = float(input('enter the length of the rod'))
    lam = smp.integrate(dr_pdt, (t, 0, L/2))
else: 
    lam = smp.integrate(dr_pdt, (t, 0, 2*smp.pi))

lam = float(lam)

print("charge density = " + str(lam))

integrand = lam * sep/sep.norm()**3 * dr_pdt
dExdt = smp.lambdify([t, x, y, z], integrand[0])
dEydt = smp.lambdify([t, x, y, z], integrand[1])
dEzdt = smp.lambdify([t, x, y, z], integrand[2])
def E(x, y, z):
    return np.array([quad(dExdt, 0, 2*np.pi, args=(x, y, z))[0],
                     quad(dEydt, 0, 2*np.pi, args=(x, y, z))[0],
                     quad(dEzdt, 0, 2*np.pi, args=(x, y, z))[0]])

def Erod(x, y, z):
    xField = quad(dExdt, 0, L/2, args=(x, y, z))[0]
    yField = quad(dEydt, 0, L/2, args=(x, y, z))[0]
    zField = quad(dEzdt, 0, L/2, args=(x, y, z))[0]

    return np.array([xField, yField, zField])

x = np.linspace(-int(xpos), int(xpos), 10)
y = np.linspace(-int(ypos), int(ypos), 10)
z = np.linspace(0, 2*np.pi, 10)
if (xpos1 == float(xpos)*t) and (ypos1 == float(ypos)*t) and (zpos1 == float(zpos)*t):
    z = np.linspace(-int(zpos),int(zpos), 10)

xv, yv, zv = np.meshgrid(x, y, z)




if L != 0:
    E_field = np.vectorize(Erod, signature='(),(),()->(n)')(xv, yv, zv)
    Ex = E_field[:,:,:,0]
    Ey = E_field[:,:,:,1]
    Ez = E_field[:,:,:,2]
else:
    E_field = np.vectorize(E, signature='(),(),()->(n)')(xv, yv, zv)
    Ex = E_field[:,:,:,0]
    Ey = E_field[:,:,:,1]
    Ez = E_field[:,:,:,2]


plt.hist(Ex.ravel(), bins=100, histtype='step',label='Ex')
plt.hist(Ey.ravel(), bins=100, histtype='step',label='Ey')
plt.hist(Ez.ravel(), bins=100, histtype='step',label='Ez')
plt.legend()
plt.xlabel('Electric Field Magnitude')
plt.ylabel('Frequency')
plt.show()

#D = float(input("Emax = "))

E_max = 150
Ex[Ex>E_max] = E_max
Ey[Ey>E_max] = E_max
Ez[Ez>E_max] = E_max

Ex[Ex<-E_max] = -E_max
Ey[Ey<-E_max] = -E_max
Ez[Ez<-E_max] = -E_max

if (xpos1 == float(xpos)*t) and (ypos1 == float(ypos)*t) and (zpos1 == float(zpos)*t):
    tt = np.linspace(0, L/2, 100)
else:
    tt = np.linspace(0, 2*np.pi, 1000)



if xpos1 == smp.cos(float(xpos)*t):
    lx = np.cos(float(xpos)*tt)
elif xpos1 == smp.sin(float(xpos)*t):
    lx = np.sin(float(xpos)*tt)
else:
    lx = float(xpos)*tt

if ypos1 == smp.cos(float(ypos)*t):
    ly = np.cos(float(ypos)*tt)
elif ypos1 == smp.sin(float(ypos)*t):
    ly = np.sin(float(ypos)*tt)
else:
    ly = float(ypos)*tt

if zpos1 == smp.cos(float(zpos)*t):
    lz = np.cos(float(zpos)*tt)
elif zpos1 == smp.sin(float(zpos)*t):
    lz = np.sin(float(zpos)*tt)
else:
    lz = float(zpos)*tt



#lx, ly, lz = np.cos(4*tt), np.sin(4*tt), tt
data = go.Cone(x=xv.ravel(), y=yv.ravel(), z=zv.ravel(),
               u=Ex.ravel(), v=Ey.ravel(), w=Ez.ravel(),
               colorscale='Inferno', colorbar=dict(title='$x^2$'),
               sizemode="scaled", sizeref=0.5)

layout = go.Layout(title=r'Plot Title',
                     scene=dict(xaxis_title=r'x',
                                yaxis_title=r'y',
                                zaxis_title=r'z',
                                aspectratio=dict(x=1, y=1, z=1),
                                camera_eye=dict(x=1.2, y=1.2, z=1.2)))

fig = go.Figure(data = data, layout=layout)
fig.add_scatter3d(x=lx, y=ly, z=lz, mode='lines',
                  line = dict(color='green', width=10))


HTML(fig.to_html(default_width=1000, default_height=600))

fig.show()