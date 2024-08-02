import warputil
import warp
from infix import *
from numpy import array,zeros
import cv2

_DEBUG_ = False
def MSG(s):
    if _DEBUG_:
        print(s)
    return
    
counter = 0

def App(f,x):
    return f(x)
def Action(u, fn):
    global xx, vv
    global counter
    MSG("Action")
    warputil.WarpOpen("warp",fn, Wres=400, Hres=400)
    warputil.WarpSetU(u)
    if xx is not None and len(u) == len(xx):
        warp.set_x(xx)
        warp.set_v(vv)
    warputil.WarpStep()
    xx = warp.get_x()
    vv = warp.get_v()
    warputil.WarpShow()
    fn_save = "./Z-imgs/tmp.jpg"
    if counter%10 == 9:
        fn_save = "./Z-imgs/semigroup-action-%03d.jpg" % counter
    warputil.WarpSave(fn_save, msg=False)    
    warputil.WarpClose(ms=1,msg=False)
    counter = counter + 1
    return fn_save
def Compose(f, g):
    def h(x):
        return f(g(x))
    return h
def T(u):
    def f(x):
        return Action(u,x)
    return f
def Pow(u, n):
    o = Infix(Compose)
    if n <= 1:
        return u
    else:
        return u |o| Pow(u, n-1)
def Schwartzschild(mass, M = 10):
    MSG("Schwartzschild")
    warputil.WarpOpen("warp",fn, M = M)
    warp.set_u_schwartzchild(mass)
    u5 = warp.get_u()
    a5 = array(u5)
    a5 = a5.reshape((M,M))
    b5 = zeros((M,M),dtype='uint8')
    for j in range(a5.shape[1]):
        for i in range(a5.shape[0]):
            b5[i,j] = a5[i,j]
    warputil.WarpClose(ms=1,msg=False)
    return u5,b5

o = Infix(Compose)
a = Infix(App)

u1 = [
    2,3,4,4,
    1,2,3,4,
    1,2,3,3,
    0,1,2,3
    ]
t1 = T(u1)

u2 = [
    1,3,3,1,
    3,5,5,3,
    3,5,5,3,
    0,3,3,1
    ]
t2 = T(u2)

u3 = [
    1,1,1,1,1,
    1,2,3,2,1,
    1,2,4,2,1,
    1,2,3,2,1,
    1,1,1,1,1
    ]
t3 = T(u3)

u4 = [
    1,1,1,1,
    1,1,1,1,
    1,1,1,1,
    1,1,1,1
    ]
t4 = T(u4)

fn = "checkboard.jpg"


mass_moon = 7.34767309e22
mass_earth = 5.972e24
mass_sun = 1.989e30
mass_sagitarius_a = 3.7e6*mass_sun
masses = [mass_moon,mass_earth,mass_sun,mass_sagitarius_a]
for mass in masses:
    xx = None
    yy = None
    print("mass = ", mass)
    u5,b5 = Schwartzschild(mass, M=8)
    print(b5)
    u5 = list(b5.flatten())
    t5 = T(u5)
    (Pow(t5,10)) |a| "graph_grid.jpg"

(Pow(t1,10) |o| t3 |o| Pow(t2,10) |o| t4 |o| t1 |o| t2) |a| fn
(Pow(t2,10) |o| t4 |o| Pow(t1,10) |o| Pow(t5,4) ) |a| fn

warputil.WarpCloseWindows(15,True)

