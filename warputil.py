import warp
from numpy import array,zeros
import cv2

_DEBUG_ = False

def MSG(s):
    if _DEBUG_:
        print "warputil:"+s

def WarpOpen(wn, fn, ms = 15, Wres = 300, Hres = 300, M = 20, k_spr = 10, dt = 0.0125):
    """
    WarpOpen(wn,fn,Wres,Hres,M,k_spr,dt).
    Open the warp session.
    wn is name of window like "warp" to display in
    fn is filename of image to warp
    Wres and Hres are the width and height resolution of warp
    M x M is number of grid rectangles that warp does
    k_spr is Hook's Spring constant for the springs in M x M grid. Eg, 1,10,100,1000,10000
    or between such that grid is not wildly out of control.
    dt is simulation time increment. Use small value to not wildly control grid.
    W and H determine resolution of warp. Smaller means faster. Bigger means
    smoother.
    """
    MSG("WarpOpen")
    cv2.namedWindow(wn)
    img = cv2.imread(fn,1)
    img = cv2.resize(img, (Wres,Hres))
    (h,w,cols) = img.shape
    L0 = list(img.flatten())
    warp.open(L0, M, w, h)
    warp.set_k_spr(k_spr)
    warp.set_dt(dt)
    return

def WarpSetU(u):
    """
    u is 2D list of values, a scalar field is the curvature or amount or
    warp expansion at a grid point
    """
    MSG("WarpSetU")
    warp.set_u(u)
    return

def WarpImage(w=500,h=500):
    """
    Get Warp image and resize to w x h.
    """
    MSG("WarpImage")
    L = warp.image()
    ww = warp.get_width()
    hh = warp.get_height()
    img = array(L, dtype='uint8')
    img = img.reshape((hh,ww,3))
    img = cv2.resize(img, (w,h))
    return img

def WarpSave(fn, w=500,h=500, msg=True):
    """
    Save Warp image resizing to w x h
    """
    MSG("WarpSave")
    img = WarpImage(w,h)
    cv2.imwrite(fn,img)
    if msg:
        print "Image written: ", fn
    return

def WarpStep(iters=5):
    """
    Perform a warp step with iters iterations. More iters is
    slower to display but show bigger changes.
    """
    MSG("WarpStep")
    for j in range(iters):
        warp.step()
    return

def WarpShow(wn="warp",ms=15, w=500, h=500):
    """
    Display the warped grid image. ms is miliseconds to wait in
    displaying image. w and h are width and height of displayed image.
    """
    MSG("WarpShow")
    img = WarpImage(w,h)
    cv2.imshow(wn,img)
    cv2.waitKey(ms)
    return

def WarpCloseWindows(ms = -1, msg = True):
    """
    Close all cv2 windows
    """
    MSG("WarpCloseWindows")
    if msg:
        print "Press space to quit and close windows"
        cv2.waitKey(ms)
        cv2.destroyAllWindows()
    return

def WarpClose(ms = -1, msg = True):
    """
    Close the warp session. ms = -1 waits else it waits ms miliseconds
    before closing existing windows.
    """
    MSG("WarpClose")
    warp.close()
    WarpCloseWindows(ms, msg)
    return
