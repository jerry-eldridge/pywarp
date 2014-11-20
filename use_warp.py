import warputil
import warp

##help(warp)
##help(warputil)

u = [
    2,3,4,4,
    1,2,3,4,
    1,2,3,3,
    0,1,2,3
    ]

warputil.WarpOpen("warp","checkboard.jpg")
warputil.WarpSetU(u)

for i in range(10):
    warputil.WarpStep()
    warputil.WarpShow()
##    fn = "./Z-imgs/imgs-%03d.jpg" % i
##    warputil.WarpSave(fn)
    print "t = ",warp.get_t()


warputil.WarpClose()
