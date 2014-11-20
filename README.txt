README.txt

This is a library to warp an image.

If warp.dll (aka warp.pyd) and warp.lib aren't created, then create
one with warp.sln in the ./src folder. Then read helper.bat and
edit Python location and check if Python location has a package
already called warp.pyd or warp.lib or warp. If no warp package
exists then you can install with helper.bat. It will rename
warp.dll to warp.pyd and copy warp.dll and warp.lib to 
(Python27)\DLLs and (Python27)\libs.

warp package is now installed and be access via

import warp

Next look at warputil.py. This is front-end to warp
accessed via

import warputil.py

You might want to add later to python or just put warputil
in current directory.

Now you can use the Warp Image Python Library (pywarp) by
running use_warp.py or warp_semigroup_action.py.

Look at warp_semigroup_action.py and use_warp.py as sample
code. Run with different scalar fields square M x M flat
lists of doubles which indicate curvature of the image.
And an input image indicated by a filename, such as a jpg
or other recognized by cv2. Look at Creating a Warp Image...
.txt file for an older version of warpmodule.cpp's help info
but still useful.

Required packages are Numpy and OpenCv2 python (copy
cv2.pyd to Python package if not already their).

The package uses a scalar field U which is M x M in length
list of doubles and a image file fn. It creates a M x M grid
of particles with springs connecting them with spring constant
k_spr. It animates in step the position x and velocity v for the
grid and forces F and warping the image to the new grid positions.

Operators were created such that converts a scalar field to a transformation
T(U) which applies to an image file T(U)(fn). Various transformations
may be composed together---they form a semigroup. (T1 |o| T2 |o| T3)(fn)
where T1 = T(u1), T2 = T(u2), and T3 = T(u3). The function application
is also written as multiplication (T1 |o| T2 |o| T3) |a| fn.

An analogy, the Ti are like operations F, R, L, B on a Rubik's cube
where fn is like the Rubik's cube. You can compose operations together
and apply that new composition function to any Rubik's cube (position) to get
a new Rubik's cube (position). Here, instead of operations, you apply
warp transformations and instead of Rubik's cube you have image files.

