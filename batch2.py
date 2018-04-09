import os
import subprocess

make_cmd = "make"
dists = [50,]
#angles = [60, 135]
#angles = [10, 30, 45]
#angles = [90, 120]
#angles = [150, 180]
#angles = [120, 150, 180]
#angles = [10, 90]
angles = [10,]
#angles = [10, 45, 90]
#angles = [10, 45, 90, 135, 180]
#angles = [10, 30, 45, 60, 90, 120, 135, 150, 180]
#angles = [5, 15, 20, 25, 35, 40, 50, 55, 65, 70, 75, 80, 85, 95, 100, 105, 110, 115]
#frecs = [0.9, 1.2, 1.4, 1.6, 1.8, 0.001]
frecs = [0.1,]
#frecs = [0.1, 0.8, 2.0]
#frecs = [0.01,]
#frecs = [0.05, 0.2, 0.5, 0.6, 0.7]
#frecs = [0.9, 3.0]
#frecs = [0.1, 0.8, 2.0, 0.4, 1.0]
#frecs = [0.1, 0.4, 0.8, 1.0 , 2.0]
#frecs = [0.1, 0.2, 0.4, 0.8, 1.0 , 1.5, 2.0]

flag_win = False
try:
    os.uname()
except AttributeError:
    flag_win = True

if flag_win:
    make_cmd = "mingw32-make"

for dist in dists:
    for frec in frecs:
        for angle in angles:
            args = [make_cmd, "--file=makefile2", "ANGLE="+str(angle), "DIST="+str(dist), "FREC="+str(frec)]
            subprocess.call(args)

    for frec in frecs:
        for angle in angles:
            subprocess.call("./{0}dist_{1}deg_{2}rec.exe".format(dist, angle, frec))
            print("{0}dist, {1}deg. {2}rec. ended.".format(dist, angle, frec))

    # subprocess.call([make_cmd, "clean"])
