import os
import subprocess

make_cmd = "make"
dists = [30,]
# angles = [10, 45, 90]
angles = [10,]

flag_win = False
try:
    os.uname()
except AttributeError:
    flag_win = True

if flag_win:
    make_cmd = "mingw32-make"

for dist in dists:
    for angle in angles:
        args = [make_cmd, "ANGLE="+str(angle), "DIST="+str(dist)]
        subprocess.call(args)
    for angle in angles:
        subprocess.call("./{0}dist_{1}deg.exe".format(dist, angle))
        print("{0}dist, {1}deg. ended.".format(dist, angle))
    # subprocess.call([make_cmd, "clean"])
