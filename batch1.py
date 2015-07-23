import os
import subprocess

make_cmd = "make"
dists = [30, 40]
angles = [10, 45]

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
		subprocess.call("./"+str(angle)+".deg")
		print(str(dist)+"dist, "+str(angle)+"deg. ended.")
	subprocess.call([make_cmd, "clean"])

