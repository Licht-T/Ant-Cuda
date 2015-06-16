import os
import subprocess
import glob

for pw in range(1,8):
	for n in range(0,500,50):
		fnamecore = "food_80deg_10e-"+str(pw)+'_'+str(n)+'normal'
		fname = fnamecore+'_sampleNo*_of_10.dat'
		files = glob.glob(fname)
		out = open("fluc_"+fnamecore+".dat", "wb")
		for f in files:
			subprocess.call(["awk","-f","calcFluc.awk", f],stdout=out)

