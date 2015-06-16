import os
import subprocess
import csv
import argparse
import glob

parser = argparse.ArgumentParser(description='Accumulate slopes.')
parser.add_argument('path', metavar='Path', type=str, help='Path to files.')
parser.add_argument('sample', metavar='Sample', type=int, help='Sample size.')
parser.add_argument('angle', metavar='Angle', type=int, help='Angle.')
parser.add_argument('foodsnum', metavar='FoodsNum', type=int, help='Foods num.')

args = parser.parse_args()
path = args.path + '/'
angle = args.angle
sample = args.sample
foodsnum = args.foodsnum

def slope(xhps,xhms,hp,hm):
	l = []
	for i in range(len(xhps)):
		l.append( (float(xhps[i])-float(xhms[i]))/(float(hp)-float(hm)) )
	return l


for pw in range(0,8):
	for n in range(0,550,50):
		slopemeanmean= 0.0
		slopeVarmean = [0.0]*foodsnum
		files = glob.glob(path+'food_' + str(angle)+"deg_10e-"+str(pw)+'_'+str(n)+'normal'+'_sampleNo'+'*'+'_of_'+str(sample)+'.dat')
		for f in files:
			sumlist = [0.0] * foodsnum
			timeList = []
			foodsDic = {}
			slopeDic = {}
			csvReader = csv.reader(open(f, 'rb'), delimiter=' ', quotechar='|')
			rowN = 0
			for row in csvReader:
				tmplist = []
				for food in range(1,foodsnum+1):
					tmplist.append(int(row[food]));
				foodsDic[int(row[0])] = tmplist;
				timeList.append(int(row[0]))
				rowN+=1
			
			for i in range(1, rowN-1):
				slopelist = slope(foodsDic[timeList[i+1]], foodsDic[timeList[i-1]], timeList[i+1],timeList[i-1])
				slopeDic[timeList[i]] = slopelist
				for j in range(0,len(slopelist)):
					sumlist[j] += slopelist[j]
			slopemean = 0.0
			for s in sumlist:
				slopemean += s
			slopemean /= ((rowN-2)*foodsnum)

			slopeVar = [0.0] * foodsnum
			for i in range(1,rowN-1):
				l = slopeDic[timeList[i]]
				for j in range(0,len(l)):
					slopeVar[j] += (l[j]-slopemean)*(l[j]-slopemean)
			for i in range(0,len(slopeVar)):
				slopeVar[i] /= (rowN-2)

			fSlope = f.replace('food_', 'slope_food_')
			csvWriter = csv.writer(open(fSlope, 'wb'), delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			for i in range(1, rowN-1):
				sl = slopeDic[timeList[i]]
				sldif = []
				for s in sl:
					sldif.append(s-slopemean)
				csvWriter.writerow( [timeList[i]] + sl + sldif + [slopemean] + slopeVar)

			for i in range(0,len(slopeVar)):
				slopeVarmean[i] += slopeVar[i]
			slopemeanmean += slopemean

		slopemeanmean/=sample
		for i in range(0,len(slopeVarmean)):
			slopeVarmean[i] /= sample
		fMean = path+'slope_food_' + str(angle)+"deg_10e-"+str(pw)+'_'+str(n)+'normal'+'.dat'
		csvWriter = csv.writer(open(fMean, 'wb'), delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		csvWriter.writerow( [slopemeanmean] + slopeVarmean)



fVarCon = path+'slope_food_' + str(angle)+'.dat'
csvWriter = csv.writer(open(fVarCon, 'wb'), delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
for pw in range(1,8):
	for n in range(0,500,50):
		fMean = path+'slope_food_' + str(angle)+"deg_10e-"+str(pw)+'_'+str(n)+'normal'+'.dat'
		csvReader = csv.reader(open(fMean, 'rb'), delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		for row in csvReader:
			csvWriter.writerow([pw,n,(float(row[1])+float(row[2]))/2])
	csvWriter.writerow([])


