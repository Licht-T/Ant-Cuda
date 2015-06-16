BEGIN{
	FS = " "
	num = 0
}
{
	time[NR]    = $1
	a11[NR] = $9
	a12[NR] = $10
	a21[NR] = $11
	a22[NR] = $12
	p1[NR] = $18
	p2[NR] = $19
	num++
}
END{
	for (i=2; i<num; i++){
		a11o=a11[i]-a11[i-1]
		a12o=a12[i]-a12[i-1]
		a21o=a21[i]-a21[i-1]
		a22o=a22[i]-a22[i-1]
		p1o=p1[i]
		p2o=p2[i]
	       	printf("%d %d %d %d %d %f %f\n",time[i],a11o,a12o,a21o,a22o,p1o,p2o)
	}
}
