import os
import argparse
p = argparse.ArgumentParser()
p.add_argument('file', help='text file')
    
args = p.parse_args()
    
file=args.file
file=str(file)

max=-1000
min=1000

xvalues=[]
yvalues=[]


with open (file,'rb') as f:

    for line in f:
        data=line.split()
        xvalues.append(data[0])
        yvalues.append(data[1])
      

    f.close()


length=len(xvalues)


for i in range(0,length-1):
    x1=float(xvalues[i])
    y1=float(yvalues[i])
    x2=float(xvalues[i+1])
    y2=float(yvalues[i+1])
    slope=(y2-y1)/(x2-x1)
    if (slope>max):
        max=slope
    if (slope<min):
        min=slope



print "max slope=" +str(max)
print "min slope="+str(min)
