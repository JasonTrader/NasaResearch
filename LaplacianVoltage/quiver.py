from plotly.tools import FigureFactory as FF
from plotly.graph_objs import *
from plotly.offline import plot
import numpy as np



x_range = raw_input().split(' ')
x_range = map(float,x_range)

y_range = raw_input().split(' ')
y_range = map(float,y_range)

xgrid = [x_range]
for i in range(len(y_range) - 1):
    xgrid.append(x_range)

ygrid = [y_range]
for i in range(len(x_range) - 1):
    ygrid.append(y_range)

#transposes ygrid
ygrid = zip(*ygrid)




gridData = [[]]
for i in range(len(y_range)):
    line = raw_input().split(' ')
    line = map(float,line)
    if i == 0:
        gridData = [line]
    else:
        gridData.append(line)

dz,dr = raw_input().split()
dz = float(dz)
dr = float(dr)

#axis = 1 gives x direction
#axis = 0 gives y direction
u = np.gradient(gridData, -500*dz, axis=1)
v = np.gradient(gridData, -500*dr, axis=0)


#print gridData


fig = FF.create_quiver(xgrid, ygrid, u, v)
plot(fig, filename="EfieldQuiver.html")
