import plotly.plotly as py
from plotly.graph_objs import *
from plotly.offline import plot

x_range = raw_input().split(' ')
for item in x_range:
    float(item)
y_range = raw_input().split(' ')
for item in y_range:
    float(item)
gridData = [[]]
for i in range(len(y_range)):
    line = raw_input().split(' ')
    for item in line:
        float(item)
    if i == 0:
        gridData = [line]
    else:
        gridData.append(line)

data = [
    Contour(
        z=gridData,
        x=x_range,
        y=y_range,
        contours=Contours(
            coloring="heatmap",

        )
    )]

fig = Figure(data=data)
plot(fig, filename="title"+".html")
