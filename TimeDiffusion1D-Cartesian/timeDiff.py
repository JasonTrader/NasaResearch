import plotly.plotly as py
from plotly.graph_objs import *
from plotly.offline import plot

x_range = raw_input().split(' ')
for item in x_range:
    float(item)

line = raw_input()
vals = line.split(' ')
for item in vals:
    float(item)
tempplot = Scatter(x=x_range,y=vals,mode='lines',line=dict(shape='linear'))
plots = [tempplot]

line = raw_input()
counter = 0
while line != 'STOP':
    vals = line.split(' ')
    for item in vals:
        float(item)
    tempplot = Scatter(x=x_range,y=vals,mode='lines',line=dict(shape='spline'))
    plots.append(tempplot)

    line = raw_input()
    counter += 1

data = plots


fig = Figure(data=data)
plot(fig, filename="results"+".html")
