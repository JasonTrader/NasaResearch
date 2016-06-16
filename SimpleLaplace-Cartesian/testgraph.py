# Get this figure: fig = py.get_figure("https://plot.ly/~jasonmtrader/0/")
# Get this figure's data: data = py.get_figure("https://plot.ly/~jasonmtrader/0/").get_data()
# Add data to this figure: py.plot(Data([Scatter(x=[1, 2], y=[2, 3])]), filename ="basic-heatmap", fileopt="extend"))
# Get z data of first trace: z1 = py.get_figure("https://plot.ly/~jasonmtrader/0/").get_data()[0]["z"]

# Get figure documentation: https://plot.ly/python/get-requests/
# Add data documentation: https://plot.ly/python/file-options/

# If you're using unicode in your file, you may need to specify the encoding.
# You can reproduce this figure in Python with the following code!

# Learn about API authentication here: https://plot.ly/python/getting-started
# Find your api_key here: https://plot.ly/settings/api

import plotly.plotly as py
from plotly.graph_objs import *
from plotly.offline import plot

title = raw_input()
xaxislabel = raw_input()
yaxislabel = raw_input()
deltax = raw_input()
lengthx = raw_input()
deltay = raw_input()
lengthy = raw_input()
num_xseg = raw_input()
num_yseg = raw_input()
gridData = [[]]
for i in range(int(num_xseg) + 1):
    dataline = raw_input().split(',')
    for j in range(int(num_yseg)):
        dataline[i] = float(dataline[i])
    gridData.append(dataline)



data = Data([
    Heatmap(
        z=gridData,
        autocolorscale=False,
        colorbar=ColorBar(
            x=1.0101337965353738,
            y=0.48881789137380194,
            bgcolor='rgba(0,0,0,0)',
            bordercolor='#444',
            borderwidth=0,
            exponentformat='B',
            len=1,
            lenmode='fraction',
            nticks=0,
            outlinecolor='#444',
            outlinewidth=1,
            showexponent='all',
            showticklabels=True,
            thickness=30,
            thicknessmode='pixels',
            tickangle='auto',
            tickfont=dict(
                family='"Open Sans", verdana, arial, sans-serif',
                size=12
            ),
            tickformat='',
            tickmode='auto',
            tickprefix='',
            ticks='',
            ticksuffix='',
            title='Click to enter colorscale title',
            titlefont=dict(
                color='#444',
                family='"Open Sans", verdana, arial, sans-serif',
                size=12
            ),
            titleside='top',
            xanchor='left',
            xpad=10,
            yanchor='middle',
            ypad=10
        ),
        colorscale=[[0, 'rgb(50, 0, 0)'], [0.35, 'rgb(153, 31, 0)'], [0.5, 'rgb(175, 50, 0)'], [0.6, 'rgb(200, 75, 0)'], [0.7, 'rgb(225, 100, 0)'], [1, 'rgb(249, 150, 0)']],
        connectgaps=False,
        dx=deltax,
        dy=deltay,
        hoverinfo='x+y+z+text',
        name='trace 0',
        opacity=1,
        reversescale=False,
        showscale=True,
        transpose=True,
        uid='ce4782',
        visible=True,
        x0=0,
        xaxis='x',
        y0=float(deltay)/2,
        yaxis='y',
        zauto=True,
        zsmooth=False,
        zsrc='jasonmtrader:1:'
    )
])
layout = Layout(
    autosize=True,
    dragmode='zoom',
    font=Font(
        color='#444',
        family='"Open Sans", verdana, arial, sans-serif',
        size=12
    ),
    height=806.312,
    hidesources=False,
    hovermode='closest',
    margin=Margin(
        r=80,
        t=100,
        autoexpand=True,
        b=80,
        l=80,
        pad=0
    ),
    paper_bgcolor='#fff',
    plot_bgcolor='#fff',
    separators='.,',
    showlegend=False,
    smith=False,
    title=title,
    titlefont=dict(
        color='#444',
        family='"Open Sans", verdana, arial, sans-serif',
        size=17
    ),
    width=1435,
    xaxis=XAxis(
        anchor='y',
        autorange=False,
        color='#444',
        domain=[0, float(deltax)],
        dtick=0.5,
        exponentformat='B',
        fixedrange=False,
        gridcolor='rgb(238, 238, 238)',
        gridwidth=1,
        hoverformat='',
        mirror=False,
        nticks=0,
        range=[0, float(lengthx) + 2*float(deltax)],
        showexponent='all',
        showgrid=True,
        showline=False,
        showticklabels=True,
        side='bottom',
        tick0=0,
        tickangle='auto',
        tickcolor='#444',
        tickfont=dict(
            color='#444',
            family='"Open Sans", verdana, arial, sans-serif',
            size=12
        ),
        tickformat='',
        ticklen=5,
        tickmode='auto',
        tickprefix='',
        ticks='outside',
        ticksuffix='',
        tickwidth=1,
        title=xaxislabel,
        titlefont=dict(
            color='#444',
            family='"Open Sans", verdana, arial, sans-serif',
            size=14
        )
    ),
    yaxis=YAxis(
        anchor='x',
        autorange=False,
        color='#444',
        domain=[0, float(deltay)],
        dtick=0.5,
        exponentformat='B',
        fixedrange=False,
        gridcolor='rgb(238, 238, 238)',
        gridwidth=1,
        hoverformat='',
        mirror=False,
        nticks=0,
        range=[0, float(lengthy) + 2*float(deltay)],
        showexponent='all',
        showgrid=True,
        showline=False,
        showticklabels=True,
        side='left',
        tick0=0,
        tickangle='auto',
        tickcolor='#444',
        tickfont=dict(
            color='#444',
            family='"Open Sans", verdana, arial, sans-serif',
            size=12
        ),
        tickformat='',
        ticklen=5,
        tickmode='auto',
        tickprefix='',
        ticks='outside',
        ticksuffix='',
        tickwidth=1,
        title=yaxislabel,
        titlefont=dict(
            color='#444',
            family='"Open Sans", verdana, arial, sans-serif',
            size=14
        )
    )
)
fig = Figure(data=data, layout=layout)
plot(fig, filename=title+".html")
