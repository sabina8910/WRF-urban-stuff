from plotly.graph_objs import *
from mpl_toolkits.basemap import Basemap
import xarray as xr
import numpy as np
from plotly.subplots import make_subplots

# These are the functions to add the country outlines to the plot:
# Function to create the scatter plot for the country outlines
def make_scatter(x,y):
    return Scatter(
        x=x,
        y=y,
        mode='lines',
        line=scatter.Line(color="black"),
        name=' '  # no name on hover
    )

# Functions converting coastline/country polygons to lon/lat traces
def polygons_to_traces(poly_paths, N_poly):
    ''' 
    pos arg 1. (poly_paths): paths to polygons
    pos arg 2. (N_poly): number of polygon to convert
    '''
    traces = []  # init. plotting list 

    for i_poly in range(N_poly):
        poly_path = poly_paths[i_poly]
        
        # get the Basemap coordinates of each segment
        coords_cc = np.array(
            [(vertex[0],vertex[1]) 
             for (vertex,code) in poly_path.iter_segments(simplify=False)]
        )
        
        # convert coordinates to lon/lat by 'inverting' the Basemap projection
        lon_cc, lat_cc = m(coords_cc[:,0],coords_cc[:,1], inverse=True)
        
        # add plot.ly plotting options
        traces.append(make_scatter(lon_cc,lat_cc))
     
    return traces

# Function generating coastline lon/lat traces
def get_coastline_traces():
    poly_paths = m.drawcoastlines().get_paths() # coastline polygon paths
    N_poly = 91  # use only the 91st biggest coastlines (i.e. no rivers)
    return polygons_to_traces(poly_paths, N_poly)

# Function generating country lon/lat traces
def get_country_traces():
    poly_paths = m.drawcountries().get_paths() # country polygon paths
    N_poly = len(poly_paths)  # use all countries
    return polygons_to_traces(poly_paths, N_poly)




m = Basemap()
ds1=xr.open_dataset("/terra/users/csag/sabina/heat_wrf/WRF/heat_exp1/WRF/regrid_temp_wrfout_d01_2016-01-01-00-06.nc")
ds2=xr.open_dataset("/terra/users/csag/sabina/heat_wrf/WRF/heat_exp1/WRF/regrid_temp_wrfout_d02_2016-01-01-00-06.nc")
ds3=xr.open_dataset("/terra/users/csag/sabina/heat_wrf/WRF/heat_exp1/WRF/regrid_temp_wrfout_d03_2016-01-01-00-06.nc")

traces_cc = get_coastline_traces()+get_country_traces()

title = "Surface Temperature"
anno_text = ""
axis_style = dict(zeroline=False,showline=False,showgrid=False,ticks='',showticklabels=True,)
layout1 = Layout(title=title,
        showlegend=False,
        hovermode="closest",        # highlight closest point on hover
        xaxis=XAxis(
            axis_style, range=[ds1.lon[0],ds1.lon[-1]]  # restrict y-axis to range of lon
        ),
        yaxis=YAxis(
            axis_style, range=[ds1.lat[0],ds1.lat[-1]]
        ),
	coloraxis={'colorbar': None},
        autosize=False,
        width=1000,
        height=1000,
    )
layout2 = Layout(title=title,
        showlegend=False,
        hovermode="closest",        # highlight closest point on hover
        xaxis2=XAxis(
            axis_style, range=[ds1.lon[0],ds1.lon[-1]]  # restrict y-axis to range of lon
        ),
        yaxis2=YAxis(
            axis_style, range=[ds1.lat[0],ds1.lat[-1]]
        ), 
	coloraxis={'colorbar': None},
        autosize=False,
        width=1000,
        height=1000,
    )
layout3 = Layout(title=title,
        showlegend=False,
        hovermode="closest",        # highlight closest point on hover
        xaxis3=XAxis(
            axis_style, range=[ds1.lon[0],ds1.lon[-1]]  # restrict y-axis to range of lon
        ),
        yaxis3=YAxis(
            axis_style, range=[ds1.lat[0],ds1.lat[-1]]
        ),
	coloraxis={'colorbar': None},
        autosize=False,
        width=1000,
        height=1000,
    )
layout4 = Layout(title=title,
        showlegend=False,
        hovermode="closest",        # highlight closest point on hover
        xaxis4=XAxis(
            axis_style, range=[ds3.lon[0],ds3.lon[-1]]  # restrict y-axis to range of lon
        ),
        yaxis4=YAxis(
            axis_style, range=[ds3.lat[0],ds3.lat[-1]]
        ),
        coloraxis={'colorbar': None},
        autosize=False,
        width=1000,
        height=1000,
    )
layout5 = Layout(title=title,
        showlegend=False,
        hovermode="closest",        # highlight closest point on hover
        xaxis5=XAxis(
            axis_style, range=[ds3.lon[0],ds3.lon[-1]]  # restrict y-axis to range of lon
        ),
        yaxis5=YAxis(
            axis_style, range=[ds3.lat[0],ds3.lat[-1]]
        ),
        coloraxis={'colorbar': None},
        autosize=False,
        width=1000,
        height=1000,
    )
layout6 = Layout(title=title,
        showlegend=False,
        hovermode="closest",        # highlight closest point on hover
        xaxis6=XAxis(
            axis_style, range=[ds3.lon[0],ds3.lon[-1]]  # restrict y-axis to range of lon
        ),
        yaxis6=YAxis(
            axis_style, range=[ds3.lat[0],ds3.lat[-1]]
        ),
        coloraxis={'colorbar': None},
        autosize=False,
        width=1000,
        height=1000,
    )

print(ds1.T_sfc)
print(ds1.lon)
print(ds1.lat)
trace1=Contour(z=ds1.T_sfc[0,:,:].values,x=ds1.lon,y=ds1.lat,zmin=280,zmax=320,ncontours=40,showscale=False)
trace2=Contour(z=ds2.T_sfc[0,:,:].values,x=ds2.lon,y=ds2.lat,zmin=280,zmax=320,ncontours=40,showscale=False)
trace3=Contour(z=ds3.T_sfc[0,:,:].values,x=ds3.lon,y=ds3.lat,zmin=280,zmax=320,ncontours=40,showscale=True,colorbar=dict(len=0.5,x=1.,y=0.769,orientation='v'))
trace4=Contour(z=ds3.T_sfc[0,:,:].values,x=ds3.lon,y=ds3.lat,zmin=290,zmax=300,ncontours=40,showscale=True,colorbar=dict(len=0.5,x=1.,y=0.231,orientation='v'))
trace5=Contour(z=ds1.T_sfc[0,:,:].values,x=ds1.lon,y=ds1.lat,zmin=290,zmax=300,ncontours=40,showscale=False)
trace6=Contour(z=ds2.T_sfc[0,:,:].values,x=ds1.lon,y=ds1.lat,zmin=290,zmax=300,ncontours=40,showscale=False)
fig = make_subplots(rows=2, cols=3,horizontal_spacing = 0.05,vertical_spacing=0.05)
fig.add_trace(trace1,row=1,col=1) #0
fig.add_trace(trace1,row=1,col=2) #1
fig.add_trace(trace1,row=1,col=3) #2
fig.add_trace(trace2,row=1,col=2) #3
fig.add_trace(trace2,row=1,col=3) #4
fig.add_trace(trace3,row=1,col=3) #5
fig.add_trace(trace5,row=2,col=1) #6
fig.add_trace(trace6,row=2,col=2) #7
fig.add_trace(trace4,row=2,col=3) #8
fig.add_traces(traces_cc,exclude_empty_subplots=False,rows=[1]*541,cols=[1]*541)
fig.add_traces(traces_cc,exclude_empty_subplots=False,rows=[1]*541,cols=[2]*541)
fig.add_traces(traces_cc,exclude_empty_subplots=False,rows=[1]*541,cols=[3]*541)
lon = np.loadtxt('lon.dat')
lat = np.loadtxt('lat.dat')
scatter1=Scatter(x=lon,y=lat,mode='lines',line=scatter.Line(color="black"))
fig.add_trace(scatter1,row=1,col=1)
fig.add_trace(scatter1,row=1,col=2)
fig.add_trace(scatter1,row=1,col=3)
fig.add_trace(scatter1,row=2,col=1)
fig.add_trace(scatter1,row=2,col=2)
fig.add_trace(scatter1,row=2,col=3)
fig.update(layout=layout1)
fig.update(layout=layout2)
fig.update(layout=layout3)
fig.update(layout=layout4)
fig.update(layout=layout5)
fig.update(layout=layout6)
frames = [dict(
               name = k,data=[Contour(z=ds1.T_sfc[k,:,:].values,x=ds1.lon,y=ds1.lat,zmin=280,zmax=320,ncontours=40),Contour(z=ds1.T_sfc[k,:,:].values,x=ds1.lon,y=ds1.lat,zmin=280,zmax=320,ncontours=40),Contour(z=ds1.T_sfc[k,:,:].values,x=ds1.lon,y=ds1.lat,zmin=280,zmax=320,ncontours=40),Contour(z=ds2.T_sfc[k,:,:].values,x=ds2.lon,y=ds2.lat,zmin=280,zmax=320,ncontours=40),Contour(z=ds2.T_sfc[k,:,:].values,x=ds2.lon,y=ds2.lat,zmin=280,zmax=320,ncontours=40),Contour(z=ds3.T_sfc[k,:,:].values,x=ds3.lon,y=ds3.lat,zmin=280,zmax=320,ncontours=40),Contour(z=ds1.T_sfc[k,:,:].values,x=ds1.lon,y=ds1.lat,zmin=290,zmax=300,ncontours=40),Contour(z=ds2.T_sfc[k,:,:].values,x=ds2.lon,y=ds2.lat,zmin=290,zmax=300,ncontours=40),Contour(z=ds3.T_sfc[k,:,:].values,x=ds3.lon,y=ds3.lat,zmin=290,zmax=302,ncontours=40)],traces=[0, 1, 2, 3, 4, 5,6,7,8] # the elements of the list [0,1,2] give info on the traces in fig.data so it knows which traces already defined should be updated  
) for k in range(0,6)]

# Create buttons - a play button to play the animation and a slider to move to the timestep required
updatemenus = [dict(type='buttons',
                    buttons=[dict(label='Play',
                                  method='animate',
                                  args=[[f'{k}' for k in range(0,6)], 
                                         dict(frame=dict(duration=500, redraw=True), 
                                              transition=dict(duration=0),
                                              easing='linear',
                                              fromcurrent=True,
                                              mode='immediate'
                                                                 )])],
                    direction= 'left', 
                    pad=dict(r= 10, t=85), 
                    showactive =True, x= 0.1, y= 0, xanchor= 'right', yanchor= 'top')
            ]

sliders = [{'yanchor': 'top',
            'xanchor': 'left', 
            'currentvalue': {'font': {'size': 16}, 'prefix': 'Frame: ', 'visible': True, 'xanchor': 'right'},
            'transition': {'duration': 500.0, 'easing': 'linear'},
            'pad': {'b': 10, 't': 50}, 
            'len': 0.9, 'x': 0.1, 'y': 0, 
            'steps': [{'args': [[k], {'frame': {'duration': 500.0, 'easing': 'linear', 'redraw': True},
                                      'transition': {'duration': 0, 'easing': 'linear'}}], 
                       'label': k, 'method': 'animate'} for k in range(0,6)       
                    ]}]
                
# Update the frames to the figure
fig.update(frames=frames),
# update the layout so that the buttons are added
fig.update_layout(updatemenus=updatemenus,
                  sliders=sliders);
#fig.update(layout=layout1)
#fig.update(layout=layout2)
#fig.update(layout=layout3)
#fig.update(layout=layout4)
#fig.update(layout=layout5)
#fig.update(layout=layout6)
fig.write_html("plot.html", auto_play=False)
