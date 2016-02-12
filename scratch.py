import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
import pres
#import simul_winds
import propigation	
#import gracon_vel2d
#import gracon_vel2d_3d
import read_rays
import radar_math
import radar_to_cart
import mathematics
import netcdf_utis
#import grad_conj_solver_plus
import grad_conj_solver_plus_plus
import grad_conj_solver_3d
import read_sounding
import met
#from matplotlib import figure
from pylab import *
from time import time as systime
from numpy import linspace, array, arctan, pi, random
#import pickle_zip
import dealias
import parse_ini
import read_radar
import cappi_v2
from matplotlib.toolkits.basemap import Basemap #map plotting things for matplotlib


ber_loc=[-12.4, 130.85] #location of Berrimah radar
gp_loc=[-12.2492,  131.0444]#location of CPOL at Gunn Point
lats=linspace(-11.5, -13.0, 100)
lons=linspace(130.0, 132.0, 100)
angs=array(propigation.make_lobe_grid(gp_loc, ber_loc, lats,lons))
f=figure()
my_contour=contour(lons, lats, angs, levels=[30])
verts=my_contour.collections[0].get_paths()[0].vertices
lon_150=[]
lat_150=[]
for xy in verts:
	lon_150.append(xy[0])
	lat_150.append(xy[1])

close(f)
dd_box=[array(lon_150).min(), array(lon_150).max(), array(lat_150).min(), array(lat_150).max()]
fat=0.05
box=[dd_box[0]-fat, dd_box[1]+fat, dd_box[2]-fat, dd_box[3]+fat]
f=figure()
mapobj= Basemap(projection='merc',lat_0=(box[2]+box[3])/2.0, lon_0=(box[0]+box[1])/2.0,llcrnrlat=box[2], llcrnrlon=box[0], urcrnrlat=box[3] , urcrnrlon=box[1], resolution='l',area_thresh=1., lat_ts=(box[2]+box[3])/2.0)
longr, latgr=meshgrid(lons,lats)
xx, yy = mapobj(longr, latgr)
mapobj.drawmapboundary()
mapobj.readshapefile('/flurry/home/scollis/shapes/cstntcd_r','coast',drawbounds=True,linewidth=0.5,color='k',antialiased=1,ax=None)
mapobj.contour(xx,yy,angs, levels=[150,30], colors=['r'])
mapobj.drawmeridians(array([130.2, 130.4, 130.6, 130.8,131.0,131.2, 131.4]), labels=[1,0,0,1])
mapobj.drawparallels(array([--12.8, -12.6, -12.4, -12.2, -12.0, -11.8, -11.6, -11.4]), labels=[1,0,0,1])
dd_box_xx, dd_box_yy=mapobj(array([dd_box[0], dd_box[1]]),array([dd_box[2], dd_box[3]]))
xy=dd_box_xx[0], dd_box_yy[0]
width, height= dd_box_xx[1]-dd_box_xx[0], dd_box_yy[1]-dd_box_yy[0]
my_patch=matplotlib.patches.Rectangle(xy, width, height, edgecolor='blue', facecolor='white')
ax=gca()
ax.add_patch(my_patch)
radar_xx, radar_yy=mapobj(array([ber_loc[1], gp_loc[1]]), array([ber_loc[0], gp_loc[0]]))
mapobj.plot(radar_xx, radar_yy, 'bo')
ax.text(radar_xx[0]+1000.0, radar_yy[0]-3000.0, 'Berrimah')
ax.text(radar_xx[1]+1000.0, radar_yy[1]-3000.0, 'Gunn Point')
savefig('/flurry/home/scollis/results/area_vis.png')
close(f)



|
 |        =================   ==============================================
 |        Property            Description
 |        =================   ==============================================
 |        alpha               float
 |        animated            [True | False]
 |        antialiased or aa   [True | False]
 |        clip_box            a matplotlib.transform.Bbox instance
 |        clip_on             [True | False]
 |        edgecolor or ec     any matplotlib color
 |        facecolor or fc     any matplotlib color
 |        figure              a matplotlib.figure.Figure instance
 |        fill                [True | False]
 |        hatch               unknown
 |        label               any string
 |        linewidth or lw     float
 |        lod                 [True | False]
 |        transform           a matplotlib.transform transformation instance
 |        visible             [True | False]
 |        zorder              any number
 |        =================   ==============================================
 |















