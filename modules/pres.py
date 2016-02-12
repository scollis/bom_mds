###################################################################
#pres.py: Presentation utilities for visulisation of radar data   #
#Near Darwin                                                      #
###################################################################
#Module for bom_mds                                               #
#Designed for use in a non-interactive mode, for job scheduling   #
#reasons                                                          #
###################################################################
#Start of BoM Python branch: Scott Collis, CAWCR, April 2008      #
###################################################################
import os
import matplotlib
# chose a non-GUI backend
matplotlib.use( 'Agg' )
from pylab import *
from os import getenv
import mathematics
from mpl_toolkits.basemap import Basemap #map plotting things for matplotlib
from numpy import linspace, unique, sqrt, zeros
import matplotlib.numerix.ma as M
import climt
from scipy import interpolate
import propigation

#Sets up a map centred around darwin useful for showing dual doppler domain of Berrimah+CPOL
def generate_darwin_plot(**kwargs):
	ber_loc=[-12.4, 130.85] #location of Berrimah radar
	gp_loc=[-12.2492,  131.0444]#location of CPOL at Gunn Point
	box=kwargs.get('box',[130.0, 132.0, -13.0, -11.5])#box for the plot
	#Set up the map and projection
	map= Basemap(projection='merc',lat_0=(box[2]+box[3])/2.0, lon_0=(box[0]+box[1])/2.0,llcrnrlat=box[2], llcrnrlon=box[0], urcrnrlat=box[3] , urcrnrlon=box[1], resolution='l',area_thresh=1., lat_ts=(box[2]+box[3])/2.0)
	#map.drawmapboundary()
	map.readshapefile(getenv('HOME')+'/bom_mds/shapes/nt_coast','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	map.drawmeridians(array([129,130,131,132,133]), labels=[1,0,0,1])
	map.drawparallels(array([-14,-13,-12,-11,-10]), labels=[1,0,0,1])
	return map

def reconstruction_plot(mapobj, lats, lons, data_cube, lvl_num, parm, u, v, angs, mask, **kwargs):
	um=M.masked_where(M.array(mask) < 0.5, M.array(u))
	vm=M.masked_where(M.array(mask) < 0.5, M.array(v))
	mag=M.array(sqrt(u**2+v**2))
	magm=M.masked_where(M.array(mask)<0.5, mag)
	fig_name=kwargs.get('fig_name', 'recon_'+parm+'_.png')
	longr, latgr=meshgrid(lons,lats)
	xx, yy = mapobj(longr, latgr)
	mapobj.drawmapboundary()
	mapobj.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	mapobj.drawmeridians(array([130,130.5, 131, 131.5, 132]), labels=[1,0,0,1])
	mapobj.drawparallels(array([-13,-12.5,-12,-11.5,-11]), labels=[1,0,0,1])
	data=data_cube[parm][:,:,lvl_num]
	if parm in ['diff']:
		data=M.masked_where(M.array(mask) < 0.5, M.array(data))
	comp_levs=linspace(-1.,1., 30)
	levs_dict={'VR':linspace(-15,15,31), 'CZ': linspace(-8,64,10), 'TEST':linspace(((abs(data)+data)/2.0).min(), data.max(), 200), 'i_comp':comp_levs,'j_comp':comp_levs,'k_comp':comp_levs,'diff':linspace(-15,15,31)}
	titles_dict={'VR': 'Radial velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST', 'i_comp':'i component','j_comp':'j component','k_comp':'k component', 'diff':'diff'}
	mapobj.contourf(xx,yy,data, levels=levs_dict[parm])
	colorbar()
	#mapobj.quiver(xx,yy,u/sqrt(u**2+v**2),v/sqrt(u**2+v**2), scale=50)
	qq=mapobj.quiver(xx,yy,um,vm, scale=kwargs.get('qscale', 200))
	quiverkey(qq, 0.1, 0.8, 10, '10 m/s', coordinates='figure')
	quiverkey(qq, 0.1, 0.75, 5, '5 m/s', coordinates='figure')
	quiverkey(qq, 0.1, 0.7, 2, '2 m/s', coordinates='figure')
	
	cobject=mapobj.contour(xx,yy,magm, colors=['k'], levels=linspace(0,20,11))
	fon = { 'fontname':'Tahoma', 'fontsize':10 }
	clabel(cobject, fmt="%d", **fon)
	mapobj.contour(xx,yy,angs, levels=[30.0, 150.0],colors=['r']) 

def reconstruction_plot_w(mapobj, lats, lons, data_cube, lvl_num, parm, u, v,w, angs, mask, **kwargs):
	um=M.masked_where(M.array(mask) < 0.5, M.array(u))
	vm=M.masked_where(M.array(mask) < 0.5, M.array(v))
	mag=M.array(sqrt(u**2+v**2))
	magm=M.masked_where(M.array(mask)<0.5, mag)
	fig_name=kwargs.get('fig_name', 'recon_'+parm+'_.png')
	longr, latgr=meshgrid(lons,lats)
	xx, yy = mapobj(longr, latgr)
	mapobj.drawmapboundary()
	mapobj.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	mapobj.drawmeridians(array([130,130.5, 131, 131.5, 132]), labels=[1,0,0,1])
	mapobj.drawparallels(array([-13,-12.5,-12,-11.5,-11]), labels=[1,0,0,1])
	data=data_cube[parm][:,:,lvl_num]
	if parm in ['diff']:
		data=M.masked_where(M.array(mask) < 0.5, M.array(data))
	comp_levs=linspace(-1.,1., 30)
	levs_dict={'VR':linspace(-15,15,31), 'CZ': linspace(-8,64,10), 'TEST':linspace(((abs(data)+data)/2.0).min(), data.max(), 200), 'i_comp':comp_levs,'j_comp':comp_levs,'k_comp':comp_levs,'diff':linspace(-15,15,31)}
	titles_dict={'VR': 'Radial velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST', 'i_comp':'i component','j_comp':'j component','k_comp':'k component', 'diff':'diff'}
	mapobj.contourf(xx,yy,data, levels=levs_dict[parm])
	colorbar()
	#mapobj.quiver(xx,yy,u/sqrt(u**2+v**2),v/sqrt(u**2+v**2), scale=50)
	qq=mapobj.quiver(xx,yy,um,vm, scale=kwargs.get('qscale', 200))
	quiverkey(qq, 0.1, 0.8, 10, '10 m/s', coordinates='figure')
	quiverkey(qq, 0.1, 0.75, 5, '5 m/s', coordinates='figure')
	quiverkey(qq, 0.1, 0.7, 2, '2 m/s', coordinates='figure')
	
	cobject=mapobj.contour(xx,yy,w, colors=['k'], levels=linspace(0,20,11))
	fon = { 'fontname':'Tahoma', 'fontsize':10 }
	clabel(cobject, fmt="%d", **fon)
	mapobj.contour(xx,yy,angs, levels=[30.0, 150.0],colors=['r']) 


def quiver_winds(mapobj, lats, lons, u,v):
	if len(lats.shape)==1:
		longr, latgr=meshgrid(lons, lats)
	else:
		longr=lons
		latgr=lats
	xx,yy=mapobj(longr, latgr)
	mapobj.quiver(xx,yy,u,v)

def quiver_contour_winds(mapobj, lats, lons, u,v, **kwargs):
	if len(lats.shape)==1:
		longr, latgr=meshgrid(lons, lats)
	else:
		longr=lons
		latgr=lats
	xx,yy=mapobj(longr, latgr)
	mag=sqrt(u**2+v**2)
	un=u/mag
	vn=v/mag
	#levs=array(unique([float(int(i)) for i in linspace(mag.min(), mag.max(), 20)]))
	bounds=kwargs.get('bounds',[mag.min(), mag.max()])
	levs=linspace(bounds[0], bounds[1], 20)
	mapobj.contourf(xx,yy,mag,levels=levs)
	colorbar()
	mapobj.quiver(xx,yy,u,v)

def contour_angs(mapobj, lats,lons,angs):
	if len(lats.shape)==1:
		longr, latgr=meshgrid(lons, lats)
	else:
		longr=lons
		latgr=lats
	#levs=[0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0]
	levs=array(unique([float(int(i)) for i in linspace(angs.min(), angs.max(), 20)]))
	xx,yy=mapobj(longr, latgr)
	mapobj.contourf(xx,yy,angs, levels=levs)
	colorbar()

def contour_comp(mapobj, lats,lons,comp):
	if len(lats.shape)==1:
		longr, latgr=meshgrid(lons, lats)
	else:
		longr=lons
		latgr=lats
	#levs=[0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0]
	#levs=array(unique([float(int(i)) for i in linspace(angs.min(), angs.max(), 20)]))
	levs=linspace(-1.0, 1.0, 11)
	xx,yy=mapobj(longr, latgr)
	mapobj.contourf(xx,yy,comp, levels=levs)
	colorbar()

def contour_vr(mapobj, lats, lons, vr):
	if len(lats.shape)==1:
		longr, latgr=meshgrid(lons, lats)
	else:
		longr=lons
		latgr=lats
	#levs=[0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0]
	#levs=array(unique([float(int(i)) for i in linspace(angs.min(), angs.max(), 20)]))
	levs=linspace(-15.0, +15.0, 31)
	xx,yy=mapobj(longr, latgr)
	mapobj.contourf(xx,yy,vr, levels=levs)
	colorbar()


def plot_ppi(sweep_dict, parm, **kwargs):
	radar_loc=kwargs.get('radar_loc', [-12.2492,  131.0444])
	fig_name=kwargs.get('fig_name', 'ppi_'+parm+'_.png')
	fig_path=kwargs.get('fig_path', getenv('HOME')+'/bom_mds/output/')
	f=figure()
	#gp_loc=[-12.2492,  131.0444]
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(radar_loc[0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=radar_loc[1]+360.0*sweep_dict['xar']/(rad_at_radar*2.0*pi)
	lats=radar_loc[0] + 360.0*sweep_dict['yar']/(Re*2.0*pi)
	def_loc_dict={'lat_0':lats.mean(), 'lon_0':lons.mean(),'llcrnrlat':lats.min(), 'llcrnrlon':lons.min(), 'urcrnrlat':lats.max() , 'urcrnrlon':lons.max(), 'lat_ts':lats.mean()}
	loc_dict=kwargs.get('loc_dict', def_loc_dict)
	map= Basemap(projection='merc', resolution='l',area_thresh=1., **loc_dict)
	xx, yy = map(lons, lats)
	#map.drawcoastlines()
	#map.drawcountries()
	map.drawmapboundary()
	map.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	map.drawmeridians(array([129,130,131,132,133]), labels=[1,0,0,1])
	map.drawparallels(array([-14,-13,-12,-11,-10]), labels=[1,0,0,1])
	levs_dict={'VR':linspace(-15,15,31), 'CZ': linspace(-8,64,10), 'PH': linspace(0,185,255), "RH": linspace(0,1.5,16), "SW":linspace(0, 5, 11), "ZD":linspace(-10,10,21), 'VE':linspace(-30,30,31), 'TI':linspace(-30,30,31), 'KD':linspace(-1.0,6.0,30)}
	titles_dict={'VR': 'Velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'PH': 'Differential Prop Phase (degrees)', 'RH':'Correlation Co-ef', 'SW':'Spectral Width (m/s)', 'ZD':'Differentail Reflectivity dBz', 'VE':'Edited Velocity (m/s)', 'TI':'Simualated winds,(m/s)', 'KD':'Specific differential Phase (Degrees/m)'}
	map.contourf(xx,yy,sweep_dict[parm], levels=levs_dict[parm])
	p=sweep_dict['date']
	dtstr='%(#1)02d-%(#2)02d-%(#3)04d %(#4)02d%(#5)02dZ ' %{"#1":p.day, "#2":p.month, "#3":p.year, "#4":p.hour, '#5':p.minute}
	title(sweep_dict['radar_name']+' '+dtstr+titles_dict[parm])
	colorbar()
	savefig(fig_path+fig_name)
	close(f)





def plot_ppi_lobes_qld(sweep_dict, parm, **kwargs):
	#[-27.669166564941406, 152.8619384765625]
	radar_loc=kwargs.get('radar_loc', [-12.2492,  131.0444])
	radar_loc2=kwargs.get('radar_loc2', [-12.2492,  131.0444])
	fig_name=kwargs.get('fig_name', 'ppi_'+parm+'_.png')
	fig_path=kwargs.get('fig_path', getenv('HOME')+'/bom_mds/output/')
	f=figure()
	#gp_loc=[-12.2492,  131.0444]
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(radar_loc[0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=radar_loc[1]+360.0*sweep_dict['xar']/(rad_at_radar*2.0*pi)
	lats=radar_loc[0] + 360.0*sweep_dict['yar']/(Re*2.0*pi)
	alats=linspace(radar_loc[0]-2.0, radar_loc[0]+2.0, 100.0)
	alons=linspace(radar_loc[1]-2.0, radar_loc[1]+2.0, 100.0)
	angs=array(propigation.make_lobe_grid(radar_loc, radar_loc2, alats,alons))
	def_loc_dict={'lat_0':lats.mean(), 'lon_0':lons.mean(),'llcrnrlat':lats.min(), 'llcrnrlon':lons.min(), 'urcrnrlat':lats.max() , 'urcrnrlon':lons.max(), 'lat_ts':lats.mean()}
	loc_dict=kwargs.get('loc_dict', def_loc_dict)
	map= Basemap(projection='merc', resolution='l',area_thresh=1., **loc_dict)
	xx, yy = map(lons, lats)
	galons, galats=meshgrid(alons, alats)
	xxl, yyl=map(galons, galats)
	#map.drawcoastlines()
	#map.drawcountries()
	map.drawmapboundary()
	map.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstqldmd_r', 'qcoast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	map.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	map.drawmeridians(array([152,152.5,153,153.5,154]), labels=[1,0,0,1])
	map.drawparallels(array([-29,-28.5,-28,-27.5,-27]), labels=[1,0,0,1])
	levs_dict={'VR':linspace(-15,15,31), 'CZ': linspace(-8,64,10), 'PH': linspace(0,185,255), "RH": linspace(0,1.5,16), "SW":linspace(0, 5, 11), "ZD":linspace(-10,10,21), 'VE':linspace(-30,30,31), 'TI':linspace(-30,30,31), 'KD':linspace(-1.0,6.0,30)}
	titles_dict={'VR': 'Velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'PH': 'Differential Prop Phase (degrees)', 'RH':'Correlation Co-ef', 'SW':'Spectral Width (m/s)', 'ZD':'Differentail Reflectivity dBz', 'VE':'Edited Velocity (m/s)', 'TI':'Simualated winds,(m/s)', 'KD':'Specific differential Phase (Degrees/m)'}
	map.contourf(xx,yy,sweep_dict[parm], levels=levs_dict[parm])
	colorbar()
	map.contour(xxl,yyl,angs, levels=[30.0, 150.0],colors=['r']) 
	p=sweep_dict['date']
	dtstr='%(#1)02d-%(#2)02d-%(#3)04d %(#4)02d%(#5)02dZ ' %{"#1":p.day, "#2":p.month, "#3":p.year, "#4":p.hour, '#5':p.minute}
	title(sweep_dict['radar_name']+' '+dtstr+titles_dict[parm])
	savefig(fig_path+fig_name)
	close(f)


def old_plot_cappi(x,y,data,sweep_dict, **kwargs):
	radar_loc=kwargs.get('radar_loc', [-12.2492,  131.0444])
	parm=kwargs.get('parm','')
	fig_name=kwargs.get('fig_name', 'cappi_'+parm+'_.png')
	f=figure()
	#gp_loc=[-12.2492,  131.0444]
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(radar_loc[0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=radar_loc[1]+360.0*x/(rad_at_radar*2.0*pi)
	lats=radar_loc[0] + 360.0*y/(Re*2.0*pi)
	def_loc_dict={'lat_0':lats.mean(), 'lon_0':lons.mean(),'llcrnrlat':lats.min(), 'llcrnrlon':lons.min(), 'urcrnrlat':lats.max() , 'urcrnrlon':lons.max(), 'lat_ts':lats.mean()}
	loc_dict=kwargs.get('loc_dict', def_loc_dict)
	map= Basemap(projection='merc', resolution='l',area_thresh=1., **loc_dict)
	longr, latgr=meshgrid(lons,lats)
	xx, yy = map(longr, latgr)
	#map.drawcoastlines()
	#map.drawcountries()
	map.drawmapboundary()
	map.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	map.drawmeridians(array([129,130,131,132,133]), labels=[1,0,0,1])
	map.drawparallels(array([-14,-13,-12,-11,-10]), labels=[1,0,0,1])
	comp_levs=linspace(-1.,1., 30)
	levs_dict={'VR':linspace(-15,15,31), 'CZ': linspace(-8,64,10), 'TEST':linspace(((abs(data)+data)/2.0).min(), data.max(), 200), 'i_comp':comp_levs,'j_comp':comp_levs,'k_comp':comp_levs}
	titles_dict={'VR': 'Velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST', 'i_comp':'i component','j_comp':'j component','k_comp':'k component'}
	map.contourf(xx,yy,data, levels=levs_dict[parm])
	p=sweep_dict['date']
	dtstr='%(#1)02d-%(#2)02d-%(#3)04d %(#4)02d%(#5)02dZ ' %{"#1":p.day, "#2":p.month, "#3":p.year, "#4":p.hour, '#5':p.minute}
	title(sweep_dict['radar_name']+' '+dtstr+titles_dict[parm])
	colorbar()
	savefig(getenv('HOME')+'/bom_mds/output/'+fig_name)
	close(f)


def plot_cappi(data,parm,lnum, **kwargs):
	x=data['xar']
	y=data['yar']
	radar_loc=data['radar_loc']
	radar_name=kwargs.get('radar_name', data['radar_name'])
	fig_name=kwargs.get('fig_name', 'cappi_'+parm+'_.png')
	f=figure()
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(radar_loc[0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=radar_loc[1]+360.0*(x-data['displacement'][0])/(rad_at_radar*2.0*pi)
	lats=radar_loc[0] + 360.0*(y-data['displacement'][1])/(Re*2.0*pi)
	def_loc_dict={'lat_0':lats.mean(), 'lon_0':lons.mean(),'llcrnrlat':lats.min(), 'llcrnrlon':lons.min(), 'urcrnrlat':lats.max() , 'urcrnrlon':lons.max(), 'lat_ts':lats.mean()}
	loc_dict=kwargs.get('loc_dict', def_loc_dict)
	#print 'Here'
	map= Basemap(projection='merc', resolution='l',area_thresh=1., **loc_dict)
	longr, latgr=meshgrid(lons,lats)
	xx, yy = map(longr, latgr)
	#print 'there'
	map.drawmapboundary()
	darwin_loc=[-12.5, 130.85]
	disp_from_darwin=mathematics.corner_to_point(darwin_loc, radar_loc)
	dist_from_darwin=sqrt(disp_from_darwin[0]**2+disp_from_darwin[1]**2)
	if dist_from_darwin < 500.*1000.: 
		map.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	else:
		map.drawcoastlines()
	map.drawmeridians(array([129,130,131,132,133]), labels=[1,0,0,1])
	map.drawparallels(array([-14,-13,-12,-11,-10]), labels=[1,0,0,1])
	comp_levs=linspace(-1.,1., 30)
	#print 'everywhere'
	levs_dict={'VR':linspace(-15,15,31), 'CZ': linspace(-8,64,10), 'PH': linspace(0,185,255), "RH": linspace(0,1.5,16), "SW":linspace(0, 5, 11), "ZD":linspace(-10,10,21), 'VE':linspace(-30,30,31), 'TI':linspace(-30,30,31)}
	#print 'or here?'
	titles_dict={'VR': 'Velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'PH': 'Differential Prop Phase (degrees)', 'RH':'Correlation Co-ef', 'SW':'Spectral Width (m/s)', 'ZD':'Differentail Reflectivity dBz', 'VE':'Edited Velocity (m/s)', 'TI':'Simualated winds,(m/s)'}
	#print parm
	if 'mask' in kwargs.keys():
		mask=kwargs['mask'][:,:,lnum]
		mdata=M.masked_where(M.array(mask) < 0.5, M.array(data[parm][:,:,lnum]))
		map.contourf(xx,yy,mdata, levels=levs_dict[parm])
	else:
		map.contourf(xx,yy,data[parm][:,:,lnum], levels=levs_dict[parm])
	colorbar()
	if 'angs' in kwargs.keys():
		map.contour(xx,yy,kwargs['angs'], levels=[30.0, 150.0],colors=['r']) 
	p=data['date']
	dtstr='%(#1)02d-%(#2)02d-%(#3)04d %(#4)02d%(#5)02dZ ' %{"#1":p.day, "#2":p.month, "#3":p.year, "#4":p.hour, '#5':p.minute}
	title(radar_name+' '+dtstr+titles_dict[parm])
	savefig(getenv('HOME')+'/bom_mds/output/'+fig_name)
	close(f)



def plot_rhi(sar,zar, data,sweep_dict, **kwargs):
	parm=kwargs.get('parm','CZ')
	fig_name=kwargs.get('fig_name', 'rhi_'+parm+'_.png')
	f=figure()
	#gp_loc=[-12.2492,  131.0444]
	levs_dict={'VR':linspace(-15,15,31), 'CZ': linspace(-8,64,10), 'TEST':linspace(((abs(data)+data)/2.0).min(), data.max(), 200)}
	titles_dict={'VR': 'Velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST'}
	contourf(sar,zar,data, levels=levs_dict[parm])
	p=sweep_dict['date']
	dtstr='%(#1)02d-%(#2)02d-%(#3)04d %(#4)02d%(#5)02dZ ' %{"#1":p.day, "#2":p.month, "#3":p.year, "#4":p.hour, '#5':p.minute}
	title(sweep_dict['radar_name']+' '+dtstr+titles_dict[parm])
	colorbar()
	savefig(getenv('HOME')+'/bom_mds/output/'+fig_name)
	close(f)

def plot_sonde(sonde_dict, fname):
	f=figure()
	climt.thermodyn.skewT(p=sonde_dict['press(hPa)'], T=sonde_dict['tdry(degs)'], Td=sonde_dict['dp(degs)'])
	savefig(getenv('HOME')+'/bom_mds/output/'+fname)
	close(f)


def plot_slice_lat(lon_sl, data_cube, lats, lons, levs, v,w, mask, par='CZ', w_mag=2.0, **kwargs):
	ksp=kwargs.get('ksp', 0.05)
	bquiver=kwargs.get('bquiver', [0.95, 0.75])
	my_lon=argsort(abs(lons-lon_sl))[0]
	print lons[my_lon]
	data=M.masked_where(M.array(mask[:,my_lon,:])<0.5, M.array(data_cube[par][:,my_lon,:]))
	wm=M.masked_where(M.array(mask[:,my_lon,:]) < 0.5, M.array(w[:,my_lon,:]))
	vm=M.masked_where(M.array(mask[:,my_lon,:]) < 0.5, M.array(v[:,my_lon,:]))
	if par in ['diff']:
		data=M.masked_where(M.array(mask) < 0.5, M.array(data))
	comp_levs=linspace(-1.,1., 30)
	levs_dict={'VR':linspace(-25,25,51), 'CZ': linspace(-8,64,10), 'TEST':linspace(((abs(data)+data)/2.0).min(), data.max(), 200), 'i_comp':comp_levs,'j_comp':comp_levs,'k_comp':comp_levs,'diff':linspace(-15,15,31), 'VE':linspace(-25,25,51)}
	titles_dict={'VR': 'Radial velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST', 'i_comp':'i component','j_comp':'j component','k_comp':'k component', 'diff':'diff', 'VE':'Edited velocity'}
	print data.shape
	print wm.shape
	print vm.shape
	print lats.shape
	print levs.shape
	xx,yy=meshgrid(lats,levs)
	#contourf(xx, yy, transpose(data), levs_dict[par])
	pcolor(xx, yy, transpose(data), vmin=levs_dict[par].min(), vmax=levs_dict[par].max())
	colorbar()
	qq=quiver(xx, yy, transpose(vm), transpose(wm*w_mag), scale=kwargs.get('qscale', 200))
	quiverkey(qq, bquiver[0], bquiver[1]+2.*ksp, 10, '10 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	print bquiver[0], bquiver[1]+2.*ksp
	quiverkey(qq, bquiver[0], bquiver[1]+ksp, 5, '5 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	quiverkey(qq, bquiver[0], bquiver[1], 2, '2 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	xlabel('Latitude (degrees)')
	ylabel('Height (Metres)')
	if 'box' in kwargs.keys():
		box=kwargs['box']
		ax=gca()
		cr=ax.axis()
		ax.axis([box[2], box[3], cr[2], cr[3]])
	return lons[my_lon]


def reconstruction_plot_w2(lats, lons, data_cube, lvl_num, parm, u, v,w, angs, mask, **kwargs):
	ber_loc=[-12.4, 130.85] #location of Berrimah radar
	gp_loc=[-12.2492,  131.0444]#location of CPOL at Gunn Point
	box=kwargs.get('box',[130.0, 132.0, -13.0, -11.5])#box for the plot
	bquiver=kwargs.get('bquiver', [0.1, 0.75])
	ksp=kwargs.get('ksp', 0.05)
	#Set up the map and projection
	mapobj= Basemap(projection='merc',lat_0=(box[2]+box[3])/2.0, lon_0=(box[0]+box[1])/2.0,llcrnrlat=box[2], llcrnrlon=box[0], urcrnrlat=box[3] , urcrnrlon=box[1], resolution='l',area_thresh=1., lat_ts=(box[2]+box[3])/2.0)
	#map.drawmapboundary()
	um=M.masked_where(M.array(mask) < 0.5, M.array(u))
	vm=M.masked_where(M.array(mask) < 0.5, M.array(v))
	mag=M.array(sqrt(u**2+v**2))
	magm=M.masked_where(M.array(mask)<0.5, mag)
	fig_name=kwargs.get('fig_name', 'recon_'+parm+'_.png')
	longr, latgr=meshgrid(lons,lats)
	xx, yy = mapobj(longr, latgr)
	mapobj.drawmapboundary()
	mapobj.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	mapobj.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstqldmd_r', 'qcoast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	#mapobj.drawmeridians(array([130,130.5, 131, 131.5, 132]), labels=[1,0,0,1])
	#mapobj.drawparallels(array([-13,-12.5,-12,-11.5,-11]), labels=[1,0,0,1])
	data=data_cube[parm][:,:,lvl_num]
	if parm in ['diff']:
		data=M.masked_where(M.array(mask) < 0.5, M.array(data))
	comp_levs=linspace(-1.,1., 30)
	levs_dict={'VR':linspace(-25,25,51), 'CZ': linspace(-8,64,10), 'TEST':linspace(((abs(data)+data)/2.0).min(), data.max(), 200), 'i_comp':comp_levs,'j_comp':comp_levs,'k_comp':comp_levs,'diff':linspace(-15,15,31), 'VE':linspace(-25,25,51)}
	titles_dict={'VR': 'Radial velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST', 'i_comp':'i component','j_comp':'j component','k_comp':'k component', 'diff':'diff', 'VE':'Edited velocity'}
	mapobj.contourf(xx,yy,data, levels=levs_dict[parm])
	colorbar()
	#mapobj.quiver(xx,yy,u/sqrt(u**2+v**2),v/sqrt(u**2+v**2), scale=50)
	qq=mapobj.quiver(xx,yy,um,vm, scale=kwargs.get('qscale', 200))
	quiverkey(qq, bquiver[0], bquiver[1]+2.*ksp, 10, '10 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	print bquiver[0], bquiver[1]+2.*ksp
	quiverkey(qq, bquiver[0], bquiver[1]+ksp, 5, '5 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	quiverkey(qq, bquiver[0], bquiver[1], 2, '2 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	cobject=mapobj.contour(xx,yy,w, colors=['k'], levels=linspace(-10,10,6))
	fon = { 'fontname':'Tahoma', 'fontsize':5 }
	clabel(cobject, fmt="%1.1f", **fon)
	mapobj.contour(xx,yy,angs, levels=[30.0, 150.0],colors=['r']) 
	mapobj.readshapefile(getenv('HOME')+'/bom_mds/shapes/nt_coast','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	mapobj.drawmeridians(linspace(0,4,41)+130., labels=[1,0,0,1], fontsize=rcParams['xtick.labelsize'])
	mapobj.drawparallels(linspace(0,4,41)-15.0, labels=[1,0,0,1], fontsize=rcParams['ytick.labelsize'])
	return mapobj

def plot_slice_lon(lat_sl, data_cube, lats, lons, levs, u,w, mask, par='CZ', w_mag=2.0, **kwargs):
	ksp=kwargs.get('ksp', 0.05)
	bquiver=kwargs.get('bquiver', [0.1, 0.2])
	my_lat=argsort(abs(lats-lat_sl))[0]
	print lats[my_lat]
	data=M.masked_where(M.array(mask[my_lat,:,:])<0.5, M.array(data_cube[par][my_lat,:,:]))
	wm=M.masked_where(M.array(mask[my_lat,:,:]) < 0.5, M.array(w[my_lat,:,:]))
	um=M.masked_where(M.array(mask[my_lat,:,:]) < 0.5, M.array(u[my_lat,:,:]))
	if par in ['diff']:
		data=M.masked_where(M.array(mask) < 0.5, M.array(data))
	comp_levs=linspace(-1.,1., 30)
	levs_dict={'VR':linspace(-25,25,51), 'CZ': linspace(-8,64,10), 'TEST':linspace(((abs(data)+data)/2.0).min(), data.max(), 200), 'i_comp':comp_levs,'j_comp':comp_levs,'k_comp':comp_levs,'diff':linspace(-15,15,31), 'VE':linspace(-25,25,51)}
	titles_dict={'VR': 'Radial velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST', 'i_comp':'i component','j_comp':'j component','k_comp':'k component', 'diff':'diff', 'VE':'Edited velocity'}
	xx,yy=meshgrid(lons,levs)
	#gca().xaxis.major.formatter.set_powerlimits((-10.,10.))
	#contourf(xx, yy, transpose(data), levs_dict[par])
	pcolor(xx, yy, transpose(data), vmin=levs_dict[par].min(), vmax=levs_dict[par].max())
	#gca().xaxis.major.formatter.set_scientific(False)
	#print myaxis.formatter.limits
	colorbar()
	print "scale ", kwargs.get('qscale', 200)
	qq=quiver(xx, yy, transpose(um), transpose(wm*w_mag), scale=kwargs.get('qscale', 200))
	quiverkey(qq, bquiver[0], bquiver[1]+2.*ksp, 10, '10 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	print bquiver[0], bquiver[1]+2.*ksp
	quiverkey(qq, bquiver[0], bquiver[1]+ksp, 5, '5 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	quiverkey(qq, bquiver[0], bquiver[1], 2, '2 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	xlabel('Longitude (degrees)')
	ylabel('Height (Metres)')
	if 'box' in kwargs.keys():
		box=kwargs['box']
		ax=gca()
		cr=ax.axis()
		ax.axis([box[0], box[1], cr[2], cr[3]])
	return lats[my_lat]
	

def cmap_discretize(cmap, N):
	cdict = cmap._segmentdata.copy()
	# N colors
	colors_i = linspace(0,1.,N)
	# N+1 indices
	indices = linspace(0,1.,N+1)
	for key in ('red','green','blue'):
		# Find the N colors
		D = array(cdict[key])
		I = interpolate.interp1d(D[:,0], D[:,1])
		colors = I(colors_i)
		# Place these colors at the correct indices.
		A = zeros((N+1,3), float)
		A[:,0] = indices
		A[1:,1] = colors
		A[:-1,2] = colors
		# Create a tuple for the dictionary.
		L = []
		for l in A:
			L.append(tuple(l))
		cdict[key] = tuple(L)
	# Return colormap object.
	return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)


def plot_slice_lon_hydro(lat_sl, data_cube, lats, lons, levs, hydro_cube, u,w, mask, par='CZ', w_mag=2.0, **kwargs):
	ksp=kwargs.get('ksp', 0.05)
	bquiver=kwargs.get('bquiver', [0.1, 0.2])
	my_lat=argsort(abs(lats-lat_sl))[0]
	my_lat_hydro=argsort(abs(hydro_cube['lats']-lat_sl))[0]
	print lats[my_lat]
	data=M.masked_where(M.array(mask[my_lat,:,:])<0.5, M.array(data_cube[par][my_lat,:,:]))
	hydro_data=M.masked_where(classify[my_lat_hydro,:,:]==-9.0, classify[my_lat_hydro,:,:])
	wm=M.masked_where(M.array(mask[my_lat,:,:]) < 0.5, M.array(w[my_lat,:,:]))
	um=M.masked_where(M.array(mask[my_lat,:,:]) < 0.5, M.array(u[my_lat,:,:]))
	if par in ['diff']:
		data=M.masked_where(M.array(mask) < 0.5, M.array(data))
	comp_levs=linspace(-1.,1., 30)
	levs_dict={'VR':linspace(-15,15,31), 'CZ': linspace(-8,64,10), 'TEST':linspace(((abs(data)+data)/2.0).min(), data.max(), 200), 'i_comp':comp_levs,'j_comp':comp_levs,'k_comp':comp_levs,'diff':linspace(-15,15,31)}
	titles_dict={'VR': 'Radial velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST', 'i_comp':'i component','j_comp':'j component','k_comp':'k component', 'diff':'diff'}
	xx,yy=meshgrid(lons,levs)
	xx_hydro, yy_hydro=meshgrid(hydro_cube['lons'], hydro_cube['zar'])
	#gca().xaxis.major.formatter.set_powerlimits((-10.,10.))
	#gca().xaxis.major.formatter.set_scientific(False)
	#print myaxis.formatter.limits
	print "my_lat_hydro:", my_lat_hydro
	print hydro_cube['lats'][my_lat_hydro]
	pcolor(xx_hydro, yy_hydro, transpose(hydro_data), vmin=-0.5, vmax=10.5, cmap=cmap_discretize(cm.jet, 11))
	myc=colorbar()
	myc.ax.text(0.0, 3.0, "HD snow")
	#yticks((0,1,2,3),("Uncl", "Dz", "Ra", "HD Sn"))
	contour(xx, yy, transpose(data), levs_dict[par])
	qq=quiver(xx, yy, transpose(um), transpose(wm*w_mag), scale=kwargs.get('qscale', 200))
	quiverkey(qq, bquiver[0], bquiver[1]+2.*ksp, 10, '10 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	print bquiver[0], bquiver[1]+2.*ksp
	quiverkey(qq, bquiver[0], bquiver[1]+ksp, 5, '5 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	quiverkey(qq, bquiver[0], bquiver[1], 2, '2 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	xlabel('Longitude (degrees)')
	ylabel('Height (Metres)')
	if 'box' in kwargs.keys():
		box=kwargs['box']
		ax=gca()
		cr=ax.axis()
		ax.axis([box[0], box[1], cr[2], cr[3]])
	return lats[my_lat]

def make_hydro_bins(hydro_cube):
	class_dict={"unclassified":0, "Drizzle":1, "Rain":2, "Dry LD snow":3, "Dry HD snow":4, "Melting Snow":5, "Dry Graupel":6, "Wet Graupel": 7, "Small Hail": 8, "Large Hail":9, "Rain + Hail":10}
	hydro_bins={"bins":class_dict.keys()}
	for class_name in class_dict.keys():
		class_lats=[]
		class_lons=[]
		class_z=[]
		cube_shape=hydro_cube['classify'].shape
		for i in range(cube_shape[0]):
			for j in range(cube_shape[1]):
				for k in range(cube_shape[2]):
					if hydro_cube['classify'][i,j,k]==class_dict[class_name]:
						class_lats.append(hydro_cube['lats'][i])
						class_lons.append(hydro_cube['lons'][j])
						class_z.append(hydro_cube['zar'][k])
		hydro_bins.update({class_name+"_lats":array(class_lats), class_name+'_lons':array(class_lons), class_name+"_z":array(class_z)})
	return hydro_bins

def make_hydro_bins_clat(hydro_cube, lat):
	i=argsort(abs(hydro_cube['lats']-lat))[0]
	class_dict={"unclassified":0, "Drizzle":1, "Rain":2, "Dry LD snow":3, "Dry HD snow":4, "Melting Snow":5, "Dry Graupel":6, "Wet Graupel": 7, "Small Hail": 8, "Large Hail":9, "Rain + Hail":10}
	
	hydro_bins={"bins":class_dict.keys()}
	for class_name in class_dict.keys():
		class_lats=[]
		class_lons=[]
		class_z=[]
		cube_shape=hydro_cube['classify'].shape
		for j in range(cube_shape[1]):
			for k in range(cube_shape[2]):
				if hydro_cube['classify'][i,j,k]==class_dict[class_name]:
					class_lats.append(hydro_cube['lats'][i])
					class_lons.append(hydro_cube['lons'][j])
					class_z.append(hydro_cube['zar'][k])
		hydro_bins.update({class_name+"_lats":array(class_lats), class_name+'_lons':array(class_lons), class_name+"_z":array(class_z)})
	return hydro_bins


def plot_slice_lon_hydro_sym(lat_sl, data_cube, lats, lons, levs, hydro_bins, u,w, mask, par='CZ', w_mag=2.0, **kwargs):
	ksp=kwargs.get('ksp', 0.05)
	bquiver=kwargs.get('bquiver', [0.1, 0.2])
	my_lat=argsort(abs(lats-lat_sl))[0]
	print lats[my_lat]
	data=M.masked_where(M.array(mask[my_lat,:,:])<0.5, M.array(data_cube[par][my_lat,:,:]))
	wm=M.masked_where(M.array(mask[my_lat,:,:]) < 0.5, M.array(w[my_lat,:,:]))
	um=M.masked_where(M.array(mask[my_lat,:,:]) < 0.5, M.array(u[my_lat,:,:]))
	if par in ['diff']:
		data=M.masked_where(M.array(mask) < 0.5, M.array(data))
	comp_levs=linspace(-1.,1., 30)
	levs_dict={'VR':linspace(-15,15,31), 'CZ': linspace(-8,64,10), 'TEST':linspace(((abs(data)+data)/2.0).min(), data.max(), 200), 'i_comp':comp_levs,'j_comp':comp_levs,'k_comp':comp_levs,'diff':linspace(-15,15,31)}
	titles_dict={'VR': 'Radial velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST', 'i_comp':'i component','j_comp':'j component','k_comp':'k component', 'diff':'diff'}
	xx,yy=meshgrid(lons,levs)
	plotting_syms={"unclassified":'k.', "Drizzle":'co', "Rain":'bo', "Dry LD snow":'cx', "Dry HD snow":'bx', "Melting Snow":'r+', "Dry Graupel":'bs', "Wet Graupel": 'rs', "Small Hail": 'cd', "Large Hail":'bD', "Rain + Hail":'g>'}
	for item in set(hydro_bins['bins'])-set(['unclassified']):
		plot(hydro_bins[item+"_lons"], hydro_bins[item+"_z"], plotting_syms[item])
	legend(set(hydro_bins['bins'])-set(['unclassified']))
	#gca().xaxis.major.formatter.set_powerlimits((-10.,10.))
	#gca().xaxis.major.formatter.set_scientific(False)
	#print myaxis.formatter.limits
	#pcolor(xx_hydro, yy_hydro, transpose(hydro_data), vmin=-0.5, vmax=10.5, cmap=cmap_discretize(cm.jet, 11))
	#colorbar()
	contour(xx, yy, transpose(data), levs_dict[par])
	qq=quiver(xx, yy, transpose(um), transpose(wm*w_mag), scale=kwargs.get('qscale', 200))
	quiverkey(qq, bquiver[0], bquiver[1]+2.*ksp, 10, '10 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	print bquiver[0], bquiver[1]+2.*ksp
	quiverkey(qq, bquiver[0], bquiver[1]+ksp, 5, '5 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	quiverkey(qq, bquiver[0], bquiver[1], 2, '2 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	xlabel('Longitude (degrees)')
	ylabel('Height (Metres)')
	#item="Rain"
	if 'box' in kwargs.keys():
		box=kwargs['box']
		ax=gca()
		cr=ax.axis()
		ax.axis([box[0], box[1], cr[2], cr[3]])
	return lats[my_lat]

def plot_slice_lon_hydro_sym_refl(lat_sl, data_cube,hydro_cube, hydro_bins, u,w, mask, par='CZ', w_mag=2.0, **kwargs):
	lats=hydro_cube['lats']
	lons=hydro_cube['lons']
	levs=hydro_cube['zar']
	ksp=kwargs.get('ksp', 0.05)
	bquiver=kwargs.get('bquiver', [0.1, 0.2])
	my_lat=argsort(abs(lats-lat_sl))[0]
	print lats[my_lat]
	data=M.masked_where(M.array(hydro_cube['CZ'][my_lat,:,:])<0.0, M.array(hydro_cube['CZ'][my_lat,:,:]))
	wm=M.masked_where(M.array(mask[my_lat,:,:]) < 0.5, M.array(w[my_lat,:,:]))
	um=M.masked_where(M.array(mask[my_lat,:,:]) < 0.5, M.array(u[my_lat,:,:]))
	if par in ['diff']:
		data=M.masked_where(M.array(mask) < 0.5, M.array(data))
	comp_levs=linspace(-1.,1., 30)
	levs_dict={'VR':linspace(-15,15,31), 'CZ': linspace(-8,64,10), 'TEST':linspace(((abs(data)+data)/2.0).min(), data.max(), 200), 'i_comp':comp_levs,'j_comp':comp_levs,'k_comp':comp_levs,'diff':linspace(-15,15,31)}
	titles_dict={'VR': 'Radial velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST', 'i_comp':'i component','j_comp':'j component','k_comp':'k component', 'diff':'diff'}
	xx,yy=meshgrid(lons,levs)
	plotting_syms={"unclassified":'k.', "Drizzle":'co', "Rain":'bo', "Dry LD snow":'cx', "Dry HD snow":'bx', "Melting Snow":'r+', "Dry Graupel":'bs', "Wet Graupel": 'rs', "Small Hail": 'cd', "Large Hail":'bD', "Rain + Hail":'g>'}
	for item in set(hydro_bins['bins'])-set(['unclassified']):
		plot(hydro_bins[item+"_lons"], hydro_bins[item+"_z"], plotting_syms[item])
	legend(set(hydro_bins['bins'])-set(['unclassified']))
	#gca().xaxis.major.formatter.set_powerlimits((-10.,10.))
	#gca().xaxis.major.formatter.set_scientific(False)
	#print myaxis.formatter.limits
	#pcolor(xx_hydro, yy_hydro, transpose(hydro_data), vmin=-0.5, vmax=10.5, cmap=cmap_discretize(cm.jet, 11))
	#colorbar()
	contour(xx, yy, transpose(data), levs_dict[par])
	print "QSCALE ", kwargs.get('qscale', 200)
	qq=quiver(xx, yy, transpose(um), transpose(wm*w_mag), scale=kwargs.get('qscale', 200))
	quiverkey(qq, bquiver[0], bquiver[1]+2.*ksp, 10, '10 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	print bquiver[0], bquiver[1]+2.*ksp
	quiverkey(qq, bquiver[0], bquiver[1]+ksp, 5, '5 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	quiverkey(qq, bquiver[0], bquiver[1], 2, '2 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	xlabel('Longitude (degrees)')
	ylabel('Height (Metres)')
	#item="Rain"
	if 'box' in kwargs.keys():
		box=kwargs['box']
		ax=gca()
		cr=ax.axis()
		ax.axis([box[0], box[1], cr[2], cr[3]])
	return lats[my_lat]



def plot_slices(lat_sl, lon_sl, lvl, data_cube, lats, lons, levs, u, v, w, angs, mask, par='CZ', w_mag=2.0, **kwargs):
	old_labelsize=rcParams['xtick.labelsize']
	rcdict={'labelsize':8}
	rc('xtick', **rcdict)
	rc('ytick', **rcdict)
	xy_kwargs=dict([(key,kwargs.get(key)) for key in set(kwargs.keys())&set(['box', 'bquiver', 'ksp','qscale']) ])
	side_args=dict([(key,kwargs.get(key)) for key in set(kwargs.keys())&set(['ksp', 'box','qscale']) ])
	lvl_num=argsort(abs(levs-lvl))[0]
	print lvl_num
	subplot(2,2,2)
	act_lon=plot_slice_lat(lon_sl, data_cube, lats, lons, levs, v,w, mask, par=par, w_mag=w_mag, **side_args)
	subplot(2,2,3)
	act_lat=plot_slice_lon(lat_sl, data_cube, lats, lons, levs, u,w, mask, par=par, w_mag=w_mag, bquiver=[0.05, 0.33], **side_args)
	ax=gcf().gca()
	mf=ax.xaxis.get_major_formatter()
	ax.xaxis.set_major_formatter(FormatStrFormatter('%4.1f'))
	subplot(2,2,1)
	#reconstruction_plot_w2(lats, lons, data_cube, lvl_num, parm, u, v,w, angs, mask, **kwargs)
	mapobj=reconstruction_plot_pcolor(lats, lons, data_cube, lvl_num, par, u[:,:,lvl_num], v[:,:,lvl_num],w[:,:,lvl_num], angs, mask[:,:,lvl_num], **xy_kwargs)
	lon_cross_lats=ones(lats.shape, dtype=float)*act_lat
	lon_cross_lons=lons
	lat_cross_lats=lats
	lat_cross_lons=ones(lons.shape, dtype=float)*act_lon
	lon_cross_x, lon_cross_y=mapobj(lon_cross_lons, lon_cross_lats)
	lat_cross_x, lat_cross_y=mapobj(lat_cross_lons, lat_cross_lats)
	mapobj.plot(lon_cross_x, lon_cross_y, 'r--')
	mapobj.plot(lat_cross_x, lat_cross_y, 'r--')
	old_rdict={'labelsize':old_labelsize}
	rc('xtick', **old_rdict)
	rc('ytick', **old_rdict)
	return act_lat, act_lon, levs[lvl_num]

def dump_radar(radar, radar_loc, **kwargs):
	subdir=kwargs.get('subdir', '/bm/gscratch/scollis/dumps/')
	date_obj=radar[0]['date']
	date_dic={'y':date_obj.year, 'm':date_obj.month, 'd':date_obj.day, 'H':date_obj.hour, 'M':date_obj.minute}
	dfname=radar[0]['radar_name']+"%(y)04d%(m)02d%(d)02d_%(H)02d%(M)02d" %date_dic
	try:
		os.mkdir(subdir+dfname)
	except OSError:
		print "directory exists" 
	gated_vars=set(['VR', 'CZ', 'RH', 'PH', 'ZD', 'SW', 'VE', 'TI'])
	av_vars=set(radar[0].keys())&gated_vars
	for parm in av_vars:
		try: 
			os.mkdir(subdir+dfname+'/'+parm)
		except OSError:
			print "directory exists" 
		for scan in radar:
			elev_str="_%(el)04.1f" %{'el':scan['Elev'][-1]}
			fname=subdir+dfname+'/'+parm+'/'+dfname+elev_str+".png"
			print fname
			plot_ppi(scan, parm, fig_path=subdir+dfname+'/'+parm+'/', fig_name=dfname+elev_str+".png", radar_loc=radar_loc)


def reconstruction_plot_pcolor(lats, lons, data_cube, lvl_num, parm, u, v,w, angs, mask, **kwargs):
	ber_loc=[-12.4, 130.85] #location of Berrimah radar
	gp_loc=[-12.2492,  131.0444]#location of CPOL at Gunn Point
	box=kwargs.get('box',[130.0, 132.0, -13.0, -11.5])#box for the plot
	bquiver=kwargs.get('bquiver', [0.1, 0.75])
	ksp=kwargs.get('ksp', 0.05)
	#Set up the map and projection
	mapobj= Basemap(projection='merc',lat_0=(box[2]+box[3])/2.0, lon_0=(box[0]+box[1])/2.0,llcrnrlat=box[2], llcrnrlon=box[0], urcrnrlat=box[3] , urcrnrlon=box[1], resolution='l',area_thresh=1., lat_ts=(box[2]+box[3])/2.0)
	#map.drawmapboundary()
	um=M.masked_where(M.array(mask) < 0.5, M.array(u))
	vm=M.masked_where(M.array(mask) < 0.5, M.array(v))
	mag=M.array(sqrt(u**2+v**2))
	magm=M.masked_where(M.array(mask)<0.5, mag)
	fig_name=kwargs.get('fig_name', 'recon_'+parm+'_.png')
	longr, latgr=meshgrid(lons,lats)
	xx, yy = mapobj(longr, latgr)
	mapobj.drawmapboundary()
	mapobj.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstntcd_r','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	mapobj.readshapefile(getenv('HOME')+'/bom_mds/shapes/cstqldmd_r', 'qcoast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	#mapobj.drawmeridians(array([130,130.5, 131, 131.5, 132]), labels=[1,0,0,1])
	#mapobj.drawparallels(array([-13,-12.5,-12,-11.5,-11]), labels=[1,0,0,1])
	data=M.masked_where(M.array(mask) < 0.5, M.array(data_cube[parm][:,:,lvl_num]))
	if parm in ['diff']:
		data=M.masked_where(M.array(mask) < 0.5, M.array(data))
	comp_levs=linspace(-1.,1., 30)
	levs_dict={'VR':linspace(-25,25,51), 'CZ': linspace(-8,64,10), 'i_comp':comp_levs,'j_comp':comp_levs,'k_comp':comp_levs,'diff':linspace(-15,15,31), 'VE':linspace(-25,25,51)}
	titles_dict={'VR': 'Radial velocity m/s', 'CZ':'Corrected Reflectivity dBz', 'TEST':'TEST', 'i_comp':'i component','j_comp':'j component','k_comp':'k component', 'diff':'diff', 'VE':'Edited velocity'}
	mapobj.pcolor(xx,yy,data, vmin=levs_dict[parm].min(), vmax=levs_dict[parm].max())
	colorbar()
	#mapobj.quiver(xx,yy,u/sqrt(u**2+v**2),v/sqrt(u**2+v**2), scale=50)
	qq=mapobj.quiver(xx,yy,um,vm, scale=kwargs.get('qscale', 200))
	quiverkey(qq, bquiver[0], bquiver[1]+2.*ksp, 10, '10 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	print bquiver[0], bquiver[1]+2.*ksp
	quiverkey(qq, bquiver[0], bquiver[1]+ksp, 5, '5 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	quiverkey(qq, bquiver[0], bquiver[1], 2, '2 m/s', coordinates='figure',fontproperties={'size':rcParams['xtick.labelsize']})
	cobject=mapobj.contour(xx,yy,w, colors=['k'], levels=linspace(-10,10,6))
	fon = { 'fontname':'Tahoma', 'fontsize':5 }
	clabel(cobject, fmt="%1.1f", **fon)
	mapobj.contour(xx,yy,angs, levels=[30.0, 150.0],colors=['r']) 
	mapobj.readshapefile(getenv('HOME')+'/bom_mds/shapes/nt_coast','coast',drawbounds=True, linewidth=0.5,color='k',antialiased=1,ax=None)
	mapobj.drawmeridians(linspace(0,4,41)+130., labels=[1,0,0,1], fontsize=rcParams['xtick.labelsize'])
	mapobj.drawparallels(linspace(0,4,41)-15.0, labels=[1,0,0,1], fontsize=rcParams['ytick.labelsize'])
	return mapobj


def put_xscale(length, units, **kwargs):
	#get keywords
	axis=kwargs.get('axis', gca())
	pos=kwargs.get('pos', [0.8,0.8])
	tof=kwargs.get('tof', 0.0)
	#scale the x and y range to graphics co-ordinates
	cr=axis.axis()
	cpos=[pos[0]*(cr[1]-cr[0]) + cr[0], pos[1]*(cr[3]-cr[2]) + cr[2]]
	#determine left and right extent
	extent=[cpos[0]-length/2.0, cpos[0]+length/2.0]
	#plot the line
	plot(extent, zeros(2, dtype=float)+cpos[1], 'k-')
	#generate text
	mytext="%(length)d %(units)s" %{'length':length, 'units':units}
	#annotate line
	text(cpos[0], cpos[1]-tof, mytext)

























