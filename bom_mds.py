###################################################################
#bom_mds: Australian Bureau of Meteorology Multi Doppler Synthesis#
###################################################################
#This top level file contains testbedding and wrapping for lower  #
#level modules. The aim of this project is to wrap fortran code   #
#from McGill University.                                          #
#McGill Authors: Stephane Laroche, Alain Protat, Christian Page,  #
#Alain Caya, Cathy Chiang.                                        #
#The core to this code is a Conjigate Gradient minimization       #
#routine: Using single-Doppler data to obtain a mesoscale         #
# environmental field. Caya A., Laroche S., Zawadzki I.,          #
# and Montmerl T.                                                 #
###################################################################
#Start of BoM Python branch: Scott Collis, CAWCR, April 2008      #
###################################################################
import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
import pres
import simul_winds
import propigation	
import gracon_vel2d
import gracon_vel2d_3d
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
from matplotlib import figure
from pylab import *
from time import time as systime
from numpy import linspace, array, arctan, pi, random
import pickle_zip
import dealias
import parse_ini
import read_radar


def test_mapping():
	f=figure()
	mapobj=pres.generate_darwin_plot()
	savefig(os.getenv('HOME')+'/bom_mds/output/test_mapping.png')
	close(f)

def test_uniform_winds():
	f=figure()
	lats=linspace(-13.0, -11.5, 20)
	lons=linspace(130., 132., 20)
	u,v=simul_winds.unif_wind(lats, lons, 10.0, 45.0)
	mapobj=pres.generate_darwin_plot()
	pres.quiver_winds(mapobj, lats, lons, u,v)
	savefig(os.getenv('HOME')+'/bom_mds/output/test_uniform.png')
	close(f)

def test_pert_winds():
	f=figure()
	lats=linspace(-13.0, -11.5, 20)
	lons=linspace(130., 132., 20)
	u,v=simul_winds.unif_wind(lats, lons, 5.0, 45.0)
	mapobj=pres.generate_darwin_plot()
	fw=0.2#degrees
	m_bump=15.0
	b1=simul_winds.speed_bump(lats, lons, [-12.0, 131.0], fw)*m_bump
	b2=-1.0*simul_winds.speed_bump(lats, lons, [-12.0, 131.2], fw)*m_bump
	pres.quiver_contour_winds(mapobj, lats, lons, u+(b1+b2),v-(b1+b2))
	savefig(os.getenv('HOME')+'/bom_mds/output/test_pert.png')
	close(f)

def test_angles():
	x=150.0*1000.0
	y=125.0*1000.0
	alt=10.0*1000.0
	#dhds(x,y,h,**kwargs)
	gradient=propigation.dhds(x,y,alt, debug=True)
	print gradient
	print arctan(gradient)*180.0/pi


def test_angle_grid2():
	lats=linspace(-13.0, -11.5, 100)
	lons=linspace(130., 132., 100)
	ber_loc=[-12.4, 130.85] #location of Berrimah radar
	gp_loc=[-12.2492,  131.0444]#location of CPOL at Gunn Point
	h=10.0*1000.0
	t0=systime()
	i_a, j_a, k_a=propigation.unit_vector_grid(lats, lons, h, ber_loc)
	print "compute took ", systime()-t0, " seconds"
	f=figure()
	mapobj=pres.generate_darwin_plot()
	pres.contour_comp(mapobj,lats, lons, k_a)
	savefig(os.getenv('HOME')+'/bom_mds/output/test_angs.png')
	close(f)

def test_radar_view():
	lats=linspace(-13.0, -11.5, 40)
	lons=linspace(130., 132., 40)
	ber_loc=[-12.4, 130.85] #location of Berrimah radar
	gp_loc=[-12.2492,  131.0444]#location of CPOL at Gunn Point
	h=2.5*1000.0
	t0=systime()
	i_a, j_a, k_a=propigation.unit_vector_grid(lats, lons, h, gp_loc)
	print "unit_vector_compute took ", systime()-t0, " seconds"
	fw=0.1#degrees
	m_bump=5.0
	#b1=simul_winds.speed_bump(lats, lons, [-12.0, 131.0], fw)*m_bump
	#b2=-1.0*simul_winds.speed_bump(lats, lons, [-12.0, 131.1], fw)*m_bump
	u,v=simul_winds.unif_wind(lats, lons, 5.0, 75.0)
	up,vp=array(simul_winds.vortex(lats, lons,  [-12.5, 131.1], fw))*m_bump
	#up=u+(b1-b2)
	#vp=v-(b1-b2)
	w=v*0.0
	v_r=i_a*(up+u)+j_a*(vp+v)+k_a*w
	f=figure()
	mapobj=pres.generate_darwin_plot()
	pres.contour_vr(mapobj,lats, lons, v_r)
	savefig(os.getenv('HOME')+'/bom_mds/output/test_radar_view_gp.png')
	close(f)
	f=figure()
	mapobj=pres.generate_darwin_plot()
	pres.quiver_contour_winds(mapobj, lats, lons, (up+u),(vp+v))
	savefig(os.getenv('HOME')+'/bom_mds/output/test_pert2.png')
	close(f)


#gv_u,gv_v,f,u_array,v_array = gracon_vel2d(gv_u,gv_v,f,u_array,v_array,i_cmpt_r1,j_cmpt_r1,i_cmpt_r2,j_cmpt_r2,vr1,vr2,nx=shape(gv_u,0),ny=shape(gv_u,1))

def test_gracon():
	#setup
	noise_level=0.0#m/s
	nx=40
	ny=40
	fw=0.1
	m_bump=10.00
	t0=systime()
	lats=linspace(-13.5, -12.0, 40)
	lons=linspace(130.5, 131.5, 40)
	ber_loc=[-12.4, 130.85] #location of Berrimah radar
	gp_loc=[-12.2492,  131.0444]#location of CPOL at Gunn Point
	h=2.5*1000.0
	print 'calculating berimah UV', systime()-t0
	i_ber, j_ber, k_ber=propigation.unit_vector_grid(lats, lons, h, ber_loc)
	print 'calculating gp UV', systime()-t0
	i_gp, j_gp, k_gp=propigation.unit_vector_grid(lats, lons, h, gp_loc)
	#make winds
	u,v=simul_winds.unif_wind(lats, lons, 10.0, 75.0)
	up,vp=array(simul_winds.vortex(lats, lons,  [-12.5, 131.1], fw))*m_bump
	#make V_r measurements
	vr_ber=i_ber*(up+u)+j_ber*(vp+v) + (random.random([nx,ny])-0.5)*(noise_level*2.0)
	vr_gp=i_gp*(up+u)+j_gp*(vp+v)+ (random.random([nx,ny])-0.5)*(noise_level*2.0)
	#try to reconstruct the wind field
	igu, igv= simul_winds.unif_wind(lats, lons, 0.0, 90.0)
	gv_u=zeros(u.shape)
	gv_v=zeros(v.shape)
	f=0.0
	print igu.mean()
	
	angs=array(propigation.make_lobe_grid(ber_loc, gp_loc, lats,lons))
	wts=zeros(angs.shape, dtype=float)+1.0
	#for i in range(angs.shape[0]):
	#	for j in range(angs.shape[1]):
	#		if (angs[i,j] < 150.0) and (angs[i,j] > 30.0): wts[i,j]=1.0
	print 'Into fortran'
	gv_u,gv_v,f,u_array,v_array = gracon_vel2d.gracon_vel2d(gv_u,gv_v,f,igu,igv,i_ber,j_ber,i_gp,j_gp,vr_ber,vr_gp,wts, nx=nx,ny=ny)
	print u_array.mean()
	print f
	bnds=[0.,20.]
	f=figure()
	mapobj=pres.generate_darwin_plot()
	pres.quiver_contour_winds(mapobj, lats, lons, (up+u),(vp+v), bounds=bnds)
	savefig(os.getenv('HOME')+'/bom_mds/output/orig_winds_clean.png')
	close(f)
	f=figure()
	mapobj=pres.generate_darwin_plot()
	pres.quiver_contour_winds(mapobj, lats, lons, (wts*u_array +0.001),(wts*v_array +0.001), bounds=bnds)
	savefig(os.getenv('HOME')+'/bom_mds/output/recon_winds_clean.png')
	close(f)
	f=figure()
	mapobj=pres.generate_darwin_plot()
	pres.quiver_contour_winds(mapobj, lats, lons, (wts*u_array - (up+u)),(wts*v_array -(vp+v)))
	savefig(os.getenv('HOME')+'/bom_mds/output/errors_clean.png')
	close(f)
	
def test_cappi():
	scanme=read_rays.construct_lassen_scan(path='/bm/gdata/scollis/gunn_pt/20060122_0357/')
	read_rays.plot_ppi(scanme[5])
	xar=linspace(-150.,150., 200)*1000.0
	yar=linspace(-150.,150., 200)*1000.0
	mycap=radar_to_cart.make_cappi(scanme, xar, yar, 2500.0, 'CZ')
	read_rays.plot_cappi(xar, yar, mycap)

def test_cappi_ber():
	ber_loc=[-12.457, 130.925]
	gp_loc=	 [-12.2492,  131.0444]
	disp=mathematics.corner_to_point(gp_loc, ber_loc)
	#scanme=read_rays.construct_uf_scan(path='/bm/gdata/scollis/berrimah/20060122_035004/')
	gp_0740=read_rays.construct_lassen_scan(path='/bm/gdata/scollis/gunn_pt/20060122_074001/')
	ber_0740=read_rays.construct_uf_scan(path='/bm/gdata/scollis/berrimah/20060122_074003/')
	read_rays.plot_ppi_vel(ber_0740[1], radar_loc=ber_loc, fig_name='ber_vppi.png')
	read_rays.plot_ppi_vel(gp_0740[1], radar_loc=gp_loc, fig_name='gp_vppi.png')
	xar=linspace(-50.,50., 200)*1000.0
	yar=linspace(-50.,50., 200)*1000.0
	cap_gp=radar_to_cart.make_cappi(gp_0740, xar, yar, 2500.0, 'VR')
	cap_ber=radar_to_cart.make_cappi(ber_0740, xar-disp[0], yar-disp[1], 2500.0, 'VR')
	read_rays.plot_cappi_vel(xar, yar, cap_ber, fig_name='ber_vcappi_2_5km.png')
	read_rays.plot_cappi_vel(xar, yar, cap_gp, fig_name='gp_vcappi_2_5km.png')

def test_newpres():
	gp_0740=read_rays.construct_lassen_scan(path='/bm/gdata/scollis/gunn_pt/20060122_074001/')
	ber_0740=read_rays.construct_uf_scan(path='/bm/gdata/scollis/berrimah/20060122_074003/')	
	ber_loc=[-12.457, 130.925]
	gp_loc=	 [-12.2492,  131.0444]
	ldict={'lat_0':gp_loc[0], 'lon_0':gp_loc[1],'llcrnrlat':-13.0, 'llcrnrlon':130.2, 'urcrnrlat':-12.0 , 'urcrnrlon':131.2, 'lat_ts':gp_loc[0]}
	az_scan=0
	pres.plot_ppi(ber_0740[az_scan], 'CZ', radar_loc=ber_loc, loc_dict=ldict, fig_name='ber_cz_ppi.png')
	pres.plot_ppi(gp_0740[az_scan], 'CZ', radar_loc=gp_loc, loc_dict=ldict, fig_name='gp_cz_ppi.png')
	pres.plot_ppi(ber_0740[az_scan], 'VR', radar_loc=ber_loc, loc_dict=ldict, fig_name='ber_vr_ppi.png')
	pres.plot_ppi(gp_0740[az_scan], 'VR', radar_loc=gp_loc, loc_dict=ldict, fig_name='gp_vr_ppi.png')
	disp=mathematics.corner_to_point(gp_loc, ber_loc)
	xar=linspace(-50.,50., 100)*1000.0
	yar=linspace(-50.,50., 100)*1000.0
	lev=1000.0
	lstr="%(lev)05d" %{'lev':lev}
	pref_dir='20062201_0740_caps/'
	cap_gp_vr=radar_to_cart.make_cappi(gp_0740, xar, yar, lev, 'VR')
	cap_ber_vr=radar_to_cart.make_cappi(ber_0740, xar-disp[0], yar-disp[1], lev, 'VR')
	pres.plot_cappi(xar,yar,cap_gp_vr,gp_0740[0], parm='VR', fig_name=pref_dir+'gp_cappi_vr_'+lstr+'.png', loc_dict=ldict, radar_loc=gp_loc)
	pres.plot_cappi(xar,yar,cap_ber_vr,ber_0740[0], parm='VR', fig_name=pref_dir+'ber_cappi_vr_'+lstr+'.png', loc_dict=ldict, radar_loc=gp_loc)
	#cap_gp_test=radar_to_cart.make_cappi_testmode(gp_0740, xar, yar, lev, 'VR')
	#plot_cappi(xar,yar,cap_gp_test,ber_0740[0], parm='TEST', fig_name='test.png', loc_dict=ldict, radar_loc=gp_loc) 
	
	gp_cube_vr=zeros([100,100,31], dtype=float)
	ber_cube_vr=zeros([100,100,31], dtype=float)
	levs=linspace(500, 10500, 31)
	xar=linspace(-50.,50., 100)*1000.0
	yar=linspace(-50.,50., 100)*1000.0
	
	for i in range(31):
		gp_cap_vr=radar_to_cart.make_cappi(gp_0740, xar, yar, levs[i], 'CZ')
		gp_cube_vr[:,:,i]=gp_cap_vr
		ber_cap_vr=radar_to_cart.make_cappi(ber_0740, xar-disp[0], yar-disp[1], levs[i], 'CZ')
		ber_cube_vr[:,:,i]=ber_cap_vr
	
	for i in range(31):
		lstr="%(lev)05d" %{'lev':levs[i]}
		pres.plot_cappi(xar,yar,gp_cube_vr[:,:,i],gp_0740[0], parm='CZ', fig_name=pref_dir+'gp_cappi_cz_'+lstr+'.png', loc_dict=ldict, radar_loc=gp_loc)
		pres.plot_cappi(xar,yar,ber_cube_vr[:,:,i],ber_0740[0], parm='CZ', fig_name=pref_dir+'ber_cappi_cz_'+lstr+'.png', loc_dict=ldict, radar_loc=gp_loc)
		print 'Done', i, ' of 31'
	

def save_cube_test(date,gpnum, bernum):
	gp_0740=read_rays.construct_lassen_scan(path='/bm/gdata/scollis/gunn_pt/'+date+gpnum+'/')
	ber_0740=read_rays.construct_uf_scan(path='/bm/gdata/scollis/berrimah/'+date+'_'+bernum+'/')	
	ber_loc=[-12.457, 130.925]
	gp_loc=	 [-12.2492,  131.0444]
	displace=mathematics.corner_to_point(gp_loc, ber_loc)
	ldict={'lat_0':gp_loc[0], 'lon_0':gp_loc[1],'llcrnrlat':-13.0, 'llcrnrlon':130.2, 'urcrnrlat':-12.0 , 'urcrnrlon':131.2, 'lat_ts':gp_loc[0]}
	levs=linspace(500,15000, 21)
	xar=linspace(-50.,50., 100)*1000.0
	yar=linspace(-50.,50., 100)*1000.0
	gp_cube=radar_to_cart.make_cube(gp_0740, xar, yar, levs)
	ber_cube=radar_to_cart.make_cube(ber_0740, xar-displace[0], yar-displace[1], levs)
	netcdf_utis.save_data_cube(ber_cube, gp_cube, '/bm/gdata/scollis/cube_data/'+date+'_'+bernum[0:4]+'_ver1.nc', gp_loc)


def simple_reconstruction(tim, lvl_str):
	#load data
	srm=array([15.0, 5.0])/sqrt(2.0)
	ber, gp=netcdf_utis.load_cube('/bm/gdata/scollis/cube_data/20060122_'+tim+'_ver1.nc')
	lvl=int(lvl_str)
	print gp['levs'][lvl]
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(gp['zero_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=gp['zero_loc'][1]+360.0*gp['xar']/(rad_at_radar*2.0*pi)
	lats=gp['zero_loc'][0] + 360.0*gp['yar']/(Re*2.0*pi)	
	ber_loc=[-12.457, 130.925]
	gp_loc=	 [-12.2492,  131.0444]
	
	angs=array(propigation.make_lobe_grid(ber_loc, gp_loc, lats,lons))
	wts_ang=zeros(angs.shape, dtype=float)
	for i in range(angs.shape[0]):
		for j in range(angs.shape[1]):
			if (angs[i,j] < 150.0) and (angs[i,j] > 30.0): wts_ang[i,j]=1.0
	
	
	#create a weighting grid
	mask_reflect=10.0#dBZ	
	mask=(gp['CZ'][:,:,lvl]/mask_reflect).round().clip(min=0., max=1.0) 
	mask_vel_ber=(ber['VR'][:,:,lvl]+100.).clip(min=0., max=1.)
	#run gracon
	print 'Into fortran'
	nx,ny=ber['CZ'][:,:,lvl].shape
	f=0.0
	gv_u=zeros(ber['CZ'][:,:,lvl].shape, dtype=float)
	gv_v=zeros(ber['CZ'][:,:,lvl].shape, dtype=float)
	igu=ones(ber['CZ'][:,:,lvl].shape, dtype=float)*srm[0]
	igv=ones(ber['CZ'][:,:,lvl].shape, dtype=float)*srm[1]
	gv_u,gv_v,f,u_array,v_array = gracon_vel2d.gracon_vel2d(gv_u,gv_v,f,igu,igv,ber['i_comp'][:,:,lvl],ber['j_comp'][:,:,lvl],gp['i_comp'][:,:,lvl],gp['j_comp'][:,:,lvl], ber['VR'][:,:,lvl],gp['VR'][:,:,lvl],mask*mask_vel_ber*wts_ang, nx=nx,ny=ny)
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(gp['zero_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=gp['zero_loc'][1]+360.0*gp['xar']/(rad_at_radar*2.0*pi)
	lats=gp['zero_loc'][0] + 360.0*gp['yar']/(Re*2.0*pi)	
	wts=mask*mask_vel_ber*wts_ang
	f=figure()
	mapobj=pres.generate_darwin_plot(box=[130.8, 131.2, -12.4, -12.0])
	pres.reconstruction_plot(mapobj, lats, lons, gp, lvl, 'CZ',u_array,v_array, angs, wts)
	#pres.quiver_contour_winds(mapobj, lats, lons, (wts*u_array).clip(min=-50, max=50),(wts*v_array).clip(min=-50, max=50))
	t1='Gunn Point CAPPI (dBZ) and reconstructed winds (m/s) at %(lev)05dm \n 22/01/06 ' %{'lev':gp['levs'][lvl]}
	title(t1+tim) 
	ff=os.getenv('HOME')+'/bom_mds/output/recons_22012006/real_%(lev)05d_' %{'lev':gp['levs'][lvl]}
	savefig(ff+tim+'_2d.png')
	close(f)	

def simple_reco(ber,gp, lvl):
	#load data
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(gp['zero_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=gp['zero_loc'][1]+360.0*gp['xar']/(rad_at_radar*2.0*pi)
	lats=gp['zero_loc'][0] + 360.0*gp['yar']/(Re*2.0*pi)	
	ber_loc=[-12.457, 130.925]
	gp_loc=	 [-12.2492,  131.0444]
	
	angs=array(propigation.make_lobe_grid(ber_loc, gp_loc, lats,lons))
	wts_ang=zeros(angs.shape, dtype=float)
	for i in range(angs.shape[0]):
		for j in range(angs.shape[1]):
			if (angs[i,j] < 150.0) and (angs[i,j] > 30.0): wts_ang[i,j]=1.0
	
	
	#create a weighting grid
	mask_reflect=10.0#dBZ	
	mask=(gp['CZ'][:,:,lvl]/mask_reflect).round().clip(min=0., max=1.0) 
	mask_vel_ber=(ber['VR'][:,:,lvl]+100.).clip(min=0., max=1.)
	#run gracon
	print 'Into fortran'
	nx,ny=ber['CZ'][:,:,lvl].shape
	f=0.0
	gv_u=zeros(ber['CZ'][:,:,lvl].shape, dtype=float)
	gv_v=zeros(ber['CZ'][:,:,lvl].shape, dtype=float)
	igu=ones(ber['CZ'][:,:,lvl].shape, dtype=float)*0.0
	igv=ones(ber['CZ'][:,:,lvl].shape, dtype=float)*0.0
	gv_u,gv_v,f,u_array,v_array = gracon_vel2d.gracon_vel2d(gv_u,gv_v,f,igu,igv,ber['i_comp'][:,:,lvl],ber['j_comp'][:,:,lvl],gp['i_comp'][:,:,lvl],gp['j_comp'][:,:,lvl], ber['VR'][:,:,lvl],gp['VR'][:,:,lvl],mask*mask_vel_ber*wts_ang, nx=nx,ny=ny)
	Re=6371.0*1000.0
	print gracon_vel2d.vel_2d_cost(gv_u*0.0,gv_v*0.0,0.0,u_array,v_array,ber['i_comp'][:,:,lvl], ber['j_comp'][:,:,lvl], gp['i_comp'][:,:,lvl], gp['j_comp'][:,:,lvl],  ber['VR'][:,:,lvl], gp['VR'][:,:,lvl],mask*mask_vel_ber*wts_ang)[2]
	return u_array, v_array, mask*mask_vel_ber*wts_ang
	
def simple_reconstruction_3d(tim, lvl_str):
	#load data
	u=0
	l=10
	ui=zeros([l],dtype=float)
	vi=zeros([l],dtype=float)
	#srm=0.0*array([1.0, 5.0])/sqrt(2.0)
	ber, gp=netcdf_utis.load_cube('/bm/gdata/scollis/cube_data/20060122_'+tim+'_ver1.nc')
	for i in range(l):
		ui[i], vi[i]=simple_reco(ber,gp,i)
	lvl=int(lvl_str)
	print gp['levs'][lvl]
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(gp['zero_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=gp['zero_loc'][1]+360.0*gp['xar']/(rad_at_radar*2.0*pi)
	lats=gp['zero_loc'][0] + 360.0*gp['yar']/(Re*2.0*pi)	
	ber_loc=[-12.457, 130.925]
	gp_loc=	 [-12.2492,  131.0444]
	
	angs=array(propigation.make_lobe_grid(ber_loc, gp_loc, lats,lons))
	wts_ang=zeros(gp['CZ'][:,:,u:l].shape, dtype=float)
	for i in range(wts_ang.shape[0]):
		for j in range(wts_ang.shape[1]):
			for k in range(wts_ang.shape[2]):
				if (angs[i,j] < 150.0) and (angs[i,j] > 30.0): wts_ang[i,j,k]=1.0
	
	
	#create a weighting grid
	mask_reflect=12.0#dBZ	
	mask=(ber['CZ'][:,:,u:l]/mask_reflect).round().clip(min=0., max=1.0) 
	mask_vel_ber=(ber['VR'][:,:,u:l]+100.).clip(min=0., max=1.)
	#run gracon
	print 'Into fortran'
	nx,ny,nz=ber['CZ'][:,:,u:l].shape
	print nx,ny,nz
	f=0.0
	gv_u=zeros(ber['CZ'][:,:,u:l].shape, dtype=float)
	gv_v=zeros(ber['CZ'][:,:,u:l].shape, dtype=float)
	igu=ones(ber['CZ'][:,:,u:l].shape, dtype=float)
	igv=ones(ber['CZ'][:,:,u:l].shape, dtype=float)
	for i in range(len(ui)):
		igu[:,:,i]=igu[:,:,i]*ui[i]
		igv[:,:,i]=igu[:,:,i]*vi[i]
	
	
	wts=mask*mask_vel_ber*wts_ang
	#gv_u,gv_v,f,u_array,v_array = gracon_vel2d_3d(gv_u,gv_v,f,u_array,v_array,i_cmpt_r1,j_cmpt_r1,i_cmpt_r2,j_cmpt_r2,vr1,vr2,weights,nx=shape(gv_u,0),ny=shape(gv_u,1),nz=shape(gv_u,2))
	gv_u,gv_v,f,u_array,v_array = gracon_vel2d_3d.gracon_vel2d_3d( gv_u, gv_v, f, igu, igv, ber['i_comp'][:,:,u:l], ber['j_comp'][:,:,u:l], gp['i_comp'][:,:,u:l], gp['j_comp'][:,:,u:l], ber['VR'][:,:,u:l], gp['VR'][:,:,u:l], mywts)#, nx=nx, ny=ny, nz=nz)
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(gp['zero_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=gp['zero_loc'][1]+360.0*gp['xar']/(rad_at_radar*2.0*pi)
	lats=gp['zero_loc'][0] + 360.0*gp['yar']/(Re*2.0*pi)	
	f=figure()
	mapobj=pres.generate_darwin_plot(box=[130.8, 131.2, -12.4, -12.0])
	diff=gp['VR'][:,:,u:l]-(u_array*gp['i_comp'][:,:,u:l]+ v_array*gp['j_comp'][:,:,u:l])
	gp.update({'diff':diff})
	pres.reconstruction_plot(mapobj, lats, lons, gp, lvl, 'diff',u_array[:,:,lvl],v_array[:,:,lvl], angs, wts[:,:,lvl])
	#pres.quiver_contour_winds(mapobj, lats, lons, (wts*u_array).clip(min=-50, max=50),(wts*v_array).clip(min=-50, max=50))
	t1='Gunn Point CAPPI (dBZ) and reconstructed winds (m/s) at %(lev)05dm \n 22/01/06 ' %{'lev':gp['levs'][lvl]}
	title(t1+tim) 
	ff=os.getenv('HOME')+'/bom_mds/output/recons_22012006/real_%(lev)05d_' %{'lev':gp['levs'][lvl]}
	savefig(ff+tim+'_2d_3d.png')
	close(f)	

def simple_reconstruction_3d_pytest(tim, lvl_str, use_guess):
	lvl=int(lvl_str)
	ber, gp=netcdf_utis.load_cube('/bm/gdata/scollis/cube_data/20060122_'+tim+'_ver_hr_big.nc')
	print gp['levs'][lvl]
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(gp['zero_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=gp['zero_loc'][1]+360.0*gp['xar']/(rad_at_radar*2.0*pi)
	lats=gp['zero_loc'][0] + 360.0*gp['yar']/(Re*2.0*pi)	
	ber_loc=[-12.457, 130.925]
	gp_loc=	 [-12.2492,  131.0444]
	if use_guess=='none':
		igu=ones(ber['CZ'].shape, dtype=float)*0.0
		igv=ones(ber['CZ'].shape, dtype=float)*0.0
	else:
		ber_ig, gp_ig=netcdf_utis.load_cube(use_guess)
		print gp_ig.keys()
		igu=gp_ig['u_array']
		igv=gp_ig['v_array']
	mywts=ones(ber['CZ'].shape, dtype=float)
	angs=array(propigation.make_lobe_grid(ber_loc, gp_loc, lats,lons))
	wts_ang=zeros(gp['CZ'][:,:,0].shape, dtype=float)
	for i in range(angs.shape[0]):
			for j in range(angs.shape[1]):
				if (angs[i,j] < 150.0) and (angs[i,j] > 30.0): wts_ang[i,j]=1.0
	for lvl_num in range(len(gp['levs'])):
		#create a weighting grid
		mask_reflect=10.0#dBZ	
		mask=(gp['CZ'][:,:,lvl_num]/mask_reflect).round().clip(min=0., max=1.0) 
		mask_vel_ber=(ber['VR'][:,:,lvl_num]+100.).clip(min=0., max=1.)
		mywts[:,:,lvl_num]=mask*mask_vel_ber*wts_ang	
	f=0.0
	gv_u=zeros(ber['CZ'].shape, dtype=float)
	gv_v=zeros(ber['CZ'].shape, dtype=float)
	wts=mask*mask_vel_ber*wts_ang
	gu,gv,f= grad_conj_solver_plus_plus.meas_cost(gv_u, gv_v, f, igu, igv, ber['i_comp'], ber['j_comp'], gp['i_comp'], gp['j_comp'],  ber['VR'], gp['VR'], mywts)
	print "Mean U gradient", gu.mean(), "gv mean", gv.mean(), "F ", f
	for i in range(len(gp['levs'])):
		print "U,V ", (igu[:,:,i]).sum()/mywts[:,:,i].sum(), (igv[:,:,i]).sum()/mywts[:,:,i].sum()
		#gv_u,gv_v,cost = vel_2d_cost(gv_u,gv_v,cost,u_array,v_array,i_cmpt_r1,j_cmpt_r1,i_cmpt_r2,j_cmpt_r2,vr1,vr2,weights,nx=shape(gv_u,0),ny=shape(gv_u,1))
      		print gracon_vel2d.vel_2d_cost(gv_u[:,:,i]*0.0,gv_v[:,:,i]*0.0,0.0,igu[:,:,i],igv[:,:,i],ber['i_comp'][:,:,i], ber['j_comp'][:,:,i], gp['i_comp'][:,:,i], gp['j_comp'][:,:,i],  ber['VR'][:,:,i], gp['VR'][:,:,i],mywts[:,:,i])[2]
	gv_u,gv_v,f,u_array,v_array = grad_conj_solver_plus_plus.gracon_vel2d_3d( gv_u, gv_v, f, igu, igv, ber['i_comp'], ber['j_comp'], gp['i_comp'], gp['j_comp'], ber['VR'], gp['VR'], mywts)#, nx=nx, ny=ny, nz=nz)
	gp.update({'u_array': u_array, 'v_array':v_array})
	netcdf_utis.save_data_cube(ber, gp, '/bm/gdata/scollis/cube_data/20060122_'+tim+'_winds_ver1.nc', gp['zero_loc'])
	plotit=True
	if plotit:
		for lvl in range(len(gp['levs'])):
			print lvl
			f=figure()
			mapobj=pres.generate_darwin_plot(box=[130.8, 131.2, -12.4, -12.0])
			diff=gp['VR']-(u_array*gp['i_comp']+ v_array*gp['j_comp'])
			gp.update({'diff':diff})
			pres.reconstruction_plot(mapobj, lats, lons, gp, lvl, 'diff',u_array[:,:,lvl],v_array[:,:,lvl], angs, mywts[:,:,lvl])
			#pres.quiver_contour_winds(mapobj, lats, lons, (wts*u_array).clip(min=-50, max=50),(wts*v_array).clip(min=-50, max=50))
			t1='Gunn Point CAPPI (dBZ) and reconstructed winds (m/s) at %(lev)05dm \n 22/01/06 ' %{'lev':gp['levs'][lvl]}
			title(t1+tim) 
			ff=os.getenv('HOME')+'/bom_mds/output/recons_22012006/real_%(lev)05d_' %{'lev':gp['levs'][lvl]}
			savefig(ff+tim+'_2d_3d.png')
			close(f)	

def std_datestr(date_obj, datetype):
	if datetype=="lassen":
		mydatestr=std_lass_datestr(date_obj)
	elif datetype=="uf":
		mydatestr=std_uf_datestr(date_obj)
	return mydatestr

def std_uf_datestr(date_obj):
	date_dict={'y':date_obj.year,'m':date_obj.month,'d':date_obj.day,'HH':date_obj.hour, 'MM':date_obj.minute}
	return "%(y)04d%(m)02d%(d)02d_%(HH)02d%(MM)02d" %date_dict

def std_lass_datestr(date_obj):
	date_dict={'y':date_obj.year,'m':date_obj.month,'d':date_obj.day,'HH':date_obj.hour, 'MM':date_obj.minute}
	return "%(y)04d%(m)02d%(d)02d%(HH)02d%(MM)02d" %date_dict

def recon(date_str, latstr, lonstr):
     use_guess='sonde'
     sonde_file='/bm/gdata/scollis/twpice/darwin.txt'
     #tim='1350'
     tim_date=num2date(datestr2num(date_str)) 
     ber, gp=netcdf_utis.load_cube('/bm/gdata/scollis/cube_data/'+std_datestr(tim_date)+'_deal.nc')
     sonde_list=read_sounding.read_sounding_within_a_day(sonde_file, tim_date)
     #launch_dates=[sonde['date_list'][0] for sonde in sonde_list]
     #launch_date_offset=[date2num(sonde['date_list'][0])- date2num(tim_date)  for sonde in sonde_list]
     #best_sonde=sonde_list[argsort(abs(array(launch_date_offset)))[0]]
     #print 'Time of radar: ', tim_date, ' Time of sonde_launch: ', best_sonde['date_list'][0], ' Time of sonde_termination: ', best_sonde['date_list'][-1]
     req=[ 'alt(m)',  'wspd(m/s)', 'tdry(degs)',  'wdir(degs)']	
     first_sonde, second_sonde=read_sounding.get_two_best_conc_sondes(date_str, req_vars=req)
     interp_sonde=read_sounding.interp_sonde_time(first_sonde, second_sonde, tim_date, gp['levs'])
     u_sonde=ones(gp['CZ'].shape, dtype=float)
     v_sonde=ones(gp['CZ'].shape, dtype=float)
     w_sonde=zeros(gp['CZ'].shape, dtype=float)
     for k in range(len(gp['levs'])):
	     u_sonde[:,:,k]=-1.0*u_sonde[:,:,k]*interp_sonde['wspd(m/s)'][k]*sin(pi*interp_sonde['wdir(degs)'][k]/180.0)
	     v_sonde[:,:,k]=-1.0*v_sonde[:,:,k]*interp_sonde['wspd(m/s)'][k]*cos(pi*interp_sonde['wdir(degs)'][k]/180.0)
     Re=6371.0*1000.0
     rad_at_radar=Re*sin(pi/2.0 -abs(gp['zero_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
     lons=gp['zero_loc'][1]+360.0*gp['xar']/(rad_at_radar*2.0*pi)
     lats=gp['zero_loc'][0] + 360.0*gp['yar']/(Re*2.0*pi)	
     ber_loc=[-12.457, 130.925]
     gp_loc= [-12.2492,  131.0444]
     
     if use_guess=='none':
	  igu=ones(ber['CZ'].shape, dtype=float)*0.0
	  igv=ones(ber['CZ'].shape, dtype=float)*0.0
	  igw=ones(ber['CZ'].shape, dtype=float)*0.0
     elif use_guess=='sonde':
	  igu=u_sonde
	  igv=v_sonde
	  igw=w_sonde
     else:
	  ber_ig, gp_ig=netcdf_utis.load_cube(use_guess)
	  print gp_ig.keys()
	  igu=gp_ig['u_array']
	  igv=gp_ig['v_array']
	  igw=ones(ber['CZ'].shape, dtype=float)*0.0
     
     #mywts=ones(ber['CZ'].shape, dtype=float)
     angs=array(propigation.make_lobe_grid(ber_loc, gp_loc, lats,lons))
     #wts_ang=zeros(gp['CZ'][:,:,0].shape, dtype=float)
     #for i in range(angs.shape[0]):
     #    for j in range(angs.shape[1]):
     #		 if (angs[i,j] < 150.0) and (angs[i,j] > 30.0): wts_ang[i,j]=1.0
     #for lvl_num in range(len(gp['levs'])):
     #	  #create a weighting grid
     #	  mask_reflect=10.0#dBZ	
     #	  mask=(gp['CZ'][:,:,lvl_num]/mask_reflect).round().clip(min=0., max=1.0) 
     #	  mask_vel_ber=(ber['VR'][:,:,lvl_num]+100.).clip(min=0., max=1.)
     #	  mywts[:,:,lvl_num]=mask*mask_vel_ber*wts_ang	
     mywts=met.make_mask(ber, gp, angs, 1.0, 80.0)
     print "Mean gp masked Velocity ", (gp['VE']*mywts).mean()
     print "min gp masked Velocity ", (gp['VE']*mywts).min()
     print "max gp masked Velocity ", (gp['VE']*mywts).max()
     print "Mean Berrimah masked Velocity ", (ber['VE']*mywts).mean()
     print "min Berrimah masked Velocity ", (ber['VE']*mywts).min()
     print "max Berrimah masked Velocity ", (ber['VE']*mywts).max()
     print "Mean gp masked CZ ", (gp['CZ']*mywts).mean()
     print "min gp masked CZ ", (gp['CZ']*mywts).min()
     print "max gp masked CZ ", (gp['CZ']*mywts).max()
     print "Mean Berrimah masked CZ ", (ber['CZ']*mywts).mean()
     print "min Berrimah masked CZ ", (ber['CZ']*mywts).min()
     print "max Berrimah masked CZ ", (ber['CZ']*mywts).max()
    
     print "Number of masked points", (mywts.shape[0]*mywts.shape[1]*mywts.shape[2])-mywts.sum()
     print "Number of unmasked points ", mywts.sum()
     print "**********************FALLSPEED INFO****************************"
     #def terminal_velocity(refl, temps, levs, display=False):
     tdry=interp_sonde['tdry(degs)']
     pressure=interp_sonde['press(hPa)']
     dummy=met.terminal_velocity(gp['CZ']*mywts, tdry, gp['levs'], display=True)
     print "**********************FALLSPEED INFO****************************"
     #print 
     
     
     f=0.0
     X=[igu,igv,igw]
     G,F,X=grad_conj_solver_3d.gracon_3d_packaged(X ,ber, gp, mywts, interp_sonde)
     u_array,v_array,w_array=X 
     lvl=2500.0
     lvl_num=argsort(abs(gp['levs']-lvl))[0]
     print "Level_num=", lvl_num
     if latstr=='max' or lonstr=='max':
	     maskedcz=gp['CZ'][:,:,lvl_num]*mywts[:,:,lvl_num]
	     i,j=mathematics.where_closest_2d(maskedcz.max(), maskedcz)
	     print i,j
	     lat=lats[j[0]]
	     lon=lons[i[0]]
	     print "Max CZ at ", lat, lon 
     else:
	     lon=float(lonstr)
	     lat=float(latstr)
     f=figure()
     #(lat_sl, lon_sl, lvl, data_cube, lats, lons, levs, u, v, w, angs, mask, par='CZ', w_mag=2.0, **kwargs)
     alat, alon, alvl=pres.plot_slices(lat, lon, lvl, gp, lats, lons, gp['levs'], u_array, v_array, w_array, angs, mywts, par='CZ', w_mag=2.0,box=[130.5, 131.5, -12.7, -12.0], bquiver=[0.05, 0.75], ksp=0.05)
     t1='Gunn Point reflectivity (dBZ) and reconstructed winds (m/s, *2.0 for w)\n sliced at %(alat)2.2fS and %(alon)3.2fE and %(alvl)d Metres on 22/01/06 at ' %{'alat':abs(alat), 'alon':alon, 'alvl':alvl}
     t2=" %(HH)02d%(MM)02dZ" %{'HH':tim_date.hour, 'MM':tim_date.minute}
     f.text( .1, .92, t1+t2) 
     inte_part=1000*(float(int(lat))-lat)
     print  {'alat':abs(alat), 'alon':alon, 'alvl':alvl}
     ff=os.getenv('HOME')+'/bom_mds/output/recons_'+std_datestr(tim_date)[0:-5]+'/slicer3_%(alat)2.02f_%(alon)3.02f_%(alvl)05d_' %{'alat':abs(alat), 'alon':alon, 'alvl':alvl}
     print ff
     savefig(ff+t2+'.png', dpi=200)
     gp.update({'u_array':u_array, 'v_array':v_array, 'w_array':w_array})
     netcdf_utis.save_data_cube(ber, gp, '/bm/gdata/scollis/cube_data/'+std_datestr(tim_date)+'_winds.nc', gp_loc)
     close(f)


def dump_pickle(pickle_name):
	locs={'Berrimah':[-12.457, 130.925], 'C-POL':[-12.2492,  131.0444], 'Berrimah_deal':[-12.457, 130.925], 'C-POL_deal':[-12.2492,  131.0444]}
	basedir='/bm/gdata/scollis/radar_shelves/'
	radar=pickle_zip.load(basedir+pickle_name)
	pres.dump_radar(radar, locs[radar[0]['radar_name']])


def test_sounding():
	fname='/bm/gdata/scollis/twpice/Darwin22Jan.txt'
	my_list=read_sounding.read_sounding(fname)
	pres.plot_sonde(my_list[0], 'Darwin22Jan0526.png')


def radar_to_winds(datestr):
	#check to see if we have the radar files
	#check to see if there are deailased files
	





if __name__ == "__main__":
	t0=systime()
	print "the uber cool test"
	print sys.argv
	#save_cube_test(sys.argv[1], sys.argv[2], sys.argv[3])
	#simple_reconstruction_3d_pytest(sys.argv[1],sys.argv[2], sys.argv[3])
	recon(sys.argv[1], sys.argv[2], sys.argv[3])
	#test_pert_winds()
	#test_sounding()
	#dump_pickle(sys.argv[1])
	print "Finished running runtime=",systime()-t0, "Seconds"

