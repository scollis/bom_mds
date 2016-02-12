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
import read_rays
import radar_math
import radar_to_cart
import mathematics
import netcdf_utis
from matplotlib import figure
from pylab import *
from time import time as systime
from numpy import linspace, array, arctan, pi, random


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

if __name__ == "__main__":
	t0=systime()
	print "the uber cool test"
	print sys.argv
	#save_cube_test(sys.argv[1], sys.argv[2], sys.argv[3])
	simple_reconstruction(sys.argv[1],sys.argv[2])
	#test_pert_winds()
	print "Finished running runtime=",systime()-t0, "Seconds"

