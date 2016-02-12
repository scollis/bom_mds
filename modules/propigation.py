###################################################################
#propigation.py: Untilities to determine propigation direction of #
#rays through a particilar CAPPI grid point                       #
###################################################################
#Called by:                                                       #
###################################################################
#Start of BoM Python branch: Scott Collis, CAWCR, April 2008      #
###################################################################
from numpy import sin, cos, arcsin, pi, sqrt, double, arctan, linspace, zeros, array, arccos
import sys
from os import getenv
sys.path.append(getenv('HOME')+'/bom_mds/modules')
from mathematics import where_fits, corner_to_point, lalo_meshgrid
from pylab import meshgrid

def ray(range_ar, ele):
	"""
	Asumes standard atmosphere, ie R=4Re/3
	"""
	Re=6371.0*1000.0
	#h=(r^2 + (4Re/3)^2 + 2r(4Re/3)sin(ele))^1/2 -4Re/3
	#s=4Re/3arcsin(rcos(ele)/(4Re/3+h))
	p_r=4.0*Re/3.0
	z=(range_ar**2 + p_r**2 + 2.0*range_ar*p_r*sin(ele*pi/180.0))**0.5 -p_r
	#arc length
	s=p_r*arcsin(range_ar*cos(ele*pi/180.)/(p_r+z))
	return s, z


def elev_of_scan(s,z):
	"""
	Given a arc distance along the earth, s, and a height
	(for example the elevation of a CAPPI scan) return
	the elevation of the radar dish to be at that height at that
	arc distance.
	
	Assumes standard atmosphere, ie R=4Re/3
	
	elev_of_scan(s,z)
	returns: Elevation angle in degrees
	"""
	Re=6371.0*1000.0
	Rr=double(4.0*Re/3.0)
	d=sqrt(2.0*(Rr**2)*(1.0-cos(s/Rr)))
	theta_e=arctan(z/d) -arcsin(d/(2.0*Rr))
	te_degs=theta_e*180.0/pi
	return te_degs

def dhds(x,y,h,debug=False, **kwargs):
	"""
	Given an x and y displacement from the radar in Metres and and height above the earth return the gradient of the beam through that layer
	"""
	max_range=kwargs.get('max_range', 150.0*1000.0)
	d_range=kwargs.get('d_range', 500.0)
	s=sqrt(x**2+y**2)
	ele_ang=elev_of_scan(s,h)
	s_ray, h_ray=ray(linspace(0.0, max_range, int(max_range/d_range)), ele_ang)
	#This catch could be replaced by a try for better speed
	if s_ray.max() < s: s_ray, h_ray=ray(linspace(0.0, 1.5*max_range, int(1.5*max_range/d_range)), ele_ang)
	lh,rh=where_fits(s, s_ray)
	ds=s_ray[rh]-s_ray[lh]
	dh=h_ray[rh]-h_ray[lh]
	if debug: 
		print "elevation angle of scan=", ele_ang
		print ds
		print dh
		print rh
		print lh 
	return dh/ds

def dhds_lalo(point, h, radar_loc):
	x,y=corner_to_point(radar_loc, point)
	grad=dhds(x,y,h)
	return grad

def unit_vector(x,y,h,debug=False, **kwargs):
	"""
	Calculate the unit vector of the ray through a layer
	"""
	max_range=kwargs.get('max_range', 250.0*1000.0)
	d_range=kwargs.get('d_range', 200.0)
	s=sqrt(x**2+y**2)
	ele_ang=elev_of_scan(s,h)
	s_ray, h_ray=ray(linspace(0.0, max_range, int(max_range/d_range)), ele_ang)
	#This catch could be replaced by a try for better speed
	if s_ray.max() < s: s_ray, h_ray=ray(linspace(0.0, 1.5*max_range, int(1.5*max_range/d_range)), ele_ang)
	lh,rh=where_fits(s, s_ray)
	ds=s_ray[rh]-s_ray[lh]
	dh=h_ray[rh]-h_ray[lh]
	angle=arctan(dh/ds)
	k_comp=sin(angle)
	r_comp=cos(angle)
	i_comp=r_comp*x/sqrt(x**2+y**2)
	j_comp=r_comp*y/sqrt(x**2+y**2)
	return i_comp, j_comp, k_comp

def unit_vector_lalo(point, h, radar_loc):
	x,y=corner_to_point(radar_loc, point)
	i,j,k=unit_vector(x,y,h)
	return i,j,k

def unit_vector_grid(lats, lons, h, radar_loc):
	longr, latgr=meshgrid(lons, lats)
	i_a=zeros(longr.shape, dtype=float)
	j_a=zeros(longr.shape, dtype=float)
	k_a=zeros(longr.shape, dtype=float)
	h=5.0*1000.0
	for i in range(longr.shape[0]):
		for j in range(longr.shape[1]):
			ic, jc, kc=unit_vector_lalo((latgr[i,j],longr[i,j]), h, radar_loc)
			#print ic,jc,kc, sqrt(ic**2+jc**2+kc**2)
			i_a[i,j]=ic
			j_a[i,j]=jc
			k_a[i,j]=kc
	return i_a, j_a, k_a


def line_angle(r1, p, r2, debug=False):
	""" 
	line_angle(r1, p, r2)
	Work out the angle (in degrees) between two lines using the cosine rule
	args: 
	r1, p, r2: 2 element tuples with the x,y coords of the three points
	"""
	pa=array(p)
	r1a=array(r1)
	r2a=array(r2)
	r1p=sqrt(((pa-r1a)**2).sum())
	r2p=sqrt(((pa-r2a)**2).sum())
	r1r2=sqrt(((r1a-r2a)**2).sum())
	#if debug: 
	#	print 'r1p %(r1p)f r2p %(r2p)f r1r2 %(r1r2)f' %{'r1p':r1p, 'r2p':r2p, 'r1r2':r1r2}
	#	print 'arg', (r1p**2 + r2p**2 - r1r2**2)/(2.0*r1p*r2p)
	theta_rads=arccos((r1p**2 + r2p**2 - r1r2**2)/(2.0*r1p*r2p))
	return theta_rads*180.0/pi


def make_lobe_grid(dop_rad1, dop_rad2, lats,lons):
	x,y=lalo_meshgrid(lats, lons)
	angs=zeros(x.shape)
	dist1=zeros(x.shape)
	dist2=zeros(x.shape)
	dop1_cart=corner_to_point([lats[0],lons[0]], dop_rad1)
	dop2_cart=corner_to_point([lats[0],lons[0]], dop_rad2)
	#print dop1_cart
	#print dop2_cart
	for i in range(len(lats)):
		for j in range(len(lons)):
			#print [x[j], y[i]]
			dist1[i,j]=sqrt(((array([x[i,j], y[i,j]])-array(dop1_cart))**2).sum())
			dist2[i,j]=sqrt(((array([x[i,j], y[i,j]])-array(dop2_cart))**2).sum())
			angs[i,j]=line_angle(dop1_cart, [x[i,j], y[i,j]], dop2_cart)
	return angs

	






