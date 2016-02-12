###################################################################
#simul_winds.py: Modules for generating wind fields for testing   #
###################################################################
#Module for bom_mds                                               #
###################################################################
#Start of BoM Python branch: Scott Collis, CAWCR, April 2008      #
###################################################################
from os import getenv
import sys
sys.path.append(getenv('HOME')+'/bom_mds/modules')
import mathematics
import propigation
from numpy import ones, sin, cos, pi, sqrt
from pylab import meshgrid, exp, log

def unif_wind(lats, lons, mag, drn):
	ucmp=mag*sin(drn*pi/180.0)*ones([len(lats),len(lons)], dtype=float)
	vcmp=mag*cos(drn*pi/180.0)*ones([len(lats),len(lons)], dtype=float)
	return ucmp, vcmp

def speed_bump(lats, lons, locn, fwhm):
	if len(lats.shape)==1:
		longr, latgr=meshgrid(lons, lats)
	else:
		longr=lons
		latgr=lats
	dist=sqrt((longr-locn[1])**2+(latgr-locn[0])**2)
	k=(2.0*fwhm**2)/(4.0*log(2.0))
	bump=exp((-1.0*dist**2)/k)
	return bump

def vortex(lats, lons, locn, fwhm):
	if len(lats.shape)==1:
		longr, latgr=meshgrid(lons, lats)
	else:
		longr=lons
		latgr=lats
	dlat=latgr-locn[0]
	dlon=longr-locn[1]
	u=speed_bump(lats,lons,locn,fwhm)*dlat/sqrt(dlat**2+dlon**2)
	v=-1.0*speed_bump(lats,lons,locn,fwhm)*dlon/sqrt(dlat**2+dlon**2)
	return u,v
