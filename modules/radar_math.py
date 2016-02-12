#Some radar propigation math
import sys
sys.path.append('/flurry/home/scollis/pylibs/lib64/python2.4/site-packages')
from Numeric import array as Nar
from numpy import *
sys.path.append('/flurry/home/scollis/python/radar/')
#from Interpolator3d import *
from time import *


def radar_coords_to_cart(rng, az, ele, debug=False):
	"""
	Asumes standard atmosphere, ie R=4Re/3
	"""
	Re=6371.0*1000.0
	#h=(r^2 + (4Re/3)^2 + 2r(4Re/3)sin(ele))^1/2 -4Re/3
	#s=4Re/3arcsin(rcos(ele)/(4Re/3+h))
	p_r=4.0*Re/3.0
	rm=rng*1000.0
	z=(rm**2 + p_r**2 + 2.0*rm*p_r*sin(ele*pi/180.0))**0.5 -p_r
	#arc length
	s=p_r*arcsin(rm*cos(ele*pi/180.)/(p_r+z))
	if debug: print "Z=", z, "s=", s
	y=s*cos(az*pi/180.0)
	x=s*sin(az*pi/180.0)
	return x,y,z



#def int_single(x, y, sweep_dict):
	



def int_radar(xar,yar,zar,par,x,y,z,npts=20):
	if sqrt((x/1000.0)**2+(y/1000.0)**2) < 150.0:
		#t0=time()
		#Determine the n closest points
		#print "Calculating distance ", time()-t0
		dar=sqrt((xar-x)**2+(yar-y)**2+(xar-z)**2)
		#print "Sorting", time()-t0
		argsorted=dar.argsort()
		#cxs=xar[argsorted[0:npts]]
		#cys=yar[argsorted[0:npts]]
		#czs=zar[argsorted[0:npts]]
		ds=dar[argsorted[0:4]]
		cps=par[argsorted[0:4]]
		#print "interpolating", time()-t0
		wts=ds[3]-ds[0:3]
		clsmean=(cps[0:3]*wts).mean()/wts.mean()
		#f = Interpolator3d(Nar(cxs), Nar(cys), Nar(czs),Nar(cps), .8)
		#dat=f.interpolate([x],[y],[z])
		#print "done", time()-t0
	else:
		clsmean=0.0
	return clsmean


def int_cone(scan,parm,x,y):
	if sqrt((x/1000.0)**2+(y/1000.0)**2) < 150.0:
		xr,yr,zr,pr=scan_list_to_ldata([scan], parm)
		xar=array(xr)
		yar=array(yr)
		zar=array(zr)
		par=array(pr)
		#t0=time()
		#Determine the n closest points
		#print "Calculating distance ", time()-t0
		dar=sqrt((xar-x)**2+(yar-y)**2)
		#print "Sorting", time()-t0
		argsorted=dar.argsort()
		#cxs=xar[argsorted[0:npts]]
		#cys=yar[argsorted[0:npts]]
		#czs=zar[argsorted[0:npts]]
		xs=xar[argsorted[0:10]]
		ys=yar[argsorted[0:10]]
		cps=par[argsorted[0:10]]
		#print "interpolating", time()-t0
		f=interpolate.interp2d(xs,ys,cps)
		clsmean=f(x,y)
		#f = Interpolator3d(Nar(cxs), Nar(cys), Nar(czs),Nar(cps), .8)
		#dat=f.interpolate([x],[y],[z])
		#print "done", time()-t0
		del f
	else:
		clsmean=0.0
	return clsmean



def make_cappi(xa,ya,za,pa,xs, ys, z, npts=2):
	xara=array(xa)
	yara=array(ya)
	zara=array(za)
	para=array(pa)
	ca=zeros([len(xs), len(ys)], dtype=float)
	for i in range(ca.shape[0]):
		print "doing row ", i, " of ", ca.shape[0]
		t=time()
		for j in range(ca.shape[1]):
			ca[i,j]=int_radar(xara,yara,zara,para, xs[i], ys[j] , z , npts=2)
		print "Row done, time taken ", time()-t, " Seconds"
	print ca.shape
	return ca





#f = Interpolator3d(Nar(X.ravel()), Nar(Y.ravel()), Nar(Z.ravel()),Nar(d.ravel()), .8)
