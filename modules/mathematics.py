#Mathematical proceedures required by various programs within the text generator 
#Scott Collis, National Meteorological and Oceanographic Centre
#Australian Bureau of Meteorology
#s.collis@bom.gov.au
#Last modified 12/6/08
import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
from pylab import meshgrid
from numpy import zeros, append, array, sqrt, arctan2, abs, pi, sin, linspace
from numpy import where as nwhere
import dif_ops
import numpy

def smooth(x,window_len=10,window='hanning'):
	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."
	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
	if window_len<3:
		return x
	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
	s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
	#print(len(s))
	if window == 'flat': #moving average
		w=ones(window_len,'d')
	else:
		w=eval('numpy.'+window+'(window_len)')
	y=numpy.convolve(w/w.sum(),s,mode='same')
	return y[window_len-1:-window_len+1]


def where_closest_2d(value, arr):
	i=linspace(0,arr.shape[0], arr.shape[0]+1)
	j=linspace(0,arr.shape[1], arr.shape[1]+1)
	ii,jj=meshgrid(i,j)
	cks=abs(arr-value)
	tfarray= cks==cks.min()
	iwhere=ii[tfarray]
	jwhere=jj[tfarray]
	return iwhere, jwhere


def where_fits(val, arr):
	"""
	finds the elements between which val fits
	needs a value and a lineraly increasing array
	args:
	val: A value
	arr: linearly increasing numpy array
	returns:
	(lh_element, rh_element) both integers 
	"""
	try:
		#Concise bit of clever code
		#which means it is difficult to explain
		#basically the LHS (arr <= val) makes
		#an array of bools which is  true
		# to the LHS of where Val fits in arr
		#nwhere (taken carefully from numpy as a different name)
		#returns the indexes of those trues and [0][-1] (shape 0,len(arr) for some unknown reason)
		#indexes the RHS-most element that is true giving the index 
		#to the left of where val fits in arr...
		#the second calc in the tupple does the same but mirrored
		mytup=(nwhere((arr <= val))[0][-1],nwhere((arr >= val))[0][0])
	except IndexError:
		#if the code above throws an IndexError the val fits at an end of arr
		#determine which end
		if val < arr.min(): mytup=(0,0)
		if val > arr.max(): mytup=(len(arr)-1, len(arr)-1)
	return mytup 



def w_mean(subject, weights):
	#compute the weighted mean of subject
	#shp=subject.shape
	integerand=subject*weights
	res=integerand.sum()/weights.sum()
	return res

#u,v to mag, direction
def uvtomd(u,v):
	#Calc magnitude
	#needs vector_direction
	mag=sqrt(u**2+v**2)
	direct=vector_direction(u,v)
	return mag, direct

def vector_direction(u,v):
	#A somewhat sophisticated way of generating a direction grid from
	#u and v grids which avoids the use of ifs so it can be run over an
	#array
	udash=u+0.000001
	flag= abs((abs(udash)-udash)/(2.0*udash))
	val=abs(arctan2(u,v))*180.0/pi
	return val+flag*(360.0-2.0*val)




def ax_radius(lat, units='radians'):
	pi=3.145
	#Determine the radius of a circle of constant longitude at a certain 
	#Latitude
	Re=6371.0*1000.0
	if units=='degrees':
		const=pi/180.0
	else:
		const=1.0
	
	R=Re*sin(pi/2.0 - abs(lat*const))
	return R



def lalo_meshgrid(lats, lons):
	#Given a 1D array of latitudes and longitudes 
	#generate an X and Y meshgrid
	#X[i,j],Y[i,j] are the X distance and Y distance (in metres) between
	#(lat[0],lon[0]) and (lat[i], lon[j])
	#The confusing scheme of accessing lat as the first parameter
	#is in keeping with the US netcdf formats and formats accepted 
	#by matplotlib.basemap (ie X.shape=(lat.shape[0],lon.shape[0]))
	pi=3.145
	x,y=meshgrid(lons,lats)
	Re=6371.0*1000.0
	#print x.shape
	#print lats.shape[0]
	#print lons.shape[0]
	for i in range(lats.shape[0]):
		Rc=ax_radius(lats[i], units='degrees')
		#print lats[i]
		#print Rc
		#print Rc/Re
		for j in range(lons.shape[0]):
			y[i,j]=((lats[i]-lats[0])/360.0)*pi*2.0*Re
			x[i,j]=((lons[j]-lons[0])/360.0)*pi*2.0*Rc
	return x,y
 


def corner_to_point(corner, point):
	pi=3.145
	Re=6371.0*1000.0
	Rc=ax_radius(point[0], units='degrees')
	#print Rc/Re
	y=((point[0]-corner[0])/360.0)*pi*2.0*Re
	x=((point[1]-corner[1])/360.0)*pi*2.0*Rc
	return x,y




def grad_inter(x,y):
	#given two 1d arrays have a linear relationship 
	#y=mx+b determine m and b
	m=(y[0]-y[1])/(x[0]-x[1])
	b=y[0]-m*x[0]
	return m,b


def extrap(x,y,nx):
	#A linear extrapolation from two points along nx
	#x: two x elements
	#y: two y elements
	m,b=dif_ops.grad_inter(x,y)
	my_extrap=m*nx+b
	return my_extrap


def dy(y):
	#Calculate dy
	#this procedure is designed to return an array of the same length as y
	#to do this it we needed to extend y by an element in each direction
	nelms=len(y)
	#There is possibly a better way to do this using zeros instead of appending...
	#seems to run fast enough...
	mdy=[]
	pdy=[]
	
	#calculate the extrapolations
	m1th=dif_ops.extrap([0,1], [y[0],y[1]], -1)
	np1th=dif_ops.extrap([nelms-2,nelms-1],[y[nelms-2],y[nelms-1]], nelms)
	
	#append them to the data
	s1=append([m1th], y)
	yext=append(s1,[np1th])
	
	#calculate the difference using N-1
	for i in range(nelms):
		j=i+2
		dy=yext[j]-yext[j-1]
		pdy.append(dy)
		
		
	#Calculate the difference using N+1
	for i in range(nelms):
		j=i+1
		dy=yext[j]-yext[j-1]
		mdy.append(dy)
		
	#Use the averate of N+1, N-1 to return dy
	return (array(mdy)+array(pdy))/2.0


def dx2d(z):
	#Run dy over a grid in the x direction
	nx=z.shape[0]
	ny=z.shape[1]
	u=zeros([nx,ny],dtype=float)
	for row in range(nx):
		thisrow=dif_ops.dy(z[row,:])
		u[row,:]=thisrow
	
	return u
	

def dy2d(z):
	#run dy over a grid in the y direction
	nx=z.shape[0]
	ny=z.shape[1]
	v=zeros([nx,ny],dtype=float)
	for col in range(ny):
		thisrow=dif_ops.dy(z[:,col])
		v[:,col]=thisrow
	
	return v
	

def curl_zcmp(u,v,x,y):
	#calculate the k component of the curl of the vector field in a x-y plane
	#curl_k=|delXv|=dv/dx-du/dy 
	t1=dx2d(v)/dx2d(x)
	t2=dy2d(u)/dy2d(y)
	return t1-t2




def grad2d(z,xg,yg):
	#Calculate the gradient vector of a scalar field
	#grad2d=(du/dx,dv/dy)
	u=dx2d(z)/dx2d(xg)
	v=dy2d(z)/dy2d(yg)
	return u,v

