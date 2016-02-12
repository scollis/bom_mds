#This module contains a variety of methods for working on grids of data
#Scott Collis NMOC/CAWCR
#April 2008

import os
def generate_paths():
	if os.name=='posix':
		base='/home/scott/python/text_formatter/'
		path_dic={'base':base, 'ncep':base+'ncep/', 'masks':base+'masks/',
		'sci':'/usr/lib/python2.5/site-packages/Scientific/IO/',
		'dll':'/usr/lib/python2.5/site-packages/Scientific/IO/',
		'modules':base+'modules/', 'shapes':base+'shapes/', 'products':base+'products/'}
	else:
		#assume windows
		base='D:\python_work\\'
		path_dic={'base':base, 'ncep':base+'ncep\\', 'masks':base+'masks\\',
		'sci':'C:\Python25\Lib\site-packages\Scientific\IO',
		'dll':'C:\Python25\Lib\site-packages\Scientific\win32',
		'modules':base+'modules\\', 'shapes':base+'shapes\\', 'products':base+'products\\'}
	
	return path_dic


paths=generate_paths()
from numpy import ndarray, abs, linspace, int, array, zeros, sqrt
from numpy import where as nwhere
from scipy import interpolate

def where_closest(value, arr):
	"""
	where_closest(value, arr)
	Given Value determine the closest index in the array arr
	Args:
	value: Value to find
	arr: a numpy array
	returns: an array of integers where it could fit 
	(unless there are two or more spots this will be of length 1)
	"""
	differ=abs(arr-value)
	cond=(differ==min(differ)) #creates an array of booleans
	index=linspace(0,len(arr), len(arr)+1) #simply an array of indexes
	mwhere=index[cond] #numpy arrays have a nice feature that accessing with an array of booleans returns an array where the boolean array is true
	return mwhere

def get_metarea(uwind10m, vwind10m, lat, lon, window_ar, debug=False):
	"""
	get_metarea(uwind10m, vwind10m, lat, lon, window_ar, debug=False)
	extract a subgrid of data from a master grid
	Args:
	lat, lon: indexing arrays
	uwind10m, vwind10m: two 2d grids of shape (len(lat),len(lon))
	window_ar:[min_lon, max_lon, min_lat, max_lat] the area of grids to be extracted
	Kwargs
	debug=False: set to true for verbose reporting 
	
	returns: (grid1, grid2, new_lats, new_lons) all numpy array objects (inherited from the input arrays)
	
	Split a grid of wind data  
	window_arr=[left lon, right lon, bottom lat, top lat]
	lats are N lons are E (eg -40 150)
	"""
	
	#calculate the indexes
	lonsindex=int(where_closest(window_ar[0], array(lon)))
	loneindex=int(where_closest(window_ar[1], array(lon)))+1
	latsindex=int(where_closest(window_ar[2], array(lat)))
	lateindex=int(where_closest(window_ar[3], array(lat)))+1
	if debug:
		print latsindex
		print lateindex
		print lonsindex
		print loneindex
		print uwind10m.shape
		print vwind10m.shape
	
	#populate the output arrays
	#print '1'
	u_metarea=uwind10m[latsindex:lateindex, lonsindex:loneindex]
	v_metarea=vwind10m[latsindex:lateindex, lonsindex:loneindex]
	#print '2'
	nlat=lat[latsindex:lateindex]
	nlon=lon[lonsindex:loneindex]
	return u_metarea, v_metarea, nlat, nlon


#def regrid(old_lat, old_lon, old_grid, new_lat, new_lon):
#	latlat,lonlon=meshgrid(old_lat,old_lon)
#	latr, lonr=map(ravel, [latlat, lonlon])
#	interp_obj=interpolate.interp2d(latr, lonr, old_grid)
#	new_grid=interp_obj(new_lat, new_lon)
#	return new_grid
#hasattr(arg, 'shape')


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



def calc_nearest_indecies(lats, lons, lat, lon):
	latsi=where_fits(lat, lats)
	lonsi=where_fits(lon, lons)
	return (latsi[0], lonsi[0]), (latsi[1], lonsi[0]), (latsi[1], lonsi[1]), (latsi[0], lonsi[1])




def lin_int_single(old_lat, old_lon, old_grid, new_lat, new_lon):
	#First check to see if a point exists at the exact point we want
	#This saves having to do the distance calculation
	md=sqrt((old_lat[1]-old_lat[0])**2 + (old_lon[1]-old_lon[0])**2)
	if (new_lat in old_lat) and (new_lon in old_lon):
		value=old_grid[old_lat.searchsorted(new_lat),old_lon.searchsorted(new_lon)]
	else:
		#if there is not the exact point then calculate the distances
		points=calc_nearest_indecies(old_lat, old_lon, new_lat, new_lon)
		#print points
		distances=array([sqrt( (new_lat-old_lat[i])**2 + (new_lon - old_lon[j])**2) for (i,j) in points])
		#print distances
		weights=md-distances
		#then calculate the weighted mean to be returned
		values=array([old_grid[point] for point in points])
		#print values
		value=(values*weights).mean()/weights.mean()
	return value





def regrid(old_lat, old_lon, old_grid, new_lat, new_lon):
	"""
	A simple 2D regridding application that does a linear interpolation
	between the four closest points
	Args:
	old_lat, old_lon: referencing arrays for old_grid
	old_grid: original data grid
	new_lat, new_lon: lon and lat arrays to reference the new grids
	"""
	#generate the blank grid to be filled
	new_grid=zeros([len(new_lat), len(new_lon)])
	#check to see that the inputs are arrays
	if not(hasattr(old_lat, 'mean')):old_lat=array(old_lat)
	if not(hasattr(old_lon, 'mean')):old_lon=array(old_lon)
	#loop over the new_lats and new_lons
	for i in range(len(new_lat)):
		for j in range(len(new_lon)):
			new_grid[i,j]=lin_int_single(old_lat, old_lon, old_grid, new_lat[i], new_lon[j])
		
	return new_grid















