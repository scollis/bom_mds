#Module to process a ray.
import numpy
import mathematics
import parse_ini
import os
__author__ = "Scott Collis"
__version__ = "1.0"


def process_scan(scan, **kwargs):
	ini_fname=kwargs.get('ini_fname', os.getenv('HOME')+'/bom_mds/bom_mds.ini')
	ini_dict=parse_ini.parse_ini(ini_fname)
	gated_vars=['VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW','KD']
	filtered_scan={}
	for item in set(scan.keys())&set(gated_vars):
		if ini_dict[item+'_bins'] !=0:
			rec_array=numpy.zeros(scan[item].shape, dtype=float)
			for az_num in range(scan[item].shape[0]):
				rec_array[az_num, :]=filter_ray(scan[item][az_num,:], scan['RH'][az_num,:], npts=ini_dict[item+'_bins'], ftype=ini_dict[item+'_window'], RH_thresh=ini_dict['RH_thresh'])
			filtered_scan.update({item:rec_array})
		else:
			filtered_scan.update({item:scan[item]})
	rec_array=numpy.zeros(scan[item].shape, dtype=float)
	for az_num in range(scan['PH'].shape[0]):
		rec_array[az_num, :]=generate_kdp(scan['PH'][az_num,:], scan['range'][az_num,:] ,scan['RH'][az_num,:] ,npts=ini_dict['KD_bins'], ftype=ini_dict['KD_window'], RH_thresh=ini_dict['RH_thresh'])
	filtered_scan.update({'KD':rec_array})
	return filtered_scan


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
		w=numpy.ones(window_len,'d')
	else:
		w=eval('numpy.'+window+'(window_len)')
	y=numpy.convolve(w/w.sum(),s,mode='same')
	return y[window_len-1:-window_len+1]


def filter_ray(data,  RH_data,npts=5, bad_data=1.31072000e+05, ftype='hanning', RH_thresh=0.5):
	#Scan ray, find subsections of the ray which are hydro echos using RH
	dest_ray=numpy.zeros([len(data)], dtype=float)+bad_data
	curr_pos=0
	starts=[]
	ends=[]
	while curr_pos < len(data):
		#see if we are in a good data echo
		if (data[curr_pos]!=bad_data) and (RH_data[curr_pos]>RH_thresh):
			starts.append(curr_pos)
			sp=curr_pos
			this_block=[]
			while (data[curr_pos]!=bad_data)  and (RH_data[curr_pos]>RH_thresh):
				#now we are scanning in a region of meteo echos
				this_block.append(data[curr_pos])
				curr_pos=curr_pos+1
				if curr_pos >= len(data): break
			ep=curr_pos 
			ends.append(ep)
			print 'Len of list:', len(this_block)
			#Do stuff with this bit of data
			if len(this_block) >npts:
				this_array=smooth(numpy.array(this_block), window_len=npts, window=ftype)
			else:
				this_array=numpy.array(this_block)
			print this_array.shape
			print sp, ep, ep-sp
			dest_ray[sp:ep]=this_array
		else:
			#we are in a non-meteo region
			curr_pos=curr_pos+1
	return dest_ray
		
def generate_kdp(data, r,  RH_data,npts=5, bad_data=1.31072000e+05, ftype='hanning', RH_thresh=0.5):
	#Scan ray, find subsections of the ray which are hydro echos using RH
	dest_ray=numpy.zeros([len(data)], dtype=float)+bad_data
	curr_pos=0
	starts=[]
	ends=[]
	while curr_pos < len(data):
		#see if we are in a good data echo
		if (data[curr_pos]!=bad_data) and (RH_data[curr_pos]>RH_thresh):
			starts.append(curr_pos)
			sp=curr_pos
			this_block=[]
			this_r=[]
			while (data[curr_pos]!=bad_data)  and (RH_data[curr_pos]>RH_thresh):
				#now we are scanning in a region of meteo echos
				this_block.append(data[curr_pos])
				this_r.append(r[curr_pos])
				curr_pos=curr_pos+1
				if curr_pos >= len(data): break
			ep=curr_pos 
			ends.append(ep)
			print 'Len of list:', len(this_block)
			#Do stuff with this bit of data
			if len(this_block) >npts:
				this_array=smooth(mathematics.dy(numpy.array(this_block))/mathematics.dy(numpy.array(this_r)), window_len=npts, window=ftype)
			elif len(this_block) > 2:
				this_array=mathematics.dy(numpy.array(this_block))/mathematics.dy(numpy.array(this_r))
			if len(this_block)>2:
				print this_array.shape
				print sp, ep, ep-sp
				dest_ray[sp:ep]=this_array
		else:
			#we are in a non-meteo region
			curr_pos=curr_pos+1
	return dest_ray
