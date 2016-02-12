import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
import pyradar
from numpy import zeros, float, linspace, array
from time import time as systime
from pylab import num2date, datestr2num
from read_rays import append_xyz, append_ijk
import mathematics
import time_utils

def screen_var_on_RH(scan, var, lev=0.4, baddata=1.31072000e+05):
	for i in range(scan[var].shape[0]):
		for j in range(scan[var].shape[1]):
			if scan['RH'][i,j] < lev:
				scan[var][i,j]=baddata
	return scan

def make_kdp(scanlist):
	new_scanlist=[]
	for scan in scanlist:
		kdp_array=mathematics.dx2d(scan['PH'])/mathematics.dx2d(scan['range'])#zeros(scan['PH'], dtype=float)
		scan.update({'KDP':kdp_array})
		screened_scan=screen_var_on_RH(scan, 'KDP')
		new_scanlist.append(screened_scan)
	return new_scanlist


def copy_DZCZ(scanlist):
	if not('CZ' in scanlist[0].keys()):
		for i in range(len(scanlist)):
			scanlist[i].update({'CZ':scanlist[i]['DZ']})
	return scanlist


def read_volume(radar_obj, **kwargs):
	if "avail_index" in kwargs.keys():
		avail_index=kwargs['avail_index']
	else:
		radar_index={'DZ':0,'VR':1,'SW':2, 'CZ':3,'ZT':4,'DR':5,'LR':6,'ZD':7,'DM':8, 'RH':9,'PH':10,'XZ':11,'CR':12,'MZ':13,'MR':14, 'ZE':15,'VE':16,'KD':17,'TI':18,'DX':19,'CH':20,'AH':21,'CV':22, 'AV':23,'SQ':24}
		avail_index={}
		for index in radar_index.keys():
			if pyradar.volume_exists(radar_obj, radar_index[index]): avail_index.update({index:radar_index[index]})
	nscans=pyradar.get_nsweeps(pyradar.get_volume(radar_obj,avail_index[avail_index.keys()[0]]))
	scanlist=[]
	for i in range(nscans):
		print "Reading scan ",i+1, " of ", nscans
		scanlist.append(append_ijk(append_xyz(read_scan(radar_obj, i, avail_index=avail_index))))
	return scanlist

def load_radar(filename, **kwargs):
	psth={}
	if 'avail_index' in kwargs.keys(): psth.append({'avail_index':kwargs['avail_index']})
	radar_obj=pyradar.get_radar(filename)
	scan_list=read_volume(radar_obj, **psth)
	pyradar.free_radar(radar_obj)
	return scan_list



def read_scan(radar_obj, sweep_num, **kwargs):
	if "avail_index" in kwargs.keys():
		avail_index=kwargs['avail_index']
	else:
		radar_index={'DZ':0,'VR':1,'SW':2, 'CZ':3,'ZT':4,'DR':5,'LR':6,'ZD':7,'DM':8, 'RH':9,'PH':10,'XZ':11,'CR':12,'MZ':13,'MR':14, 'ZE':15,'VE':16,'KD':17,'TI':18,'DX':19,'CH':20,'AH':21,'CV':22, 'AV':23,'SQ':24}
		avail_index={}
		for index in radar_index.keys():
			if pyradar.volume_exists(radar_obj, radar_index[index]): avail_index.update({index:radar_index[index]})
	
	test_index=avail_index[avail_index.keys()[0]]
	myvol=pyradar.get_volume(radar_obj, test_index)
	mysweep=pyradar.get_sweep(myvol, sweep_num)
	myray=pyradar.get_ray(mysweep, 0)
	ngates=pyradar.get_ngates(myray)
	nrays=pyradar.get_nrays(mysweep)
	rng=linspace(pyradar.get_first_gate(myray), pyradar.get_first_gate(myray)+pyradar.get_gate_size(myray)*(pyradar.get_ngates(myray)-1.0), pyradar.get_ngates(myray))/1000.0
	date_dict={'y':pyradar.get_yr(myray),'m':pyradar.get_month(myray),'d':pyradar.get_day(myray),'HH':pyradar.get_hour(myray), 'MM':pyradar.get_min(myray)}
	my_str="%(y)04d%(m)02d%(d)02d %(HH)02d%(MM)02d" %date_dict
	print ngates
	print nrays
	my_dict={}
	
	for index in avail_index.keys():
		darr=zeros([nrays,ngates], dtype=float)
		rarr=zeros([nrays,ngates], dtype=float)
		elevs=zeros([nrays], dtype=float)
		azimuths=zeros([nrays], dtype=float)
		myvol=pyradar.get_volume(radar_obj, avail_index[index])
		mysweep=pyradar.get_sweep(myvol, sweep_num)
		for i in range(nrays):
			elevs[i]=pyradar.get_elev(myray)
			azimuths[i]=pyradar.get_azimuth(myray)
			rarr[i,:]=rng
			for j in range(ngates):
				myray=pyradar.get_ray(mysweep, i)
				darr[i,j]=pyradar.get_value(myray, j)
		my_dict.update({index:darr})
	my_dict.update({'Azmth':azimuths})
	my_dict.update({'range':rarr})
	my_dict.update({'Elev':elevs})
	my_dict.update({'date':num2date(datestr2num(my_str))})
	my_dict.update({'radar_name':pyradar.get_radarname(radar_obj)})
	my_dict.update({'radar_loc':[pyradar.get_lat(radar_obj),pyradar.get_lon(radar_obj)]})
	return my_dict
	



def flip_rhi(rhi_scan, loud=False):
	gated_vars=list(set(rhi_scan[0].keys()) & set(['range', 'xar', 'yar', 'zar', 'VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW', 'KD', 'i_comp', 'j_comp', 'k_comp', 'u_array', 'v_array', 'w_array']) )
	range_indep_vars=['Azmth', 'Elev']
	n_azimuth=len(rhi_scan)
	n_gates=rhi_scan[0][gated_vars[0]].shape[1]
	n_elev=rhi_scan[0][gated_vars[0]].shape[0]
	PPI_scan=[]
	for elev in range(n_elev):
		PPI_scan.append({})
		if loud: print "Compositing Elevation ", elev+1, " Of ", n_elev 
		for item in range_indep_vars:
			blanker=zeros([n_azimuth], dtype=float)
			for az_num in range(n_azimuth):
				blanker[az_num]=rhi_scan[az_num][item][elev]
			PPI_scan[elev].update({item:blanker})
		for item in gated_vars:
			blanker=zeros([n_azimuth, n_gates], dtype=float)
			for az_num in range(n_azimuth):
				for g_num in range(n_gates):
					blanker[az_num, g_num]=rhi_scan[az_num][item][elev, g_num]
					PPI_scan[elev].update({item:blanker})
	times=time_utils.date_linspace(rhi_scan[0]['date'], rhi_scan[-1]['date'], len(PPI_scan))
	for i in range(len(PPI_scan)): 
		PPI_scan[i].update({'date':times[i]})
		PPI_scan[i].update({'radar_loc':rhi_scan[0]['radar_loc']})
		PPI_scan[i].update({'radar_name':rhi_scan[0]['radar_name']})
	return PPI_scan



