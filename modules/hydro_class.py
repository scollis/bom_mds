import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
from pylab import num2date, date2num, datestr2num
import parse_ini


def hydro_std_datetime(date_obj, **kwargs):
	path=kwargs.get('path', '/bm/gscratch/scollis/ascii_hydro/')
	prefix=kwargs.get('prefix', 'cpol_hydroclass_')
	postfix=kwargs.get('postfix', '.ascii')
	fdict={'y':date_obj.year,'m':date_obj.month,'d':date_obj.day,'HH':date_obj.hour, 'MM':date_obj.minute}
	dtstr="%(y)04d%(m)02d%(d)02d_%(HH)02d%(MM)02d" %fdict
	return path+prefix+dtstr+postfix



def read_class(date_obj, **kwargs):
	
	hydro_filename='/data/cpol_hyrdo_class/'+dt+'_'+tm+'.ascii'
	hydro_file=open(hydro_filename, 'r')
	hydro_ascii=hydro_file.readlines()
	hydro_file.close()
	hydro_dateobj=num2date(datestr2num(hydro_ascii[0]))
	data_list=[float(item) for item in hydro_ascii[1].split()]
	data_dict={'zero_loc':[data_list[0], data_list[1]], 'xrng':[data_list[2], data_list[3]], 'yrng':[data_list[5], data_list[6]], 'zrng':[data_list[8], data_list[9]], 'res':[data_list[4], data_list[7], data_list[10]]}
	
	float_list=[]
	for item in hydro_ascii[2:len(hydro_ascii)]:
		for number in item.split():
			float_list.append(float(number))
	
	classifyT=zeros(data_dict['res'], dtype='float')
	reflT=zeros(data_dict['res'], dtype='float')
	#i+nx*(j-1)+nx*ny*(k-1)
	classify=zeros(data_dict['res'], dtype='float')
	refl=zeros(data_dict['res'], dtype='float')
	
	for i in range(data_dict['res'][0]):
		for j in range(data_dict['res'][1]):
			for k in range(data_dict['res'][2]):
				classifyT[i,j,k]=float_list[int(i+j*data_dict['res'][0]+data_dict['res'][1]*data_dict['res'][0]*k)*2 +1]
				reflT[i,j,k]=float_list[int(i+j*data_dict['res'][0]+data_dict['res'][1]*data_dict['res'][0]*k)*2]
	
	for k in range(data_dict['res'][2]):
		classify[:,:,k]=transpose(classifyT[:,:,k])
		refl[:,:,k]=transpose(reflT[:,:,k])
	
	xar=linspace(data_dict['xrng'][0], data_dict['xrng'][1],  data_dict['res'][0])*1000.0
	yar=linspace(data_dict['yrng'][0], data_dict['yrng'][1],  data_dict['res'][1])*1000.0
	zar=linspace(data_dict['zrng'][0], data_dict['zrng'][1],  data_dict['res'][2])*1000.0
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(data_dict['zero_loc'][0]*pi/180.0))#ax_radius(float(lat_cpol), units='degrees')
	lons=data_dict['zero_loc'][1]+360.0*xar/(rad_at_radar*2.0*pi)
	lats=data_dict['zero_loc'][0] + 360.0*yar/(Re*2.0*pi)	
		
	data_dict.update({'xar':xar, 'yar':yar, 'zar':zar,'CZ':refl, 'classify':classify, 'lats':lats, 'lons':lons})
	
	
	
	datestr=dt+' '+tm
	kwargs={}
	ini_fname=kwargs.get('ini_fname', os.getenv('HOME')+'/bom_mds/bom_mds.ini')
	ini_dict=parse_ini.parse_ini(ini_fname)
	
	dateobj=num2date(datestr2num(datestr))
	tim_date=num2date(datestr2num(datestr))
	radar1, radar2=netcdf_utis.load_cube('/data/cube_data/'+std_datestr(tim_date)+'_winds.nc')
	Re=6371.0*1000.0
	rad_at_radar=Re*sin(pi/2.0 -abs(radar1['radar_loc'][0]*pi/180.0))
	lons=radar1['radar_loc'][1]+360.0*radar1['xar']/(rad_at_radar*2.0*pi)
	lats=radar1['radar_loc'][0] + 360.0*radar1['yar']/(Re*2.0*pi)	
	ber_loc=[-12.457, 130.925]
	gp_loc=[-12.2492,  131.0444]




	##angs=array(propigation.make_lobe_grid(radar2['radar_loc'], radar1['radar_loc'], lats,lons))
	##mywts=met.make_mask(radar2, radar1, angs, 1.0, 80.0)
	
	#hydro_bins=make_hydro_bins_clat(data_dict, lat)
	
	#f=figure()
	#plot_slice_lon_hydro_sym(lat, radar1, lats, lons, radar1['levs'], hydro_bins, radar1['u_array'], radar1['w_array'],mywts, par='CZ', w_mag=2.0, box=[130.8, 131.1, 10,10])
	#savefig('/flurry/home/scollis/bom_mds/output/class_test_sym'+dt+tm+'.png')
	#close(f)
	
	#f=figure()
	#plot_slice_lon_hydro(lat, radar1, lats, lons, radar1['levs'], data_dict, radar1['u_array'], radar1['w_array'],mywts, par='CZ', w_mag=2.0, box=[130.8, 131.1, 10,10])
	#savefig('/flurry/home/scollis/bom_mds/output/class_test_pcolor'+dt+tm+'.png')
	#close(f)
	

