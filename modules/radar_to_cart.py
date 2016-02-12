###################################################################
#radar_to_cart.py: given a list of scans genrate a constant height#
#XY grid of radar properties                                      # 
###################################################################
#Called by:                                                       #
###################################################################
#Scott Collis, CAWCR, May 2008                                    #
#Last Modified: 22/5/08                                           #
###################################################################
import sys
from os import getenv
sys.path.append(getenv('HOME')+'/bom_mds/modules')
import propigation
import mathematics
from numpy import sqrt, array, argsort, abs, zeros, float, pi
from pylab import meshgrid, num2date, date2num

def cond_lin_int(v1, v2, wghts, baddata=1.31072000e+05):
	ans=baddata
	if not(v1==baddata  or v2 ==baddata):
		ans=v1*wghts[0]+v2*wghts[1]
	elif v1 ==baddata and not(v2 ==baddata):
		ans=v2
	elif v2 ==baddata and not(v1==baddata):
		ans=v1
	return ans


def get_value_at_xyz(scan_list, x,y,z, var, **kwargs):
	"""
	scan_list: a list of scans 
	x,y,z: position to find value (interpolation position) in metres
	var: parameter to interpolate (eg VR, DZ,...)
	"""
	#Generate a array of elevations for the scan list
	#We use the minus onth (last) element in the scan as the radar has settled by then
	elevs=array([scan['Elev'][-1]  for scan in scan_list])
	#pick the above and below elevations
	s=sqrt(x**2+y**2)
	az=mathematics.vector_direction(x,y)
	baddata=kwargs.get('baddata',1.31072000e+05)
	my_elevation=propigation.elev_of_scan(s,z)
	elev_numbers=mathematics.where_fits(my_elevation, elevs)
	#for each elevation pick the ray on either side of the point
	if not(elev_numbers[0] == elev_numbers[1]):
		#we are neither at the top or bottom
		top_scan=scan_list[elev_numbers[1]]
		bottom_scan=scan_list[elev_numbers[0]]
		#Scan above the point ----------TOP ELEVATION-----------
		#work out what gates to use (ie range in S direction)
		sar=sqrt(top_scan['xar'][0,:]**2 + top_scan['yar'][0,:]**2)
		gates=mathematics.where_fits(s, sar)
		ds=sar[gates[1]]-sar[gates[0]]
		if ds==0.0:
			top_gates_comp=[0.5,0.5]
		else:
			top_gates_comp=[ abs(s-sar[gates[1]])/ds, abs(s-sar[gates[0]])/ds] 
		scan_order=top_scan['Azmth'].argsort()
		fits_order=mathematics.where_fits(az, top_scan['Azmth'][scan_order])
		ray_numbers=[scan_order[fits_order[0]], scan_order[fits_order[1]]]
		d_angle=top_scan['Azmth'][ray_numbers[1]]-top_scan['Azmth'][ray_numbers[0]]
		top_ang_comp=[abs(az-top_scan['Azmth'][ray_numbers[1]])/d_angle, abs(az-top_scan['Azmth'][ray_numbers[0]])/d_angle]
		if fits_order[0] == fits_order[1]:
			#we are near 360
			azd=az
			ray_numbers=[scan_order[-1], scan_order[0]]
			d_angle=top_scan['Azmth'][ray_numbers[1]]-(top_scan['Azmth'][ray_numbers[0]]-360.0)
			if az > d_angle: azd=az-360.0
			top_ang_comp=[abs(azd-top_scan['Azmth'][ray_numbers[1]])/d_angle,abs(azd-(top_scan['Azmth'][ray_numbers[0]]-360.0))/d_angle ]
		tp1=cond_lin_int(top_scan[var][ray_numbers[0], gates[0]], top_scan[var][ray_numbers[0], gates[1]], top_gates_comp)
		tp2=cond_lin_int(top_scan[var][ray_numbers[1], gates[0]], top_scan[var][ray_numbers[1], gates[1]], top_gates_comp)
		top_point=cond_lin_int(tp1, tp2, top_ang_comp)
		top_z=(top_scan['zar'][ray_numbers[0], gates[0]]*top_gates_comp[0]+ top_scan['zar'][ray_numbers[0], gates[1]]*top_gates_comp[1])
		#BOTTTTTTTTOOMMMMMMMMMMMM!
		sar=sqrt(bottom_scan['xar'][0,:]**2 + bottom_scan['yar'][0,:]**2)
		gates=mathematics.where_fits(s, sar)
		ds=sar[gates[1]]-sar[gates[0]]
		if ds==0.:
			bottom_gates_comp=[0.5,0.5]
		else:
			bottom_gates_comp=[ abs(s-sar[gates[1]])/ds, abs(s-sar[gates[0]])/ds] 
		scan_order=bottom_scan['Azmth'].argsort()
		fits_order=mathematics.where_fits(az, bottom_scan['Azmth'][scan_order])
		ray_numbers=[scan_order[fits_order[0]], scan_order[fits_order[1]]]
		d_angle=bottom_scan['Azmth'][ray_numbers[1]]-bottom_scan['Azmth'][ray_numbers[0]]
		bottom_ang_comp=[abs(az-bottom_scan['Azmth'][ray_numbers[1]])/d_angle, abs(az-bottom_scan['Azmth'][ray_numbers[0]])/d_angle]
		if fits_order[0] == fits_order[1]:
			#we are near 360
			azd=az
			ray_numbers=[scan_order[-1], scan_order[0]]
			d_angle=bottom_scan['Azmth'][ray_numbers[1]]-(bottom_scan['Azmth'][ray_numbers[0]]-360.0)
			if az > d_angle: azd=az-360.0
			bottom_ang_comp=[abs(azd-bottom_scan['Azmth'][ray_numbers[1]])/d_angle,abs(azd-(bottom_scan['Azmth'][ray_numbers[0]]-360.0))/d_angle ]
		bp1=cond_lin_int(bottom_scan[var][ray_numbers[0], gates[0]], bottom_scan[var][ray_numbers[0], gates[1]], bottom_gates_comp)
		bp2=cond_lin_int(bottom_scan[var][ray_numbers[1], gates[0]], bottom_scan[var][ray_numbers[1], gates[1]], bottom_gates_comp)
		bottom_point=cond_lin_int(bp1, bp2, bottom_ang_comp)
		bottom_z=(bottom_scan['zar'][ray_numbers[0], gates[0]]*bottom_gates_comp[0]+ bottom_scan['zar'][ray_numbers[0], gates[1]]*bottom_gates_comp[1])
		#calculate components
		d_z=top_z-bottom_z
		#print d_z
		z_comps=[abs(z-top_z)/d_z, abs(z-bottom_z)/d_z]
		my_point=cond_lin_int(bottom_point, top_point, z_comps)
		#bottom_point*z_comps[0]+top_point*z_comps[1]
	else:
		my_point=0.0
	if "az_range" in kwargs.keys():
		a1=kwargs['az_range'][0]
		a2=kwargs['az_range'][0]
		if a1 > a2: #we go through 360
			if not((az >= a1) or (az <=a2)): 
				my_point=baddata
		else:
			if not((az >= a1) or (az <= a2)):
				my_point=baddata
	return my_point

def get_values_at_xyz(scan_list, x,y,z,**kwargs):
	"""
	scan_list: a list of scans 
	x,y,z: position to find value (interpolation position) in metres
	var: parameter to interpolate (eg VR, DZ,...)
	"""
	#Generate a array of elevations for the scan list
	#We use the minus onth (last) element in the scan as the radar has settled by then
	avail_vars=set(scan_list[0].keys()) & set(['VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW', 'i_comp', 'j_comp', 'k_comp']) 
	baddata=kwargs.get( 'baddata',1.31072000e+05)
	#tp1={}
	#tp2={}
	#bp1={}
	#bp2={}
	bottom_point={}
	top_point={}
	my_point={}
	elevs=array([scan['Elev'][-1]  for scan in scan_list])
	#pick the above and below elevations
	s=sqrt(x**2+y**2)
	az=mathematics.vector_direction(x,y)
	my_elevation=propigation.elev_of_scan(s,z)
	elev_numbers=mathematics.where_fits(my_elevation, elevs)
	#for each elevation pick the ray on either side of the point
	if not(elev_numbers[0] == elev_numbers[1]):
		#we are neither at the top or bottom
		top_scan=scan_list[elev_numbers[1]]
		bottom_scan=scan_list[elev_numbers[0]]
		#Scan above the point ----------TOP ELEVATION-----------
		#work out what gates to use (ie range in S direction)
		sar=sqrt(top_scan['xar'][0,:]**2 + top_scan['yar'][0,:]**2)
		gates=mathematics.where_fits(s, sar)
		ds=sar[gates[1]]-sar[gates[0]]
		if ds==0.0:
			top_gates_comp=[0.5,0.5]
		else:
			top_gates_comp=[ abs(s-sar[gates[1]])/ds, abs(s-sar[gates[0]])/ds] 
		scan_order=top_scan['Azmth'].argsort()
		fits_order=mathematics.where_fits(az, top_scan['Azmth'][scan_order])
		ray_numbers=[scan_order[fits_order[0]], scan_order[fits_order[1]]]
		d_angle=top_scan['Azmth'][ray_numbers[1]]-top_scan['Azmth'][ray_numbers[0]]
		top_ang_comp=[abs(az-top_scan['Azmth'][ray_numbers[1]])/d_angle, abs(az-top_scan['Azmth'][ray_numbers[0]])/d_angle]
		if fits_order[0] == fits_order[1]:
			#we are near 360
			azd=az
			ray_numbers=[scan_order[-1], scan_order[0]]
			d_angle=top_scan['Azmth'][ray_numbers[1]]-(top_scan['Azmth'][ray_numbers[0]]-360.0)
			if az > d_angle: azd=az-360.0
			top_ang_comp=[abs(azd-top_scan['Azmth'][ray_numbers[1]])/d_angle,abs(azd-(top_scan['Azmth'][ray_numbers[0]]-360.0))/d_angle ]
		for var in avail_vars:
			tp1=cond_lin_int(top_scan[var][ray_numbers[0], gates[0]], top_scan[var][ray_numbers[0], gates[1]], top_gates_comp)
			tp2=cond_lin_int(top_scan[var][ray_numbers[1], gates[0]], top_scan[var][ray_numbers[1], gates[1]], top_gates_comp)
			tp=cond_lin_int(tp1, tp2, top_ang_comp)
			top_point.update({var:tp})
		
		top_z=(top_scan['zar'][ray_numbers[0], gates[0]]*top_gates_comp[0]+ top_scan['zar'][ray_numbers[0], gates[1]]*top_gates_comp[1])
		#BOTTTTTTTTOOMMMMMMMMMMMM!
		sar=sqrt(bottom_scan['xar'][0,:]**2 + bottom_scan['yar'][0,:]**2)
		gates=mathematics.where_fits(s, sar)
		ds=sar[gates[1]]-sar[gates[0]]
		if ds==0.:
			bottom_gates_comp=[0.5,0.5]
		else:
			bottom_gates_comp=[ abs(s-sar[gates[1]])/ds, abs(s-sar[gates[0]])/ds] 
		scan_order=bottom_scan['Azmth'].argsort()
		fits_order=mathematics.where_fits(az, bottom_scan['Azmth'][scan_order])
		ray_numbers=[scan_order[fits_order[0]], scan_order[fits_order[1]]]
		d_angle=bottom_scan['Azmth'][ray_numbers[1]]-bottom_scan['Azmth'][ray_numbers[0]]
		bottom_ang_comp=[abs(az-bottom_scan['Azmth'][ray_numbers[1]])/d_angle, abs(az-bottom_scan['Azmth'][ray_numbers[0]])/d_angle]
		if fits_order[0] == fits_order[1]:
			#we are near 360
			azd=az
			ray_numbers=[scan_order[-1], scan_order[0]]
			d_angle=bottom_scan['Azmth'][ray_numbers[1]]-(bottom_scan['Azmth'][ray_numbers[0]]-360.0)
			if az > d_angle: azd=az-360.0
			bottom_ang_comp=[abs(azd-bottom_scan['Azmth'][ray_numbers[1]])/d_angle,abs(azd-(bottom_scan['Azmth'][ray_numbers[0]]-360.0))/d_angle ]
		for var in avail_vars:
			bp1=cond_lin_int(bottom_scan[var][ray_numbers[0], gates[0]], bottom_scan[var][ray_numbers[0], gates[1]], bottom_gates_comp)
			bp2=cond_lin_int(bottom_scan[var][ray_numbers[1], gates[0]], bottom_scan[var][ray_numbers[1], gates[1]], bottom_gates_comp)
			bp=cond_lin_int(bp1, bp2, bottom_ang_comp)
			bottom_point.update({var:bp})
		
		bottom_z=(bottom_scan['zar'][ray_numbers[0], gates[0]]*bottom_gates_comp[0]+ bottom_scan['zar'][ray_numbers[0], gates[1]]*bottom_gates_comp[1])
		#calculate components
		d_z=top_z-bottom_z
		#print d_z
		z_comps=[abs(z-top_z)/d_z, abs(z-bottom_z)/d_z]
		for var in avail_vars:
			mp=cond_lin_int(bottom_point[var], top_point[var], z_comps)
			my_point.update({var:mp})
		#bottom_point*z_comps[0]+top_point*z_comps[1]
		approx_time=date2num(bottom_scan['date'])*z_comps[0] + date2num(top_scan['date'])*z_comps[1]
		my_point.update({'datenum':approx_time})
	else:
		my_point=dict([(var,baddata) for var in avail_vars])
		my_point.update({'datenum':baddata})
	if "az_range" in kwargs.keys():
		#print az
		a1=kwargs['az_range'][0]
		a2=kwargs['az_range'][1]
		if a1 > a2: #we go through 360
			if not((az >= a1) or (az <=a2)): 
				my_point=dict([(var,baddata) for var in avail_vars])
				my_point.update({'datenum':baddata})
		else:
			if not((az >= a1) and (az <= a2)):
				#print "Bad bad", az
				my_point=dict([(var,baddata) for var in avail_vars])
				my_point.update({'datenum':baddata})
	return my_point


#for i in linspace(700,800,10):
#	print get_value_at_xyz(gp_0740, 25.0*1000.0,25.0*1000.0,i, 'VR')

#print get_value_at_xyz(scanme, 0.0,64.0*1000.0,2500.0, 'CZ')
#print get_value_at_xyz(gp_0740, 25.0*1000.0,25.0*1000.0,900.0, var)
#mycap=make_cappi(scanme, xar, yar, 2500.0, 'CZ') 

def get_value_at_xyz_testmode(scan_list, x,y,z, var):
	"""
	scan_list: a list of scans 
	x,y,z: position to find value (interpolation position) in metres
	var: parameter to interpolate (eg VR, DZ,...)
	"""
	#Generate a array of elevations for the scan list
	#We use the minus onth (last) element in the scan as the radar has settled by then
	elevs=array([scan['Elev'][-1]  for scan in scan_list])
	#pick the above and below elevations
	s=sqrt(x**2+y**2)
	az=mathematics.vector_direction(x,y)
	my_elevation=propigation.elev_of_scan(s,z)
	elev_numbers=mathematics.where_fits(my_elevation, elevs)
	#for each elevation pick the ray on either side of the point
	if not(elev_numbers[0] == elev_numbers[1]):
		#we are neither at the top or bottom
		top_scan=scan_list[elev_numbers[1]]
		bottom_scan=scan_list[elev_numbers[0]]
		#Scan above the point ----------TOP ELEVATION-----------
		#work out what gates to use (ie range in S direction)
		sar=sqrt(top_scan['xar'][0,:]**2 + top_scan['yar'][0,:]**2)
		gates=mathematics.where_fits(s, sar)
		ds=sar[gates[1]]-sar[gates[0]]
		if ds==0.0:
			top_gates_comp=[0.5,0.5]
		else:
			top_gates_comp=[ abs(s-sar[gates[1]])/ds, abs(s-sar[gates[0]])/ds] 
		scan_order=top_scan['Azmth'].argsort()
		fits_order=mathematics.where_fits(az, top_scan['Azmth'][scan_order])
		ray_numbers=[scan_order[fits_order[0]], scan_order[fits_order[1]]]
		d_angle=top_scan['Azmth'][ray_numbers[1]]-top_scan['Azmth'][ray_numbers[0]]
		top_ang_comp=[abs(az-top_scan['Azmth'][ray_numbers[1]])/d_angle, abs(az-top_scan['Azmth'][ray_numbers[0]])/d_angle]
		if fits_order[0] == fits_order[1]:
			#we are near 360
			azd=az
			ray_numbers=[scan_order[-1], scan_order[0]]
			d_angle=top_scan['Azmth'][ray_numbers[1]]-(top_scan['Azmth'][ray_numbers[0]]-360.0)
			if az > d_angle: azd=az-360.0
			top_ang_comp=[abs(azd-top_scan['Azmth'][ray_numbers[1]])/d_angle,abs(azd-(top_scan['Azmth'][ray_numbers[0]]-360.0))/d_angle ]
		tp1=cond_lin_int(elevs[elev_numbers[1]], elevs[elev_numbers[1]], top_gates_comp)
		tp2=cond_lin_int(elevs[elev_numbers[1]], elevs[elev_numbers[1]], top_gates_comp)
		top_point=cond_lin_int(tp1, tp2, top_ang_comp)
		top_z=(top_scan['zar'][ray_numbers[0], gates[0]]*top_gates_comp[0]+ top_scan['zar'][ray_numbers[0], gates[1]]*top_gates_comp[1])
		#BOTTTTTTTTOOMMMMMMMMMMMM!
		sar=sqrt(bottom_scan['xar'][0,:]**2 + bottom_scan['yar'][0,:]**2)
		gates=mathematics.where_fits(s, sar)
		ds=sar[gates[1]]-sar[gates[0]]
		if ds==0.:
			bottom_gates_comp=[0.5,0.5]
		else:
			bottom_gates_comp=[ abs(s-sar[gates[1]])/ds, abs(s-sar[gates[0]])/ds] 
		scan_order=bottom_scan['Azmth'].argsort()
		fits_order=mathematics.where_fits(az, bottom_scan['Azmth'][scan_order])
		ray_numbers=[scan_order[fits_order[0]], scan_order[fits_order[1]]]
		d_angle=bottom_scan['Azmth'][ray_numbers[1]]-bottom_scan['Azmth'][ray_numbers[0]]
		bottom_ang_comp=[abs(az-bottom_scan['Azmth'][ray_numbers[1]])/d_angle, abs(az-bottom_scan['Azmth'][ray_numbers[0]])/d_angle]
		if fits_order[0] == fits_order[1]:
			#we are near 360
			azd=az
			ray_numbers=[scan_order[-1], scan_order[0]]
			d_angle=bottom_scan['Azmth'][ray_numbers[1]]-(bottom_scan['Azmth'][ray_numbers[0]]-360.0)
			if az > d_angle: azd=az-360.0
			bottom_ang_comp=[abs(azd-bottom_scan['Azmth'][ray_numbers[1]])/d_angle,abs(azd-(bottom_scan['Azmth'][ray_numbers[0]]-360.0))/d_angle ]
		bp1=cond_lin_int(elevs[elev_numbers[0]],elevs[elev_numbers[0]], bottom_gates_comp)
		bp2=cond_lin_int(elevs[elev_numbers[0]],elevs[elev_numbers[0]], bottom_gates_comp)
		bottom_point=cond_lin_int(bp1, bp2, bottom_ang_comp)
		bottom_z=(bottom_scan['zar'][ray_numbers[0], gates[0]]*bottom_gates_comp[0]+ bottom_scan['zar'][ray_numbers[0], gates[1]]*bottom_gates_comp[1])
		#calculate components
		d_z=top_z-bottom_z
		#print d_z
		z_comps=[abs(z-top_z)/d_z, abs(z-bottom_z)/d_z]
		my_point=cond_lin_int(bottom_point, top_point, z_comps)
		#bottom_point*z_comps[0]+top_point*z_comps[1]
	else:
		my_point=0.0
	return my_point

#print get_value_at_xyz_testmode(gp_0740, 25.0*1000.0,25.0*1000.0,900.0, var)
#mycap=make_cappi(scanme, xar, yar, 2500.0, 'CZ') 





def make_cappi(scan_list, xar, yar, z, var, **kwargs):
	max_range=kwargs.get('max_range', 150.0*1000.0)
	baddata=kwargs.get('baddata', 1.31072000e+05)
	xx,yy=meshgrid(xar,yar)
	cappi=zeros(xx.shape, dtype=float)
	for i in range(xx.shape[0]):
		print 'Doing another row'
		for j in range(yy.shape[1]):
			if sqrt(xx[i,j]**2+yy[i,j]**2) < max_range:
				cappi[i,j]=get_value_at_xyz(scan_list, xx[i,j],yy[i,j],z, var)
			else:
				cappi[i,j]=baddata
	return cappi

def make_cappi_dict(scan_list, xar, yar, z, **kwargs):
	max_range=kwargs.get('max_range', 150.0*1000.0)
	baddata=kwargs.get('baddata', 1.31072000e+05)
	avail_vars=get_values_at_xyz(scan_list, 100.0, 100.0, 100.0).keys()
	xx,yy=meshgrid(xar,yar)
	cappi=dict([(var,zeros(xx.shape, dtype=float)) for var in avail_vars])
	for i in range(xx.shape[0]):
		#print 'Doing another row', yy[i,j]
		for j in range(yy.shape[1]):
			if sqrt(xx[i,j]**2+yy[i,j]**2) < max_range:
				all_vars=get_values_at_xyz(scan_list, xx[i,j],yy[i,j],z, **kwargs)
				for var in avail_vars:
					cappi[var][i,j]=all_vars[var]
			else:
				for var in avail_vars:
					cappi[var][i,j]=baddata
	return cappi

def make_cube(scan_list,xar, yar, levs, **kwargs):
	displacement=kwargs.get('displacement', [0.0, 0.0])
	avail_vars=get_values_at_xyz(scan_list, 100.0, 100.0, 100.0).keys()
	avail_vars.append('datenum')
	cube_dict=dict([(var,zeros([len(xar), len(yar), len(levs)], dtype=float)) for var in avail_vars])
	for levnum in range(len(levs)):
		print "doing level ", levnum+1, " of ", len(levs)
		this_cap=make_cappi_dict(scan_list, xar-displacement[0], yar-displacement[1], levs[levnum], **kwargs)
		for var in avail_vars:
			cube_dict[var][:,:,levnum]=this_cap[var]
	cube_dict.update({'xar':xar, 'yar':yar,'levs':levs})
	cube_dict.update({'radar_loc': scan_list[0]['radar_loc']})
	cube_dict.update({'displacement': displacement})
	cube_dict.update({'date':scan_list[0]['date']})
	cube_dict.update({'end_date':scan_list[-1]['date']})
	cube_dict.update({'radar_name':scan_list[0]['radar_name']})
	return cube_dict


def make_cappi_testmode(scan_list, xar, yar, z, var, **kwargs):
	max_range=kwargs.get('max_range', 150.0*1000.0)
	baddata=kwargs.get('baddata', 1.31072000e+05)
	xx,yy=meshgrid(xar,yar)
	cappi=zeros(xx.shape, dtype=float)
	for i in range(xx.shape[0]):
		print 'Doing another row'
		for j in range(yy.shape[1]):
			if sqrt(xx[i,j]**2+yy[i,j]**2) < max_range:
				cappi[i,j]=get_value_at_xyz_testmode(scan_list, xx[i,j],yy[i,j],z, var)
			else:
				cappi[i,j]=baddata
	return cappi

def make_rhi(scan_list, sar, zar, az, var, **kwargs):
	max_range=kwargs.get('max_range', 150.0*1000.0)
	baddata=kwargs.get('baddata', 1.31072000e+05)
	ss,zz=meshgrid(sar,zar)
	rhi=zeros(ss.shape, dtype=float)
	for i in range(ss.shape[0]):
		print 'Doing another row'
		for j in range(zz.shape[1]):
			xx=ss[i,j]*sin(az*pi/180.)
			yy=ss[i,j]*cos(az*pi/180.)
			#print xx, yy, zz[i,j]
			if sqrt(xx**2+yy**2) < max_range:
				rhi[i,j]=get_value_at_xyz(scan_list, xx,yy,zz[i,j], var)
			else:
				rhi[i,j]=baddata
	return rhi


























