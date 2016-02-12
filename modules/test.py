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
from numpy import sqrt, array, argsort, abs

def get_value_at_xyz(scan_list, x,y,z, var):
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
		top_scan=scan_list[elev_numbers[0]]
		bottom_scan=scan_list[elev_numbers[1]]
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
		top_ang_comp=[abs(az-top_scan['Azmth'][ray_numbers[0]])/d_angle, abs(az-top_scan['Azmth'][ray_numbers[1]])/d_angle]
		if fits_order[0] == fits_order[1]:
			#we are near 360
			azd=az
			ray_numbers=[scan_order[-1], scan_order[0]]
			d_angle=top_scan['Azmth'][ray_numbers[1]]-(top_scan['Azmth'][ray_numbers[0]]-360.0)
			if az > d_angle: azd=az-360.0
			top_angle_comp=[abs(azd-top_scan['Azmth'][ray_numbers[1]])/d_angle,abs(azd-(top_scan['Azmth'][ray_numbers[0]]-360.0))/d_angle ]
		top_point=top_ang_comp[0]*(top_scan[var][ray_numbers[0], gates[0]]*top_gates_comp[0]+ top_scan[var][ray_numbers[0], gates[1]]*top_gates_comp[1])+ top_ang_comp[1]*(top_scan[var][ray_numbers[1], gates[0]]*top_gates_comp[0]+ top_scan[var][ray_numbers[1], gates[1]]*top_gates_comp[1])
		
		top_z=(top_scan['zar'][ray_numbers[0], gates[0]]*top_gates_comp[0]+ top_scan['zar'][ray_numbers[0], gates[1]]*top_gates_comp[1])
		
		#Scan below the point ----------BOTTOM ELEVATION-----------
		#work out what gates to use (ie range in S direction)
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
		bottom_ang_comp=[abs(az-bottom_scan['Azmth'][ray_numbers[0]])/d_angle, abs(az-bottom_scan['Azmth'][ray_numbers[1]])/d_angle]
		if fits_order[0] == fits_order[1]:
			#we are near 360
			azd=az
			ray_numbers=[scan_order[-1], scan_order[0]]
			d_angle=bottom_scan['Azmth'][ray_numbers[1]]-(bottom_scan['Azmth'][ray_numbers[0]]-360.0)
			if az > d_angle: azd=az-360.0
			bottom_angle_comp=[abs(azd-bottom_scan['Azmth'][ray_numbers[1]])/d_angle,abs(azd-(bottom_scan['Azmth'][ray_numbers[0]]-360.0))/d_angle ]
		bottom_point=bottom_ang_comp[0]*(bottom_scan[var][ray_numbers[0], gates[0]]*bottom_gates_comp[0]+ bottom_scan[var][ray_numbers[0], gates[1]]*bottom_gates_comp[1])+ bottom_ang_comp[1]*(bottom_scan[var][ray_numbers[1], gates[0]]*bottom_gates_comp[0]+ bottom_scan[var][ray_numbers[1], gates[1]]*bottom_gates_comp[1])
		bottom_z=(bottom_scan['zar'][ray_numbers[0], gates[0]]*bottom_gates_comp[0]+ bottom_scan['zar'][ray_numbers[0], gates[1]]*bottom_gates_comp[1])
		
		#calculate components
		d_z=top_z-bottom_z
		z_comps=[abs(z-top_z)/d_z, abs(z-bottom_z)/d_z]
		my_point=bottom_point*z_comps[0]+top_point*z_comps[1]
	else:
		my_point=0.0
	return my_point
