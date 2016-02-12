import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
import mathematics
import propigation
from numpy import zeros, concatenate, abs,sqrt, array, where, unique, sin, pi, transpose


def meshcube(x,y,z):
	xcube=zeros([len(x),len(y),len(z)], dtype=float)
	ycube=zeros([len(x),len(y),len(z)], dtype=float)
	zcube=zeros([len(x),len(y),len(z)], dtype=float)
	for i in range(len(x)):
		for j in range(len(y)):
			for k in range(len(z)):
				xcube[i,j,k]=x[i]
				ycube[i,j,k]=y[j]
				
				
				zcube[i,j,k]=z[k]
	return xcube, ycube, zcube

def radarcoords(cx,cy,cz):
	cs=sqrt(cx**2+cy**2)
	cele=propigation.elev_of_scan(cs,cz)
	caz=mathematics.vector_direction(cx,cy)
	return cs, caz, cele


def get_ray_bundle(scanlist, az, ele, az_inf, ele_inf ):
	#get a bundle of rays within a certain azimuthal and elevation angle 
	gated_vars=set(scanlist[0].keys()) & set(['VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW', 'i_comp', 'j_comp', 'k_comp', 'xar', 'yar', 'zar'])
	elevs=array([scanlist[i]['Elev'][0] for i in range(len(scanlist))])
	az_in_infl=[]
	ele_in_infl=[]
	rays=[]
	if not((ele > elevs.max()) or (ele < elevs.min())):
		#we are bounded by radar scans
		ele_inf_list=where(abs(elevs-ele) < ele_inf)[0]
		for ele_num in ele_inf_list:
			azs1=where(abs(scanlist[ele_num]['Azmth']-az) < az_inf)[0]
			azs2=where(abs(scanlist[ele_num]['Azmth']+360.0-az) < az_inf)[0]
			azs3=where(abs(scanlist[ele_num]['Azmth']-360.0-az) < az_inf)[0]
			azs=concatenate((azs1,azs2,azs3), axis=1)
			for az_num in azs:
				az_in_infl.append(scanlist[ele_num]['Azmth'][az_num])
				ele_in_infl.append(elevs[ele_num])
				rays.append(dict([(var, scanlist[ele_num][var][az_num,:]) for var in gated_vars]))
	return rays, az_in_infl, ele_in_infl


def get_ray_slice(scanlist, az, az_inf):
	#get a bundle of rays within a certain azimuthal and elevation angle 
	gated_vars=set(scanlist[0].keys()) & set(['VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW', 'i_comp', 'j_comp', 'k_comp', 'xar', 'yar', 'zar'])
	elevs=array([scanlist[i]['Elev'][0] for i in range(len(scanlist))])
	az_in_infl=[]
	ele_in_infl=[]
	rays=[]
	for ele_num in range(len(elevs)):
		azs1=where(abs(scanlist[ele_num]['Azmth']-az) < az_inf)[0]
		azs2=where(abs(scanlist[ele_num]['Azmth']+360.0-az) < az_inf)[0]
		azs3=where(abs(scanlist[ele_num]['Azmth']-360.0-az) < az_inf)[0]
		azs=concatenate((azs1,azs2,azs3), axis=1)
		for az_num in azs:
			az_in_infl.append(scanlist[ele_num]['Azmth'][az_num])
			ele_in_infl.append(elevs[ele_num])
			rays.append(dict([(var, scanlist[ele_num][var][az_num,:]) for var in gated_vars]))
	return rays, az_in_infl, ele_in_infl



def get_ray_s(ray, s, badval=1.31072000e+05):
	sar=sqrt(ray['xar']**2 + ray['yar']**2)
	ht=ray['zar']
	elements=list(mathematics.where_fits(s, sar))
	point={}
	if elements[0] ==len(sar)-1:
		point.update(dict([(var, badval) for var in set([ 'i_comp', 'j_comp', 'k_comp', 'VE','VR', 'CZ', 'RH', 'PH', 'ZD', 'SW'])& set(ray.keys())]))
		#point.update(dict([(var, badval) for var in set([ 'VE',])& set(ray.keys())]))
	else:
		#print elements
		if sar[elements[1]]-sar[elements[0]]!=0:
			w=1.0-abs((sar[elements]-s)/(sar[elements[1]]-sar[elements[0]]))
		else:
			w=[0.5,0.5]
		#print w
		my_ht=(ht[elements]*w).sum()
		if not(badval in ray['VE'][elements]):
			point.update(dict([(var,(ray[var][elements]*w).sum()) for var in [ 'i_comp', 'j_comp', 'k_comp', 'VE']]))
			#point.update(dict([(var,(ray[var][elements]*w).sum()) for var in [ 'VE']]))
		else:
			point.update(dict([(var,badval) for var in [ 'i_comp', 'j_comp', 'k_comp', 'VE']]))
			#point.update(dict([(var,badval) for var in [ 'VE']]))
		for var in set(['VR', 'CZ', 'RH', 'PH', 'ZD', 'SW'])& set(ray.keys()):
			if not(badval in ray[var][elements]):
				point.update({var:(ray[var][elements]*w).sum()})
			else:
				point.update({var:badval})
	return point, my_ht


def get_ray_ht(ray, ht, badval=1.31072000e+05):
	sar=sqrt(ray['xar']**2 + ray['yar']**2)
	elements=list(mathematics.where_fits(ht, ray['zar']))
	point={}
	if elements[0] ==len(ray['zar'])-1:
		point.update(dict([(var, badval) for var in set([ 'i_comp', 'j_comp', 'k_comp', 'VE','VR', 'CZ', 'RH', 'PH', 'ZD', 'SW'])& set(ray.keys())]))
		#point.update(dict([(var, badval) for var in set([ 'VE',])& set(ray.keys())]))
	else:
		#print elements
		if ray['zar'][elements[1]]-ray['zar'][elements[0]] !=0:
			w=1.0-abs((ray['zar'][elements]-ht)/(ray['zar'][elements[1]]-ray['zar'][elements[0]]))
		else:
			w=[0.5,0.5]
		#print w
		my_s=(sar[elements]*w).sum()
		if not(badval in ray['VE'][elements]):
			point.update(dict([(var,(ray[var][elements]*w).sum()) for var in [ 'i_comp', 'j_comp', 'k_comp', 'VE']]))
			#point.update(dict([(var,(ray[var][elements]*w).sum()) for var in [ 'VE']]))
		else:
			point.update(dict([(var,badval) for var in [ 'i_comp', 'j_comp', 'k_comp', 'VE']]))
			#point.update(dict([(var,badval) for var in [ 'VE']]))
		for var in set(['VR', 'CZ', 'RH', 'PH', 'ZD', 'SW'])& set(ray.keys()):
			if not(badval in ray[var][elements]):
				point.update({var:(ray[var][elements]*w).sum()})
			else:
				point.update({var:badval})
	return point, my_s


def make_stack_horiz(bundle,baz, bele, ele_array,az, ht_array, badval=1.31072000e+05):
	abele=array(bele)
	p=0
	#print len(ht_array)
	data_dict=dict([(var,[]) for var in get_ray_ht(bundle[0], ht_array[0]).keys()])
	cbaz=zeros([len(baz)], dtype=float)
	for i in range(len(baz)): #deal with the 360 degree problem
		poss=array([baz[i], (baz[i]+360.0), (baz[i]-360.0)])
		cbaz[i]=poss[where(abs(poss-az) ==abs(poss-az).min())[0]] 
	for i in range(len(ht_array)):
		if ((ele_array[i] > abele.min()) and (ele_array[i] < abele.max())):
			uniq_eles=unique(bele)
			my_ele=uniq_eles[list(mathematics.where_fits(ele_array[i], uniq_eles))]
			az1=cbaz[where(bele==my_ele[0])[0]][list(mathematics.where_fits(az, cbaz[where(bele==my_ele[0])[0]]))]
			az2=cbaz[where(bele==my_ele[1])[0]][list(mathematics.where_fits(az, cbaz[where(bele==my_ele[1])[0]]))]
			rays=[]
			rays.append(where(cbaz==az1[0])[0][where(abele[where(cbaz==az1[0])[0]]==my_ele[0])[0]][0])
			rays.append(where(cbaz==az1[1])[0][where(abele[where(cbaz==az1[1])[0]]==my_ele[0])[0]][0])
			rays.append(where(cbaz==az2[0])[0][where(abele[where(cbaz==az2[0])[0]]==my_ele[1])[0]][0])
			rays.append(where(cbaz==az2[1])[0][where(abele[where(cbaz==az2[1])[0]]==my_ele[1])[0]][0])
			small_bundle=array(bundle)[rays]
			small_az=cbaz[rays]
			small_ele=abele[rays]
			dat=[]
			for spray in small_bundle:
				dat.append(get_ray_ht(spray, ht_array[i]))
			dat=array(dat)
			#print rays
			#print small_ele
			#print small_az
			#print dat
			if small_az[1]-small_az[0] !=0.0:
				wdown=1.0-abs((small_az[[0,1]]-az)/(small_az[1]-small_az[0]))
			else:
				wdown=[0.5,0.5]
			if small_az[3]-small_az[2] != 0.0:
				wup=1.0-abs((small_az[[2,3]]-az)/(small_az[3]-small_az[2]))
			else:
				wup=[0.5,0.5]
			if small_ele[2]-small_ele[0]!=0.0:
				w=1.0-abs((small_ele[[0,2]]-ele_array[i])/(small_ele[2]-small_ele[0]))
			else:
				w=[0.5,0.5]
			for var in dat[0].keys():
					#print 'here', [
				if not(badval in [dat[0][var], dat[1][var]]):
					#data_dict[var].append((wdown*dat[0:1][var]).sum())
					vdown=(wdown*array([dat[0][var], dat[1][var]])).sum()
				else:
					vdown=badval
				if not(badval in [dat[2][var], dat[3][var]]):
					#data_dict[var].append((wdown*dat[0:1][var]).sum())
					vup=(wdown*array([dat[2][var], dat[3][var]])).sum()
				else:
					vup=badval
				if not(badval in [vdown, vup]):
					data_dict[var].append((w*array([vdown,vup])).sum())
					#print "good val",p
					#p=p+1
				else:
					data_dict[var].append(badval)
					#print 'badval', p
					#p=p+1
		else:
			for var in data_dict.keys():
				data_dict[var].append(badval)
				#print "out of range", p
				#p=p+1
					
	#print data_dict['VE']
	return data_dict

def make_stack_all(bundle,baz, bele, ele_array,s,az, ht_array, max_elev, badval=1.31072000e+05, order=0.5):
	abele=array(bele)
	#print len(ht_array)
	data_dict=dict([(var,[]) for var in get_ray_ht(bundle[0], ht_array[0])[0].keys()])
	cbaz=zeros([len(baz)], dtype=float)
	for i in range(len(baz)): #deal with the 360 degree problem
		poss=array([baz[i], (baz[i]+360.0), (baz[i]-360.0)])
		cbaz[i]=poss[where(abs(poss-az) ==abs(poss-az).min())[0]] 
	for i in range(len(ht_array)):
		if ((ele_array[i] > abele.min()) and (ele_array[i] < abele.max())):
			uniq_eles=unique(bele)
			my_ele=uniq_eles[list(mathematics.where_fits(ele_array[i], uniq_eles))]
			az1=cbaz[where(bele==my_ele[0])[0]][list(mathematics.where_fits(az, cbaz[where(bele==my_ele[0])[0]]))]
			az2=cbaz[where(bele==my_ele[1])[0]][list(mathematics.where_fits(az, cbaz[where(bele==my_ele[1])[0]]))]
			rays=[]
			rays.append(where(cbaz==az1[0])[0][where(abele[where(cbaz==az1[0])[0]]==my_ele[0])[0]][0])
			rays.append(where(cbaz==az1[1])[0][where(abele[where(cbaz==az1[1])[0]]==my_ele[0])[0]][0])
			rays.append(where(cbaz==az2[0])[0][where(abele[where(cbaz==az2[0])[0]]==my_ele[1])[0]][0])
			rays.append(where(cbaz==az2[1])[0][where(abele[where(cbaz==az2[1])[0]]==my_ele[1])[0]][0])
			small_bundle=array(bundle)[rays]
			small_az=cbaz[rays]
			small_ele=abele[rays]
			small_z=[]
			small_s=[]
			dat=[]
			vdat=[]
			for spray in small_bundle:
				ht_v, my_s=get_ray_ht(spray, ht_array[i])
				s_v, my_ht=get_ray_s(spray, s)
				dat.append(ht_v)
				small_s.append(my_s)
				vdat.append(s_v)
				small_z.append(my_ht)
			dat=array(dat)
			vdat=array(vdat)
			small_s=array(small_s)
			small_z=array(small_z)
			#print rays
			#print small_ele
			#print small_az
			#print dat
			if small_az[1]-small_az[0] !=0.0:
				wdown=1.0-abs((small_az[[0,1]]-az)/(small_az[1]-small_az[0]))
			else:
				wdown=[0.5,0.5]
			if small_az[3]-small_az[2] != 0.0:
				wup=1.0-abs((small_az[[2,3]]-az)/(small_az[3]-small_az[2]))
			else:
				wup=[0.5,0.5]
			if small_s[2]-small_s[0]!=0.0:
				w=1.0-abs((small_s[[0,2]]-s)/(small_s[2]-small_s[0]))
			else:
				w=[0.5,0.5]
			if small_z[2]-small_z[0]!=0.0:
				wz=1.0-abs((small_z[[0,2]]-ht_array[i])/(small_z[2]-small_z[0]))
			else:
				wz=[0.5,0.5]
			hwt=h_wts(ele_array[i], max_elev, order)
			#print hwt
			for var in dat[0].keys():
					#print 'here', [
				if not(badval in [dat[0][var], dat[1][var]]):
					#data_dict[var].append((wdown*dat[0:1][var]).sum())
					vdown=(wdown*array([dat[0][var], dat[1][var]])).sum()
				else:
					vdown=badval
				if not(badval in [dat[2][var], dat[3][var]]):
					#data_dict[var].append((wdown*dat[0:1][var]).sum())
					vup=(wdown*array([dat[2][var], dat[3][var]])).sum()
				else:
					vup=badval
				if not(badval in [vdown, vup]):
					horiz_int=(w*array([vdown,vup])).sum()
					#print "good val",p
					#p=p+1
				else:
					horiz_int=badval
					#print 'badval', p
					#p=p+1
				#VERTICAL
				if not(badval in [vdat[0][var], vdat[1][var]]):
					#data_dict[var].append((wdown*dat[0:1][var]).sum())
					vdown=(wdown*array([vdat[0][var], vdat[1][var]])).sum()
				else:
					vdown=badval
				if not(badval in [vdat[2][var], vdat[3][var]]):
					#data_dict[var].append((wdown*dat[0:1][var]).sum())
					vup=(wdown*array([vdat[2][var], vdat[3][var]])).sum()
				else:
					vup=badval
				if not(badval in [vdown, vup]):
					vert_int=(wz*array([vdown,vup])).sum()
					#print "good val",p
					#p=p+1
				else:
					vert_int=badval
					#print 'badval', p
					#p=p+1
				#data_dict[var].append(horiz_int)
				hbad=horiz_int==badval
				vbad=vert_int==badval
				if not(hbad) and not(vbad):
					data_dict[var].append(hwt*horiz_int+(1.0-hwt)*vert_int)
				elif hbad and not(vbad):
					data_dict[var].append(vert_int)
				elif vbad and not(hbad):
					data_dict[var].append(horiz_int)
				else:
					data_dict[var].append(badval)
		else:
			for var in data_dict.keys():
				data_dict[var].append(badval)
				#print "out of range", p
				#p=p+1
	#print data_dict['VE']
	return data_dict


def make_cube_horiz(scan_list,xar, yar, levs, **kwargs):
	displacement=kwargs.get('displacement', [0.0, 0.0])
	badval=kwargs.get( 'badval',1.31072000e+05)
	bundle_size=kwargs.get('bundle_size', 1.6)
	good_thresh=kwargs.get('good_thresh', 0.3)
	#avail_vars=get_values_at_xyz(scan_list, 100.0, 100.0, 100.0).keys()
	#avail_vars.append('datenum')
	avail_vars=set(scan_list[0].keys()) & set(['VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW', 'i_comp', 'j_comp', 'k_comp'])
	cube_dict=dict([(var,zeros([len(xar), len(yar), len(levs)], dtype=float)) for var in avail_vars])
	cx,cy,cz=meshcube(xar,yar,levs)
	cs,caz,cele=radarcoords(cx,cy,cz)
	for i in range(cx.shape[0]):
		print "doing row ",i," of ", cx.shape[0]
		for j in range(cy.shape[1]):
			#print j
			bundle,azmths,elevs=get_ray_slice(scan_list, caz[i,j,0], bundle_size)
			mystack=make_stack_horiz(bundle,array(azmths),array(elevs),cele[i,j,:], caz[i,j,0], cz[i,j,:])# make_stack_horiz(bundle,baz, bele, ele_array,az, ht_array, badval=1.31072000e+05):
			#print mystack
			#print  caz[i,j,0]
			for var in mystack.keys():
				#print var
				#print array(mystack[var]).shape
				cube_dict[var][i,j,:]=array(mystack[var])
	cube_dict.update({'xar':xar, 'yar':yar,'levs':levs})
	cube_dict.update({'radar_loc': scan_list[0]['radar_loc']})
	cube_dict.update({'displacement': displacement})
	cube_dict.update({'date':scan_list[0]['date']})
	cube_dict.update({'end_date':scan_list[-1]['date']})
	cube_dict.update({'radar_name':scan_list[0]['radar_name']})
	return cube_dict


def make_cube_all(scan_list,xar, yar, levs, **kwargs):
	displacement=kwargs.get('displacement', [0.0, 0.0])
	badval=kwargs.get( 'badval',1.31072000e+05)
	bundle_size=kwargs.get('bundle_size', 1.6)
	good_thresh=kwargs.get('good_thresh', 0.3)
	order=kwargs.get('order', 0.5)
	max_el=array([scan['Elev'][0] for scan in scan_list]).max()
	#make_stack_all(bundle,baz, bele, ele_array,s,az, ht_array, max_elev, badval=1.31072000e+05, order=0.5):
	#avail_vars=get_values_at_xyz(scan_list, 100.0, 100.0, 100.0).keys()
	#avail_vars.append('datenum')
	avail_vars=set(scan_list[0].keys()) & set(['VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW', 'i_comp', 'j_comp', 'k_comp'])
	cube_dict=dict([(var,zeros([len(xar), len(yar), len(levs)], dtype=float)) for var in avail_vars])
	cx,cy,cz=meshcube(xar-displacement[0],yar-displacement[1],levs)
	cs,caz,cele=radarcoords(cx,cy,cz)
	for i in range(cx.shape[0]):
		print "doing row ",i," of ", cx.shape[0]
		for j in range(cy.shape[1]):
			#print j
			bundle,azmths,elevs=get_ray_slice(scan_list, caz[i,j,0], bundle_size)
			mystack=make_stack_all(bundle,array(azmths),array(elevs),cele[i,j,:], cs[i,j,0], caz[i,j,0], cz[i,j,:], max_el, order=order)# make_stack_horiz(bundle,baz, bele, ele_array,az, ht_array, badval=1.31072000e+05):
			#print mystack
			#print  caz[i,j,0]
			for var in mystack.keys():
				#print var
				#print array(mystack[var]).shape
				cube_dict[var][j,i,:]=array(mystack[var]) #uber dodge i,j reveresal as some where we transposed... WATCH OUT!
	cube_dict.update({'xar':xar, 'yar':yar,'levs':levs})
	cube_dict.update({'radar_loc': scan_list[0]['radar_loc']})
	cube_dict.update({'displacement': displacement})
	cube_dict.update({'date':scan_list[0]['date']})
	cube_dict.update({'end_date':scan_list[-1]['date']})
	cube_dict.update({'radar_name':scan_list[0]['radar_name']})
	return cube_dict



def h_wts(elev, max_elev, n):
	return 0.0*((sin(pi*elev/max_elev -pi/2.0)+1.0)/2.0)**n


def blend(vert_cube, horiz_cube, max_el,**kwargs):
	badval=kwargs.get( 'badval',1.31072000e+05)
	loud=kwargs.get('loud', False)
	order=kwargs.get('order', 0.5)
	new_dict={}
	gated_vars=(set(['VE', 'VR', 'CZ', 'RH', 'PH', 'ZD', 'SW'])&set(vert_cube.keys()))&set(horiz_cube.keys())
	copy_keys=set(vert_cube.keys())-gated_vars
	if loud:
		print "Gated Vars", gated_vars
		print "Copy keys", copy_keys
	cx,cy,cz=meshcube(vert_cube['xar'],vert_cube['yar'],vert_cube['levs'])
	cs,caz,cele=radarcoords(cx,cy,cz)
	wts_h=h_wts(cele, max_el, order)
	wts_v=1.0-wts_h
	for key in copy_keys:
		new_dict.update({key:vert_cube[key]})
	for key in gated_vars:
		new_cube=zeros(wts_h.shape, dtype=float)
		for k in range(vert_cube[key].shape[2]):
			h_slice=transpose(horiz_cube[key][:,:,k])
			v_slice=vert_cube[key][:,:,k]
			for i in range(h_slice.shape[0]):
				for j in range(h_slice.shape[1]):
					hbad=h_slice[i,j]==badval
					vbad=v_slice[i,j]==badval
					if not(hbad or vbad):
						new_cube[i,j,k]=wts_v[i,j,k]*v_slice[i,j] +wts_h[i,j,k]*h_slice[i,j]
					elif hbad and not(vbad):
						new_cube[i,j,k]=v_slice[i,j]
					elif vbad and not(hbad):
						new_cube[i,j,k]=h_slice[i,j]
					elif vbad and hbad:
						new_cube[i,j,k]=badval
		new_dict.update({key:new_cube})
	return new_dict


