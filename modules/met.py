#Need to comment
import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
from numpy import *
import mathematics
from time import time as systime

def make_mask(radar1, radar2, angs, min_refl, max_refl, baddata=1.31072000e+05,vel='VE'):
	mask=ones(radar1['CZ'].shape, dtype=float)
	for i in range(radar1['CZ'].shape[0]):
		for j in range(radar1['CZ'].shape[1]):
			ang_test= not((angs[i,j] < 150.0) and (angs[i,j] > 20.0))
			for k in range(radar1['CZ'].shape[2]):
				upper_limit=(radar1['CZ'][i,j,k] > max_refl) or (radar2['CZ'][i,j,k] > max_refl)
				lower_limit=(radar1['CZ'][i,j,k] < min_refl) or (radar2['CZ'][i,j,k] < min_refl)
				vrmasked=(radar1[vel][i,j,k] == baddata) or (radar2[vel][i,j,k] ==baddata)
				vrupper=(radar1[vel][i,j,k] > 100.) or (radar2[vel][i,j,k] >100.)
				if upper_limit or lower_limit or vrmasked or ang_test or vrupper: mask[i,j,k]=0.0
	return mask

def make_mask_bad2(radar1, radar2, angs, min_refl, max_refl, baddata=1.31072000e+05,vel='VE'):
	mask=ones(radar1['CZ'].shape, dtype=float)
	for i in range(radar1['CZ'].shape[0]):
		for j in range(radar1['CZ'].shape[1]):
			ang_test= not((angs[i,j] < 150.0) and (angs[i,j] > 20.0))
			for k in range(radar1['CZ'].shape[2]):
				upper_limit=(radar1['CZ'][i,j,k] > max_refl)
				lower_limit=(radar1['CZ'][i,j,k] < min_refl)
				vrmasked=(radar1[vel][i,j,k] == baddata) or (radar2[vel][i,j,k] ==baddata)
				vrupper=(radar1[vel][i,j,k] > 100.) or (radar2[vel][i,j,k] >100.)
				if upper_limit or lower_limit or vrmasked or ang_test or vrupper: mask[i,j,k]=0.0
	return mask

def make_mask_bad1(radar1, radar2, angs, min_refl, max_refl, baddata=1.31072000e+05,vel='VE', acc_angs=[20.0, 150.0]):
	mask=ones(radar1['CZ'].shape, dtype=float)
	for i in range(radar1['CZ'].shape[0]):
		for j in range(radar1['CZ'].shape[1]):
			ang_test= not((angs[i,j] < acc_angs[1]) and (angs[i,j] > acc_angs[0]))
			for k in range(radar1['CZ'].shape[2]):
				upper_limit=(radar2['CZ'][i,j,k] > max_refl)
				lower_limit= (radar2['CZ'][i,j,k] < min_refl)
				vrmasked=(radar1[vel][i,j,k] == baddata) or (radar2[vel][i,j,k] ==baddata)
				vrupper=(radar1[vel][i,j,k] > 100.) or (radar2[vel][i,j,k] >100.)
				if upper_limit or lower_limit or vrmasked or ang_test or vrupper: mask[i,j,k]=0.0
	return mask

def make_mask_RH(radar1, radar2, angs, min_refl, max_refl, baddata=1.31072000e+05,vel='VE', rh=0.5):
	if 'RH' in radar1.keys():
		pol_radar=radar1
	elif  'RH' in radar2.keys():
		pol_radar=radar2
	else:
		raise ValueError, 'One of the CAPPIs needs to have RH'
	mask=ones(radar1['CZ'].shape, dtype=float)
	for i in range(radar1['CZ'].shape[0]):
		for j in range(radar1['CZ'].shape[1]):
			ang_test= not((angs[i,j] < 150.0) and (angs[i,j] > 20.0))
			for k in range(radar1['CZ'].shape[2]):
				upper_limit=(radar1['CZ'][i,j,k] > max_refl) or (radar2['CZ'][i,j,k] > max_refl)
				lower_limit=(radar1['CZ'][i,j,k] < min_refl) or (radar2['CZ'][i,j,k] < min_refl)
				vrmasked=(radar1[vel][i,j,k] == baddata) or (radar2[vel][i,j,k] ==baddata)
				vrupper=(radar1[vel][i,j,k] > 100.) or (radar2[vel][i,j,k] >100.)
				not_hydrometeor=pol_radar['RH'][i,j,k]  < rh
				if not_hydrometeor or upper_limit or lower_limit or vrmasked or ang_test or vrupper: mask[i,j,k]=0.0
	return mask


def terminal_velocity(refl, temps, levs, display=False):
	bi=0.25
	H=10000.0#I assume this is metres
	ai=4.836
	Ni=1000.0*exp(-0.639)
	vt_param2=1.38*ai*(pi*900.0*Ni)**(-1.0*bi/3.0)
	#r0=refl
	VT_PARAM=14.08
	term_vel=zeros(refl.shape, dtype=float)
	mask_reflect=1.0#dBZ
	for i in range(len(levs)):
		mask=(refl[:,:,i]/mask_reflect).round().clip(min=0., max=1.0)
		r0=10.0**((refl[:,:,i]-43.1)/17.5 -3.0)
		if temps[i] >= 0.0:
			cz=-1.0*VT_PARAM*exp(levs[i]/(2.0*H))
			vt_slice=cz*(r0**0.125)
		else:
			cz=-1.0*vt_param2*exp(levs[i]/(2.0*H))
			vt_slice=cz*(r0**(bi/3.0))
		term_vel[:,:,i]=vt_slice*mask
		if display: print 'Height=', levs[i], ' T=',temps[i], ' Av Fallspeed=', (vt_slice*mask).sum()/mask.sum(), " Max fallspeed=", (vt_slice*mask).min()
	return term_vel


#print ((interp_sonde['tdry(degs)'][i]+273.0)*interp_sonde['press(hPa)'][i-1])/((interp_sonde['tdry(degs)'][i-1]+273.0)*interp_sonde['press(hPa)'][i])
def w_from_continuity(u_array, v_array, dx, dy, levs, pres, temps_c, submask):
	#calculate divergence
	#t0=systime()
	w_array_up=u_array*0.0
	w_array_down=u_array*0.0
	divr=u_array*0.0
	nl=u_array.shape[2]
	#print systime()-t0
	for k in range(nl):
		divr[:,:,k]=(mathematics.dx2d(u_array[:,:,k])/dx + mathematics.dy2d(v_array[:,:,k])/dy)*submask[:,:,k]
	#upward integeral, 
	#initial surface
	#w_array[:,:,0]=-1.0*((p0*(interp_sonde['tdry(degs)'][0]+273.0))/((t0+273.0)*interp_sonde['press(hPa)'][0]))*
	#print systime()-t0
	for k in range(nl-1):
		w_array_up[:,:,k+1]=(((temps_c[k+1]+273.0)*pres[k])/((temps_c[k]+273.0)*pres[k+1]))*(w_array_up[:,:,k]-divr[:,:,k]*(levs[k+1]-levs[k]))
		k_bot=nl-(k+1)
		#print (((temps_c[k+1]+273.0)*pres[k])/((temps_c[k]+273.0)*pres[k+1]))
		w_array_down[:,:,k_bot-1]=(((temps_c[k_bot-1]+273.0)*pres[k_bot])/((temps_c[k_bot]+273.0)*pres[k_bot-1]))*(w_array_down[:,:,k_bot]+divr[:,:,k_bot]*(levs[k_bot-1]-levs[k_bot]))
		#print  (((temps_c[k_bot-1]+273.0)*pres[k_bot])/((temps_c[k_bot]+273.0)*pres[k_bot-1]))
	#print systime()-t0
	return w_array_up, w_array_down

def make_submask(mask):
	submask=mask*0.0
	nx,ny,nz=mask.shape
	for i in range(nx):
		for j in range(ny):
			for k in range(nz):
				if not((i==0) or (i==(nx-1)) or (j==0) or (j==(ny-1))):
					submask[i,j,k]=min([mask[i,j,k], mask[i-1,j-1,k], mask[i+1, j+1, k], mask[i-1, j+1, k], mask[i,j+1,k], mask[i+1, j, k], mask[i, j-1, k], mask[i-1, j, k]])
	return submask 

#for k in range(nl-1):
#	print k
#	print k+1
#	k_bot=nl-(k+1)
#	print k_bot
#	print k_bot-1
	