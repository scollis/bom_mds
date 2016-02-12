import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
from numpy import float, array, zeros, ones, sqrt, min, arange, append, copy
import met
import smooth_cost
from time import time as systime

def meas_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, vr1, vr2, weights):
       vr1_sim=i_cmpt_r1*u_array+j_cmpt_r1*v_array
       vr2_sim=i_cmpt_r2*u_array+j_cmpt_r2*v_array
       cost=0.5*(weights*((vr1-vr1_sim)**2 + (vr2-vr2_sim)**2)).sum()
       gv_u=gv_u + weights*(i_cmpt_r1*(vr1_sim-vr1) + i_cmpt_r2*(vr2_sim-vr2))
       gv_v=gv_v+ weights*(j_cmpt_r1*(vr1_sim-vr1) + j_cmpt_r2*(vr2_sim-vr2))
       return gv_u,gv_v,cost


def gracon_vel2d_3d(gv_u,gv_v,f,u_array,v_array,i_cmpt_r1,j_cmpt_r1,i_cmpt_r2,j_cmpt_r2,vr1,vr2,weights):
	#!dimensions
	submask=met.make_submask(weights)
	nx, ny, nz=u_array.shape
	X=zeros([2*nx*ny*nz], dtype=float)
	G=zeros([2*nx*ny*nz], dtype=float)
	XK1=zeros([2*nx*ny*nz], dtype=float)
	D=zeros([2*nx*ny*nz], dtype=float)
	GK1=zeros([2*nx*ny*nz], dtype=float)
	Fal=0.0
	i_ter=1
	grad=1.0
	nerr=0
	dimx=2*nx*ny*nz
	F=0.0
	acc=0.0
	itr=100
	#!Initial guess packaging
	#do i=1,nx
	print "Doing transfer", nx, ny, nz
	X=array([u_array,v_array])
	G=array([u_array,v_array])
	gradi=sqrt((G**2).sum())
	print "entering loop"
	while ((i_ter <= itr+1) and (grad > acc)):
		if (i_ter==1):
			Fm1=F
			#!----------------------COST FUNCTION--------------------
			#!pack variables
			#do i=1,nx
			#HERE
			u_array=X[0]
			v_array=X[1]
			gv_u=G[0]*0.0
			gv_v=G[1]*0.0
			F=0.0
			#call meas_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1,&
			#i_cmpt_r2, j_cmpt_r2, nx, ny,nz, vr1, vr2, weights)
			gv_u,gv_v,F = meas_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, vr1, vr2, weights)
			#!unpack variables
			#do i=1,nx
			X=array([u_array,v_array])
			G=array([gv_u,gv_v])
			gradi=sqrt((G**2).sum())
			#!----------------------COST FUNCTION--------------------      
		#!subroutine rms2(iterc, f, gradi)
		print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
		#call rms2(iter,F,gradi)        
		#!     -------- To stop in the first iteration -------------------------
		if (itr==0):
			grad = -1.0
			print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
			return gv_u, gv_v, F, u_array, v_array#need to put vars in
		#!     -----------------------------------------------------------------
		#do i=1,dimx
		#for i in range(dimx):
		XK1=X
		#XK1(i)=X(i)
		#!     ===============================================================
		#!     ------------------------Steepest descent-----------------------
		if (i_ter==1):
			gg=0.0
			#do i=1,dimx
			#for i in range(dimx):
			gg=gg+(G**2).sum()#G[i]*G[i]
			D=-G
			GK1=G
			gd=-gg
		else:
			#!    ----------------------Conjugate Gradient-----------------------
			gg=0.0
			upk=0.0
			dok=0.0
			gg=gg+(G**2).sum()
			upk=upk+(G*(G-GK1)).sum()
			dok=dok+(D*(G-GK1)).sum()
			gd=0.0
			#do i=1,dimx
			#for i in range(dimx):
			GK1=G
			if (dok==0.0): dok = 1.0e-10
         		if (upk==0.0): upk = 1.0e-10
         		#D(i)=-G(i)+(upk/dok)*D(i)
			D=-G+(upk/dok)*D
			#!<<<alain
			#!dd=dd+D(i)*D(i)
			#!>>>alain
         		gd=gd+(G*D).sum()
			#!            print*,'Angle between search direction and gradient:',180.0/(4.0*atan(1.0))*acos(-gd/sqrt(dd*gg))
			#!     --------------------------------------------------------------
		#endif
		#!     ===============================================================
		#!     ----------------------- step length ---------------------------
		if (i_ter==1):
			delF=0.5*F
		elif (i_ter!=1):
			 delF=Fm1-F
		if (F!=0.0): rdF=delF/F
		if (rdF < 1.0e-5):
			print'variation of F too weak'
			return gv_u, gv_v, F, u_array, v_array#need to put vars in
		if (gd > 0.0): gd=-gd
		axpect = -2*delF/gd
		al = min(2.0,axpect)
		Fm1=F
		ncomp=0
		subfinish = 0
		#do while (subfinish.eq.0)
		while (subfinish==0):
			#print "in subloop"
			ncomp=ncomp+1
			nerr=0
			subfinish = 1
			if (ncomp>200):
				print'too many subiterations'
				return gv_u, gv_v, F, u_array, v_array#need vars
			#endif
			#do i=1,dimx
			#for i in range(dimx):
			#X(i)=XK1(i) + al*D(i)
			X=XK1 + al*D
			#enddo
			#do i=1,nx
			u_array=X[0]
			v_array=X[1]
			gv_u=G[0]*0.0
			gv_v=G[1]*0.0
			F=0.0
			gv_u,gv_v,F = meas_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, vr1, vr2, weights)
			X=array([u_array,v_array])
			G=array([gv_u,gv_v])
			gradi=sqrt((G**2).sum())
			#!----------------------COST FUNCTION--------------------
			if (nerr==1):
				al=al/2
				subfinish = 0
			#endif
			if (subfinish==1):
				gald=0.0
				#do i=1,dimx
				#for i in range(dimx):
				gald=gald+(G*D).sum()
				#enddo
				if ((gald < 0.0) and (Fal < F)):
					al = 4*al
					subfinish = 0
				#endif
				if (subfinish==1):
					z = 3.0*(F - Fal)/al + gd + gald
					w = sqrt(z**2 - gd*gald)
					amin = al*(1 - (gald + w - z)/(gald - gd + 2*w))
					#do i=1,dimx
					#for i in range(dimx):
					#X(i)=XK1(i) + amin*D(i)
					X=XK1 + amin*D
					#enddo
					#!----------------------COST FUNCTION--------------------
					#!pack variables
					#do i=1,nx
					u_array=X[0]
					v_array=X[1]
					gv_u=G[0]*0.0
					gv_v=G[1]*0.0
					F=0.0
					gv_u,gv_v,F = meas_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, vr1, vr2, weights)
					#!unpack variables
					#do i=1,nx
					X=array([u_array,v_array])
					G=array([gv_u,gv_v])
					gradi=sqrt((G**2).sum())
					#!----------------------COST FUNCTION--------------------         
					if (F > Fm1) or (nerr==1):
						al=al/2.5
						subfinish = 0
					#endif
					if (subfinish==1):
						i_ter=i_ter+1
						grad=gg**0.5
					#endif
				#endif
			#endif
		#enddo
	#enddo
	#do i=1,nx
	u_array=X[0]
	v_array=X[1]
	gv_u=G[0]
	gv_v=G[1]
	print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
	return gv_u, gv_v, F, u_array, v_array
	#end subroutine gracon_vel2d_3d
	
def gracon_3d(gv_u,gv_v,gv_w, f, u_array, v_array, w_array ,i_cmpt_r1, j_cmpt_r1, k_cmpt_r1, i_cmpt_r2, j_cmpt_r2, k_cmpt_r2, vr1, vr2, weights, term_velocity):
	#!dimensions
	nx, ny, nz=u_array.shape
	Fal=0.0
	i_ter=1
	grad=1.0
	nerr=0
	dimx=2*nx*ny*nz
	F=0.0
	acc=0.0
	itr=100
	#!Initial guess packaging
	#do i=1,nx
	print "Doing transfer", nx, ny, nz
	X=array([u_array,v_array, w_array])
	G=array([u_array,v_array, w_array])*0.0
	gradi=sqrt((G**2).sum())
	print "entering loop"
	while ((i_ter <= itr+1) and (grad > acc)):
		if (i_ter==1):
			Fm1=F
			#!----------------------COST FUNCTION--------------------
			#!pack variables
			#do i=1,nx
			#HERE
			u_array=X[0]
			v_array=X[1]
			w_array=X[2]
			gv_u=G[0]*0.0
			gv_v=G[1]*0.0
			gv_w=G[2]*0.0
			
			F=0.0
			#call meas_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1,&
			#i_cmpt_r2, j_cmpt_r2, nx, ny,nz, vr1, vr2, weights)
			gv_u,gv_v, gv_w, F = meas_cost_3d(gv_u, gv_v, gv_w, F, u_array, v_array, w_array, i_cmpt_r1, j_cmpt_r1, k_cmpt_r1, i_cmpt_r2, j_cmpt_r2,k_cmpt_r2, vr1, vr2, weights, term_velocity)
			#!unpack variables
			#do i=1,nx
			X=array([u_array,v_array, w_array])
			G=array([gv_u,gv_v, gv_w])
			gradi=sqrt((G**2).sum())
			#!----------------------COST FUNCTION--------------------      
		#!subroutine rms2(iterc, f, gradi)
		print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
		#call rms2(iter,F,gradi)        
		#!     -------- To stop in the first iteration -------------------------
		if (itr==0):
			grad = -1.0
			print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
			return gv_u, gv_v, F, u_array, v_array#need to put vars in
		#!     -----------------------------------------------------------------
		#do i=1,dimx
		#for i in range(dimx):
		XK1=X
		#XK1(i)=X(i)
		#!     ===============================================================
		#!     ------------------------Steepest descent-----------------------
		if (i_ter==1):
			gg=0.0
			#do i=1,dimx
			#for i in range(dimx):
			gg=gg+(G**2).sum()#G[i]*G[i]
			D=-G
			GK1=G
			gd=-gg
		else:
			#!    ----------------------Conjugate Gradient-----------------------
			gg=0.0
			upk=0.0
			dok=0.0
			gg=gg+(G**2).sum()
			upk=upk+(G*(G-GK1)).sum()
			dok=dok+(D*(G-GK1)).sum()
			gd=0.0
			#do i=1,dimx
			#for i in range(dimx):
			GK1=G
			if (dok==0.0): dok = 1.0e-10
         		if (upk==0.0): upk = 1.0e-10
         		#D(i)=-G(i)+(upk/dok)*D(i)
			D=-G+(upk/dok)*D
			#!<<<alain
			#!dd=dd+D(i)*D(i)
			#!>>>alain
         		gd=gd+(G*D).sum()
			#!            print*,'Angle between search direction and gradient:',180.0/(4.0*atan(1.0))*acos(-gd/sqrt(dd*gg))
			#!     --------------------------------------------------------------
		#endif
		#!     ===============================================================
		#!     ----------------------- step length ---------------------------
		if (i_ter==1):
			delF=0.5*F
		elif (i_ter!=1):
			 delF=Fm1-F
		if (F!=0.0): rdF=delF/F
		if (rdF < 1.0e-5):
			print'variation of F too weak'
			return gv_u, gv_v, gv_w, F, u_array, v_array, w_array#need to put vars in
		if (gd > 0.0): gd=-gd
		axpect = -2*delF/gd
		al = min(2.0,axpect)
		Fm1=F
		ncomp=0
		subfinish = 0
		#do while (subfinish.eq.0)
		while (subfinish==0):
			#print "in subloop"
			ncomp=ncomp+1
			nerr=0
			subfinish = 1
			if (ncomp>200):
				print'too many subiterations'
				return gv_u, gv_v, gv_w, F, u_array, v_array, w_array#need vars
			#endif
			#do i=1,dimx
			#for i in range(dimx):
			#X(i)=XK1(i) + al*D(i)
			X=XK1 + al*D
			#enddo
			#do i=1,nx
			u_array=X[0]
			v_array=X[1]
			w_array=X[2]
			gv_u=G[0]*0.0
			gv_v=G[1]*0.0
			gv_w=G[2]*0.0
			F=0.0
			gv_u,gv_v,gv_w, F = meas_cost_3d(gv_u, gv_v, gv_w, F, u_array, v_array, w_array, i_cmpt_r1, j_cmpt_r1, k_cmpt_r1, i_cmpt_r2, j_cmpt_r2, k_cmpt_r2, vr1, vr2, weights, term_velocity)
			X=array([u_array,v_array,w_array])
			G=array([gv_u,gv_v, gv_w])
			gradi=sqrt((G**2).sum())
			#!----------------------COST FUNCTION--------------------
			if (nerr==1):
				al=al/2
				subfinish = 0
			#endif
			if (subfinish==1):
				gald=0.0
				#do i=1,dimx
				#for i in range(dimx):
				gald=gald+(G*D).sum()
				#enddo
				if ((gald < 0.0) and (Fal < F)):
					al = 4*al
					subfinish = 0
				#endif
				if (subfinish==1):
					z = 3.0*(F - Fal)/al + gd + gald
					w = sqrt(z**2 - gd*gald)
					amin = al*(1 - (gald + w - z)/(gald - gd + 2*w))
					#do i=1,dimx
					#for i in range(dimx):
					#X(i)=XK1(i) + amin*D(i)
					X=XK1 + amin*D
					#enddo
					#!----------------------COST FUNCTION--------------------
					#!pack variables
					#do i=1,nx
					u_array=X[0]
					v_array=X[1]
					w_array=X[2]
					gv_u=G[0]*0.0
					gv_v=G[1]*0.0
					gv_w=G[2]*0.0
					F=0.0
					gv_u,gv_v,gv_w,F = meas_cost_3d(gv_u, gv_v,gv_w, F, u_array, v_array, w_array, i_cmpt_r1, j_cmpt_r1, k_cmpt_r1, i_cmpt_r2, j_cmpt_r2, k_cmpt_r1, vr1, vr2, weights, term_velocity)
					#!unpack variables
					#do i=1,nx
					X=array([u_array,v_array,w_array])
					G=array([gv_u,gv_v,gv_w])
					gradi=sqrt((G**2).sum())
					#!----------------------COST FUNCTION--------------------         
					if (F > Fm1) or (nerr==1):
						al=al/2.5
						subfinish = 0
					#endif
					if (subfinish==1):
						i_ter=i_ter+1
						grad=gg**0.5
					#endif
				#endif
			#endif
		#enddo
	#enddo
	#do i=1,nx
	u_array=X[0]
	v_array=X[1]
	w_array=X[2]
	gv_u=G[0]
	gv_v=G[1]
	gv_w=G[3]
	print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
	return gv_u, gv_v,gv_w, F, u_array, v_array, w_array

def gracon_3d_packaged(X ,radar1, radar2, weights, sounding, **kwargs):
	#!dimensions
	submask=met.make_submask(weights)
	nx, ny, nz=radar1['CZ'].shape
	Fal=0.0
	i_ter=1
	grad=1.0
	nerr=0
	dimx=2*nx*ny*nz
	F=0.0
	acc=0.0
	itr=kwargs.get('itr',100)
	#!Initial guess packaging
	#do i=1,nx
	print "Doing transfer", nx, ny, nz
	#X=array([u_array,v_array, w_array])
	#gradi=sqrt((G**2).sum())
	print "entering loop"
	while ((i_ter <= itr+1) and (grad > acc)):
		if (i_ter==1):
			Fm1=F
			#!----------------------COST FUNCTION--------------------
			G,F,X=cost_function(X ,radar1, radar2, weights, sounding, submask, loud=True)
			gradi=sqrt((G**2).sum())
			#!----------------------COST FUNCTION--------------------      
		#!subroutine rms2(iterc, f, gradi)
		print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
		c_X=copy(X)
		c_F=F
		#call rms2(iter,F,gradi)        
		#!     -------- To stop in the first iteration -------------------------
		if (itr==0):
			grad = -1.0
			print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
			return G,F,X#need to put vars in
		#!     -----------------------------------------------------------------
		#do i=1,dimx
		#for i in range(dimx):
		XK1=X
		#XK1(i)=X(i)
		#!     ===============================================================
		#!     ------------------------Steepest descent-----------------------
		if (i_ter==1):
			gg=0.0
			#do i=1,dimx
			#for i in range(dimx):
			gg=gg+(G**2).sum()#G[i]*G[i]
			D=-G
			GK1=G
			gd=-gg
		else:
			#!    ----------------------Conjugate Gradient-----------------------
			gg=0.0
			upk=0.0
			dok=0.0
			gg=gg+(G**2).sum()
			upk=upk+(G*(G-GK1)).sum()
			dok=dok+(D*(G-GK1)).sum()
			gd=0.0
			#do i=1,dimx
			#for i in range(dimx):
			GK1=G
			if (dok==0.0): dok = 1.0e-10
         		if (upk==0.0): upk = 1.0e-10
         		#D(i)=-G(i)+(upk/dok)*D(i)
			D=-G+(upk/dok)*D
			#!<<<alain
			#!dd=dd+D(i)*D(i)
			#!>>>alain
         		gd=gd+(G*D).sum()
			#!            print*,'Angle between search direction and gradient:',180.0/(4.0*atan(1.0))*acos(-gd/sqrt(dd*gg))
			#!     --------------------------------------------------------------
		#endif
		#!     ===============================================================
		#!     ----------------------- step length ---------------------------
		if (i_ter==1):
			delF=0.5*F
		elif (i_ter!=1):
			 delF=Fm1-F
		if (F!=0.0): rdF=delF/F
		if (rdF < 1.0e-5):
			print'variation of F too weak'
			return G,F,X#need to put vars in
		if (gd > 0.0): gd=-gd
		axpect = -2.0*delF/gd
		al = min(2.0,axpect)
		Fm1=F
		ncomp=0
		subfinish = 0
		#do while (subfinish.eq.0)
		while (subfinish==0):
			#print "in subloop"
			ncomp=ncomp+1
			nerr=0
			subfinish = 1
			if (ncomp>200):
				print 'too many subiterations'
				print "Current search point: Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
				print "Begining of Subloop Cost: ", c_F
				if c_F < F:
					print "Using start of subloop"
					X=c_X
					F=c_F
				return G,F,X
			#endif
			#do i=1,dimx
			#for i in range(dimx):
			#X(i)=XK1(i) + al*D(i)
			X=XK1 + al*D
			#enddo
			#do i=1,nx
			G,F,X=cost_function(X ,radar1, radar2, weights, sounding, submask)
			gradi=sqrt((G**2).sum())
			#!----------------------COST FUNCTION--------------------
			if (nerr==1):
				al=al/2.
				subfinish = 0
			#endif
			if (subfinish==1):
				gald=0.0
				#do i=1,dimx
				#for i in range(dimx):
				gald=gald+(G*D).sum()
				#enddo
				if ((gald < 0.0) and (Fal < F)):
					al = 4.0*al
					subfinish = 0
				#endif
				if (subfinish==1):
					z = 3.0*(F - Fal)/al + gd + gald
					w = sqrt(z**2 - gd*gald)
					amin = al*(1.0 - (gald + w - z)/(gald - gd + 2.0*w))
					#do i=1,dimx
					#for i in range(dimx):
					#X(i)=XK1(i) + amin*D(i)
					X=XK1 + amin*D
					#enddo
					#!----------------------COST FUNCTION--------------------
					#!pack variables
					#do i=1,nx
					G,F,X=cost_function(X ,radar1, radar2, weights, sounding, submask)
					gradi=sqrt((G**2).sum())
					#!----------------------COST FUNCTION--------------------         
					if (F > Fm1) or (nerr==1):
						al=al/2.5
						subfinish = 0
					#endif
					if (subfinish==1):
						i_ter=i_ter+1
						grad=gg**0.5
					#endif
				#endif
			#endif
		#enddo
	#enddo
	#do i=1,nx
	print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
	return G,F,X



def cost_function(X, radar1, radar2, weights, sounding, submask,loud=False):
	vel='VE'
	#t0=systime()
	costs=zeros(3,dtype=float)
	#unpack X
	gv_u=zeros(radar1['CZ'].shape, dtype=float)
	gv_v=zeros(radar1['CZ'].shape, dtype=float)
	gv_w=zeros(radar1['CZ'].shape, dtype=float)
	u_array=X[0]
	v_array=X[1]
	w_array=X[2]
	
	#unpack radar data and calculate info
	#Soudings
	tdry=sounding['tdry(degs)']
	pressure=sounding['press(hPa)']
	#radar1
	levs=radar1['levs']
	dx=radar1['xar'][1]-radar1['xar'][0]
	dy=radar1['yar'][1]-radar1['yar'][0]
	i_cmpt_r1=radar1['i_comp']
	j_cmpt_r1=radar1['j_comp']
	k_cmpt_r1=radar1['k_comp']
	vr1=radar1[vel]
	
	#radar2
	i_cmpt_r2=radar2['i_comp']
	j_cmpt_r2=radar2['j_comp']
	k_cmpt_r2=radar2['k_comp']
	vr2=radar2[vel]
	
	reflectivity=radar2['CZ']
	#print "Packing finished ",systime()-t0 
	#t1=systime()
	#terminal velocity calculation
	term_velocity=met.terminal_velocity(reflectivity*weights, tdry, levs)
	#print "terminal velocities calculated ", systime()-t1
	#t2=systime() 
	#calculate 
	#gv_w, costs[0], w_array=continuity_cost(gv_w, u_array, v_array, dx, dy, levs, pressure, tdry, weights, submask)
	gv_w, costs[0]=continuity_cost2(gv_w, u_array, v_array,w_array, dx, dy, levs, pressure, tdry, weights, submask)
	#print "Vertical velocities done ", systime()-t2
	#t3=systime()
	gv_u, gv_v, gv_w, costs[1]=meas_cost_3d(gv_u, gv_v, gv_w, 0.0, u_array, v_array, w_array, i_cmpt_r1, j_cmpt_r1, k_cmpt_r1, i_cmpt_r2, j_cmpt_r2, k_cmpt_r2, vr1, vr2, weights, term_velocity)
	#print "VR cost ", systime()-t3
	#t4=systime()
	#gpsy,f = smoothing(gpsy,f,psy,aph,dixy,alt,nx=shape(gpsy,0),ny=shape(gpsy,1),nz=shape(gpsy,2))
	smooth_wt=0.05
	gv_u, costs[2]=smooth_cost.smoothing(gv_u, costs[2], u_array, smooth_wt, dx, levs)
	gv_v, costs[2]=smooth_cost.smoothing(gv_v, costs[2], v_array, smooth_wt, dx, levs)
	#gv_w_pc, costs[2]=smooth_cost.smoothing(gv_w, costs[2], w_array, smooth_wt, dx, levs)
	#print "smoothness done ", systime()-t4
	#package gradients
	G=[gv_u, gv_v, gv_w]
	X=[u_array, v_array, w_array]
	F=costs.sum()
	if loud: print "Costs:  Continuity: ", costs[0], " Vr measurement error: ", costs[1], " smoothness cost: ", costs[2]
	#print "costs  function done", systime()-t0
	return  array(G), F, array(X)

def return_cost(X, radar1, radar2, weights, sounding, submask,loud=False):
	vel='VE'
	#t0=systime()
	costs=zeros(3,dtype=float)
	#unpack X
	gv_u=zeros(radar1['CZ'].shape, dtype=float)
	gv_v=zeros(radar1['CZ'].shape, dtype=float)
	gv_w=zeros(radar1['CZ'].shape, dtype=float)
	u_array=X[0]
	v_array=X[1]
	w_array=X[2]
	
	#unpack radar data and calculate info
	#Soudings
	tdry=sounding['tdry(degs)']
	pressure=sounding['press(hPa)']
	#radar1
	levs=radar1['levs']
	dx=radar1['xar'][1]-radar1['xar'][0]
	dy=radar1['yar'][1]-radar1['yar'][0]
	i_cmpt_r1=radar1['i_comp']
	j_cmpt_r1=radar1['j_comp']
	k_cmpt_r1=radar1['k_comp']
	vr1=radar1[vel]
	
	#radar2
	i_cmpt_r2=radar2['i_comp']
	j_cmpt_r2=radar2['j_comp']
	k_cmpt_r2=radar2['k_comp']
	vr2=radar2[vel]
	
	reflectivity=radar2['CZ']
	#print "Packing finished ",systime()-t0 
	#t1=systime()
	#terminal velocity calculation
	term_velocity=met.terminal_velocity(reflectivity*weights, tdry, levs)
	#print "terminal velocities calculated ", systime()-t1
	#t2=systime() 
	#calculate 
	#gv_w, costs[0], w_array=continuity_cost(gv_w, u_array, v_array, dx, dy, levs, pressure, tdry, weights, submask)
	gv_w, costs[0]=continuity_cost2(gv_w, u_array, v_array,w_array, dx, dy, levs, pressure, tdry, weights, submask)
	#print "Vertical velocities done ", systime()-t2
	#t3=systime()
	gv_u, gv_v, gv_w, costs[1]=meas_cost_3d(gv_u, gv_v, gv_w, 0.0, u_array, v_array, w_array, i_cmpt_r1, j_cmpt_r1, k_cmpt_r1, i_cmpt_r2, j_cmpt_r2, k_cmpt_r2, vr1, vr2, weights, term_velocity)
	#print "VR cost ", systime()-t3
	#t4=systime()
	#gpsy,f = smoothing(gpsy,f,psy,aph,dixy,alt,nx=shape(gpsy,0),ny=shape(gpsy,1),nz=shape(gpsy,2))
	smooth_wt=0.05
	gv_u, costs[2]=smooth_cost.smoothing(gv_u, costs[2], u_array, smooth_wt, dx, levs)
	gv_v, costs[2]=smooth_cost.smoothing(gv_v, costs[2], v_array, smooth_wt, dx, levs)
	#gv_w_pc, costs[2]=smooth_cost.smoothing(gv_w, costs[2], w_array, smooth_wt, dx, levs)
	#print "smoothness done ", systime()-t4
	#package gradients
	G=[gv_u, gv_v, gv_w]
	X=[u_array, v_array, w_array]
	F=costs.sum()
	if loud: print "Costs:  Continuity: ", costs[0], " Vr measurement error: ", costs[1], " smoothness cost: ", costs[2]
	#print "costs  function done", systime()-t0
	return  costs


def continuity_cost2(gv_w, u_array, v_array, w_array, dx, dy, levs, press, temps_c, weights, submask):
     w_array_up, w_array_down=met.w_from_continuity(u_array, v_array, dx, dy, levs, press, temps_c,submask)
     w_from_cont=zeros(w_array.shape, dtype=float)
     upw=(levs[-1]-levs)/levs[-1]
     dow=(levs)/levs[-1]
     for k in range(len(levs)):
	  w_from_cont[:,:,k]=upw[k]*w_array_up[:,:,k] + dow[k]*w_array_down[:,:,k]
     cont_cost=0.5*((w_array-w_from_cont)**2).sum()
     gv_w=gv_w+(w_array-w_from_cont)
     return gv_w, cont_cost
     
     


def meas_cost_3d(gv_u, gv_v, gv_w, F, u_array, v_array, w_array, i_cmpt_r1, j_cmpt_r1, k_cmpt_r1, i_cmpt_r2, j_cmpt_r2, k_cmpt_r2, vr1, vr2, weights, term_velocity):
       vr1_sim=i_cmpt_r1*u_array+j_cmpt_r1*v_array + k_cmpt_r1*(w_array + term_velocity)
       vr2_sim=i_cmpt_r2*u_array+j_cmpt_r2*v_array + k_cmpt_r2*(w_array + term_velocity)
       cost=0.5*(weights*((vr1-vr1_sim)**2 + (vr2-vr2_sim)**2)).sum()
       gv_u=gv_u + weights*(i_cmpt_r1*(vr1_sim-vr1) + i_cmpt_r2*(vr2_sim-vr2))
       gv_v=gv_v+ weights*(j_cmpt_r1*(vr1_sim-vr1) + j_cmpt_r2*(vr2_sim-vr2))
       gv_w=gv_w+ weights*(k_cmpt_r1*(vr1_sim-vr1) + k_cmpt_r2*(vr2_sim-vr2))
       return gv_u,gv_v, gv_w, cost
	