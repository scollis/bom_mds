import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
from numpy import float, array, zeros, ones, sqrt, min, arange, append

def meas_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, vr1, vr2, weights):
       vr1_sim=i_cmpt_r1*u_array+j_cmpt_r1*v_array
       vr2_sim=i_cmpt_r2*u_array+j_cmpt_r2*v_array
       cost=0.5*(weights*((vr1-vr1_sim)**2 + (vr2-vr2_sim)**2)).sum()
       gv_u=gv_u + weights*(i_cmpt_r1*(vr1_sim-vr1) + i_cmpt_r2*(vr2_sim-vr2))
       gv_v=gv_v+ weights*(j_cmpt_r1*(vr1_sim-vr1) + j_cmpt_r2*(vr2_sim-vr2))
       return gv_u,gv_v,cost


def gracon_vel2d_3d(gv_u,gv_v,f,u_array,v_array,i_cmpt_r1,j_cmpt_r1,i_cmpt_r2,j_cmpt_r2,vr1,vr2,weights):
	#!dimensions
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
