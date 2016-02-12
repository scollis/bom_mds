import os
import sys
sys.path.append(os.getenv('HOME')+'/bom_mds/modules')
sys.path.append(os.getenv('HOME')+'/bom_mds/fortran')
from numpy import float, array, zeros, ones, sqrt, min, arange, append
from gracon_vel2d_3d import vel_2d_3d_cost
def gracon_vel2d_3d(gv_u,gv_v,f,u_array,v_array,i_cmpt_r1,j_cmpt_r1,i_cmpt_r2,j_cmpt_r2,vr1,vr2,weights):
	#!dimensions
	nx, ny, nz=u_array.shape
	#!inputs
	#!gradients (initially zeros), first guess velocities
	#real(kind=8), dimension(nx, ny, nz):: gv_u, gv_v, u_array, v_array
	#!Observations
	#real(kind=8), dimension(nx, ny, nz):: vr1, vr2, weights
	#!unit vectors for radar 1 and radar 2
	#real(kind=8), dimension(nx, ny, nz):: i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2
	#!F=cost (initially zero)
	#real(kind=8):: F, Fm1 
	#!local variables
	#!Packing and unpacking variables
	#real(kind=8), dimension(2*nx*ny*nz):: X, G
	X=zeros([2*nx*ny*nz], dtype=float)
	G=zeros([2*nx*ny*nz], dtype=float)
	#!loopers and determiners
	#integer:: iter, nerr, dimx, i, j,k, subfinish, ncomp, acc, itr
	#NEED TO REPLACE ITER
	#!informative
	#real(kind=8):: grad, gradi
	#!variables used in the gradient conjigate proccess
	#real(kind=8), dimension(2*nx*ny*nz):: XK1, D, GK1
	XK1=zeros([2*nx*ny*nz], dtype=float)
	D=zeros([2*nx*ny*nz], dtype=float)
	GK1=zeros([2*nx*ny*nz], dtype=float)
	Fal=0.0
	#real(kind=8):: gg, gd, upk, dok, delF, rdF, axpect, al, Fal, gald, z,w, amin
	#!Python bindings
	#!f2py intent(in) gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, nx, ny,nz, vr1, vr2, weights
	#!f2py intent(out) gv_u, gv_v, F, u_array, v_array
	#!the gradient conjigate solver to minimise the cost function
	#!this iteratively calls the 2d cost function calculator
	#!the actually step solver is from McGill university (S.Laroche etal)
	#!DEPENDS: 
	#!Scott Collis CAWCR May 2008
	#!Pack vars
	#!initial cost from cost function
	#!run the conjigate solver
	#!unpack and return
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
	for i in arange(nx):
		#do j=1,ny
		for j in arange(ny):
			#do k=1,nz
			for k in arange(nz):
				#!package
				#X(i+nx*(j-1)+nx*ny*(k-1))=u_array(i,j,k)
				X[i+nx*j+nx*ny*k]=u_array[i,j,k]
         			#X(i+nx*(j-1)+nx*ny*(k-1) + nx*ny*nz)=v_array(i,j,k)
				X[i+nx*j+nx*ny*k + nx*ny*nz]=v_array[i,j,k]
	#do while ((iter.le.itr+1).and.(grad.gt.acc))
	print "entering loop"
	while ((i_ter <= itr+1) and (grad > acc)):
		if (i_ter==1):
			Fm1=F
			#!----------------------COST FUNCTION--------------------
			#!pack variables
			#do i=1,nx
			for i in range(nx):
				#do j=1,ny
				for j in range(ny):
					#do k=1,nz
					for k in range(nz):
						#!package
						#u_array(i,j,k)=X(i+nx*(j-1)+nx*ny*(k-1))
						u_array[i,j,k]=X[i+nx*j+nx*ny*k]
						#v_array(i,j,k)=X(i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz)
						v_array[i,j,k]=X[i+nx*j+nx*ny*k+nx*ny*nz]
						#gv_u(i,j,k)=0.
						#gv_v(i,j,k)=0.
						gv_u[i,j,k]=0.
						gv_v[i,j,k]=0.
						F=0.
			#call vel_2d_3d_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1,&
			#i_cmpt_r2, j_cmpt_r2, nx, ny,nz, vr1, vr2, weights)
			gv_u,gv_v,F = vel_2d_3d_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, vr1, vr2, weights)
			print "done first cost"
			#!unpack variables
			gradi=0.0
			#do i=1,nx
			for i in range(nx):
				#do j=1,ny
				for j in range(ny):
					#do k=1,nz
					for k in range(nz):
						#!unpackage
						#X(i+nx*(j-1)+nx*ny*(k-1))=u_array(i,j,k)
						#X(i+nx*(j-1)+nx*ny*(k-1) + nx*ny*nz)=v_array(i,j,k)
						#G(i+nx*(j-1)+nx*ny*(k-1))=gv_u(i,j,k)
						#gradi=gradi+G(i+nx*(j-1)+nx*ny*(k-1))**2
						#G(i+nx*(j-1) + nx*ny*(k-1)+ nx*ny*nz)=gv_v(i,j,k)
						#gradi=gradi+G(i+nx*(j-1) + nx*ny*(k-1) + nx*ny*nz)**2
						X[i+nx*j+nx*ny*k]=u_array[i,j,k]
						X[i+nx*j+nx*ny*k + nx*ny*nz]=v_array[i,j,k]
						G[i+nx*j+nx*ny*k]=gv_u[i,j,k]
						gradi=gradi+G[i+nx*j+nx*ny*k]**2
						G[i+nx*j + nx*ny*k+ nx*ny*nz]=gv_v[i,j,k]
						gradi=gradi+G[i+nx*j + nx*ny*k + nx*ny*nz]**2
			gradi=gradi**0.5
			#!----------------------COST FUNCTION--------------------      
		#!subroutine rms2(iterc, f, gradi)
		print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
		#call rms2(iter,F,gradi)        
		#!     -------- To stop in the first iteration -------------------------
		if (itr==0):
			grad = -1.0
			print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
			return #need to put vars in
		#!     -----------------------------------------------------------------
		#do i=1,dimx
		for i in range(dimx):
			XK1[i]=X[i]
			#XK1(i)=X(i)
		#!     ===============================================================
		#!     ------------------------Steepest descent-----------------------
		if (i_ter==1):
			gg=0.0
			#do i=1,dimx
			for i in range(dimx):
				gg=gg+G[i]*G[i]
				D[i]=-G[i]
				GK1[i]=G[i]
			gd=-gg
		else:
			#!    ----------------------Conjugate Gradient-----------------------
			gg=0.0
			upk=0.0
			dok=0.0
			#do i=1,dimx
			for i in range(dimx):
				#gg=gg+G(i)*G(i)
				#upk=upk+G(i)*(G(i)-GK1(i))
				#dok=dok+D(i)*(G(i)-GK1(i))
				#!===Num.Recipies.               dok=dok+GK1(i)*GK1(i)
				gg=gg+G[i]*G[i]
				upk=upk+G[i]*(G[i]-GK1[i])
				dok=dok+D[i]*(G[i]-GK1[i])
			gd=0.0
			#do i=1,dimx
			for i in range(dimx):
				GK1[i]=G[i]
				if (dok==0.0): dok = 1.0e-10
         			if (upk==0.0): upk = 1.0e-10
         			#D(i)=-G(i)+(upk/dok)*D(i)
				D[i]=-G[i]+(upk/dok)*D[i]
				#!<<<alain
				#!dd=dd+D(i)*D(i)
				#!>>>alain
         			gd=gd+G[i]*D[i]
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
			return #need to put vars in
		if (gd > 0.0): gd=-gd
		axpect = -2*delF/gd
		al = min(2.0,axpect)
		Fm1=F
		ncomp=0
		subfinish = 0
		#do while (subfinish.eq.0)
		while (subfinish==0):
			print "in subloop"
			ncomp=ncomp+1
			nerr=0
			subfinish = 1
			if (ncomp>200):
				print'too many subiterations'
				return #need vars
			#endif
			#do i=1,dimx
			#for i in range(dimx):
			#X(i)=XK1(i) + al*D(i)
			X=XK1 + al*D
			#enddo
			#!----------------------COST FUNCTION--------------------
#      do i=1,nx
#       	do j=1,ny
#         !package
#             do k=1,nz
#               u_array(i,j,k)=X(i+nx*(j-1)+nx*ny*(k-1))
#               v_array(i,j,k)=X(i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz)
#               gv_u(i,j,k)=0.
#               gv_v(i,j,k)=0.
#               F=0.
#           enddo
#         enddo
#      enddo
#      call vel_2d_3d_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1,&
#                          j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, nx, ny,nz, vr1, vr2, weights)
#      !unpack variables
#      gradi=0.0
#      do i=1,nx
#         do j=1,ny
#            do k=1,nz
#               !unpackage
#               X(i+nx*(j-1)+nx*ny*(k-1))=u_array(i,j,k)
#               X(i+nx*(j-1)+nx*ny*(k-1) + nx*ny*nz)=v_array(i,j,k)
#               G(i+nx*(j-1)+nx*ny*(k-1))=gv_u(i,j,k)
#               G(i+nx*(j-1) +nx*ny*(k-1) + nx*ny*nz)=gv_v(i,j,k)
#               gradi=gradi+G(i+nx*(j-1)+nx*ny*(k-1))**2 + G(i+nx*(j-1) +nx*ny*(k-1) + nx*ny*nz)**2
#            enddo
#         enddo
#      enddo
#!----------------------COST FUNCTION--------------------
			#!pack variables
			#do i=1,nx
			for i in range(nx):
				#do j=1,ny
				for j in range(ny):
					#do k=1,nz
					for k in range(nz):
						#!package
						#u_array(i,j,k)=X(i+nx*(j-1)+nx*ny*(k-1))
						u_array[i,j,k]=X[i+nx*j+nx*ny*k]
						#v_array(i,j,k)=X(i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz)
						v_array[i,j,k]=X[i+nx*j+nx*ny*k+nx*ny*nz]
						#gv_u(i,j,k)=0.
						#gv_v(i,j,k)=0.
						gv_u[i,j,k]=0.
						gv_v[i,j,k]=0.
						F=0.
			#call vel_2d_3d_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1,&
			#i_cmpt_r2, j_cmpt_r2, nx, ny,nz, vr1, vr2, weights)
			gv_u,gv_v,F = vel_2d_3d_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, vr1, vr2, weights)
			#!unpack variables
			gradi=0.0
			#do i=1,nx
			for i in range(nx):
				#do j=1,ny
				for j in range(ny):
					#do k=1,nz
					for k in range(nz):
						#!unpackage
						#X(i+nx*(j-1)+nx*ny*(k-1))=u_array(i,j,k)
						#X(i+nx*(j-1)+nx*ny*(k-1) + nx*ny*nz)=v_array(i,j,k)
						#G(i+nx*(j-1)+nx*ny*(k-1))=gv_u(i,j,k)
						#gradi=gradi+G(i+nx*(j-1)+nx*ny*(k-1))**2
						#G(i+nx*(j-1) + nx*ny*(k-1)+ nx*ny*nz)=gv_v(i,j,k)
						#gradi=gradi+G(i+nx*(j-1) + nx*ny*(k-1) + nx*ny*nz)**2
						X[i+nx*j+nx*ny*k]=u_array[i,j,k]
						X[i+nx*j+nx*ny*k + nx*ny*nz]=v_array[i,j,k]
						G[i+nx*j+nx*ny*k]=gv_u[i,j,k]
						gradi=gradi+G[i+nx*j+nx*ny*k]**2
						G[i+nx*j + nx*ny*k+ nx*ny*nz]=gv_v[i,j,k]
						gradi=gradi+G[i+nx*j + nx*ny*k + nx*ny*nz]**2
			gradi=gradi**0.5
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
					for i in range(nx):
						#do j=1,ny
						for j in range(ny):
							#do k=1,nz
							for k in range(nz):
								#!package
								#u_array(i,j,k)=X(i+nx*(j-1)+nx*ny*(k-1))
								u_array[i,j,k]=X[i+nx*j+nx*ny*k]
								#v_array(i,j,k)=X(i+nx*(j-1)+nx*ny*(k-1)+nx*ny*nz)
								v_array[i,j,k]=X[i+nx*j+nx*ny*k+nx*ny*nz]
								#gv_u(i,j,k)=0.
								#gv_v(i,j,k)=0.
								gv_u[i,j,k]=0.
								gv_v[i,j,k]=0.
								F=0.
					#call vel_2d_3d_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1,&
					#i_cmpt_r2, j_cmpt_r2, nx, ny,nz, vr1, vr2, weights)
					gv_u,gv_v,F = vel_2d_3d_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1,	i_cmpt_r2, j_cmpt_r2, vr1, vr2, weights)
					#!unpack variables
					gradi=0.0
					#do i=1,nx
					#for i in range(nx):
						#do j=1,ny
					#	for j in range(ny):
							#do k=1,nz
					#		for k in range(nz):
								#!unpackage
								#X(i+nx*(j-1)+nx*ny*(k-1))=u_array(i,j,k)
								#X(i+nx*(j-1)+nx*ny*(k-1) + nx*ny*nz)=v_array(i,j,k)
								#G(i+nx*(j-1)+nx*ny*(k-1))=gv_u(i,j,k)
								#gradi=gradi+G(i+nx*(j-1)+nx*ny*(k-1))**2
								#G(i+nx*(j-1) + nx*ny*(k-1)+ nx*ny*nz)=gv_v(i,j,k)
								#gradi=gradi+G(i+nx*(j-1) + nx*ny*(k-1) + nx*ny*nz)**2
					#			X[i+nx*j+nx*ny*k]=u_array[i,j,k]
					#			X[i+nx*j+nx*ny*k + nx*ny*nz]=v_array[i,j,k]
					#			G[i+nx*j+nx*ny*k]=gv_u[i,j,k]
					#			gradi=gradi+G[i+nx*j+nx*ny*k]**2
					#			G[i+nx*j + nx*ny*k+ nx*ny*nz]=gv_v[i,j,k]
					#			gradi=gradi+G[i+nx*j + nx*ny*k + nx*ny*nz]**2
					G=append(gv_u.flatten(), gv_v.flatten())
					gradi=sqrt((G**2).sum())
					#gradi=gradi**0.5
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
	for i in range(nx):
		#do j=1,ny
		for j in range(ny):
			#do k=1,nz
			for k in range(nz):
				#!package
				u_array[i,j,k]=X[i+nx*j+nx*ny*k]
				v_array[i,j,k]=X[i+nx*j + nx*ny*k + nx*ny*nz]
				gv_u[i,j,k]=G[i+nx*j+nx*ny*k]
				gv_v[i,j,k]=G[i+nx*j + nx*ny*k + nx*ny*nz]
			#enddo
		#enddo
	
	print "Iter:", i_ter, " Cost: ", F, " Gradi:", gradi
	return
	#end subroutine gracon_vel2d_3d
