subroutine gracon_vel3d(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, nx, ny, vr1, vr2, weights)
implicit none

!variables in and out

!dimensions
integer:: nx, ny, nl
!inputs
!gradients (initially zeros), first guess velocities
real(kind=8), dimension(nx, ny, nl):: gv_u, gv_v, u_array, v_array
!Observations
real(kind=8), dimension(nx,ny):: vr1, vr2, weights
!unit vectors for radar 1 and radar 2
real(kind=8), dimension(nx,ny):: i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2

!F=cost (initially zero)
real(kind=8):: F, Fm1 

!local variables
!Packing and unpacking variables
real(kind=8), dimension(2*nx*ny):: X, G
!loopers and determiners
integer:: iter, nerr, dimx, i, j, subfinish, ncomp, acc, itr
!informative
real(kind=8):: grad, gradi
!variables used in the gradient conjigate proccess
real(kind=8), dimension(2*nx*ny):: XK1, D, GK1
real(kind=8):: gg, gd, upk, dok, delF, rdF, axpect, al, Fal, gald, z,w, amin

!Python bindings
!f2py intent(in) gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, nx, ny, vr1, vr2, weights
!f2py intent(out) gv_u, gv_v, F, u_array, v_array


!the gradient conjigate solver to minimise the cost function
!this iteratively calls the 2d cost function calculator
!the actually step solver is from McGill university (S.Laroche etal)
!DEPENDS: 
!Scott Collis CAWCR May 2008

!Pack vars

!initial cost from cost function

!run the conjigate solver

!unpack and return
iter=1
grad=1.0
nerr=0
dimx=2*nx*ny
F=0.0
acc=0.0
itr=100


do while ((iter.le.itr+1).and.(grad.gt.acc))
   if (iter.eq.1) then
      Fm1=F
!----------------------COST FUNCTION--------------------
      !pack variables
      do i=1,nx
         do j=1,ny
         !package
            u_array(i,j)=X(i+nx*(j-1))
            v_array(i,j)=X(i+nx*(j-1)+nx*ny)
            gv_u(i,j)=0.
            gv_v(i,j)=0.
            F=0.
         enddo
      enddo
      call vel_2d_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, nx, ny, vr1, vr2, weights)
      !unpack variables
      gradi=0.0
      do i=1,nx
         do j=1,ny
            !unpackage
            X(i+nx*(j-1))=u_array(i,j)
            X(i+nx*(j-1) + nx*ny)=v_array(i,j)
            G(i+nx*(j-1))=gv_u(i,j)
            gradi=gradi+G(i+nx*(j-1))**2
            G(i+nx*(j-1) + nx*ny)=gv_v(i,j)
            gradi=gradi+G(i+nx*(j-1) + nx*ny)**2
         enddo
      enddo
      gradi=gradi**0.5
!----------------------COST FUNCTION--------------------      
   endif
!subroutine rms2(iterc, f, gradi)
   call rms2(iter,F,gradi)        
!     -------- To stop in the first iteration -------------------------
   if (itr.eq.0) then 
      grad = -1.0 
      call rms2(iter-1,F,gradi)
      return
   endif
!     -----------------------------------------------------------------
   do i=1,dimx
      XK1(i)=X(i)
   enddo
!     ===============================================================
!     ------------------------Steepest descent-----------------------
   if (iter.eq.1) then
      gg=0.0
      do i=1,dimx
         gg=gg+G(i)*G(i)
         D(i)=-G(i)
         GK1(i)=G(i)
      enddo
      gd=-gg
   else
!    ----------------------Conjugate Gradient-----------------------
      gg=0.0
      upk=0.0
      dok=0.0
      do i=1,dimx
         gg=gg+G(i)*G(i)
         upk=upk+G(i)*(G(i)-GK1(i))
         dok=dok+D(i)*(G(i)-GK1(i))
!===Num.Recipies.               dok=dok+GK1(i)*GK1(i)
      enddo
      gd=0.0
      do i=1,dimx
         GK1(i)=G(i)
         if (dok.eq.0.0) then
            dok = 1.0e-10
         endif
         if (upk.eq.0.0) then
            upk = 1.0e-10
         endif
         D(i)=-G(i)+(upk/dok)*D(i)
!<<<alain
!               dd=dd+D(i)*D(i)
!>>>alain
         gd=gd+G(i)*D(i)
      enddo
!            print*,'Angle between search direction and gradient:',180.0/(4.0*atan(1.0))*acos(-gd/sqrt(dd*gg))
!     --------------------------------------------------------------
   endif
!     ===============================================================
!     ----------------------- step length ---------------------------

   if (iter.eq.1) then
      delF=0.5*F
   elseif (iter.ne.1) then
      delF=Fm1-F
   endif
   if (F.ne.0.0) rdF=delF/F
   if (rdf.lt.1.0e-5) then
      print*,'variation of F too weak'
      return
   endif
   if (gd.gt.0.0) gd=-gd
   axpect = -2*delF/gd
   al = min(2.0,axpect)
   Fm1=F
   ncomp=0
   subfinish = 0
   do while (subfinish.eq.0)
      ncomp=ncomp+1
      nerr=0
      subfinish = 1
      if (ncomp.gt.200) then
         print*,'too many subiterations'
         return
      endif
      do i=1,dimx
         X(i)=XK1(i) + al*D(i)
      enddo
!----------------------COST FUNCTION--------------------
      do i=1,nx
         do j=1,ny
         !package
            u_array(i,j)=X(i+nx*(j-1))
            v_array(i,j)=X(i+nx*(j-1) + nx*ny)
            gv_u(i,j)=0.
            gv_v(i,j)=0.
            F=0.
         enddo
      enddo
      call vel_2d_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, nx, ny, vr1, vr2, weights)
      !unpack variables
      gradi=0.0
      do i=1,nx
         do j=1,ny
            !unpackage
            X(i+nx*(j-1))=u_array(i,j)
            X(i+nx*(j-1) + nx*ny)=v_array(i,j)
            G(i+nx*(j-1))=gv_u(i,j)
            G(i+nx*(j-1) + nx*ny)=gv_v(i,j)
            gradi=gradi+G(i+nx*(j-1))**2 + G(i+nx*(j-1) + nx*ny)**2
         enddo
      enddo
      gradi=gradi**0.5
!----------------------COST FUNCTION--------------------
      if (nerr.eq.1) then
         al=al/2
         subfinish = 0
      endif
      if (subfinish.eq.1) then
         gald=0.0
         do i=1,dimx
            gald=gald+G(i)*D(i)
         enddo
         if ((gald.lt.0.0).and.(Fal.lt.F)) then
            al = 4*al
            subfinish = 0
         endif
         if (subfinish.eq.1) then
            z = 3*(F - Fal)/al + gd + gald
            w = sqrt(z**2 - gd*gald)
            amin = al*(1 - (gald + w - z)/(gald - gd + 2*w))
            do i=1,dimx
               X(i)=XK1(i) + amin*D(i)
            enddo
!----------------------COST FUNCTION--------------------
            do i=1,nx
               do j=1,ny
                  !package
                  u_array(i,j)=X(i+nx*(j-1))
                  v_array(i,j)=X(i+nx*(j-1) + nx*ny)
                  gv_u(i,j)=0.
                  gv_v(i,j)=0.
                  F=0.
               enddo
            enddo
            call vel_2d_cost(gv_u, gv_v, F, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, nx, ny, vr1, vr2, weights)
            !unpack variables
            gradi=0.0
            do i=1,nx
               do j=1,ny
                  !unpackage
                  X(i+nx*(j-1))=u_array(i,j)
                  X(i+nx*(j-1) + nx*ny)=v_array(i,j)
                  G(i+nx*(j-1))=gv_u(i,j)
                  G(i+nx*(j-1) + nx*ny)=gv_v(i,j)
                  gradi=gradi+G(i+nx*(j-1))**2 + G(i+nx*(j-1) + nx*ny)**2
               enddo
            enddo
            gradi=gradi**0.5
!----------------------COST FUNCTION--------------------         
            if (F.gt.Fm1.or.nerr.eq.1) then
               al=al/2.5
               subfinish = 0
            endif
            if (subfinish.eq.1) then
               iter=iter+1
               grad=gg**0.5
            endif
         endif
      endif
   enddo
enddo
do i=1,nx
   do j=1,ny
      !package
      u_array(i,j)=X(i+nx*(j-1))
      v_array(i,j)=X(i+nx*(j-1) + nx*ny)
      gv_u(i,j)=G(i+j)
      gv_v(i,j)=G(i+j+nx+ny)
   enddo
enddo
call rms2(iter-1,F,gradi)
return
end subroutine gracon_vel2d
