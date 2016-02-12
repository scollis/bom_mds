subroutine vel_2d_3d_cost(gv_u, gv_v, cost, u_array, v_array, i_cmpt_r1, j_cmpt_r1,&
 i_cmpt_r2, j_cmpt_r2, nx, ny, nz, vr1, vr2, weights)
implicit none

!determine the cost function and gradient of the cost function for a very simple (initially simulated) wind field
!this program will be nested inside a gradient conjigate solver, which in turn will be nested in a Python program
!using f2py  
!Scott Collis, CAWCR, May 2008

!variables fed in and out of the program
!Dimensions
integer:: nx,ny,nz
!control variables
real(kind=8),dimension(nx,ny,nz):: u_array, v_array
!ray unit vectors for projecting back to v_r for radar 1 and 2, not for feeding out, 
!THIS is what distinguishes our two radars 
real(kind=8), dimension(nx,ny,nz):: i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, weights
!Measurements, radial velocities as measured by the radars
real(kind=8), dimension(nx,ny,nz):: vr1, vr2
!Gradients in the cost function with respect to the controls
real(kind=8), dimension(nx,ny,nz):: gv_u, gv_v
!The cost function to minimise
real(kind=8):: cost

!Local variables which do not get moved in or out
!counting vars
integer::i,j,k
!stepwise calculated radial velocities for U and V as seen by the two radars
real(kind=8):: vr1_sim, vr2_sim, p1, p2

!Python bindings
!f2py intent(in) gv_u, gv_v, cost, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, nx, ny, nz, vr1, vr2, weights
!f2py intent(out) gv_u, gv_v, cost

!Now we actually set to work 
!print*, 'Running velocity cost'
do i=1,nx
  do j=1,ny
    do k=1,nz
       vr1_sim=i_cmpt_r1(i,j,k)*u_array(i,j,k)+j_cmpt_r1(i,j,k)*v_array(i,j,k)
       vr2_sim=i_cmpt_r2(i,j,k)*u_array(i,j,k)+j_cmpt_r2(i,j,k)*v_array(i,j,k)
       cost=cost+0.5*weights(i,j,k)*((vr1(i,j,k)-vr1_sim)**2 + (vr2(i,j,k)-vr2_sim)**2)
       !p1=(i_cmpt_r1(i,j,k)*j_cmpt_r1(i,j,k)*v_array(i,j,k) + u_array(i,j,k)*i_cmpt_r1(i,j,k)**2 - vr1(i,j,k)*i_cmpt_r1(i,j,k))
       !p2=(i_cmpt_r2(i,j,k)*j_cmpt_r2(i,j,k)*v_array(i,j,k) + u_array(i,j,k)*i_cmpt_r2(i,j,k)**2 - vr2(i,j,k)*i_cmpt_r2(i,j,k))
       gv_u(i,j,k)=gv_u(i,j,k) + weights(i,j,k)*(i_cmpt_r1(i,j,k)*(vr1_sim-vr1(i,j,k)) + i_cmpt_r2(i,j,k)*(vr2_sim-vr2(i,j,k)))
       !p1=(i_cmpt_r1(i,j,k)*j_cmpt_r1(i,j,k)*u_array(i,j,k) + v_array(i,j,k)*j_cmpt_r1(i,j,k)**2 - vr1(i,j,k)*j_cmpt_r1(i,j,k))
       !p2=(i_cmpt_r2(i,j,k)*j_cmpt_r2(i,j,k)*u_array(i,j,k) + v_array(i,j,k)*j_cmpt_r2(i,j,k)**2 - vr2(i,j,k)*j_cmpt_r2(i,j,k))
       gv_v(i,j,k)=gv_v(i,j,k)+ weights(i,j,k)*(j_cmpt_r1(i,j,k)*(vr1_sim-vr1(i,j,k)) + j_cmpt_r2(i,j,k)*(vr2_sim-vr2(i,j,k)))
    enddo
  enddo
enddo
!print*,'cost is ', cost 
return
end subroutine vel_2d_3d_cost