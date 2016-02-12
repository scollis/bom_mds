subroutine vel_2d_cost(gv_u, gv_v, cost, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, nx, ny, vr1, vr2)
implicit none

!determine the cost function and gradient of the cost function for a very simple (initially simulated) wind field
!this program will be nested inside a gradient conjigate solver, which in turn will be nested in a Python program
!using f2py  
!Scott Collis, CAWCR, May 2008

!variables fed in and out of the program
!Dimensions
integer:: nx,ny
!control variables
real(kind=8),dimension(nx,ny):: u_array, v_array
!ray unit vectors for projecting back to v_r for radar 1 and 2, not for feeding out, 
!THIS is what distinguishes our two radars 
real(kind=8), dimension(nx,ny):: i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2
!Measurements, radial velocities as measured by the radars
real(kind=8), dimension(nx,ny):: vr1, vr2
!Gradients in the cost function with respect to the controls
real(kind=8), dimension(nx,ny):: gv_u, gv_v
!The cost function to minimise
real(kind=8):: cost

!Local variables which do not get moved in or out
!counting vars
integer::i,j
!stepwise calculated radial velocities for U and V as seen by the two radars
real(kind=8):: vr1_sim, vr2_sim, p1, p2

!Python bindings
!f2py intent(in) gv_u, gv_v, cost, u_array, v_array, i_cmpt_r1, j_cmpt_r1, i_cmpt_r2, j_cmpt_r2, nx, ny, vr1, vr2
!f2py intent(out) gv_u, gv_v, cost

!Now we actually set to work 
print*, 'Running velocity cost'
do i=1, nx
  do j=1, ny
    vr1_sim=i_cmpt_r1(i,j)*u_array(i,j)+j_cmpt_r1(i,j)*v_array(i,j)
    vr2_sim=i_cmpt_r2(i,j)*u_array(i,j)+j_cmpt_r2(i,j)*v_array(i,j)
    cost=cost+0.5*((vr1(i,j)-vr1_sim)**2 + (vr2(i,j)-vr2_sim)**2)
    p1=(i_cmpt_r1(i,j)*j_cmpt_r1(i,j)*v_array(i,j) + u_array(i,j)*i_cmpt_r1(i,j)**2 - vr1(i,j)*i_cmpt_r1(i,j))
    p2=(i_cmpt_r2(i,j)*j_cmpt_r2(i,j)*v_array(i,j) + u_array(i,j)*i_cmpt_r2(i,j)**2 - vr2(i,j)*i_cmpt_r2(i,j))
    gv_u(i,j)=gv_u(i,j) + p1 + p2
    p1=(i_cmpt_r1(i,j)*j_cmpt_r1(i,j)*u_array(i,j) + v_array(i,j)*j_cmpt_r1(i,j)**2 - vr1(i,j)*j_cmpt_r1(i,j))
    p2=(i_cmpt_r2(i,j)*j_cmpt_r2(i,j)*u_array(i,j) + v_array(i,j)*j_cmpt_r2(i,j)**2 - vr2(i,j)*j_cmpt_r2(i,j))
    gv_v(i,j)=gv_v(i,j)+ p1 + p2
   enddo
enddo
print*,'cost is ', cost 
return
end subroutine vel_2d_cost