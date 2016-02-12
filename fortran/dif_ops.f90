subroutine grad_inter(x,y,m,b)
implicit none

!Scott Collis, CAWCR, May 2008
real(kind=8), dimension(2):: x,y
real(kind=8):: m,b

!Python bindings
!f2py intent(in) x,y
!f2py intent(out) m,b
m=(y(1)-y(2))/(x(1)-x(2))
b=y(1)-m*x(1)
return
end subroutine grad_inter

subroutine extrap(x, y, nx, my_extrap)
implicit none
real(kind=8), dimension(2):: x,y
real(kind=8):: m,b,my_extrap
real(kind=8)::nx
!f2py intent(in) x,y,nx
!f2py intent(out) my_extrap
call grad_inter(x,y,m,b)
!print *, m
!print *, b
!print *, nx
my_extrap=m*nx+b
!print *, my_extrap
return
end subroutine extrap


subroutine dy(y,ny, diffy)
implicit none

real(kind=8), dimension(ny+2):: yext
real(kind=8), dimension(ny):: y, diffy
real(kind=8):: m1th, np1th
real(kind=8), dimension(2):: st_pos, end_pos, vals_start, vals_end
integer:: ny,i
real(kind=8):: nyr
!Python bindings
!f2py intent(in) y, ny
!f2py intent(out) diffy
nyr=ny*1.0
st_pos(1)=1.0
st_pos(2)=2.0
end_pos(1)=nyr-1.0
end_pos(2)=nyr
vals_start(1)=y(1)
vals_start(2)=y(2)
vals_end(1)=y(ny-1)
vals_end(2)=y(ny)
!print *, nyr

call extrap(st_pos, vals_start, 0.0, m1th)
call extrap(end_pos, vals_end, nyr+1.0, np1th)

yext(1)=m1th
yext(ny+2)=np1th
do i=2,ny+1
   yext(i)=y(i-1)
enddo
!print *, yext
do i=1,ny
   diffy(i)=((yext(i+2)-yext(i+1))+(yext(i+1)-yext(i)))/2.0
enddo
return
end subroutine dy



!def dy(y):
!	#Calculate dy
!	#this procedure is designed to return an array of the same length as y
!	#to do this it we needed to extend y by an element in each direction
!	nelms=len(y)
!	#There is possibly a better way to do this using zeros instead of appending...
!	#seems to run fast enough...
!	mdy=[]
!	pdy=[]
!	
!	#calculate the extrapolations
!	m1th=dif_ops.extrap([0,1], [y[0],y[1]], -1)
!	np1th=dif_ops.extrap([nelms-2,nelms-1],[y[nelms-2],y[nelms-1]], nelms)
!	
!	#append them to the data
!	s1=append([m1th], y)
!	yext=append(s1,[np1th])
!	
!	#calculate the difference using N-1
!	for i in range(nelms):
!		j=i+2
!		dy=yext[j]-yext[j-1]
!		pdy.append(dy)
!		
!		
!	#Calculate the difference using N+1
!	for i in range(nelms):
!		j=i+1
!		dy=yext[j]-yext[j-1]
!		mdy.append(dy)
!		
!	#Use the averate of N+1, N-1 to return dy
!	return (array(mdy)+array(pdy))/2.0
!
!
!def grad_inter(x,y):
!	#given two 1d arrays have a linear relationship 
!	#y=mx+b determine m and b
!	m=(y[0]-y[1])/(x[0]-x[1])
!	b=y[0]-m*x[0]
!	return m,b
!
!
!def extrap(x,y,nx):
!	#A linear extrapolation from two points along nx
!	#x: two x elements
!	#y: two y elements
!	m,b=grad_inter(x,y)
!	my_extrap=m*nx+b
!	return my_extrap
