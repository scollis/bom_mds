      subroutine con_v3_mc2_com(gvz, F, vx, vy, vz, poidw, dixy, alt_t,
     $           mag, magzb, magzt, nx, ny, nz, nt, nc)
      implicit none

c     Weighted upward and downward integrations of the continuity equation

      integer mag,magzb,magzt,nx,ny,nz, nt, nc
      real F
      real dixy, alt_t(nz)
      real vx(nx,ny,nz,nt-1), vy(nx,ny,nz,nt-1), vz(nx,ny,nz,nt-1)
      real gvz(nx,ny,nz,nt-1), poidw(nx,ny,nz,nt-1)

      integer i,j,k
      real ax,fac,sigw, psup
      real vz2(nx,ny,nz), gvz2(nx,ny,nz,nt-1)
      real div(nx,ny,nz)

c
c*** Initialize buffers
c
      do k=magzb+1,nz-magzt
         do j=mag+1,ny-mag-1
            do i=mag+1,nx-mag-1
               vz2(i,j,k) = 0.0
               vz(i,j,k,nc) = 0.0
               gvz2(i,j,k,nc)=0.0
            enddo
         enddo
      enddo
c
c** Calculate divergence
c
      do k=magzb+1,nz-magzt
         do j=mag+1,ny-mag-1
            do i=mag+1,nx-mag-1
               div(i,j,k) = vx(i+1,j,k,nc) - vx(i,j,k,nc) +
     $                      vy(i,j+1,k,nc) - vy(i,j,k,nc)
            enddo
         enddo
      enddo
c     
c     *** UPWARD INTEGRATION ***
c     
      k=magzb+1
      ax = alt_t(k)/dixy
      do j=mag+1,ny-mag-1
         do i=mag+1,nx-mag-1
            vz(i,j,k,nc)= -ax*div(i,j,k)
         enddo
      enddo

      do k=magzb+2,nz-magzt
         fac = 1.0 + (alt_t(k)-alt_t(k-1)) * 1.0e-4
         ax = (alt_t(k)-alt_t(k-1))/dixy
         do j=mag+1,ny-mag-1
            do i=mag+1,nx-mag-1
               vz(i,j,k,nc)=vz(i,j,k-1,nc)*fac - ax*div(i,j,k-1)
            enddo
         enddo
      enddo

c     
c     *** DOWNWARD INTEGRATION ***
c     
      k=nz-magzt
      fac = 1.0 + (alt_t(k)-alt_t(k-1)) * 1.0e-4
      ax = (alt_t(k)-alt_t(k-1))/dixy
      do j=mag+1,ny-mag-1
         do i=mag+1,nx-mag-1
            vz2(i,j,k) = ax*div(i,j,k)/fac
         enddo
      enddo

      do k=nz-magzt-1,magzb+1,-1
         fac = 1.0 + (alt_t(k+1)-alt_t(k)) * 1.0e-4
         ax = (alt_t(k+1)-alt_t(k))/dixy
         do j=mag+1,ny-mag-1
            do i=mag+1,nx-mag-1
               vz2(i,j,k)=( vz2(i,j,k+1) + ax*div(i,j,k) ) / fac
            enddo
         enddo
      enddo

      sigw = 3.0                ! m/s
      psup = (sigw)**(-2.0)

      do k=magzb+1,nz-magzt
         do j=mag+1,ny-mag-1
            do i=mag+1,nx-mag-1
               F = F + 0.5*psup*(vz(i,j,k,nc)-vz2(i,j,k))**2
               gvz2(i,j,k,nc)=gvz2(i,j,k,nc)-psup*(vz(i,j,k,nc)-vz2(i,j,k))
               gvz(i,j,k,nc)=gvz(i,j,k,nc)+psup*(vz(i,j,k,nc)-vz2(i,j,k))
            enddo
         enddo
      enddo      

      return
      end
