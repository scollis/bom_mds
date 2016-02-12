      subroutine smoothing(GPSY, F, psy, aph, dixy, alt, 
     $                     nx, ny, nz)
      implicit none

c     Pascal Guillemette Tue Feb 23 1999
c     Cathy Chiang Jan 15 2003
c     modified for python Scott Collis June 2008

c     Parameters
      integer nx,ny,nz
      real gpsy(nx,ny,nz)
      real psy(nx,ny,nz)
      real F, aph, aphz

      real dixy, alt(nz)

c     Local variables
      integer i,j,k
      real dxz
      real buf,buf1
      integer tot
cf2py intent(in) GPSY, F, psy, aph, dixy, alt, nx, ny, nz
cf2py intent(out) GPSY, F
      aphz = aph

c     ---------------- Smoothness constraint -------------

c     *** DXX ***
      do k=1,nz
	 do j=1,ny
	    do i=2,nx-1
               if (psy(i-1,j,k).ne.0.0.and.psy(i+1,j,k).ne.0.0.and.
     $              psy(i,j,k).ne.0.0) then

                  buf=psy(i-1,j,k)+psy(i+1,j,k)-2*psy(i,j,k)

                  F = F + 0.5*aph*(buf**2)

                  buf = aph*buf

                  gpsy(i-1,j,k) = gpsy(i-1,j,k) + buf
                  gpsy(i+1,j,k) = gpsy(i+1,j,k) + buf
                  gpsy(i,j,k)   = gpsy(i,j,k)   - 2*buf

               endif
	    enddo
	 enddo
      enddo

c     *** DYY ***
      do k=1,nz
	 do j=2,ny-1
	    do i=1,nx
               if (psy(i,j-1,k).ne.0.0.and.psy(i,j+1,k).ne.0.0.and.
     $              psy(i,j,k).ne.0.0) then

                  buf=psy(i,j-1,k)+psy(i,j+1,k)-2*psy(i,j,k)

                  F =F + 0.5*aph*(buf**2)

                  buf = aph*buf
 
                  gpsy(i,j-1,k) = gpsy(i,j-1,k) + buf
                  gpsy(i,j+1,k) = gpsy(i,j+1,k) + buf
                  gpsy(i,j,k)   = gpsy(i,j,k)   - 2*buf

               endif
	    enddo
	 enddo
      enddo

c     *** DXY ***
      do k=1,nz
	 do j=1,ny-1
	    do i=1,nx-1
               if (psy(i+1,j+1,k).ne.0.0.and.psy(i+1,j,k).ne.0.0.and.
     $              psy(i,j+1,k).ne.0.0.and.psy(i,j,k).ne.0.0) then

                  buf=psy(i+1,j+1,k)+psy(i,j,k)-
     $                 psy(i+1,j,k)-psy(i,j+1,k)

                  F =F + 0.5*aph*(buf**2)

                  buf = aph*buf
 
                  gpsy(i+1,j+1,k) = gpsy(i+1,j+1,k) + buf
                  gpsy(i,j,k) = gpsy(i,j,k) + buf
                  gpsy(i+1,j,k) = gpsy(i+1,j,k) - buf
                  gpsy(i,j+1,k) = gpsy(i,j+1,k) - buf

               endif
	    enddo
	 enddo
      enddo

      if (aphz.ne.0.0) then

c     *** DZZ ***
      do k=2,nz-1
         buf1=2.0*(dixy**2)/(alt(k+1)-alt(k-1))
	 do j=1,ny
	    do i=1,nx
               if (psy(i,j,k-1).ne.0.0.and.psy(i,j,k+1).ne.0.0.and.
     $              psy(i,j,k).ne.0.0) then
                   
                  buf=buf1*((psy(i,j,k+1)-psy(i,j,k))/(alt(k+1)-alt(k))-
     $                 (psy(i,j,k)-psy(i,j,k-1))/(alt(k)-alt(k-1)) )

                  F =F + 0.5*aphz*(buf**2)

                  buf = aphz*buf*buf1

                  gpsy(i,j,k-1)= gpsy(i,j,k-1) + buf/(alt(k)-alt(k-1))
                  gpsy(i,j,k+1)= gpsy(i,j,k+1) + buf/(alt(k+1)-alt(k))
                  gpsy(i,j,k)  = gpsy(i,j,k)   - buf/(alt(k+1)-alt(k))
                  gpsy(i,j,k)  = gpsy(i,j,k)   - buf/(alt(k)-alt(k-1))

               endif
	    enddo
	 enddo
      enddo

c     *** DXZ ***
      do k=1,nz-1
         dxz = (alt(k+1)-alt(k))/dixy
	 do j=1,ny
	    do i=1,nx-1
               if (psy(i+1,j,k+1).ne.0.0.and.psy(i+1,j,k).ne.0.0.and.
     $              psy(i,j,k+1).ne.0.0.and.psy(i,j,k).ne.0.0) then

                  buf=(psy(i+1,j,k+1)+psy(i,j,k)-
     $                 psy(i+1,j,k)-psy(i,j,k+1))/dxz

                  F = F + 0.5*aphz*(buf**2)

                  buf = aphz*buf/dxz

                  gpsy(i+1,j,k+1) = gpsy(i+1,j,k+1) + buf
                  gpsy(i,j,k) = gpsy(i,j,k) + buf
                  gpsy(i+1,j,k) = gpsy(i+1,j,k) - buf
                  gpsy(i,j,k+1) = gpsy(i,j,k+1) - buf

               endif
	    enddo
	 enddo
      enddo


c     *** DYZ ***
      do k=1,nz-1
         dxz = (alt(k+1)-alt(k))/dixy
	 do j=1,ny-1
	    do i=1,nx
               if (psy(i,j+1,k+1).ne.0.0.and.psy(i,j+1,k).ne.0.0.and.
     $              psy(i,j,k+1).ne.0.0.and.psy(i,j,k).ne.0.0) then

                  buf=(psy(i,j+1,k+1)+psy(i,j,k)-
     $                 psy(i,j+1,k)-psy(i,j,k+1))/dxz

                  F = F + 0.5*aphz*(buf**2)

                  buf = aphz*buf/dxz

                  gpsy(i,j+1,k+1) = gpsy(i,j+1,k+1) + buf
                  gpsy(i,j,k) = gpsy(i,j,k) + buf
                  gpsy(i,j+1,k) = gpsy(i,j+1,k) - buf
                  gpsy(i,j,k+1) = gpsy(i,j,k+1) - buf

               endif
	    enddo
	 enddo
      enddo

      endif                     ! aphz

      return
      end
