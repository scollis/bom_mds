      subroutine vel3_v2_int_mc2(GVXv, GVYv, GVZv, F, vx, vy, vz,
     $     cra, cre, sra, sre, vt_i, vr1_i, poidv, poidv2, poidv3, poidw,
     $     mag, magzb, magzt, nx, ny, nz, nt, psup, nc, k, cover, u0s, v0s)
      implicit none

      integer nx,ny,nz,nt,nc,k, cover
      integer mag,magzb,magzt
      real u0s, v0s
      real gvxv(nx,ny,nz,nt-1),gvyv(nx,ny,nz,nt-1),gvzv(nx,ny,nz,nt-1)
      real F,psup
      real vx(nx,ny,nz,nt-1),vy(nx,ny,nz,nt-1),vz(nx,ny,nz,nt-1)
      real cra(nx,ny,nt-1),cre(nx,ny,nz,nt-1)
      real sra(nx,ny,nt-1),sre(nx,ny,nz,nt-1)
      real vt_i(nx,ny,nz,nt-1), vr1_i(nx,ny,nz,nt-1)
      real poidv(nx,ny,nz,nt-1),poidv2(nx,ny,nz,nt-1),poidv3(nx,ny,nz,nt-1)
      real poidw(nx,ny,nz,nt-1)

c     Local variables

      integer i,j
      real vrm,ce,se,cr,sr,uqm,vqm,wqm,buf

c     Hypothesis : radial velocity is on q point, momentum level
c     The constraint is applied on q points, momentum levels

      do j=mag+1,ny-mag-1
         do i=mag+1,nx-mag-1

            if ( poidv(i,j,k,nc).ne.0.0
     $           .and.(poidv2(i,j,k,nc).ne.0.0.or.
     $           poidv3(i,j,k,nc).ne.0.0.or.cover.eq.1)
     $           ) then

               sr = sra(i,j,nc)
               cr = cra(i,j,nc)
               ce = cre(i,j,k,nc)
               se = sre(i,j,k,nc)

               uqm = (vx(i,j,k,nc) + vx(i+1,j,k,nc))/2.0 + u0s
               vqm = (vy(i,j,k,nc) + vy(i,j+1,k,nc))/2.0 + v0s
               if (k.ne.nz) then
                  wqm = (vz(i,j,k,nc) + vz(i,j,k+1,nc))/2.0
               else
                  wqm = vz(i,j,k,nc)/2.0
               endif

               vrm = uqm*cr*ce+vqm*sr*ce+(wqm+vt_i(i,j,k,nc))*se

               F = F + 0.5*psup*poidv(i,j,k,nc)*
     $              (vrm-vr1_i(i,j,k,nc))**2

               buf = psup*poidv(i,j,k,nc)*
     $              (vrm-vr1_i(i,j,k,nc))/2.0

               gvxv(i,j,k,nc) = gvxv(i,j,k,nc) + buf*cr*ce
               gvxv(i+1,j,k,nc) = gvxv(i+1,j,k,nc) + buf*cr*ce

               gvyv(i,j,k,nc) = gvyv(i,j,k,nc) + buf*sr*ce
               gvyv(i,j+1,k,nc) = gvyv(i,j+1,k,nc) + buf*sr*ce

c            if (poidw(i,j,k,nc).ne.0.0) then
               if (k.ne.nz) then
                  gvzv(i,j,k+1,nc) = gvzv(i,j,k+1,nc) + buf*se
               endif
               gvzv(i,j,k,nc) = gvzv(i,j,k,nc) + buf*se
c            endif

            endif
         enddo
      enddo

      return
      end
