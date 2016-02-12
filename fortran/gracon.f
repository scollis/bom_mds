      subroutine gracon(VX, VY, VZ, PRES, TEMP, RAIN, MC,
     $     u0, au, bu, v0, av, bv, GVX, GVY, GVZ, GVP, GVT,
     $     GVR, GVMC, gvu0, gvau, gvbu, gvv0, gvav, gvbv, G, X, dimx,
     $     cra, cre, sra, sre, cra2, cre2, sra2, sre2, cra3, cre3, sra3,
     $     sre3, vt_i, vr1_i, vr2_i, vr3_i, ref_i, qr, refrac_i,
     $     poidv, poidv2, poidv3, poidr, poidw, prt, prq,
     $     vxb, vyb, vzb, presb, tempb, rainb, mcb, vrp, poidsp, profx, profy, 
     $     vtpr, nombre_rcvr, aph, nmo, acc, itr, mag, magzb, magzt,
     $     dixy, alt, alt_t, dtemp, nxr, nyr, nzr, nx, ny, nz, nt, test,
     $     sol, cover, vort, rx1, ry1, F, stage, u0s, v0s, sig, nsig,
     $     kpass, hwl, vxv, vyv, vzv, presv, tempv, rainv, mcv)
      implicit none

c
c ** AProtat : a mode 10 has been added (12/12/97) for thermo retrieval
c ** AProtat : a mode 11 has been added (10/03/98) for retrieval of temporal derivatives
c ** ACaya : a mode 12 has been added (21/04/98) for thermo retrieval with MC2 equations
c
      integer nx, ny, nz, nt, nsig, dimx, nxr, nyr, nzr
      integer magzb, magzt, mag, test, cover, vort, stage
      integer kpass
      real hwl
      real u0(nz,nt-1),au(nz,nt-1),bu(nz,nt-1)
      real v0(nz,nt-1),av(nz,nt-1),bv(nz,nt-1)
      real gvu0(nz,nt-1),gvau(nz,nt-1),gvbu(nz,nt-1)
      real gvv0(nz,nt-1),gvav(nz,nt-1),gvbv(nz,nt-1)
      real vx(nx,ny,nz,nt-1), vy(nx,ny,nz,nt-1), vz(nx,ny,nz,nt-1)
      real vxv(nx,ny,nz,nt-1), vyv(nx,ny,nz,nt-1), vzv(nx,ny,nz,nt-1)
      real vxb(nx,ny,nz,nt-1), vyb(nx,ny,nz,nt-1), vzb(nx,ny,nz,nt-1)
      real pres(nx,ny,nz,nt-1), temp(nx,ny,nz,nt-1), rain(nx,ny,nz,nt-1)
      real presv(nx,ny,nz,nt-1), tempv(nx,ny,nz,nt-1), rainv(nx,ny,nz,nt-1)
      real presb(nx,ny,nz,nt-1), tempb(nx,ny,nz,nt-1), rainb(nx,ny,nz,nt-1)
      real mc(nx,ny,nz,nt-1), mcb(nx,ny,nz,nt-1), mcv(nx,ny,nz,nt-1)
      real vr1_i(nx,ny,nz,nt-1), vr2_i(nx,ny,nz,nt-1), vr3_i(nx,ny,nz,nt-1)
      real vt_i(nx,ny,nz,nt-1), ref_i(nx,ny,nz,nt-1), refrac_i(nx,ny,(nt-1))
      real gvx(nx,ny,nz,nt-1), gvy(nx,ny,nz,nt-1), gvz(nx,ny,nz,nt-1)
      real gvp(nx,ny,nz,nt-1), gvt(nx,ny,nz,nt-1), gvr(nx,ny,nz,nt-1)
      real gvmc(nx,ny,nz,nt-1)
      real poidw(nx,ny,nz,nt-1),sig(nsig,nz)
      real poidv(nx,ny,nz,nt-1), poidv2(nx,ny,nz,nt-1), poidv3(nx,ny,nz,nt-1)
      real poidr(nx,ny,nz,nt-1), qr(nx,ny,nz,nt)
      real prt(nz), prq(nz)
      real cra(nx,ny,nt-1),cre(nx,ny,nz,nt-1)
      real sra(nx,ny,nt-1),sre(nx,ny,nz,nt-1)
      real cra2(nx,ny,nt-1),cre2(nx,ny,nz,nt-1)
      real sra2(nx,ny,nt-1),sre2(nx,ny,nz,nt-1)
      real cra3(nx,ny,nt-1),cre3(nx,ny,nz,nt-1)
      real sra3(nx,ny,nt-1),sre3(nx,ny,nz,nt-1)
      integer itr, nmo
      real acc, aph, sol
      real X(dimx), G(dimx)
      real XK1(dimx), GK1(dimx)
      real D(dimx)
      real gradi
      integer nombre_rcvr
      real dtemp
      real dixy, alt(nz), alt_t(nz)
      real F, Fal

      integer iter,nerr,ncomp
      integer i, subfinish

      real grad,gg,gd,dd,upk,dok,rdf,delF,Fm1
      real al,axpect,z,w,amin,gald
      real f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,Ji(14)
      real vrp(nx,ny,nz,(nt-1)), poidsp(nx,ny,nz,(nt-1))
      real profx, profy, vtpr(nx,ny,nz,(nt-1))
      real rx1, ry1, u0s, v0s

      iter=1
      grad=1.0
      nerr=0
      F=0.0
      f1 = 0.0
      f2 = 0.0
      f3 = 0.0
      f4 = 0.0
      f5 = 0.0
      f6 = 0.0
      f7 = 0.0
      f8 = 0.0
      f9 = 0.0
      f10 = 0.0
      f11 = 0.0
      do i=1,14
         Ji(i)=0.0
      enddo

      do while ((iter.le.itr+1).and.(grad.gt.acc))
         
         if (iter.eq.1) then
            Fm1=F
            if (nmo.eq.3) then
               call funct3(u0, au, bu, v0, av, bv,
     $              gvu0, gvau, gvbu, gvv0, gvav, gvbv,
     $              G, F, GRADI, nerr, vxb, vyb, X,
     $              cra, cre, sra,
     $              sre, cra2, cre2, sra2, sre2, cra3, cre3, sra3, sre3,
     $              vt_i, vr1_i, vr2_i, vr3_i, ref_i,
     $              poidv, poidv2, poidv3, nombre_rcvr,
     $              mag, magzb, magzt, dixy, alt, alt_t, dimx,
     $              nx, ny, nz, nt, test, f1, f2, f3, rx1, ry1, stage, sig, nsig)
            elseif (nmo.eq.5) then
               call funct5(G, F, gradi, nerr, u0s, v0s, gvu0, gvv0,
     $              X, qr, poidr, magzb, magzt, dixy, dtemp, dimx,
     $              nx, ny, nz, nt, alt)
            elseif (nmo.eq.11) then
               call funct11(vx, vy, vz, gvx, gvy, gvz,
     $              G, F, gradi, nerr, vxb, vyb, vzb, X, cra, cre,
     $              sra, sre, cra2, cre2, sra2, sre2, cra3, cre3,
     $              sra3, sre3, vt_i, vr1_i, vr2_i, vr3_i,
     $              poidw, poidr, poidv, poidv2, poidv3,
     $              nombre_rcvr, mag, magzb, magzt,
     $              dixy, alt, alt_t, dtemp, dimx, nx, ny, nz, nt,
     $              test, Ji, sol, cover, stage, sig, nsig)
            elseif (nmo.eq.12) then
               call funct12(vx, vy, vz, pres, temp, rain,
     $              gvp, gvt, G, F, gradi, nerr, prt, prq, X,
     $              aph, mag, magzb, magzt, dixy, alt, alt_t, dtemp,
     $              dimx, nxr, nyr, nzr, nx, ny, nz, nt, test,
     $              sol, cover)
            elseif (nmo.eq.13) then
               call funct13(vx, vy, vz, pres, temp, rain, mc,
     $              gvx, gvy, gvz, gvp, gvt, gvr, gvmc,
     $              G, F, gradi, nerr,
     $              prt, prq, vxb, vyb, vzb, presb, tempb, rainb, mcb, X, cra, cre,
     $              sra, sre, cra2, cre2, sra2, sre2, cra3, cre3,
     $              sra3, sre3, vt_i, vr1_i, vr2_i, vr3_i,
     $              ref_i, refrac_i, poidw, poidr, poidv, poidv2, poidv3,
     $              nombre_rcvr, u0s, v0s, mag, magzb, magzt,
     $              dixy, alt, alt_t, dtemp, dimx, nx, ny, nz, nt,
     $              test, Ji, sol, cover, stage, sig, nsig, kpass, hwl,
     $              vxv, vyv, vzv, presv, tempv, rainv, mcv)
            endif
         endif
         
         call rms2(iter,F,gradi,Ji,nmo,stage)
         
c     -------- To stop in the first iteration -------------------------
         if (itr.eq.0) then 
            grad = -1.0 
            call rms2(iter-1,F,gradi,Ji,nmo,stage)
            return
         endif
c     -----------------------------------------------------------------
         do i=1,dimx
            XK1(i)=X(i)
         enddo

c     ===============================================================
c     ------------------------Steepest descent-----------------------
         if (iter.eq.1) then
            gg=0.0
            do i=1,dimx
               gg=gg+G(i)*G(i)
               D(i)=-G(i)
               GK1(i)=G(i)
            enddo
            gd=-gg
            
         else

c     ----------------------Conjugate Gradient-----------------------
            gg=0.0
            upk=0.0
            dok=0.0
            do i=1,dimx
               gg=gg+G(i)*G(i)
               upk=upk+G(i)*(G(i)-GK1(i))
               dok=dok+D(i)*(G(i)-GK1(i))
c===Num.Recipies.               dok=dok+GK1(i)*GK1(i)
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
c<<<alain
c               dd=dd+D(i)*D(i)
c>>>alain
               gd=gd+G(i)*D(i)
            enddo

c            print*,'Angle between search direction and gradient:',180.0/(4.0*atan(1.0))*acos(-gd/sqrt(dd*gg))
            
c     --------------------------------------------------------------
            
         endif

c     ===============================================================
c     ----------------------- step length ---------------------------

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
            if (ncomp.gt.20) then
               print*,'too many subiterations'
               return
            endif
            
            do i=1,dimx
               X(i)=XK1(i) + al*D(i)
            enddo
            
            if (nmo.eq.3) then
               call funct3(u0, au, bu, v0, av, bv,
     $              gvu0, gvau, gvbu, gvv0, gvav, gvbv,
     $              G, Fal, GRADI, nerr, vxb, vyb, X,
     $              cra, cre, sra,
     $              sre, cra2, cre2, sra2, sre2, cra3, cre3, sra3, sre3,
     $              vt_i, vr1_i, vr2_i, vr3_i, ref_i,
     $              poidv, poidv2, poidv3, nombre_rcvr,
     $              mag, magzb, magzt, dixy, alt, alt_t, dimx,
     $              nx, ny, nz, nt, test, f1, f2, f3, rx1, ry1, stage, sig, nsig)
            elseif (nmo.eq.5) then
               call funct5(G, Fal, gradi, nerr, u0s, v0s, gvu0, gvv0,
     $              X, qr, poidr, magzb, magzt, dixy, dtemp, dimx,
     $              nx, ny, nz, nt)
            elseif (nmo.eq.11) then
               call funct11(vx, vy, vz, gvx, gvy, gvz,
     $              G, F, gradi, nerr, vxb, vyb, vzb, X, cra, cre,
     $              sra, sre, cra2, cre2, sra2, sre2, cra3, cre3,
     $              sra3, sre3, vt_i, vr1_i, vr2_i, vr3_i,
     $              poidw, poidr, poidv, poidv2, poidv3,
     $              nombre_rcvr, mag, magzb, magzt,
     $              dixy, alt, alt_t, dtemp, dimx, nx, ny, nz, nt,
     $              test, Ji, sol, cover, stage, sig, nsig)
            elseif (nmo.eq.12) then
               call funct12(vx, vy, vz, pres, temp, rain,
     $              gvp, gvt, G, Fal, gradi, nerr, prt, prq, X,
     $              aph, mag, magzb, magzt, dixy, alt, alt_t, dtemp,
     $              dimx, nxr, nyr, nzr, nx, ny, nz, nt, test,
     $              sol, cover)
            elseif (nmo.eq.13) then
               call funct13(vx, vy, vz, pres, temp, rain, mc,
     $              gvx, gvy, gvz, gvp, gvt, gvr, gvmc,
     $              G,
     $              Fal, gradi, nerr, prt, prq,
     $              vxb, vyb, vzb, presb, tempb, rainb, mcb, X, cra, cre,
     $              sra, sre, cra2, cre2, sra2, sre2, cra3, cre3,
     $              sra3, sre3, vt_i, vr1_i, vr2_i, vr3_i,
     $              ref_i, refrac_i, poidw, poidr, poidv, poidv2, poidv3,
     $              nombre_rcvr, u0s, v0s, mag, magzb, magzt,
     $              dixy, alt, alt_t, dtemp, dimx, nx, ny, nz, nt,
     $              test, Ji, sol, cover, stage, sig, nsig, kpass, hwl,
     $              vxv, vyv, vzv, presv, tempv, rainv, mcv)
            endif

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
                  
                  if (nmo.eq.3) then
                     call funct3(u0, au, bu, v0, av, bv,
     $                    gvu0, gvau, gvbu, gvv0, gvav, gvbv,
     $                    G, F, GRADI, nerr, vxb, vyb, X,
     $                    cra, cre, sra,
     $                    sre, cra2, cre2, sra2, sre2, cra3, cre3, sra3, sre3,
     $                    vt_i, vr1_i, vr2_i, vr3_i, ref_i,
     $                    poidv, poidv2, poidv3, nombre_rcvr,
     $                    mag, magzb, magzt, dixy, alt, alt_t, dimx, nx, ny,
     $                    nz, nt, test, f1, f2, f3, rx1, ry1, stage, sig, nsig)
                  elseif (nmo.eq.5) then
                     call funct5(G, F, gradi, nerr, u0s, v0s, gvu0, gvv0, X,
     $                    qr, poidr, magzb, magzt, dixy, dtemp, dimx,
     $                    nx, ny, nz, nt)
                  elseif (nmo.eq.11) then
                     call funct11(vx, vy, vz, gvx, gvy, gvz,
     $                    G, F, gradi, nerr, vxb, vyb, vzb, X,
     $                    cra, cre, sra, sre, cra2, cre2, sra2, sre2,
     $                    cra3, cre3, sra3, sre3, vt_i, vr1_i, vr2_i,
     $                    vr3_i, poidw, poidr, poidv, poidv2, poidv3,
     $                    nombre_rcvr, mag, magzb, magzt,
     $                    dixy, alt, alt_t, dtemp, dimx,
     $                    nx, ny, nz, nt, test,
     $                    Ji, sol, cover, stage, sig, nsig)
                  elseif (nmo.eq.12) then
                     call funct12(vx, vy, vz, pres, temp,
     $                    rain, gvp, gvt, G, F, gradi, nerr, prt, prq, X,
     $                    aph, mag, magzb, magzt, dixy,
     $                    alt, alt_t, dtemp, dimx, nxr, nyr, nzr, nx, ny, nz,
     $                    nt, test, sol, cover)
                  elseif (nmo.eq.13) then
                     call funct13(vx, vy, vz, pres, temp,
     $                    rain, mc, gvx, gvy, gvz, gvp, gvt, gvr, gvmc,
     $                    G, F, gradi, nerr,
     $                    prt, prq, vxb, vyb, vzb, presb, tempb, rainb, mcb, X,
     $                    cra, cre, sra, sre, cra2, cre2, sra2, sre2,
     $                    cra3, cre3, sra3, sre3, vt_i, vr1_i, vr2_i,
     $                    vr3_i, ref_i, refrac_i,
     $                    poidw, poidr, poidv, poidv2, poidv3,
     $                    nombre_rcvr, u0s, v0s, mag, magzb, magzt,
     $                    dixy, alt, alt_t, dtemp, dimx,
     $                    nx, ny, nz, nt, test,
     $                    Ji, sol, cover, stage, sig, nsig, kpass, hwl,
     $                    vxv, vyv, vzv, presv, tempv, rainv, mcv)
                  endif
                  
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
      
      call rms2(iter-1,F,gradi,Ji,nmo,stage)

      return
      end
