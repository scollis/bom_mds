      subroutine rms2(iterc, f, gradi, Ji, nmo, pass)
      implicit none

      integer ITERC,nmo, pass
      real F,gradi
      real Ji(14)

c      write(4,*) ITERC-1,F,gradi
      if (mod(ITERC-1,1).eq.0) then
         write(*,103) ITERC-1,F,gradi
      endif

 103  format(' iteration ',I6,5X,'cout =',1PE12.6,5X
     $     ,'gradient =',1PE12.6)
      
      return
      end
      
