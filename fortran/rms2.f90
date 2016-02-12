subroutine rms2(iterc, f, gradi)

      implicit none

      integer:: ITERC
      real(kind=8):: F,gradi

      if (mod(ITERC-1,1).eq.0) then
         write(*,103) ITERC-1,F,gradi
      endif

 103  format(' iteration ',I6,5X,'Cost =',1PE13.6,5X,'gradient =',1PE12.6)
      
      return
end subroutine rms2
      
