! NOT IN USE 
        subroutine checking_actual_config(N)
        use mparameters
        use mgraft
        use mporesystem
        use mrands
        use mparameters_chain
        use mparameters_monomer
        use posmk

        implicit none
        integer ncha

        integer N,k
        integer ncha_current
        common /comncha/ ncha_current
        integer i,j1,j2
        real*8  dist2
        print *,'checking .. ',ncha_current
        do j1=1,N-1
        do j2=j1+1,N            
         dist2 = 0.D0
         do k=1,3
          dist2 = dist2 + (current(j1,k)-current(j2,k))**2
         enddo
         if (dist2.lt.d2) then
          print *,j1,j2
          print *,dist2,' < ',d2
          stop 'BUG'
         end if 
        enddo
        enddo
        return
        end
