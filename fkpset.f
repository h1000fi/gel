* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* * *
c     The routine kpreco is the preconditioner setup routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:

      subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)

      use mparameters
      use mkinsol
      use mparameters_chain
      use mparameters_monomer

      implicit none

      integer ier
      integer*8 neq, i
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vtemp1(*), vtemp2(*)
      integer n

      common /psize/ neq

      n = neq/(2+N_poorsol)

      do i = 1, n
         pp(i) = 1.0 / (1.0+(1.0-udata(i))*exp(1.0-udata(i)))
      enddo

      do i = n+1,2*n
         pp(i) = 1.0
      enddo

      do i = 2*n+1, neq
         pp(i) = 1.0 / (1.0+(1.0-udata(i))*exp(1.0-udata(i)))
      enddo



      ier = 0

      return
      end
