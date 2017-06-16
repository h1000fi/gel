c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
* * *
c     The routine kpsol is the preconditioner solve routine. It must have
c     that specific name be used in order that the c code can find and link
c     to it.  The argument list must also be as illustrated below:

      subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)

      use mparameters
      use mkinsol
      use mparameters_chain
      use mparameters_monomer

      implicit none

      integer ier
      integer*8 neq, i
      double precision udata(*), uscale(*), fdata(*), fscale(*)
      double precision vv(*), ftem(*)

      common /psize/ neq

      do  i = 1, neq
         vv(i) = vv(i) * pp(i)
      enddo
      ier = 0

      return
      end
