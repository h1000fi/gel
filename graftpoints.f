      subroutine graftpoints
      use mparameters
      use mncells
      use mporesystem
      use mgraft       
      use mparameters_chain
      use mparameters_monomer

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h'
      integer midR, iZ, iR, ii, zmax, zmin
       real*8 zpos, rpos, slope, ord, rcurve
      integer rpos_int

      real*8 a_param, b_param, c_param

        logical   outside
        external  outside

      real*8 x(3)

      ii = rank+1

      select case (poretype)

!-------------------- Simplest case: long pore ---------------------
      case (1)


       midR = int(dfloat(dimR)/2)

       posgraft(1) = dfloat(dimR)/2*delta
       posgraft(2) = dfloat(dimZ)/2*delta



!-------------------- conic pore  ---------------------
c      case (4)
c
c       slope =
c     & (dfloat(BdimR)-dfloat(TdimR))/(dfloat(dimZ)-2*dfloat(RdimZ))
c       ord = TdimR*delta-RdimZ*delta*slope
c       print*, ord, slope
c
c       zpos = zposition(ii) 
c       rcurve = slope*zpos + ord ! position on the curve
c
c       rpos_int = int(rcurve/delta) + 1  ! lattice where the position on the curve is
c
c       rpos = (dfloat(rpos_int)-0.5)*delta ! center of that lattice
c
c       if(rcurve.lt.rpos)rpos = rpos - delta ! is lattice site part of the membrane? 
c 
c       rpos = rpos - delta/2 + 1e-5 ! border of the lattice site
c
cc       posgraft(1) = rpos 
c       posgraft(2) = zpos




!-------------------- hour glass half circle  ---------------------
       case (10)

       zmin = int((dimZ - CdimZ)/2)
       zmax = zmin + CdimZ

       zpos = zposition(ii) 
       rpos = dimR*delta - delta/2.0

       x(1) = rpos
       x(2) = 0.0
       x(3) = zpos

       do while (outside(x))
       rpos = rpos - delta/2.0
       x(1) = rpos
       x(2) = 0.0
       x(3) = zpos
       enddo

       posgraft(1) = rpos
       posgraft(2) = zpos


!------------------------ Hour glass like  ----------------------

       case (11)

       zmin = int((dimZ - CdimZ)/2)
       zmax = zmin + CdimZ

       zpos = zposition(ii) 
       rpos = dimR*delta - delta/2.0
       
       x(1) = rpos
       x(2) = 0.0
       x(3) = zpos

       do while (outside(x))
       rpos = rpos - delta/2.0
       x(1) = rpos
       x(2) = 0.0
       x(3) = zpos
       enddo 

       posgraft(1) = rpos
       posgraft(2) = zpos

!------------------------ cylinder  ----------------------

       case (12)

       zmin = int((dimZ - CdimZ)/2)
       zmax = zmin + CdimZ

       zpos = zposition(ii) 
       rpos = (dfloat(CdimR)-0.5)*delta

       posgraft(1) = rpos
       posgraft(2) = zpos

!-------------------- hour glass parabola an protein  ---------------------
c       case (5)
c
c       a_param = curvature
c       b_param = -a_param*dimz*delta
c       c_param = CdimR*delta
c     & -a_param*(RdimZ*delta)**2 - b_param*RdimZ*delta
c
c       zpos = zposition(ii) 
c       rcurve = a_param*zpos**2 + b_param*zpos + c_param
c
c       rpos_int = int(rcurve/delta) + 1  ! lattice where the position on the curve is
c
c       rpos = (dfloat(rpos_int)-0.5)*delta ! center of that lattice
c
c       if(rcurve.lt.rpos)rpos = rpos - delta ! is lattice site part of the membrane?
c
c       rpos = rpos - delta/2 + 1e-5 ! border of the lattice site
c
c       posgraft(1) = rpos
c       posgraft(2) = zpos


!-------------------- hour glass parabola  ---------------------
c      case (3)

c       a_param = curvature
c       b_param = -a_param*dimz*delta
c       c_param = CdimR*delta
c     & -a_param*(RdimZ*delta)**2 - b_param*RdimZ*delta
c
c       zpos = zposition(ii) 
c       rcurve = a_param*zpos**2 + b_param*zpos + c_param
c
c       rpos_int = int(rcurve/delta) + 1  ! lattice where the position on the curve is
c
c       rpos = (dfloat(rpos_int)-0.5)*delta ! center of that lattice
c
c       if(rcurve.lt.rpos)rpos = rpos - delta ! is lattice site part of the membrane? 
c
c       rpos = rpos - delta/2 + 1e-5 ! border of the lattice site
c
cc       posgraft(1) = rpos
c       posgraft(2) = zpos
c
       endselect
c
      if (N_chains.ne.size) then
      print*, 'The number of chains,', N_chains, ', and the
     & the number of processors',size,
     & 'should be the same'
      call MPI_FINALIZE(ierr) ! finaliza MPI
       stop
       endif

      print*, 'Chain #', rank+1, ' grafted at x=', posgraft(1), ', y=',
     & posgraft(2)

      end 
