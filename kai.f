!#####################################################################
! This program calculate the kai parameters for poor-solvent 
! using MC integration for a 2D nanopore
!
!#####################################################################

      subroutine kai

      use mvariables
      use mparameters
      use mkai
      use mparameters_chain
      use mparameters_monomer
      use mrands

      implicit none

      include 'mpif.h'
      include 'MPI.h'

!      real*8  suma(dimR, dimR, -3:3)  ! Xu(a,b, c) used to be -1, 1
                                                     ! interaction of a
                                                     ! segment at R = a
                                                     ! and z = 0, with a
                                                     ! segment at R = b
                                                     ! and Z = c
      real*8 R,theta,z
      real*8 thetamin, thetamax, rmax, rmin, zmax, zmin, r0
      real*8 cutoff

      real*8 rn
      integer i, ii, jj
      real*8 radio ! total radius
      real*8 x1,x2,y1, y2, z1, z2, vect
      integer Rj, Zj, j
      character*40 filename

      real*8 rounding

      if(readkai.ne.1) then
         radio = delta*dimR
         if(rank.eq.0)print*,'Kai calculation, readkai =', readkai
         write(filename,'(A18, I4.4, A4)')
     & 'kais-', dimR, '.kai'
         open(unit=111, file=filename)
         Xu = 0.0               ! vector Xu
         seed = readseed
c         MCsteps = 1

         cutoff = paircut ! 1.5*delta
         rounding = 1e-6 ! used to prevent errors due to rounding

         do ii = 1, dimR       ! loop over segment positions
            do i = 1, MCsteps
       
            zmin = -cutoff+rounding
            zmax = cutoff-rounding
            
            r0 = (dfloat(ii) - 0.5)*delta

            if(r0.gt.cutoff) then
            rmin = r0-cutoff+rounding
            rmax = min(r0+cutoff, radio)-rounding
            thetamax = asin(cutoff/r0)
            thetamin = -thetamax
            else
            rmin = 0+rounding
            rmax = min(r0+cutoff, radio)-rounding
            thetamin = 0
            thetamax = 2*pi
            endif

            Z = rands(seed)*(zmax-zmin)+zmin ! Z between zmax and zmin 
            theta = rands(seed)*(thetamax-thetamin)+thetamin
            R = rands(seed)*(rmax-rmin)+rmin

!     segment coordinates (x1,y1,z1) and point to integrate coordinates
!     (x2,y2,z2)
               x1 = r0 ! theta segment = 0, z segment =
               y1 = 0.0
               z1 = 0.0
               x2 = R*cos(theta)
               y2 = R*sin(theta)
               z2 = z
               vect = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2) ! vector diff
               Rj = int(R/delta)+1 ! R cell
               Zj = nint(z/delta) ! from -.5 to .5 => layer 0
                                ! from .5 to 1.5 => layer 1
                                ! from -.5 to -1.5 => layer -1
                                ! z is the distance to the segment,
                                ! while R
                                ! is the position in the lattice
               if(Zj.eq.4) then !changed 2 to 4
               print*, z, z/delta,i,ii
               stop
               endif
               if(vect.gt.(cutoff))goto 15 ! outside cut-off sphere
!               if(vect.lt.lseg) goto 15
               do jj = 1, N_poorsol
                  if(jj.eq.1) then
                  Xu(jj,ii,Rj,Zj)=Xu(jj,ii,Rj,Zj)+1.0*R !((lseg/vect)**6)*R ! incluye el jacobiano R(segmento)
                  else
                    if(vect.lt.lseg) then !goto 15
                      Xu(jj,ii,Rj,Zj)=Xu(jj,ii,Rj,Zj)+1.0*R
                    else
                      Xu(jj,ii,Rj,Zj)=Xu(jj,ii,Rj,Zj)+((lseg/vect)**6)*R
                    endif
                  endif
               enddo !jj
 15         enddo ! MCsteps
            do Rj = 1, dimR
               do Zj = -3,3 ! changed -1,1 to -3,3
                  do jj = 1, N_poorsol
                       Xu(jj, ii, Rj, Zj) = Xu(jj, ii, Rj, Zj)/MCsteps*
     &                 (zmax-zmin)*(thetamax-thetamin)*(rmax-rmin)
                  enddo !jj
                  write(111,*)ii,Rj,Zj,Xu(:,ii,Rj,Zj)
                       !write(111,*,advance='no')Xu(jj,ii,Rj,Zj)
               enddo !Zj
            enddo !Rj
         enddo                  ! ii
         close(111)
      else                      ! readkai = 1
        write(filename,'(A5, I4.4, A4)')'kais-', dimR, '.kai'
         open(unit=111, file=filename)
         print*,'Read Kais'
         do j = 1, dimR*dimR*5
            read(111,*)ii,Rj,Zj,Xu(:,ii,Rj,Zj)
         enddo
      endif
 
      end
