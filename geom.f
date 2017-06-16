       subroutine geom
       use mparameters
       use mncells
       use mporesystem
       use mparameters_chain
       use mparameters_monomer

       implicit none
       integer midR, iR, iZ, midZ

       integer*4 temp(0:dimR+1, 0:dimZ+1)
       real*8 zpos, rpos, slope, ord, rcurve, rcurve2

       real*8 b_param, a_param, c_param
       real*8 volprot

       integer Zmin, Zmax

      integer PposZ

      common /PposZ/ PposZ
      


C codes matriz = ! NCELLS+1: bulk iz
                 ! NCELLS+2: bulk derecho
                 ! -1:symmetry
                 ! -2: wall
                 !  0: canal
                 ! -3: no definido


       matriz = -3 ! -3 => undefined

       ncells = 0
       temp = 0

!!!!!!!!!!! pore !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      select case (poretype)

! -----------------------  Simplest case: long pore --------------------------
      case (1)
       midR = int(dfloat(dimR)/2)



       do iZ = 1, dimZ
       do iR = 1, midR
  
       matriz(iR, iZ) = 0 ! canal
       ncells = ncells + 1

       enddo
       enddo

!------------------------ Hour glass, using a parabola ---------------------- 
       case (3)
       a_param = curvature
       b_param = -a_param*dimz*delta
       c_param = CdimR*delta
     & -a_param*(RdimZ*delta)**2 - b_param*RdimZ*delta

       do iZ = 1, dimZ
       do iR = 1, dimR

       zpos = (dfloat(iZ)-0.5)*delta
       rpos = (dfloat(iR)-0.5)*delta

       rcurve = a_param*zpos**2 + b_param*zpos + c_param

       if((rcurve.gt.rpos).or.(iZ.le.RdimZ).or.(iZ.gt.(dimZ-RdimZ))) 
     & then
       matriz(iR, iZ) = 0 ! canal
       ncells = ncells + 1
       endif

       enddo
       enddo

! ---------------------------------- conic pore ----------------------
       case (4)

       do iZ = 1, dimZ
       do iR = 1, dimR

       zpos = (dfloat(iZ)-0.5)*delta
       rpos = (dfloat(iR)-0.5)*delta

       slope =
     & (dfloat(BdimR)-dfloat(TdimR))/(dfloat(dimZ)-2*dfloat(RdimZ))
       ord = dfloat(TdimR)*delta-RdimZ*delta*slope

      

       rcurve = slope*zpos + ord

       if((rcurve.gt.rpos).or.(iZ.le.RdimZ).or.(iZ.gt.(dimZ-RdimZ)))
     & then
       matriz(iR, iZ) = 0 ! canal
       ncells = ncells + 1
       endif

       enddo
       enddo

! ------------------------------ hour glass (parabola) and protein --------
c       case (5)
c
c       matriz = 0
c
c       temp = 0 
c       do iZ = 1, dimZ
c       temp(0,iZ) = -1
c       enddo
c       ncells = dimZ*dimR
c
c       a_param = curvature
c       b_param = -a_param*dimz*delta
c       c_param = CdimR*delta
c     & -a_param*(RdimZ*delta)**2 - b_param*RdimZ*delta
c
c       do iZ = 1, dimZ
c       do iR = 1, dimR
c
c       zpos = (dfloat(iZ)-0.5)*delta
c       rpos = (dfloat(iR)-0.5)*delta
c
c       rcurve = a_param*zpos**2 + b_param*zpos + c_param
c
c       if((rcurve.le.rpos).and.(iZ.gt.RdimZ).and.(iZ.le.(dimZ-RdimZ)))
c     & then ! es pared
c
c       matriz(iR, iZ) = -3 
c       ncells = ncells - 1
c
c      endif
c
c       if((zpos.gt.(dfloat(PposZ)-dfloat(PdimR))*delta)
c     & .and.(zpos.lt.(dfloat(PposZ)+dfloat(PdimR))*delta)) then
c
c       rcurve2 = sqrt((float(PdimR)*delta)**2
cc     &- (zpos-dfloat(PposZ)*delta)**2) 
c
c      if(rcurve2.gt.rpos) then ! es particula
c
c      if(matriz(iR, iZ).eq.-3) then
c      print*, 'Colision entre pared y particula'
c      endif
c
c      protein(iR,iZ) = prot_vol
c      temp(iR,iZ) = 1
c
c      endif
c      endif
c      
c      enddo
c      enddo
c
c      do iR = 1, dimR
c      do iZ = 1, dimZ
c
c      if(temp(iR,iZ).eq.1) then
c
c      if((temp(iR+1,iZ).eq.0).or.(temp(iR-1,iZ).eq.0).or.
c     & (temp(iR,iZ+1).eq.0).or.(temp(iR,iZ-1).eq.0)) then
c
c      volprot = volprot + delta**3*2*pi*(dfloat(iR)-0.5)
c      proteinq(iR,iZ) = prot_q
c
c      endif
c      endif
c     
c      enddo
c      enddo
c
c      if(volprot.ne.0)proteinq = proteinq / volprot
c
c       case (40)
c
c       do iZ = 1, dimZ
c       do iR = 1, dimR
c
c       zpos = (dfloat(iZ)-0.5)*delta
c       rpos = (dfloat(iR)-0.5)*delta
c
c       slope =
c     & (dfloat(BdimR)-dfloat(TdimR))/(dfloat(dimZ)-2*dfloat(RdimZ))
c       ord = dfloat(TdimR)*delta-RdimZ*delta*slope
c
c       rcurve = slope*zpos + ord
c
c       if((rcurve.gt.rpos).or.(iZ.le.RdimZ).or.(iZ.gt.(dimZ-RdimZ)))
c     & then
c       matriz(iR, iZ) = 0 ! canal
c       ncells = ncells + 1
c       endif
c
c      enddo
c      enddo

!------------------------ Hour glass, using half circle  ---------------------- 
       case (10)
     
       ncells = dimZ*dimR
       matriz = 0

       zmin = int((dimZ - CdimZ)/2)
       zmax = zmin + CdimZ

       do iZ = zmin+1, zmax
       do iR = 1, dimR

       rcurve = sqrt((dfloat(CdimZ)/2*delta)**2 
     & - ((dfloat(iZ-zmin)-0.5)*delta - dfloat(CdimZ)/2*delta)**2)
       rcurve = CdimR*delta + dfloat(CdimZ)/2*delta - rcurve 

       zpos = (dfloat(iZ)-0.5)*delta
       rpos = (dfloat(iR)-0.5)*delta

       if((rcurve.lt.rpos))
     & then
       matriz(iR, iZ) = -3 
       ncells = ncells - 1
       endif

       enddo
       enddo

!------------------------ Hour glass like  ---------------------- 

         case (11)
  
         ncells = dimZ*dimR
         matriz = 0
 
         zmin = (dimZ - CdimZ)/2
         zmax = zmin + CdimZ
 
         do iZ = zmin+1, zmin+CdimZ/2
         do iR = CdimR+1, dimR
 
         slope = (dfloat(CdimZ-CdimZmin)/2)/(dfloat(CdimZ)/2)
         slope = -1/slope
         ord = (dfloat(CdimR)+dfloat(CdimZ)/2)*delta-zmin*delta*slope
 
         zpos = (dfloat(iZ)-0.5)*delta
         rpos = (dfloat(iR)-0.5)*delta

         rcurve = slope*zpos + ord

         if((rcurve.lt.rpos))
     &   then
         matriz(iR, iZ) = -3
         ncells = ncells - 1
         endif
 
         enddo
         enddo

         do iZ = zmin+CdimZ/2+1, zmax
         do iR = CdimR+1, dimR

         slope = (dfloat(CdimZ-CdimZmin)/2)/(dfloat(CdimZ)/2)
         slope = 1/slope
         ord = (dfloat(CdimR)+dfloat(CdimZ)/2)*delta-zmax*delta*slope

         zpos = (dfloat(iZ)-0.5)*delta
         rpos = (dfloat(iR)-0.5)*delta
         rcurve = slope*zpos + ord
         if((rcurve.lt.rpos))
     &   then
         matriz(iR, iZ) = -3
         ncells = ncells - 1
         endif

         enddo
         enddo

!------------------------------ Short Cycilnder --------------------

         case (12)

         ncells = dimZ*dimR
         matriz = 0

         zmin = int((dimZ - CdimZ)/2)
         zmax = zmin + CdimZ

         do iZ = zmin+1, zmax
         do iR = CdimR+1, dimR

         matriz(iR, iZ) = -3
         ncells = ncells - 1

         enddo
         enddo

       endselect


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reservoirs

       
       midZ = int(dfloat(dimZ)/2)
       do iZ = 0, midZ
       matriz(dimR+1, iZ) = ncells+1
       enddo

       do iZ = midZ, dimZ+1
       matriz(dimR+1, iZ) = ncells+2
       enddo

       do iR = 1, dimR+1
       matriz(iR, 0) = ncells+1
       matriz(iR, dimZ+1) = ncells+2
       enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! symmetry

       do iZ = 0, dimZ+1
       matriz(0, iZ) = -1
       enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! search for borders
       temp = matriz

       do iR = 1, dimR
       do iZ = 1, dimZ

          if(temp(iR, iZ).eq.-3) then
               if((temp(iR, iZ+1).eq.0).or.(temp(iR, iZ-1).eq.0).or.
     &   (temp(iR+1, iZ).eq.0).or.(temp(iR-1, iZ).eq.0)) then

                matriz(iR, iZ) = -2 ! pared
           endif
       endif

       enddo
       enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! save to disk


c      do iR = 0,dimR+1
c      do iZ = 0,dimZ+1
c      write(999, *)iR, iZ, matriz(iR, iZ)
c      enddo
c      enddo
      end
