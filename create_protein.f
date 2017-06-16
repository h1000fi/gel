        subroutine create_protein

        use mprotein
        use mvariables
        use mncells
        use mlookup
        use mparameters
        use mncells
        use mporesystem

        implicit none

      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries


        integer iC, iR, iZ
        real*8 volprot, volprot2
        real*8 posz, posr, rcurve
        integer temp(0:dimR+1, 0:dimZ+1)
        real*8 protcutoff
        integer flag
        real*8 totalvol

        real*8 PdimRr
        real*8 PdimHr
        real*8 Hr
        real*8 Pposr
        real*8 center
        real*8 temp1
        integer cc, ccc
        character*5  title

        Hr = dfloat(hrange)*delta
        PdimRr = dfloat(PdimR)*delta
        Pposr = dfloat(Ppos)*delta
        PdimHr = dfloat(PdimH)*delta
       
        temp = 0

        volprot = 0.0
        volprot2 = 0.0
        protein = 0.0
        proteinC = 0.0 
        proteinq = 0.0
        proteinqC = 0.0 
        proteinhC = 0.0 
        proteinh = 0.0 
 
! assigns volume

        totalvol = 0.0

       do iZ = 1, dimZ
       do iR = 1, dimR
 
       posz = (dfloat(iZ)-0.5)*delta
       posr = (dfloat(iR)-0.5)*delta

       rcurve = 0.0 

      if((Ptype.eq.1).or.(Ptype.eq.3)) then ! sphere
      if((posz.gt.(Pposr-PdimRr))
     & .and.(posz.lt.(Pposr+PdimRr))) then
       rcurve = sqrt((PdimRr)**2
     &- (posz-Pposr)**2)
      endif
      endif

      if(Ptype.eq.4) then ! ellipsoid
      if((posz.gt.(Pposr-PdimHr))
     & .and.(posz.lt.(Pposr+PdimHr))) then
       rcurve = PdimRr/PdimHr*sqrt((PdimHr)**2
     &- (posz-Pposr)**2)
      endif
      endif

      if(Ptype.eq.2) then ! rounded cylinder

      temp1 = PdimRr+PdimHr/2

       if((posz.gt.(Pposr-temp1)) ! lower cap
     & .and.(posz.lt.(Pposr-PdimHr/2.0))) then
       center =  Pposr-PdimHr/2.0
       rcurve = sqrt((PdimRr)**2
     &- (posz-center)**2)
       endif


       if((posz.lt.(Pposr+temp1)) ! upper cap
     & .and.(posz.gt.(Pposr+PdimHr/2.0))) then
       center = Pposr+PdimHr/2.0
       rcurve = sqrt((PdimRr)**2
     &- (posz-center)**2)
       endif

       if((posz.le.(Pposr+PdimHr/2.0)) ! cylinder
     & .and.(posz.ge.(Pposr-PdimHr/2.0))) then
       rcurve = PdimRr 
       endif
       endif

      if(Ptype.eq.5) then ! cylinder
      temp1 = PdimRr+PdimHr/2
       if((posz.le.(Pposr+PdimHr/2.0)) ! cylinder
     & .and.(posz.ge.(Pposr-PdimHr/2.0))) then
       rcurve = PdimRr 
       endif
       endif

 

      if(rcurve.gt.posr) then ! is protein
      if(matriz(iR, iZ).eq.-3) then
      print*, 'There is a collision 
     & between the protein and the pore walls!'
      endif

      protein(iR,iZ) = prot_vol
      temp(iR,iZ) = 1
      totalvol = totalvol 
     & + (dfloat(iR)-0.5)*delta*2*pi*delta**2
      endif

      enddo
      enddo

      print*, 'Total prot vol', totalvol

! assigns charge to the surface

      if((Ptype.eq.1).or.(Ptype.eq.2).or.
     & (Ptype.eq.4).or.(Ptype.eq.5)) then
      do iR = 1, dimR
      do iZ = 1, dimZ
 
      flag = 0

      if(temp(iR,iZ).eq.1) then

      if((temp(iR+1,iZ).eq.0).and.(matriz(iR+1,iZ).ne.-1))flag=1
      if((temp(iR-1,iZ).eq.0).and.(matriz(iR-1,iZ).ne.-1))flag=1
      if((temp(iR,iZ+1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1
      if((temp(iR,iZ-1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1

      if(flag.eq.1) then
      volprot = volprot + delta**3*2.0*pi*(dfloat(iR)-0.5)
      proteinq(iR,iZ) = prot_q
      endif
      endif

      enddo
      enddo

      if(volprot.ne.0.0)proteinq = proteinq / volprot
      endif


      if(Ptype.eq.3) then
      do iR = 1, dimR
      do iZ = 1, dimZ

      flag = 0
      if(temp(iR,iZ).eq.1) then
      if((temp(iR+1,iZ).eq.0).and.(matriz(iR+1,iZ).ne.-1))flag=1
      if((temp(iR-1,iZ).eq.0).and.(matriz(iR-1,iZ).ne.-1))flag=1
      if((temp(iR,iZ+1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1
      if((temp(iR,iZ-1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1
      if(flag.eq.1) then
      posz = (dfloat(iZ)-0.5)*delta
      if(posz.gt.Pposr) then
      volprot = volprot + delta**3*2.0*pi*(dfloat(iR)-0.5)
      else
      volprot2 = volprot2 + delta**3*2.0*pi*(dfloat(iR)-0.5)
      endif
      endif
      endif
      enddo
      enddo


      do iR = 1, dimR
      do iZ = 1, dimZ
      flag = 0
      if(temp(iR,iZ).eq.1) then
      if((temp(iR+1,iZ).eq.0).and.(matriz(iR+1,iZ).ne.-1))flag=1
      if((temp(iR-1,iZ).eq.0).and.(matriz(iR-1,iZ).ne.-1))flag=1
      if((temp(iR,iZ+1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1
      if((temp(iR,iZ-1).eq.0).and.(matriz(iR,iZ-1).ne.-1))flag=1
      if(flag.eq.1) then
       posz = (dfloat(iZ)-0.5)*delta
      if(posz.gt.Pposr) then
      proteinq(iR,iZ) = prot_q/volprot
      else
      proteinq(iR,iZ) = prot_q2/volprot2
      endif
      endif
      endif
      enddo
      enddo

      end if



! assigns hydrophobic interaction matrix

        do iZ = 1, dimZ
        do iR = 1, dimR

       posz = (dfloat(iZ)-0.5)*delta
       posr = (dfloat(iR)-0.5)*delta

       rcurve = 0.0

      if((Ptype.eq.1).or.(Ptype.eq.3)) then ! sphere
       if((posz.gt.(Pposr-PdimRr-Hr))
     & .and.(posz.lt.(Pposr+PdimRr+Hr))) then
       rcurve = sqrt(((PdimRr+Hr))**2
     &- (posz-Pposr)**2)
      endif
      endif

      if(Ptype.eq.2) then ! rounded cylinder

      temp1 = PdimRr+PdimHr/2

       if((posz.gt.(Pposr-temp1-Hr)) ! lower cap
     & .and.(posz.lt.(Pposr-PdimHr/2.0))) then
       center =  Pposr-PdimHr/2.0
       rcurve = sqrt((PdimRr+Hr)**2
     &- (posz-center)**2)
       endif


       if((posz.lt.(Pposr+temp1+Hr)) ! upper cap
     & .and.(posz.gt.(Pposr+PdimHr/2.0))) then
       center = Pposr+PdimHr/2.0
       rcurve = sqrt((PdimRr+Hr)**2
     &- (posz-center)**2)
       endif

       if((posz.le.(Pposr+PdimHr/2.0)) ! cylinder
     & .and.(posz.ge.(Pposr-PdimHr/2.0))) then
      rcurve = PdimRr+Hr
       endif
       endif

      if(Ptype.eq.5) then ! cylinder
      temp1 = PdimRr+PdimHr/2
       if((posz.le.(Pposr+PdimHr/2.0)) ! cylinder
     & .and.(posz.ge.(Pposr-PdimHr/2.0))) then
      rcurve = PdimRr+Hr
       endif
       endif


      if(Ptype.eq.4) then ! ellipsoid
      if((posz.gt.(Pposr-PdimHr))
     & .and.(posz.lt.(Pposr+PdimHr))) then
       rcurve = PdimRr/PdimHr*sqrt((PdimHr)**2
     &- (posz-Pposr)**2)
      endif
      endif


      if(rcurve.gt.posr) then ! is protein
      if((matriz(iR, iZ).ge.1).and.(protein(iR, iZ).eq.0.0)) then ! if the point is inside the system and is not protein 
      proteinh(iR, iZ) = 1.0
      endif
      endif

      enddo
      enddo

! create proteinC and proteinqC

      do iC = 1, ncells
      proteinC(iC) = protein(indexa(iC,1),indexa(iC,2))
      proteinqC(iC) = proteinq(indexa(iC,1),indexa(iC,2))
      proteinhC(iC) = proteinh(indexa(iC,1),indexa(iC,2))
      enddo

      if(savetodisk_flag.eq.12) then

      if(rank.eq.0) then
      print*, 'Saving protein properties and exiting'
      cc = 1
      ccc = 1

      title = 'protC'
      call savetodisk(proteinC, title, cc ,ccc)

      title = 'proqC'
      call savetodisk(proteinqC, title, cc ,ccc)

      title = 'prohC'
      call savetodisk(proteinhC, title, cc ,ccc)
      endif
       
      call MPI_FINALIZE(ierr) ! finalize MPI
      stop
      endif 

      end

