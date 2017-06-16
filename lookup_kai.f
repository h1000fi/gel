      subroutine lookup_kai

      use mlookup
      use mparameters
      use mkai
      use mncells
      use mporesystem
      use mparameters_chain
      use mparameters_monomer
        

      implicit none
      integer ii, iC, iR, iZ,  iiR, iiZ, imin, imax

! define Xu

      nXu = 0
      Xulist_cell = 0
      Xulist_value = 0

      do ii = 1, N_poorsol
      do iC = 1, ncells ! loop over cells

      iR = indexa(iC,1) ! position of cell in the system
      iZ = indexa(iC,2)

      do iiR = 1, dimR ! loop over R-neighbors

      imin = -min(iZ,3)
      imax = min((dimZ-iZ),3)

      do iiZ = imin,imax ! loop over Z-neighbors
      if((matriz(iiR,iiZ+iZ).gt.0).and.(Xu(ii,iR,iiR,iiZ).ne.0d0)) then ! inside the system and different from zero

      nXu(ii, iC) = nXu(ii, iC) + 1
      if(nXu(ii, iC).gt.size(Xulist_cell,2)) then
      print*, size(Xulist_cell,2)
      print*, 'Error in lookup_kai, increase Xulist_cell dimensions'
      stop
      endif

      Xulist_cell(ii, iC, nXu(ii, iC)) = matriz(iiR, iiZ+iZ)
      Xulist_value(ii, iC, nXu(ii, iC)) = Xu(ii, iR, iiR, iiZ)
      endif

      enddo ! iiZ
      enddo ! iiR

      enddo ! iC
      enddo ! ii
      return
      end
