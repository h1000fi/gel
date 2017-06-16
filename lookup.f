      subroutine lookup

      use mlookup
      use mparameters
      use mporesystem
      use mncells
      use mparameters_chain
      use mparameters_monomer

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h'
      integer iR, iZ, indice, iC

c ----------------- contruye matriz -------------------------


C codes matriz = ! NCELLS+1: bulk iz
                 ! NCELLS+2: bulk derecho
                 ! -1:symmetry
                 ! -2: wall
                 !  0: canal
                 ! -3: no definido


c------------------ crea lookup list --------------------------

! indexa matriz


      indice = 0

      do iR = 1,dimR
      do iZ = 1,dimZ

      if(matriz(iR, iZ).eq.0) then
      indice = indice + 1
       matriz(iR, iZ) = indice
      indexa(indice,1) = iR
      indexa(indice,2) = iZ
      endif
      enddo
      enddo

! save to disk

c      do iR = 0,dimR+1
c      do iZ = 0,dimZ+1
c      write(998, *)iR, iZ, matriz(iR, iZ)
c      enddo
c      enddo

! Ahora define las matrices...


      do iR = 1,dimR
      do iZ = 1,dimZ

      if(matriz(iR,iZ).gt.0) then ! punto dentro del sistema

! Radio +1

      select case (matriz(iR+1,iZ))
 
      case (1 : )
      rp(matriz(iR,iZ))=matriz(iR+1,iZ)
      case (-2:-1)
      rp(matriz(iR,iZ))=matriz(iR,iZ) ! condicion de simetria
      case (-3)
      if(rank.eq.0)print*, "Error en la matriz!", iR, iZ
      
      end select

! Radio -1

      select case (matriz(iR-1,iZ))

      case (1 : )
      rm(matriz(iR,iZ))=matriz(iR-1,iZ)
      case (-2:-1)
      rm(matriz(iR,iZ))=matriz(iR,iZ) ! condicion de simetria
      case (-3)
      if(rank.eq.0)print*, "Error en la matriz!", iR, iZ

      end select

! Zeta +1

      select case (matriz(iR,iZ+1))

      case (1 : )
      zp(matriz(iR,iZ))=matriz(iR,iZ+1)
      case (-2:-1)
      zp(matriz(iR,iZ))=matriz(iR,iZ) ! condicion de simetria
      case (-3)
      if(rank.eq.0)print*, "Error en la matriz!", iR, iZ

      end select

! Zeta - 1

      select case (matriz(iR,iZ-1))

      case (1 : )
      zm(matriz(iR,iZ))=matriz(iR,iZ-1)
      case (-2:-1)
      zm(matriz(iR,iZ))=matriz(iR,iZ) ! condicion de simetria
      case (-3)
      if(rank.eq.0)print*, "Error en la matriz!", iR, iZ

      end select
      endif

      enddo
      enddo
  
      return
      end
