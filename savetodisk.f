      subroutine savetodisk(array2, title, counter, counter2)

      use mlookup
      use mparameters
      use mvariables
      use mncells
      use mparameters_chain
      use mparameters_monomer

       implicit none
       real*4 single
       integer ix, iy, iR, i, jx, jy, jz, iZ, iC
       integer maxT, iT
       real*8 posx, posy, posz
       character*30 filename, tempc
       character*6 titlez
       character*5 title
       integer counter, counter2
       real*8 arrayC(ncells)
       real*8 array2(dimR*dimZ)

       real*8 array(dimR, dimZ)


! Variables

! Crea array

      array = -10 ! hay que extraer las celdas con valor -1000

      do iC = 1, ncells
      array(indexa(iC,1),indexa(iC,2))=array2(iC)
      enddo

! Archivo paraview 3D

      if((savetodisk_type.eq.0).or.(savetodisk_type.eq.2)) then

      maxT = 1

      write(filename,'(A5, A1, A2, F3.1, A4)')title,'.',
     &  'pH', pHbulk, '.vtk'
      open (unit=45, file=filename)
      write(45,'(A)')'# vtk DataFile Version 2.0'
      write(45,'(A)')title
      write(45,'(A)')'ASCII'
      write(45,'(A)')'DATASET STRUCTURED_GRID '
      write(45,'(A, I5, A1, I5, A1, I5)')
     & 'DIMENSIONS', dimZ+1, ' ', maxT+1, ' ',dimR+1
      write(45,'(A, I8, A)')'POINTS ',(dimZ+1)*(maxT+1)
     % *(dimR+1),' float'

      do iR = 0, dimR
        do iT = 0, maxT
          do iz = 0, dimZ

      posx = sin(dfloat(iT)/maxT*2.0*3.14159)*(iR)*delta
      posy = cos(dfloat(iT)/maxT*2.0*3.14159)*(iR)*delta
      posz = iz*delta

            write(45, *)
     & posx,'   ', ! sistema de coord x y
     & posy
     &, '   ', posz ! grafico en el viejo sistema x,y no en v,u
          enddo
        enddo
      enddo

      write(45,'(A, I8)')'CELL_DATA ', dimR*dimZ*maxT
      tempc = 'SCALARS ' // title // ' float 1'
      write(45,'(A)')tempc
      write(45,'(A)')'LOOKUP_TABLE default'

       do iR = 1, dimR
        do iT = 1, maxT
          do iZ = 1, dimZ

             single = array(iR, iZ) ! Lo necesito en single presicion

            write(45,*)single
          enddo
        enddo
      enddo
      close (45)

      endif

      if((savetodisk_type.eq.1).or.(savetodisk_type.eq.2)) then

      write(filename,'(A5, A1, A2, F3.1, A4)')title,'.',
     &  'pH', pHbulk, '.dat'

      open(unit=45, file=filename)

      do iR=1,dimR
         do iZ=1,dimZ
            write(45,*)(iR-0.5)*delta, (iZ-0.5)*delta, array(iR, iZ)
         enddo
      enddo

      close(45)

      endif

      return
      end
