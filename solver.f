      subroutine solver(x1,xg1,ier)

      use mlookup
      use mparameters ! parametros     
      use mncells
      use mkai
      use mporesystem
      use mvariables
      use mparameters_chain
      use mparameters_monomer
      use mprotein
      use mrands

      implicit none

      include 'mpif.h'
      include 'MPI.h' ! MPI libraries

C--------------------------------------------

! Kinsol

      integer *4 ier ! Kinsol error flag
      integer *8 neq ! Kinsol number of equations

! IMSL

      external fcnelect

C---------------------------------------------
      common /psize/ neq ! Kinsol
      external fcn


      real*8 x1((2+N_poorsol)*ncells)
      real*8 xg1((2+N_poorsol)*ncells)
      real*8 xflag((2+N_poorsol)*ncells)
      real*8 xflag2((2+N_poorsol)*ncells)
      real*8 f((2+N_poorsol)*ncells)

! JEFE

        if(rank.eq.0) then ! solo el jefe llama al solver

          iter = 0
          print*, 'Enter Solver ', ncells*(2+N_poorsol), ' equations.'
            

          if(infile.ne.5) then       
             call call_kinsol(x1, xg1, ier)
          else
             call fkfun(x1, f, ier) ! data analisis
          endif
          
          flagsolver = 0
          
          CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, 
     &    MPI_COMM_WORLD,err)
        endif

! Subordinados

        if(rank.ne.0) then

          do

            flagsolver = 0
            source = 0

            call MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, 
     &           MPI_COMM_WORLD,err)

            if(flagsolver.eq.1) then
               call MPI_BCAST(x1, (2+N_poorsol)*ncells, 
     &              MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,err)

               call fkfun(x1, f, ier) ! todavia no hay solucion => fkfun 
            endif ! flagsolver

            if(flagsolver.eq.0)goto 1010 ! Detiene el programa para este nodo

          enddo

        endif

 1010   continue
          
! Retrives ier and nomr
! This way slaves know if solver process in master has converged
!

         ! master

        if (rank.eq.0) then

            norma_tosend = norma
            ier_tosend = ier 

            call MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,  
     &      MPI_COMM_WORLD,err)

            call MPI_BCAST(ier_tosend, 1, MPI_INTEGER, 0,
     &      MPI_COMM_WORLD,err)

        endif

        ! slaves

        if (rank.ne.0) then

            call MPI_BCAST(norma_tosend, 1, MPI_DOUBLE_PRECISION,0,  
     &      MPI_COMM_WORLD,err)

            call MPI_BCAST(ier_tosend, 1, MPI_INTEGER, 0,
     &      MPI_COMM_WORLD,err)

            norma = norma_tosend
            ier = ier_tosend

        endif
      end
