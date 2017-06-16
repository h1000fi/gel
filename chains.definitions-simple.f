         subroutine chains_definitions

         use mparameters
         use mparameters_chain
         implicit none
         include 'mpif.h'
         include 'MPI.h'
         real*8 z_center
         integer i

         if(rank.eq.0) then
         print*, 'Using chain definitions for the NPC'
         print*, '20 different chains in the system'
         endif

         N_chains = 4
         ALLOCATE (long(N_chains))
         ALLOCATE (chainsperdelta(N_chains))
         ALLOCATE (zposition(N_chains))
         ALLOCATE (segtype(N_chains, maxlong)) 

         z_center   = delta*dimZ/2.0  ! NPC center z[nm]

         i=0
         chainsperdelta(:) = polyCov
!         long(:) = polyLen

!         i=i+1
!         call read_seq(i) 
!         zposition(i)      = zspace + 3.0*deltaZp + z_center   ! in nm           

!         i=i+1
!         call read_seq(i)
!         zposition(i)      = zspace + 2.0*deltaZp + z_center

         i=i+1
         call read_seq(i)
         zposition(i)      = zspace + 1.0*deltaZp + z_center

         i=i+1
         call read_seq(i)
         zposition(i)      = zspace + 0.0*deltaZp + z_center

         i=i+1
         call read_seq(i)
         zposition(i)      = -zspace -0.0*deltaZp + z_center

         i=i+1
         call read_seq(i)
         zposition(i)      = -zspace -1.0*deltaZp + z_center

!         i=i+1
!         call read_seq(i)
!         zposition(i)      = -zspace -2.0*deltaZp + z_center

!         i=i+1
!         call read_seq(i)
!         zposition(i)      = -zspace -3.0*deltaZp + z_center

        if (i.ne.N_chains) then 
        if(rank.eq.0)print*, 'BAD NUMBER OF CHAINS for seq_type 13'
        stop
        end if                                  ! MK END ADDED

        end

C*************************************************************
      subroutine read_seq(i)

      use mparameters_chain

      implicit none
      include 'mpif.h'
      include 'MPI.h'

      character*100 filename 
      integer i, nseq, j

          write(filename,'(A15, I3.3, A4)') 'polymer_2color_',i,'.txt'  

      open(unit=2110+i,file=filename)
          
          !print*, filename, i, 'opened'
          
          read(2110+i, *), nseq
c          print*, 'nseq', nseq
          
!          if (nseq.eq.long(i)) then 
                long(i)=nseq
                do j=1,nseq !long(i)
                        read(2110+i, *), segtype(i,j) 
                enddo
!          else
!                print*, 'Error in sequence file', 2110+i
!                call MPI_FINALIZE(ierr) ! finaliza MPI
!                stop
!          endif

          close(2110+i) 
          
          return 
      end



