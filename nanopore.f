!#############################################################
!
! NPC program v3.14

      program nanopore

      use mncells
      use mparameters
      use mparameters_chain
      use mparameters_monomer

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries

      integer counter
      counter = 1

c--------------------------------------------------
!
!  Inits MPI
!

! MPI

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

       if(rank.eq.0)print*, 'Nuclear Pore Complex Program v3.7'
c--------------------------------------------------

       call globals       ! read global variables from file and allocate arrays
       if(rank.eq.0)print*, 'Globals OK'
       if(rank.eq.0)print*, ' '

c--------------------------------------------------

      call monomer_definitions ! set different types of monomer
      if(rank.eq.0)print*, 'Monomer definitions OK'
       if(rank.eq.0)print*, ' '
      call chains_definitions   ! set different types of chains and load protein sequences
      if(rank.eq.0)print*, 'Chains definitions OK'
       if(rank.eq.0)print*, ' '

c--------------------------------------------------

      call array_alloc ! array memory allocation
      if(rank.eq.0)print*, 'Memoray allocation OK'
       if(rank.eq.0)print*, ' '

c--------------------------------------------------

      call kai                 ! calculate or load kai
      if(rank.eq.0)print*, 'Kai OK'
       if(rank.eq.0)print*, ' '

      call geom                ! define pore shape
      if(rank.eq.0)print*, 'Geom OK'
       if(rank.eq.0)print*, ' '

      call lookup              ! define conectivity in matriz
      if(rank.eq.0)print*, 'Lookup OK'
       if(rank.eq.0)print*, ' '

      call create_protein
      if(rank.eq.0)print*, 'Create OK'
       if(rank.eq.0)print*, ' '

      call graftpoints         ! define graftpoint positions
      if(rank.eq.0)print*, 'Graftpoints OK'
       if(rank.eq.0)print*, ' '

      call creador             ! creates chains
      if(rank.eq.0)print*, ' '
      if(rank.eq.0)print*, 'Creador OK'

      call lookup_kai          ! sets up kai matrixes in lattice
      if(rank.eq.0)print*, 'Lookup_kai OK'
       if(rank.eq.0)print*, ' '

      call mainsolver(counter)          ! perform calculation and save results to disk

      call MPI_FINALIZE(ierr) ! finalize MPI

      end

