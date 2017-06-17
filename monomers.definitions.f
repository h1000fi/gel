      subroutine monomer_definitions

      use mparameters_monomer
      use mparameters   

      implicit none

      include 'mpif.h'
      include 'MPI.h'

      integer i, j

      N_poorsol = 2 ! number of different kais
      N_monomer = 4 

      ALLOCATE (st_matrix(N_poorsol, N_poorsol)) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
      ALLOCATE (zpol(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      ALLOCATE (hydroph(N_monomer)) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
      ALLOCATE (pKa(N_monomer), Ka(N_monomer), K0(N_monomer))
      ALLOCATE (henergy(N_poorsol))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      st_matrix(:,:)= 0.0

      do i = 1, N_poorsol
        do j = 1, N_poorsol
          st_matrix(i,j) = 0.0
        enddo
      enddo

      st_matrix(1,1)= 0.0 !1.1 !2.70
      st_matrix(2,2)= kaiHO !0.0 !1.0

      st_matrix(1,2)= kaiHO
      st_matrix(2,1)= kaiHO

      henergy(1) =  pairst !pairst !st_matrix(1,1)
      henergy(2) =  0 !kaiHI !st_matrix(1,2)

! Segment type 1, pairing monomers

      zpol(1) = 0
      hydroph(1) = 1
      pKa(1) = 1 ! set any number if zpol = 0...

! Segment type 2, spacers

      zpol(2) = 0
      hydroph(2) = 2
      pKa(2) = 1 ! set any number if zpol = 0...

! Segment type 3, charged group

      zpol(3) = zpolHO1 ! 1
      hydroph(3) = 0
      pKa(3) = pKaHO1 ! 5.0 ! set any number if zpol = 0....

! Segment type 4, charged group

      zpol(4) = zpolHO2
      hydroph(4) = 0
      pKa(4) = pKaHO2

      end

