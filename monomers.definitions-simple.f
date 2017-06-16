      subroutine monomer_definitions
     
      use mparameters
      use mparameters_monomer

      implicit none

      include 'mpif.h'
      include 'MPI.h'

      N_poorsol = 2 ! number of different kais
      N_monomer = 4 

      ALLOCATE (st_matrix(N_poorsol, N_poorsol)) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
      ALLOCATE (zpol(N_monomer))    ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      ALLOCATE (hydroph(N_monomer)) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
      ALLOCATE (pKa(N_monomer), Ka(N_monomer), K0(N_monomer))
      ALLOCATE (henergy(N_poorsol))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      st_matrix(1,1)= kaiHO !1.1 !2.70
      st_matrix(2,2)= kaiHI !0.0 !1.0
      st_matrix(2,1)= 0 !sqrt(st_matrix(1,1)*st_matrix(2,2))
      st_matrix(1,2)=st_matrix(2,1)

      henergy(1) =  pairst !st_matrix(1,1)
      henergy(2) =  st_matrix(1,2)

! Segment type 1, neutral , always

      zpol(1) = 0
      hydroph(1) = 1
      pKa(1) = 1 ! set any number if zpol = 0...

! Segment type 1, neutral , always

      zpol(2) = 0
      hydroph(2) = 2
      pKa(2) = 1 ! set any number if zpol = 0...

! Segment type 2, hydrophobic or charged

      zpol(3) = zpolHO1 ! 1
      hydroph(3) = 0
      pKa(3) = pKaHO1 ! 5.0 ! set any number if zpol = 0....

! Segment type 3, hydrophobic or charged

      zpol(4) = zpolHO2
      hydroph(4) = 0
      pKa(4) = pKaHO2

      end

