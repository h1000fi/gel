      integer N_chains ! number of different chains types
      parameter (N_chains = 20)

      integer cadenastype
      parameter (cadenastype = 1) ! 1: old cadenas routine, 2: MK's new cadenas routine

      integer segtype(N_chains, 
