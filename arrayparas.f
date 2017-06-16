      module mprotein
      real*8, allocatable :: protein(:,:)
      real*8, allocatable :: proteinq(:,:)
      real*8, allocatable :: proteinC(:)
      real*8, allocatable :: proteinqC(:)
      real*8, allocatable :: proteinh(:, :)
      real*8, allocatable :: proteinhC(:)
      integer Ppos ! protein position in delta units
      integer PdimR ! protein radius in delta units      
      integer PdimH ! protein length in delta units, only for type 2
      integer Ptype
      real*8 prot_vol ! internal protein volume fraction
      real*8 prot_q   ! total protein charge
      real*8 prot_q2   
      real*8 hst
      real*8 hst2
      integer hrange
      real*8 pKaP, KaP, K0P
      integer weakP

      endmodule mprotein

      module posmk
      real*8, allocatable :: current(:,:)
      integer*2, allocatable :: firstcell(:,:,:)
      integer*2, allocatable :: nextbead(:)
      endmodule posmk

      module mparameters_chain
      integer N_chains ! number of chains, should be equal to size in the actual implementation
      integer, allocatable :: long(:) ! length of chain i, see chains.definitions.f
      real*8, allocatable :: chainsperdelta(:) ! number of equivalent chains grafted at position i
      real*8, allocatable :: zposition(:)      ! z position of chain i
      integer, allocatable :: segtype(:,:)
      real*8, allocatable :: endtoend(:)
      real*8 endtoend_av
      end module mparameters_chain

      module mparameters_monomer
      integer N_poorsol ! number of different kais
      integer N_monomer ! number of different monomer types
      real*8, allocatable :: st_matrix(:,:) ! interaction between monomer types in fraction of st, scaled by st-scale during running....
      integer, allocatable :: zpol(:)  ! charge of monomer segment: 1: base, -1: acid, 0:neutral
      integer, allocatable :: hydroph(:) ! 0: hydrophilic, 1 < x < N_poorsol, type of poor solvent
      real*8, allocatable ::  pKa(:), Ka(:), K0(:)
      real*8, allocatable :: henergy(:) 
      integer monomerz
      real*8 monomerpKa
      endmodule mparameters_monomer

      module mkinsol
      double precision, allocatable :: pp(:)
      endmodule 

      module mlookup
      integer*4, allocatable :: rp(:), rm(:), zp(:), zm(:)
      integer*4, allocatable :: indexa(:,:) ! 1 = R, 2 = Z
      endmodule mlookup

      module mchains
      integer*4, allocatable :: fs(:)
      real*8, allocatable :: pbias(:,:)
      integer*1, allocatable :: displ(:,:)
      endmodule mchains
       
      module mgraft
      real*8 posgraft(2) ! grafting point positions coord, GP
      endmodule mgraft

      module mkai
      real*8 st
      real*8, allocatable :: Xu(:,:,:,:)
      integer, allocatable :: nXu(:,:), Xulist_cell(:,:,:)
      real*8, allocatable :: Xulist_value(:,:,:)
      endmodule mkai

      module mncells
      integer ncells
      endmodule mncells

      module mporesystem
      integer*4, allocatable :: matriz(:, :)
      endmodule mporesystem

      module mrands
      real*8 rands
      external rands
      integer seed
      endmodule mrands

      module msegme
      real*8, allocatable :: in1(:,:)  ! Posicion del segmento
      integer, allocatable :: celda(:)
      integer flag
      endmodule msegme

      module mvariables
! saveindex
      integer saveindex

! Volumen
      real*8 vsol
      real*8 vpol
      real*8 vsalt
! Volumen fraction
      real*8, allocatable :: avpol(:, :, :)
      real*8, allocatable :: qtot(:) 
      real*8, allocatable :: psi2(:) 
      real*8, allocatable :: xsol(:) 
      real*8, allocatable :: xtotal2(:, :) 
      real*8, allocatable :: xpos(:) ! pos ion
      real*8, allocatable :: xneg(:) ! neg ioni
      real*8, allocatable :: xHplus(:) ! H+
      real*8, allocatable :: xOHmin(:) ! OH-
      real*8, allocatable :: aveP(:) ! average HO pair within the cutoff
      real*8, allocatable :: aveC(:) ! average absolute charge within the cutoff
      real*8, allocatable :: qtot_amp(:) ! absolute charge within the cutoff
      real*8, allocatable :: Spair(:) ! pair entropy compensation
      real*8, allocatable :: Fpair(:) ! pair free energy
      real*8, allocatable :: Rpair(:) ! pair rate as a reaction
      real*8, allocatable :: Fpair_tot(:) ! total pair free energy as a reaction
! Bulk
      real*8 xsolbulk ! volume fraction solvent in bulk
      real*8 xposbulk, xposbulk2, xnegbulk, xsalt,csalt
      real*8 expmupos, expmuneg
      real*8 chainperdelta
      real*8 sigma, sigmadelta,  mupol1, mupol2, mupold
      real*8 sigma1, sigma2
      integer*1 sigmaflag
! Charge
      real*8 zpos, zneg
      real*8 sigmaq
      real*8 lb
      real*8 constq
      real*8 betae
! Weak pol
      real*8, allocatable :: fdis(:, :)
      real*8, allocatable :: fdisP(:)
      real*8 Kw, pKw, pHbulk, expmuHplus, expmuOHmin
      real*8 xHplusbulk, cHplus, pOHbulk, xOHminbulk, cOHmin
      integer strongpoly
      real*8 shift
      real*8 norma
      integer iter
! chains
      real*8 lnpro, pro, q, sumprolnpro
! solver
      integer savetodisk_flag
! others 
      integer infile 
      endmodule mvariables

      module mparameters
      real*8 sts(100)
      integer nst

      integer maxlong
      integer poretype
      integer cadenastype
      integer dimR, dimZ         ! system total size is dimR x dimZ
      real*8 delta
      real*8 lseg
      integer polyLen
      integer zpolHO1
      real*8 pKaHO1
      integer zpolHO2
      real*8 pKaHO2
      real*8 kaiHO
      real*8 kaiHI
      real*8 paircut
      real*8 pairst
      real*8 pHstart,pHend,looper
      real*8 pHstep
      real*8 kC
      real*8 pairsize
      real*8 hring
      real*8 oval
      integer nkp
      real*8 kp
      real*8 kps(100) 
      real*8 polyCov
      real*8 deltaZp
      real*8 zspace
      real*8 tol
      real*8 kBias
      real*8 Rbias
      real*8 Zbias
      integer cuantas
      integer, allocatable :: newcuantas(:)
      real*8 Na
      parameter (Na=6.02d23)    ! Avogadro's number
      real*8 pi
      parameter (pi=3.14159)    ! pi
      integer CdimR
      integer CdimZ
      integer CdimZmin
      integer BdimR
      integer TdimR
      integer RdimZ
      real*8 Curvature
      integer savetodisk_type
      integer readkai
      integer save_every
      real*8 dseg
      integer nearbonds
      integer mcube
      integer ncha_max
      real*8 b2,d2                    ! here calculated once from lseg and dseg (keep unchanged)
      integer wantedCPUsecs
      integer calq
      real*8 qprob0
      integer MCsteps
      integer temp_long
      integer readseed
      endmodule mparameters
