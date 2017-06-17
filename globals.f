      subroutine globals

      use mparameters
      use mparameters_monomer
      use mvariables
      use mprotein

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries


! Input related variables
      character (len=100)  buffer,label
      integer pos
      integer, parameter :: fh = 15
      integer ios 
      integer line, linemax
      integer i

      character(len=50) :: filename = 'DEFINITIONS.txt'


! Control file variables

      line = 0
      ios = 0

      open(fh, file=filename)

      if(rank.eq.0)print*, 'Reading parameters from ', filename

! ios is negative  if an end of record condition is encountered or if
! an endfile condition was detected.  It is positive  if an error was
! detected.  ios is zero otherwise.
     
     
      do while (ios == 0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios == 0) then
        line = line + 1

! Find the first instance of whitespace.  Split label and data.
 

        pos = scan(buffer, ' ')
 
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)

       case ('hst')
           read(buffer, *, iostat=ios) hst
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

       case ('hst2')
           read(buffer, *, iostat=ios) hst2
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

       case ('Ptype')
           read(buffer, *, iostat=ios) Ptype
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

       case ('hrange')
           read(buffer, *, iostat=ios) hrange
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

       case ('saveindex')
           read(buffer, *, iostat=ios) saveindex
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('poretype')
           read(buffer, *, iostat=ios) poretype
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('cadenastype')
           read(buffer, *, iostat=ios) cadenastype
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('dimR')
           read(buffer, *, iostat=ios) dimR
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('dimZ')
           read(buffer, *, iostat=ios) dimZ
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('delta')
           read(buffer, *, iostat=ios) delta
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('lseg')
           read(buffer, *, iostat=ios) lseg
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('polyLen')
           read(buffer, *, iostat=ios) polyLen
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('zpolHO1')
           read(buffer, *, iostat=ios) zpolHO1
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('pKaHO1')
           read(buffer, *, iostat=ios) pKaHO1
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('zpolHO2')
           read(buffer, *, iostat=ios) zpolHO2
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('pKaHO2')
           read(buffer, *, iostat=ios) pKaHO2
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('kaiHO')
           read(buffer, *, iostat=ios) kaiHO
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('kaiHI')
           read(buffer, *, iostat=ios) kaiHI
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('paircut')
           read(buffer, *, iostat=ios) paircut
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('pairst')
           read(buffer, *, iostat=ios) pairst
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('kC')
           read(buffer, *, iostat=ios) kC
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('pairsize')
           read(buffer, *, iostat=ios) pairsize
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('avepair')
           read(buffer, *, iostat=ios) avepair
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('hguess')
           read(buffer, *, iostat=ios) hguess
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('hring')
           read(buffer, *, iostat=ios) hring
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('oval')
           read(buffer, *, iostat=ios) oval
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('nkp')
           read(buffer, *, iostat=ios) nkp
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)
           do i = 1, nkp
           read(fh,*)kps(i)
           enddo

        case ('polyCov')
           read(buffer, *, iostat=ios) polyCov
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('deltaZp')
           read(buffer, *, iostat=ios) deltaZp
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('zspace')
           read(buffer, *, iostat=ios) zspace
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

!        case ('tol')
!           read(buffer, *, iostat=ios) tol
!         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('kBias')
           read(buffer, *, iostat=ios) kBias
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('Rbias')
           read(buffer, *, iostat=ios) Rbias
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('Zbias')
           read(buffer, *, iostat=ios) Zbias
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('cuantas')
           read(buffer, *, iostat=ios) cuantas
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)
         print*, label, buffer, 'OK'

        case ('CdimR')
           read(buffer, *, iostat=ios) CdimR
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)
         print*, label, buffer, 'OK'

        case ('CdimZ')
           read(buffer, *, iostat=ios) CdimZ
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('CdimZmin')
           read(buffer, *, iostat=ios) CdimZmin
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('BdimR')
           read(buffer, *, iostat=ios) BdimR
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('TdimR')
           read(buffer, *, iostat=ios) TdimR
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('RdimZ')
           read(buffer, *, iostat=ios) RdimZ
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('Curvature')
           read(buffer, *, iostat=ios) curvature
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('savetodisk_type')
           read(buffer, *, iostat=ios) savetodisk_type
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('save_every')
           read(buffer, *, iostat=ios) save_every
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('dseg')
           read(buffer, *, iostat=ios) dseg
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('nearbonds')
           read(buffer, *, iostat=ios) nearbonds
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('mcube')
           read(buffer, *, iostat=ios) mcube
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('ncha_max')
           read(buffer, *, iostat=ios) ncha_max
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('wantedCPUsecs')
           read(buffer, *, iostat=ios) wantedCPUsecs
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('calq')
           read(buffer, *, iostat=ios) calq
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('qprob')
           read(buffer, *, iostat=ios) qprob0
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('csalt')
           read(buffer, *, iostat=ios) csalt
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('pHbulk')
           read(buffer, *, iostat=ios) pHbulk
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('infile')
           read(buffer, *, iostat=ios) infile
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('st')
           read(buffer, *, iostat=ios) nst
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)
           read(fh, *), (sts(i), i=1, nst)
         if(rank.eq.0)print*, 'sts: ', sts

        case ('savetodisk_flag')
           read(buffer, *, iostat=ios) savetodisk_flag
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('sigma')
           read(buffer, *, iostat=ios) sigma
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('maxlong')
           read(buffer, *, iostat=ios) maxlong
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('MCsteps')
           read(buffer, *, iostat=ios) MCsteps
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

        case ('readkai')
           read(buffer, *, iostat=ios) readkai
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('Ppos')
          read(buffer, *, iostat=ios) Ppos
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('PdimR')
          read(buffer, *, iostat=ios) PdimR
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('PdimH')
          read(buffer, *, iostat=ios) PdimH
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('prot_vol')
          read(buffer, *, iostat=ios) prot_vol
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('prot_q')
          read(buffer, *, iostat=ios) prot_q
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('prot_q2')
          read(buffer, *, iostat=ios) prot_q2
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('weakP')
          read(buffer, *, iostat=ios) weakP
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('pKaP')
          read(buffer, *, iostat=ios) pKaP
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('readseed')
          read(buffer, *, iostat=ios) readseed
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('monomerz')
          read(buffer, *, iostat=ios) monomerz
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)

         case ('monomerpKa')
          read(buffer, *, iostat=ios) monomerpKa
         if(rank.eq.0)print*,'Set ',trim(label),' = ',trim(buffer)



        end select

        endif 
        enddo     

        ncha_max=cuantas/1000+200
       if(rank.eq.0)print*,'Set ncha_max = ', ncha_max

        b2 = lseg*lseg
       if(rank.eq.0)print*,'Set b2 = ', b2
        d2 = dseg*dseg
       if(rank.eq.0)print*,'Set d2 = ', d2

       end

       subroutine array_alloc

      use mparameters
      use mchains
      use mlookup
      use mkai
      use mporesystem
      use msegme
      use mvariables
      use mkinsol
      use mparameters_chain
      use mparameters_monomer
      use posmk
      use mprotein

      implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! ALLOCATE ARRAYS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

! chains

        ALLOCATE (fs(cuantas))
        ALLOCATE (pbias(N_chains,cuantas))        
        ALLOCATE (displ(cuantas, int(maxlong/2)))
        ALLOCATE (endtoend(cuantas))

! lookup

      ALLOCATE (rp(dimR*dimZ))
      ALLOCATE (rm(dimR*dimZ))
      ALLOCATE (zp(dimR*dimZ))
      ALLOCATE (zm(dimR*dimZ))
      ALLOCATE (indexa(dimR*dimZ,2)) ! 1 = R, 2 = Z

! kai
     
      ALLOCATE (Xu(N_poorsol, dimR, dimR, -3:3))
      ALLOCATE (nXu(N_poorsol, dimR*dimZ))
      ALLOCATE (Xulist_cell(N_poorsol, dimR*dimZ, 50)) !used to be 10
      ALLOCATE (Xulist_value(N_poorsol, dimR*dimZ, 50)) !used to be 10

! pore.system

      ALLOCATE (matriz(0:dimR+1, 0:dimZ+1))

! segme

      ALLOCATE (in1(maxlong,3))
      ALLOCATE (celda(maxlong))

! variables

      ALLOCATE (avpol(N_monomer, N_chains, dimR*dimZ))
      ALLOCATE (qtot(dimR*dimZ))
      ALLOCATE (psi2(dimR*dimZ))
      ALLOCATE (xsol(dimR*dimZ))
      ALLOCATE (xtotal2(N_poorsol, dimR*dimZ+2))
      ALLOCATE (xpos(dimR*dimZ))
      ALLOCATE (xneg(dimR*dimZ))
      ALLOCATE (xHplus(dimR*dimZ))
      ALLOCATE (xOHmin(dimR*dimZ))
      ALLOCATE (aveP(dimR*dimZ))
      ALLOCATE (aveC(dimR*dimZ))
      ALLOCATE (qtot_amp(dimR*dimZ))
      ALLOCATE (Spair(dimR*dimZ))
      ALLOCATE (Fpair(dimR*dimZ))
      ALLOCATE (Rpair(dimR*dimZ))
      ALLOCATE (Fpair_tot(dimR*dimZ))

      ALLOCATE (fdis(N_monomer, dimR*dimZ))
      ALLOCATE (fdisP(dimR*dimZ))
      ALLOCATE (newcuantas(N_chains))

! kinsol

      ALLOCATE (pp((2+N_poorsol)*dimR*dimZ))

! posmk

      ALLOCATE (firstcell(-mcube:mcube,-mcube:mcube,-mcube:mcube))
      ALLOCATE (current(maxlong, 3))
      ALLOCATE (nextbead(maxlong))

! protein

      ALLOCATE (protein(dimR, dimZ))
      ALLOCATE (proteinq(dimR, dimZ))
      ALLOCATE (proteinh(dimR, dimZ))
      ALLOCATE (proteinC(dimR*dimZ))
      ALLOCATE (proteinqC(dimR*dimZ))
      ALLOCATE (proteinhC(dimR*dimZ))

      end




