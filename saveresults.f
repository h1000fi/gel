      subroutine saveresults(cc, ccc) ! save results

      use mlookup
      use mparameters ! parametros     
      use mncells
      use mkai
      use mporesystem
      use mvariables
      use mparameters_chain
      use mprotein
      use mparameters_monomer

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries

      integer ccc, i, ii, jj, j, cc, im, is,iC,iR
      real*8 avpol_all(dimR*dimZ)
      real*8 avpol_temp(dimR*dimZ)
      real*8 fdis_temp(dimR*dimZ)
      character*5  title
      character*24 filename 
      real*8 temp
      real*8 temp1

C----------------------------------------------------------
C  OUTPUT!
C----------------------------------------------------------

      print*, 'SAVETODISK', savetodisk_flag

      if(rank.eq.0) then ! only master can save to disk

!!!!!!!!!!!!!!!!!!! saves files  !!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Polymer, sum

      if(savetodisk_flag.ne.10) then

      avpol_all(:) = 0
      do im = 1, N_monomer
      do ii = 1, N_chains
      avpol_all(:) = avpol_all(:) + avpol(im, ii, :)
      enddo
      enddo

      title = 'avpol'
      call savetodisk(avpol_all, title, cc ,ccc)

      endif

      if((savetodisk_flag.eq.3).or.savetodisk_flag.eq.0) then

! Fraction hydrophobic - hydrophylic

      do j = 1, N_poorsol

      avpol_temp(:) = 0.0
      do im = 1, N_monomer

      if(hydroph(im).eq.j) then 

         do ii = 1, N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
         enddo

      endif

      enddo ! im

      write(title,'(A4, I1.1)')'hpho',j
      call savetodisk(avpol_temp, title, cc ,ccc)

      enddo ! j
       

      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(hydroph(im).eq.0) then

         do ii = 1, N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
         enddo

      endif 
      enddo ! im

      write(title,'(A5)')'hphil'
      call savetodisk(avpol_temp, title, cc ,ccc)

! pos, neg, neutral

      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(zpol(im).eq.0) then

         do ii = 1, N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
         enddo

      endif
      enddo ! im

      write(title,'(A5)')'neutr'
      call savetodisk(avpol_temp, title, cc ,ccc)

!

      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(zpol(im).eq.1) then

         do ii = 1, N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
         enddo

      endif
      enddo ! im

      write(title,'(A5)')'posit'
      call savetodisk(avpol_temp, title, cc ,ccc)

!

      avpol_temp(:) = 0.0
      do im = 1, N_monomer
      if(zpol(im).eq.-1) then

         do ii = 1, N_chains
         avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)
         enddo

      endif
      enddo ! im

      write(title,'(A5)')'negat'
      call savetodisk(avpol_temp, title, cc ,ccc)

      endif
      if((savetodisk_flag.eq.2).or.(savetodisk_flag.eq.0)) then

! all polymers


      do ii = 1, N_chains

      avpol_temp(:) = 0.0

      do im = 1, N_monomer

      avpol_temp(:) = avpol_temp(:)+avpol(im, ii, :)

      enddo

      write(title,'(A3, I2.2)')'cha',ii
      call savetodisk(avpol_temp, title, cc ,ccc)

      enddo

      endif



      if((savetodisk_flag.eq.0).or.(savetodisk_flag.eq.5)) then


! Total charge

      title = 'qtodo'
      call savetodisk(qtot, title, cc ,ccc)

! Solvente

c      title = 'avsol'
c      call savetodisk(xsol, title, cc, ccc)

! Cationes

      title = 'avpos'
      call savetodisk(xpos, title, cc, ccc)

! Aniones

      title = 'avneg'
      call savetodisk(xneg, title, cc, ccc)

! H+

      title = 'avHpl'
      call savetodisk(xHplus, title, cc, ccc)

! OH-

      title = 'avOHm'
      call savetodisk(xOHmin, title, cc,ccc)

! fdis
      fdis_temp(:) = 0

      do im = 1, N_monomer

      fdis_temp(:) = fdis(im,:)

      write(title,'(A3, I1.1)')'fdis',im
      call savetodisk(fdis_temp, title, cc ,ccc)

      enddo


! fdisP

      title = 'fdisP_'
      call savetodisk(fdisP, title, cc ,ccc)

! Pair rate

      title = 'Rpair'
      call savetodisk(Rpair, title, cc, ccc)

! Pair energy

      title = 'Fpair'
      call savetodisk(Fpair, title, cc, ccc)

! abs charge

      title = 'qtot_amp'
      call savetodisk(qtot_amp, title, cc, ccc)


      endif ! savetodisk flag


       if(((savetodisk_flag.eq.0).or.(savetodisk_flag.eq.5)
     & .or.(savetodisk_flag.eq.11))) then
! protein

      title = 'protC'
      call savetodisk(proteinC, title, cc ,ccc)

! Potencial electrostatico

      title = 'poten'
      call savetodisk(psi2, title, cc, ccc)
      endif

!!!!!!!!!!!!!!!!! Informacion del sistema !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         write(filename,'(A7, I3.3, A1, I3.3, A4)')'system.',ccc,'.',cc,
     &     '.dat'
         open (unit=510, file=filename)

         write(510,*)'st          = ',st ! residual size of iteration vector
         write(510,*)'fnorm       = ',norma ! residual size of iteration vector
         write(510,*)'delta       = ',delta
         write(510,*)'vsol        = ',vsol
         write(510,*)'vsol        = ',vpol
         write(510,*)'vsalt       = ',vsalt*vsol
         write(510,*)'csalt       = ',csalt
         write(510,*)'pHbulk      = ',pHbulk
         write(510,*)'pKw         = ',pKw
         write(510,*)'zpos        = ',zpos
         write(510,*)'zneg        = ',zneg
         write(510,*)'cuantas     = ',cuantas
         write(510,*)'newcuantas     = ',newcuantas
         write(510,*)'iterations  = ',iter
         write(510,*)'CHAINS PER DELTA  = ',chainsperdelta
         write(510,*)'system length / nm = ', dimZ*delta
         write(510,*)'sites with grafted chains = ', CdimZ

!!!!!!!!!!!!!!!!!!!!!! CHECK TOTAL NUMBER OF CHAINS !!!!!!!!!!!!!!!!!
        temp = 0.0
        
        do iC = 1, ncells
        do im = 1, N_monomer
        do ii = 1, N_chains

        iR = indexa(iC, 1)
        temp = temp + avpol(im,ii,iC)
     & *(dfloat(iR)-0.5)*delta*2*pi*delta**2

        enddo
        enddo
        enddo

        temp = temp/vpol/vsol
        write(510,*)'Total number of segments in system', temp

        temp = 0.0
        do ii = 1, N_chains
        temp = temp + long(ii)*chainsperdelta(ii)
        enddo
        write(510,*)'Total number of segments in system should be',
     &  temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!! PARTICLE DISOC !!!!!!!!!!!!!!!!!
        if(weakP.eq.1) then

        temp = 0.0
        temp1 = 0.0

        do iC = 1, ncells
        iR = indexa(iC, 1)
        temp = temp + proteinqC(iC)
     & *(dfloat(iR)-0.5)*delta*2*pi*delta**2
        enddo

        write(510,*)'Total number of cargo ionizable charges', temp


        do iC = 1, ncells
         
        iR = indexa(iC, 1)
        temp1 = temp1 + proteinqC(iC)*fdisP(iC)
     & *(dfloat(iR)-0.5)*delta*2*pi*delta**2

        enddo

        write(510,*)'Total number of cargo ionized charges', temp1

        write(510,*)'Fraction of cargo ionized charges', temp1/temp
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




         close(510)
 
         endif ! rank
        
      end
