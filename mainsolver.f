      subroutine mainsolver(counter)

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

      integer i,j, k, ic, ii,ir,iz, im, ikp
      integer counter
      real*8 stok

!      real*8 pHstep,pHstart,pHend,looper
      real*8 avpol_all(dimR*dimZ)

      real*8 charge_all
      real*8 temp
      real*8 temp1

      character basura
      character*24 filename 
      character*5  title

      real*8 error
      real*8 sumpol, fmedio     

C-----  solving variables  -----------

      real*8 algo

      integer flag
       
      real*8 errel, fnorm
      integer n, itmax

      integer cc,ccc,cccc


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

! Variables

      shift = 1d0

      lb = 0.714 ! bjerrum lenght in nm

      zpos = 1.0
      zneg = -1.0
      
      vsol = 0.030     
      vsalt=((4.0/3.0)*pi*(0.27)**3)/vsol  ! volume salt in units of vsol 0.2=radius salt  
      vpol= 0.095/vsol! ((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume polymer segment in units of vsol 

      constq=delta*delta*4.0*pi*lb/vsol   ! multiplicative factor in poisson eq  

      pKw = 14
      Kw = 1.0e-14

      error = 1e-6 ! para comparar con la norma...

      betae = 38.94             ! beta * e


      n=ncells
      
      errel=1d-6
      itmax=200

C-------------------------------------------------------
C MAIN LOOP
C-------------------------------------------------------

      flag = 0
      ccc = 1
      cc = 1

      do while(ccc.le.nst) ! loop in st

 555  if(rank.eq.0)print*, 'ccc', ccc, 'of', nst

       st = sts(ccc)


 257  cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
      xHplusbulk = (cHplus*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol
      pOHbulk= pKw -pHbulk
      cOHmin = 10**(-pOHbulk)   ! concentration OH- in bulk
      xOHminbulk = (cOHmin*Na/(1.0d24))*(vsol)  ! volume fraction H+ in bulk vH+=vsol  

      xsalt=(csalt*Na/(1.0d24))*(vsalt*vsol)   ! volume fraction salt,csalt in mol/l 

      if(pHbulk.le.7) then  ! pH<= 7
            xposbulk=xsalt/zpos
            xnegbulk=
     &  - xsalt/zneg + (xHplusbulk-xOHminbulk)*vsalt ! NaCl+ HCl  
      else                  ! pH >7 
            xposbulk=xsalt/zpos +(xOHminbulk-xHplusbulk)*vsalt ! NaCl+ NaOH   
            xnegbulk=-xsalt/zneg 
      endif

         xsolbulk=1.0 -xHplusbulk -xOHminbulk -
     &        xnegbulk -xposbulk 

         do im = 1, N_monomer
         Ka(im)=10**(-pKa(im))
         select case (zpol(im))
         case (-1) ! acid
         K0(im) = (Ka(im)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Ka
         case (1) ! base
         K0(im) = ((Kw/Ka(im))*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Kb 
         end select
         enddo

         KaP=10**(-pKaP)
         if (prot_q < 0) then
         K0P = (KaP*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Ka
         else ! base
         K0P = ((Kw/KaP)*vsol/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant, Kb 
         endif
 


         expmupos = xposbulk /xsolbulk**vsalt
         expmuneg = xnegbulk /xsolbulk**vsalt
         expmuHplus = xHplusbulk /xsolbulk   ! vsol = vHplus 
         expmuOHmin = xOHminbulk /xsolbulk   ! vsol = vOHmin 

         if(rank.eq.0) then
         print*, 'Bulk composition / volume fracion'
         print*, 'Cations ', xposbulk
         print*, 'Anions  ', xnegbulk
         print*, 'H+      ', xHplusbulk
         print*, 'OH-     ', xOHminbulk
         print*, 'Solvent ', xsolbulk
         print*, 'Charge / vsol units',
     & xposbulk/vsalt*zpos+xnegbulk/vsalt*zneg+xHplusbulk-xOHminbulk


         endif

c-------------------------------------------------------------------------------
c
c INITIAL GUESS 
c
c-------------------------------------------------------------------------------

      if(infile.eq.2) then
      if(ccc.ne.1) then
      do i = 1, (2+N_poorsol)*ncells ! dejar aca, asegura guess correcto para infile = 2 (cc > 1)
      xg1(i) = xflag(i)     
      x1(i) = xflag(i)
      enddo
      endif
      endif
c
      if(infile.eq.0) then
         do i=1,n
           xg1(i)=xsolbulk
           x1(i)=xsolbulk
         enddo
         do i=n+1, n*2
           xg1(i)=0.0d0
           x1(i)=0.0d0
         enddo
         do i=(1+N_poorsol)*ncells+1, n*(2+N_poorsol)
           xg1(i)=0.0
           x1(i)=0.0
         enddo
      endif
 
      if(infile.eq.1) then
         open(unit=45, file='in.txt')
         do i=1, (2+N_poorsol)*ncells
            read(45,*) xg1(i)
            x1(i) = xg1(i)
         enddo
         close(45)
      endif 

      if(infile.eq.5) then
         write(filename,'(A4, I3.3, A4)')'out.', ccc, '.dat'
         print*, 'reading ', filename
         open(unit=45, file=filename)
         do i = 1, (2+N_poorsol)*ncells
            read(45, *) x1(i)
         enddo
         xg1 = x1
         close(45)
      endif

C--------------------------------------------------------------
C               +++ SOLVER +++ 
C--------------------------------------------------------------

!      pairst = 0.0

!      st = sts(1)
      kp = 1.0d10+kps(1)
      do ikp = 1, nkp
       do while (kp.ne.kps(i))
        kp = kps(ikp)
        if(rank.eq.0)write(stdout,*)'Switch to kp = ', kp

        call presolver
        call solver(x1,xg1,ier)
          
! Retrives xsol y psi (NOT COMMON!)

        do iC = 1, n

            xsol(iC)=x1(iC)
            psi2(iC)=x1(iC+n)
  
            do i = 1, N_poorsol
               xtotal2(i, iC)=x1(iC+(i+1)*n)
            enddo

        enddo

! Solution converged, store xflag

        do i = 1, (2+N_poorsol)*n
           xflag(i) = x1(i) ! xflag will be used as input for the next iteration
        end do

        if(infile.ne.5)infile = 2 ! no vuelve a leer infile
           stok = st

        if(flag.eq.1) then  ! recorvers from an error
           if(rank.eq.0) then
              print*, 'Recovers from an error'
              print*, 'OK', stok
           endif
           flag = 0
           goto 555
        endif

! Analysis

        call calc_free_energy(dfloat(ccc), dfloat(cc))
        call saveresults(ikp, cc)

        if(rank.eq.0) then

           print*, 'Bulk composition / volume fracion'
           print*, 'Cations ', xposbulk
           print*, 'Anions  ', xnegbulk
           print*, 'H+      ', xHplusbulk
           print*, 'OH-     ', xOHminbulk
           print*, 'Solvent ', xsolbulk

           charge_all = 0
           do iC = 1, n
             do im = 1, N_monomer
                do ii = 1, N_chains
                   charge_all = charge_all + avpol(im, ii, iC)
     &             *fdis(im, iC)/vpol*(dfloat(indexa(iC,1))-0.5)*2*pi
     &             *delta**3/vsol*zpol(im)
                enddo
             enddo
           enddo
           print*, 'polymer charge', charge_all

           avpol_all(:) = 0
           do im = 1, N_monomer
              do ii = 1, N_chains
                 avpol_all(:) = avpol_all(:) + avpol(im, ii, :)
              enddo
           enddo

           if(infile.ne.5) then
              write(filename,'(A3, I3.3, A4)')'in.', ikp , '.dat'
              open(unit=45, file=filename)
              do i = 1, (2+N_poorsol)*ncells
                 write(45, *) x1(i)
              enddo
              close(45)
           endif

           title = 'steps'
           call savetodisk(avpol_all, title, ikp ,ccc)

           do i=1, (2+N_poorsol)*ncells
              xg1(i) = x1(i)
           enddo

        endif

      enddo ! kp
      enddo


      if(rank.eq.0) then ! only master can save to disk

! saves infile

         print*, 'finished!!!!'

         if(infile.ne.5) then
            write(filename,'(A4, I3.3, A4)')'out.', ccc, '.dat'
            open(unit=45, file=filename)
            do i = 1, (2+N_poorsol)*ncells
               write(45, *) x1(i)
            enddo
            close(45)
         endif

         cc = saveindex 
!         call saveresults(ccc, cc) ! save results

      endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      call calc_free_energy(dfloat(ccc), dfloat(cc)) ! calculate free energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ccc = ccc + 1

      enddo ! ccc

      end
