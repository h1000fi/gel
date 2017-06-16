      subroutine fkfun(x,f,ier)

      use mlookup
      use mchains
      use mparameters
      use mncells
      use mkai
      use mporesystem
      use mvariables
      use mparameters_chain
      use mparameters_monomer
      use mprotein

      implicit none
 
      include 'mpif.h' ! MPI libraries
      include 'MPI.h' ! MPI libraries

      double precision Factorcurv

      integer ntot
      real*8 x((2+N_poorsol)*ncells)
      real*8 f((2+N_poorsol)*ncells)

      real*8  protemp, protemp1, avePnorm, xave, smoothstep

      integer i,j, k, ix, iy, iC, ii, aR, temp, aZ, jj, im, imc, iic, jp

      real*8 temp2
      real*8 hdistance, hfactor
      real*8 apair, bpair, cpair

      real*8 psi(ncells+2) ! psi se define asi para evitar problemas al construir las fs
      real*8 xtotal(N_poorsol, ncells+2) ! psi se define asi para evitar problemas al construir las fs

      real*8 xh(ncells) ! psi se define asi para evitar problemas al construir las fs
 
! Kinsol

      integer*8 neq
      integer*4 ier

! MPI

      integer spp ! sitios por procesador

      integer cuantas_p(2)

! solving auxiliary variables

      real*8 avpol_temp(N_monomer, N_chains, ncells)
      real*8 avpol_tosend(N_monomer, N_chains, ncells)

      real*8 xpot(N_monomer,ncells)
       
     
! bit operations

      integer*1 displ_one(0:7)
      integer*4 pC, pZ, pR

      integer*1 displacement(maxlong+1,2)
      integer*1 binary(int(maxlong/2))

C-----------------------------------------------------

! Jefe

      if(rank.eq.0) then ! llama a subordinados y pasa vector x
       flagsolver = 1

       CALL MPI_BCAST(flagsolver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,err)

       CALL MPI_BCAST(x, (2+N_poorsol)*ncells , MPI_DOUBLE_PRECISION,0
     &  , MPI_COMM_WORLD,err)

      endif


! Recupera xh y psi

      ntot = ncells

            do iC = 1,ntot

            xh(iC)=x(iC)
            psi(iC)=x(iC+ntot)

            do i=1, N_poorsol
            xtotal(i,iC)=x(iC+(1+i)*ntot) ! v-frac of the good solvent polymers
            enddo

            enddo
  
      qtot_amp = 0.0      

! Condiciones de borde psi 
 
! This implictly considers SIGMAQ = 0
! 


! bulk
      psi(ncells+1) = 0.0
      psi(ncells+2) = 0.0 

! las condiciones en la pared y condicion de simetria esta en el lookuptable

! Fracciones de volumen inicial y fdis

           avpol=0
           avpol_tosend = 0.0
           fdis = 0.0
           fdisP = 0.0

           do iC = 1, ncells
     
           xpos(iC) = expmupos*(xh(iC)**vsalt)
     &     *dexp(-psi(iC)*zpos) ! ion plus volume fraction

           xneg(iC) = expmuneg*(xh(iC)**vsalt)
     &     *dexp(-psi(iC)*zneg) ! ion neg volume fraction
     
           xHplus(iC) = expmuHplus*(xh(iC))
     &     *dexp(-psi(iC))           ! H+ volume fraction
     
           xOHmin(iC) = expmuOHmin*(xh(iC))
     &     *dexp(+psi(iC))           ! OH-  volume fraction



         if (prot_q < 0) then
               fdisP(iC) =
     &              1.0 /(1.0 + xHplus(iC)/(K0P*xh(iC)) )
         else ! base
               fdisP(iC) =
     &              1.0 /(1.0 + xOHmin(iC)/(K0P*xh(iC)) )
         endif

 
         do im =1,N_monomer
            if (zpol(im).eq.1) then !BASE
               fdis(im,iC) =
     &              1.0 /(1.0 + xOHmin(iC)/(K0(im)*xh(iC)) )
            else if (zpol(im).eq.-1) then !ACID
               fdis(im,iC) =
     &              1.0 /(1.0 + xHplus(iC)/(K0(im)*xh(iC)) )
            endif
         enddo

            enddo


! Calculo de xtotal para poor solvent

          do i = 1, N_poorsol
          xtotal(i, ncells+1) = 0.0
          xtotal(i, ncells+2) = 0.0
          enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculo de xpot

      do im =1, N_monomer
      do iC = 1, ntot

      select case (hguess)

      case (1) !wall

      hdistance = (dfloat(indexa(iC,1))*delta-hring)**2
     &  +(oval*(dfloat(indexa(iC,2))-0.5*dimZ)*delta)**2

      case (2) !center

      hdistance = (dfloat(indexa(iC,1))*delta)**2
     &  +(oval*(dfloat(indexa(iC,2))-0.5*dimZ)*delta)**2

      case (3) !wall-center

      hdistance = (min(dfloat(indexa(iC,1))*delta,
     & abs(dfloat(indexa(iC,1))*delta-hring)))**2
     &  +(oval*(dfloat(indexa(iC,2))-0.5*dimZ)*delta)**2

      endselect

      hfactor = dexp(-(kp**2)*hdistance)

      xpot(im,iC) = dlog(xh(iC))*vpol

      if(zpol(im).ne.0) then
      xpot(im,iC) = xpot(im,iC)
     & -psi(iC)*zpol(im)
     & -dlog(fdis(im, iC))
      endif

      if(hydroph(im).ne.0) then 

               xpot(im,iC) = xpot(im,iC)+henergy(hydroph(im))
     & *proteinhC(iC)*hst

               do ii = 1, N_poorsol ! loop over different poor sv types
               do jj = 1, nXu(ii, iC) ! loop over kai neighbors

                        xpot(im,iC) = xpot(im,iC) +
     &                       (st_matrix(hydroph(im),ii) ! st_matrix(x, y) : interaction of hydrophobic segments of type x with those of type y ( should be diagonal)
     &                       *hfactor*st/(vsol*vpol)*           ! st in kT/monomer          
     &                       Xulist_value(2, iC, jj)*
     &                       xtotal(ii, Xulist_cell(ii, iC, jj)))

               enddo ! jj
               enddo ! ii

               jp = 1
               if(hydroph(im).eq.jp) then

                 aveP(iC) = 0.0
                 avePnorm = 0.0

                 do jj = 1, nXu(jp, iC) ! loop over kai neighbors
                 avePnorm = avePnorm + Xulist_value(jp, iC, jj)
     &                      *(dfloat(indexa(iC,1))-0.5)
                 aveP(iC) = aveP(iC) +
     &                      (Xulist_value(jp, iC, jj)*
     &                      (dfloat(indexa(iC,1))-0.5)*
     &                      (xtotal(jp, Xulist_cell(jp, iC, jj))))
                 enddo
                 if(avepair.eq.1) then
                  aveP(iC) = aveP(iC)/avePnorm
                 else
                  aveP(iC) = xtotal(iC)  ! +decouple*xtotal(2,iC)
                 endif

                 if(aveP(iC).gt.0) then
                 apair = dexp(pairst*hfactor)-0.5
                 bpair=2*dexp(pairst*hfactor)-1+1/pairsize/aveP(iC)
                 cpair = dexp(pairst*hfactor)
                 Rpair(iC) =
     &                (bpair-sqrt(bpair**2-4*apair*cpair))/apair/2
                 Fpair(iC) = dlog(1-Rpair(iC))
     &               -dlog(1-pairsize*aveP(iC)*(1-0.5*Rpair(iC)))
                 xpot(im,iC)=xpot(im,iC)-Fpair(iC)
                 endif
               endif

      endif ! hydrophob
      enddo ! iC
      enddo ! im
      ii = rank+1

      avpol_temp = 0.0

      q = 0.0
      sumprolnpro = 0.0
      endtoend_av = 0.0

      do i=1,newcuantas(ii)

         lnpro=0.0


!! DECODER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         pC = fs(i)                                        !
         pR = indexa(pC,1)                                 !
         pZ = indexa(pC,2)                                 !
         binary(:) = displ(i,:)                            !
         call decode(displacement, binary, long(ii))                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do j=1,long(ii)

         pR = pR + displacement(j,1) 
         pZ = pZ + displacement(j,2) 
         pC = matriz(pR, pZ)

         lnpro = lnpro + xpot(segtype(ii, j),pC) !+ dlog(pbias(i))

         enddo ! j
            
            pro = dexp(lnpro)*pbias(ii,i)
           ! lnpro = dlog(pro)

            q=q+pro
            sumprolnpro = sumprolnpro + pro*lnpro

            endtoend_av = endtoend_av + pro*endtoend(i)

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         pC = fs(i)                                     !
         pR = indexa(pC,1)                              !
         pZ = indexa(pC,2)                              !
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do j=1,long(ii)

         pR = pR + displacement(j,1)
         pZ = pZ + displacement(j,2)
         pC = matriz(pR, pZ)

         avpol_temp(segtype(ii,j), ii, pC) 
     & = avpol_temp(segtype(ii,j), ii, pC)
     & + pro*chainsperdelta(ii)*vsol
     & *vpol*Factorcurv(pR)

         enddo ! j
         enddo !i

          avpol_temp=avpol_temp/q
          endtoend_av = endtoend_av/q
          avpol_tosend = avpol_temp

          avpol_temp = 0
        
c------------------ MPI ----------------------------------------------


         call MPI_REDUCE(avpol_tosend, avpol_temp
     & , ncells*N_monomer*N_chains,
     &   MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

      CALL MPI_BARRIER(MPI_COMM_WORLD,err)

      do ii = 1, N_chains
      do iC = 1,ntot
      do im = 1, N_monomer
      avpol(im,ii,iC) = avpol_temp(im,ii,iC)
      enddo ! im
      enddo ! iC
      enddo ! ii

      if(rank.ne.0)goto 3333
!!!!!!!!!!! IMPORTANT, SLAVES STOP HERE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C----------------------------------------------------------------------------------------------
C   Construct f vector
C----------------------------------------------------------------------------------------------

! Qtot

      do iC=1,ntot

         qtot(iC) = 
     &   (zpos*xpos(iC)+zneg*xneg(iC))/vsalt 
     &   + xHplus(iC)-xOHmin(iC)

         qtot_amp(iC) =
     &   (abs(zpos*xpos(iC))+abs(zneg*xneg(iC)))/vsalt
     &   + abs(xHplus(iC))+abs(xOHmin(iC))

         if(weakP.eq.1) then
         qtot(iC) = qtot(iC)  + proteinqC(iC)*vsol*fdisP(iC) 
         qtot_amp(iC)=qtot_amp(iC)+abs(proteinqC(iC))*vsol
     &   *fdisP(iC)
         else
         qtot(iC) = qtot(iC)  + proteinqC(iC)*vsol
         qtot_amp(iC)=qtot_amp(iC)+abs(proteinqC(iC))*vsol
         endif

         do im = 1, N_monomer
         do ii = 1, N_chains

         qtot(iC) = qtot(iC) + avpol(im, ii, iC)
     &   *fdis(im, iC)*zpol(im)/vpol
         qtot_amp(iC) = qtot_amp(iC) + avpol(im, ii, iC)
     &   *fdis(im, iC)*abs(zpol(im))/vpol

         enddo ! ii
         enddo ! im
 
      enddo



! Volume fraction

      do iC=1, ntot
         
               f(iC)=
     &   xh(iC) + xneg(iC) 
     &  + xpos(iC) + xHplus(iC) + xOHmin(iC)
     &  + proteinC(iC) - 1.0d0
   
      do im = 1, N_monomer
      do ii = 1, N_chains

      f(iC) = f(iC) + avpol(im, ii, iC)

      enddo ! ii
      enddo ! im

      enddo

! Poisson eq.

            do iC=1,ntot

! Cilindro (derivada al centro), ver notas     
     
               f(iC+ntot)=
     & psi(rp(iC)) -2*psi(iC) + psi(rm(iC)) +
     & (0.5/(dfloat(indexa(iC,1))-0.5))*(psi(rp(iC))-psi(rm(iC))) + ! termino de curvatura
     & psi(zp(iC)) -2*psi(iC) + psi(zm(iC)) + ! derivada en z
     & qtot(iC)*constq
     
               f(iC+ntot)=f(iC+ntot)/(-2.0) ! mejora kinsol...
      enddo

! poor solvent

      do ii = 1, N_poorsol
      do iC = 1,ntot

      f(iC+(1+ii)*ntot) = xtotal(ii,iC)

      do im = 1, N_monomer
      if(hydroph(im).eq.ii) then
      do jj = 1, N_chains
        f(iC+(1+ii)*ntot) = f(iC+(1+ii)*ntot) - avpol(im,jj,iC)
      enddo ! jj
      endif
      enddo ! im
      enddo ! iC
      enddo ! ii

      iter = iter + 1

      norma = 0.0

      do i = 1, (2+N_poorsol)*ntot
      norma = norma +(f(i))**2    
      enddo
      
      print*, iter, norma, q
      open(unit=6660, file='iter.dat', access='append')
      write(6660,*)iter, norma, q
      close (6660)
      ier = 0.0

! saves infile

      if(mod(iter, save_every).eq.0) then
      print*, 'Save resume file'
      open(unit=45, file = 'out.temp.dat')
      do i = 1, (2+N_poorsol)*ncells
      write(45, *)x(i)
      enddo
      close(45)
      endif

 3333 ier = 0.0

      return
      end
