      subroutine presolver

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

      integer im

C--------------------------------------------

      zpos = 1.0
      zneg = -1.0
      
      vsol = 0.030     
      vsalt=((4.0/3.0)*pi*(0.27)**3)/vsol  ! volume salt in units of vsol 0.2=radius salt  
      vpol= 0.095/vsol! ((4.0/3.0)*pi*(0.2)**3)/vsol  ! volume polymer segment in units of vsol 

      constq=delta*delta*4.0*pi*lb/vsol   ! multiplicative factor in poisson eq  

      pKw = 14
      Kw = 1.0e-14

      betae = 38.94             ! beta * e


c-------------------------------------------------------
c
c  CASE DEPENDENT VARIABLES
c

      cHplus = 10**(-pHbulk)    ! concentration H+ in bulk
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

      end
