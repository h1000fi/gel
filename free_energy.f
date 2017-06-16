      subroutine calc_free_energy(counter, counter2)
      
      use mparameters                                             
      use mparameters_chain
      use mparameters_monomer
      use mvariables
      use mlookup
      use mgraft
      use mncells
      use mkai
      use mprotein 

      implicit none                                 
      include 'mpif.h'                              
      include 'MPI.h'                               
                                                    
      double precision Factorcurv                   
                                                    
      real*8 Free_energy, F_Mix_s, F_Mix_pos        
      real*8 F_Mix_neg, F_Mix_Hplus                 
      real*8 Free_energy2, sumpi, sumrho, sumel, sum, mupol, pilat
      real*8 sumfcargo
      real*8 F_Mix_OHmin, F_Conf, F_Conf_temp                     
      real*8  F_Conf2, F_Conf_temp2                               
     & , F_Eq, F_Eq_P, F_vdW, F_eps, F_electro, F_pair                 
                                                                  
      real*8 counter, counter2                                    
                                                                  
      integer iC, ii, iiC, i, jj, im, iii                              
      real*8 avpol_monom(N_monomer, dimR*dimZ)    
      real*8 avpol_hydroph(N_poorsol, dimR*dimZ)    
      real*8 q0(N_chains), sumprolnpro0(N_chains)
      real*8 q0_tosend(N_chains), sumprolnpro0_tosend(N_chains)

      real*8 proteinN(ncells)

      real*8 endtoend_av_all(N_chains)
      real*8 endtoend_av_tosend(N_chains)

! pass chain statistic to master process and quit 

       if(rank.eq.0)print*, 'Starting free energy calculation...'

! proteinN has the number density of protein charged groups

       do iC = 1, ncells
       proteinN(iC) = abs(proteinqC(iC)) 
       enddo

! avpol preliminary calculations

      avpol_monom = 0.0
      do im = 1, N_monomer
      do ii = 1, N_chains
      avpol_monom(im, :) = avpol_monom(im, :) +  avpol(im,ii,:)
      enddo
      enddo


c      avpol_hydroph(:, :) = 0.0
c
c      do j = 1, N_poorsol
c      do im = 1, N_monomer
c
c      if(hydroph(im).eq.j) then
c         do ii = 1, N_chains
c         avpol_hydroph(j, :) = avpol_hydroph(j, :)+ avpol(j, ii, :)
c         enddo
c      endif
c
c      enddo ! im
c      enddo ! j

! open files

       if(pHbulk.eq.pHstart) then
       open(unit=301, file='Free_energy.dat')                       
         write(301,*) 'pHbulk ','Free_energy ','Free_energy2 ',
     &  'F_Mix_s ','F_Mix_pos ','F_Mix_neg ','F_Mix_Hplus ',
     &  'F_Mix_OHmin ','F_Conf ','F_Eq ','F_Eq_P ','F_vdW ','F_eps ',
     &  'F_pair ','F_electro '
       close(301)
       endif                         
                                                                  
ccc---------------------------------------------------------------------
cc        Calculate Free Energy                                        
ccC-------------------------------------------------------------------- 
                                                                       
      Free_Energy = 0.0                                                
      Free_Energy2 = 0.0                                               
                                                                       
c! 1. Solvent translational entropy
                                                                       
      F_Mix_s = 0.0                                                    
                                                                       
      do iC = 1, ncells                                                
      F_Mix_s = F_Mix_s + xsol(iC)*(dlog(xsol(iC))-1.0)*               
     & (dfloat(indexa(iC,1))-0.5)                                      
      F_Mix_s = F_Mix_s - xsolbulk*(dlog(xsolbulk)-1.0)*               
     & (dfloat(indexa(iC,1))-0.5)                                      
      enddo                                                            
                                                                       
      F_Mix_s = F_Mix_s * delta**3/vsol*2*pi                           
      Free_Energy = Free_Energy + F_Mix_s                              
                                                                       
                                                                       
c! 2. Pos ion translational entropy
                                     
                                                                       
      F_Mix_pos = 0.0                                                  
                                                                       
      do iC = 1, ncells                                                
                                                                       
      F_Mix_pos = F_Mix_pos + xpos(iC)*(dlog(xpos(iC)/vsalt)-1.0       
     &  - dlog(expmupos) + dlog(vsalt))*                               
     & (dfloat(indexa(iC,1))-0.5)                                      
                                                                       
      F_Mix_pos = F_Mix_pos - xposbulk*(dlog(xposbulk/vsalt)-1.0       
     &  - dlog(expmupos) + dlog(vsalt))*                               
     & (dfloat(indexa(iC,1))-0.5)                                      
                                                                       
      enddo                                                            
      F_Mix_pos = F_Mix_pos * delta**3/vsol/vsalt*2*pi                 
      Free_Energy = Free_Energy + F_Mix_pos                            
                                                                       
c!3. Neg ion translational entropy
                                                                       
      F_MIX_neg = 0.0                                                  
                                                                       
      do iC = 1, ncells                                                
      F_Mix_neg = F_Mix_neg + xneg(iC)*(dlog(xneg(iC)/vsalt)-1.0       
     & - dlog(expmuneg) + dlog(vsalt))*                                
     & (dfloat(indexa(iC,1))-0.5)                                      
                                                                       
      F_Mix_neg = F_Mix_neg - xnegbulk*(dlog(xnegbulk/vsalt)-1.0       
     & - dlog(expmuneg) + dlog(vsalt))*                                
     & (dfloat(indexa(iC,1))-0.5)                                      
                                                                       
      enddo                                                            
      F_Mix_neg = F_Mix_neg * delta**3/vsol/vsalt*2*pi                 
      Free_Energy = Free_Energy + F_Mix_neg                            
                                                                       
                                                                       
c! 4. H+ translational entropy
                                                                       
      F_Mix_Hplus = 0.0                                                
                                                                       
      do iC = 1, ncells                                                
      F_Mix_Hplus = F_Mix_Hplus + xHplus(iC)*(dlog(xHplus(iC))-1.0     
     & -dlog(expmuHplus))*                                             
     & (dfloat(indexa(iC,1))-0.5)                                      
                                                                       
      F_Mix_Hplus = F_Mix_Hplus - xHplusbulk*(dlog(xHplusbulk)-1.0     
     & -dlog(expmuHplus))*                                             
     & (dfloat(indexa(iC,1))-0.5)                                      
                                                                       
      enddo                                                            
      F_Mix_Hplus = F_Mix_Hplus * delta**3/vsol*2*pi                   
      Free_Energy = Free_Energy + F_Mix_Hplus                          
                                                                       
                                                                       
c! 5. OH- translational entropy 
                                                                       
      F_Mix_OHmin = 0.0                                                
                                                                       
      do iC = 1, ncells                                                
      F_Mix_OHmin = F_Mix_OHmin + xOHmin(iC)*(dlog(xOHmin(iC))-1.0     
     & -dlog(expmuOHmin))*                                             
     & (dfloat(indexa(iC,1))-0.5)                                      
                                                                       
      F_Mix_OHmin = F_Mix_OHmin - xOHminbulk*(dlog(xOHminbulk)-1.0     
     & -dlog(expmuOHmin))*                                             
     & (dfloat(indexa(iC,1))-0.5)                                      
                                                                       
      enddo                                                            
      F_Mix_OHmin = F_Mix_OHmin * delta**3/vsol*2*pi                   
      Free_Energy = Free_Energy + F_Mix_OHmin                          
                                                                       
                                                                       
c! 6. Polymer conformational entropy                                         
      
      F_conf = 0.0
      q0 = 0.0
      q0_tosend = 0.0
      sumprolnpro0_tosend = 0.0
     
      ii = rank+1
 
      q0_tosend(ii) = q
      sumprolnpro0_tosend(ii) = sumprolnpro 
                                                            
         call MPI_REDUCE(q0_tosend, q0
     & , N_chains,
     &   MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

         call MPI_REDUCE(sumprolnpro0_tosend, sumprolnpro0
     & , N_chains,
     &   MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

      do ii = 1, N_chains
      F_Conf = F_Conf + (sumprolnpro0(ii)/q0(ii) - dlog(q0(ii)))
     &  * chainsperdelta(ii)
      enddo                                                           

      Free_Energy = Free_Energy + F_Conf
                                                                       
c! 7. Chemical Equilibrium                                              

      F_Eq = 0.0                                                       

      do im = 1, N_monomer                                              
      do iC  = 1, ncells                                               
                         
      if(zpol(im).ne.0) then
                                              
      F_Eq = F_Eq + fdis(im, iC)*dlog(fdis(im, iC))
     & *avpol_monom(im, iC)/vpol     
     & *(dfloat(indexa(iC,1))-0.5)                                     
                                                                       
      F_Eq = F_Eq + (1.0-fdis(im, iC))                                 
     & *dlog(1.0-fdis(im, iC))*avpol_monom(im, iC)/vpol                
     & *(dfloat(indexa(iC,1))-0.5)                                     
                                                                       
      F_Eq = F_Eq + (1.0-fdis(im, iC))*
     & dlog(K0(im))*avpol_monom(im, iC)/vpol     
     & *(dfloat(indexa(iC,1))-0.5)                                     
            

         select case (zpol(im))
         case (-1) ! acid
      F_Eq = F_Eq + (1.0-fdis(im, iC))                                 
     & *(-dlog(expmuHplus))*avpol_monom(im, iC)/vpol                   
     & *(dfloat(indexa(iC,1))-0.5)
         case (1) ! base
      F_Eq = F_Eq + (1.0-fdis(im, iC))                                 
     & *(-dlog(expmuOHmin))*avpol_monom(im, iC)/vpol                   
     & *(dfloat(indexa(iC,1))-0.5)
         end select
              
      endif
      enddo                                                            
      enddo                                                            
                                                                       
      F_eq = F_eq *delta**3/vsol*2*pi                                  
                                                                       
      Free_Energy = Free_Energy + F_Eq                                 
                                                                       
c! 8.vdW ! Ojo, los  son negativos => atraccion                         

      F_VdW = 0.0

      do ii = 1, N_poorsol
      do iii = 1, N_poorsol
        do iC = 1, ncells
               do iiC = 1, nXu(ii, iC) ! loop over kai neighbors     

             F_vdW = F_vdW - 0.5000*delta**3*xtotal2(ii,iC)*
     &       xtotal2(iii,Xulist_cell(ii, iC, iiC))*
     &       Xulist_value(ii, iC, iiC)*st_matrix(ii, iii)*st
     &       /((vpol*vsol)**2)
     &       *(dfloat(indexa(iC,1))-0.5)*2*pi

               enddo  ! iiC                                            
        enddo ! iC          
      enddo ! iii                                             
      enddo ! ii

      Free_Energy = Free_Energy + F_vdW                            
                                                                       
c! 9. Electrostati -- no charge on surfaces                            
                                                                       
      F_electro = 0.0                                                  
                                                                       
      do iC  = 1, ncells                                               
                                                                       
      F_electro = F_electro + delta**3*psi2(iC)*qtot(iC)/2.0/vsol      
     &               *(dfloat(indexa(iC,1))-0.5)*2*pi                  
                                                                       
      enddo                                                            
                                                                       
      Free_Energy = Free_Energy + F_electro                            

c! 10. Protein-sup

      F_eps = 0.0

      do iC = 1, ncells
      do ii = 1, N_poorsol
      F_eps = F_eps - xtotal2(ii,iC)*henergy(ii)*proteinhC(iC)*
     & hst*(dfloat(indexa(iC,1))-0.5)*2.0*pi
     & *(delta**3)/vpol/vsol
      enddo ! ii
      enddo ! iC

      Free_Energy = Free_Energy + F_eps

c! 11. Protein-sup

      F_pair = 0.0

      do iC = 1, ncells
      F_pair = F_pair - xtotal2(1,iC)*Fpair(iC)*
     & (dfloat(indexa(iC,1))-0.5)*2.0*pi
     & *(delta**3)/(vpol*vsol)
      enddo ! iC

      Free_Energy = Free_Energy + F_pair
!

c! 12. Chemical Equilibrium Particle                                             

      F_Eq_P = 0.0                                                       

      if(weakP.eq.1) then

      do iC  = 1, ncells                                               
                         
      F_Eq_P = F_Eq_P + fdisP(iC)*dlog(fdisP(iC))
     & *proteinN(iC)*vsol     
     & *(dfloat(indexa(iC,1))-0.5)                                     
                                                                       
      F_Eq_P = F_Eq_P + (1.0-fdisP(iC))                                     
     & *dlog(1.0-fdisP(iC))*proteinN(iC)*vsol                      
     & *(dfloat(indexa(iC,1))-0.5)                                     
                                                                       
      F_Eq_P = F_Eq_P + (1.0-fdisP(iC))*
     & dlog(K0P)*proteinN(iC)*vsol
     & *(dfloat(indexa(iC,1))-0.5)                                     
            

         if (prot_q < 0) then
      F_Eq_P = F_Eq_P + (1.0-fdisP(iC))                                     
     & *(-dlog(expmuHplus))*proteinN(iC)*vsol                    
     & *(dfloat(indexa(iC,1))-0.5)
         else
      F_Eq_P = F_Eq_P + (1.0-fdisP(iC))                                     
     & *(-dlog(expmuOHmin))*proteinN(iC)*vsol                
     & *(dfloat(indexa(iC,1))-0.5)
         endif
              
      enddo                                                            
      endif
                                                                       
      F_eq_P = F_eq_P *delta**3/vsol*2*pi                                  
                                                                       
      Free_Energy = Free_Energy + F_Eq_P                                 


      if(rank.eq.0)print*, Free_energy, F_electro
      if(rank.eq.0)print*, 'Free Energy, method I: ', Free_Energy
 
c! Method II                                                          
                                                                       
      Free_Energy2 = 0.0                                               
                                                                       
         sumpi = 0.0                                                   
         sumrho=0.0                                                    
         sumel=0.0                                                     
         sumfcargo = 0.0                                                                       

         do iC=1,ncells                                                
                                                                       
            sumpi =                                                    
     & sumpi+(1.0-proteinC(iC))*dlog(xsol(iC))                         
     & *(dfloat(indexa(iC,1))-0.5)*2*pi                                
 
c            sumpi =
c     & sumpi+(1.0)*dlog(xsol(iC))
c     & *(dfloat(indexa(iC,1))-0.5)*2*pi
                                                                      
            sumpi = sumpi-dlog(xsolbulk)*(dfloat(indexa(iC,1))-0.5)*2*pi
            sumrho = sumrho + ( - xsol(iC) -xHplus(iC) -xOHmin(iC)      
     &               -(xpos(iC)+xneg(iC))/vsalt)! sum over  rho_i i=+,-,si
     &               *(dfloat(indexa(iC,1))-0.5)*2*pi                     
                                                                          
            sumrho = sumrho - ( - xsolbulk -xHplusbulk -xOHminbulk        
     &               -(xposbulk+xnegbulk)/vsalt)! sum over  rho_i i=+,-,si
     &               *(dfloat(indexa(iC,1))-0.5)*2*pi                     
                                

            sumel = sumel - qtot(iC)*psi2(iC)/2.0
     &               *(dfloat(indexa(iC,1))-0.5)*2*pi

            sumel = sumel + proteinqC(iC)*psi2(iC)*vsol                   
     &               *(dfloat(indexa(iC,1))-0.5)*2*pi   

      if(weakP.eq.1) then

            sumfcargo = sumfcargo + dlog(fdisP(iC))
     & *proteinN(iC)*vsol
     & *(dfloat(indexa(iC,1))-0.5)*2*pi

      endif

         enddo                                                      
                                                                    
         sumpi = (delta**3/vsol)*sumpi                              
         sumrho = (delta**3/vsol)*sumrho                            
         sumel = (delta**3/vsol)*sumel                               
         sumfcargo = (delta**3/vsol)*sumfcargo                        
 
         Free_Energy2 = sumpi + sumrho + sumel + sumfcargo         
         Free_Energy2 = Free_Energy2 - F_vdW
                          
         do ii = 1, N_chains                                         
         Free_Energy2 = Free_Energy2 -
     & chainsperdelta(ii)*dlog(q0(ii)/shift)
         enddo
                                                          
      if(rank.eq.0)print*, 'Free Energy, method II: ', Free_Energy2
         
         if(rank.eq.0) then                                            
         open(301,file='Free_energy.dat',access='append',
     & status='old',action='write')
         write(301,*) pHbulk,Free_energy,Free_energy2,F_Mix_s,F_Mix_pos,
     &   F_Mix_neg,F_Mix_Hplus,F_Mix_OHmin,F_Conf,F_Eq,F_Eq_P,F_vdW,
     &   F_eps,F_pair,F_electro           
         close(301)                      
         endif                           

c! Save end-to-end distances         

      ii = rank+1
      endtoend_av_tosend = 0.0
      endtoend_av_tosend(ii) = endtoend_av

         call MPI_REDUCE(endtoend_av_tosend, endtoend_av_all
     & , N_chains,
     &   MPI_DOUBLE_PRECISION, MPI_SUM,0, MPI_COMM_WORLD, err)

       if(rank.eq.0) then
       open(unit=502, file='endtoend.dat')                          
       do ii = 1, N_chains
       write(502,*)ii, endtoend_av_all(ii)
       enddo                      
       endif
                                      
 1515    continue
                                                                          
         end                                                              
                                                                          
                                                                          

