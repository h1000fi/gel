      subroutine creador ! crea cadenaes

      use mlookup
      use mchains
      use mparameters
      use mncells
      use mgraft
      use mparameters_chain
      use mparameters_monomer
      use mrands
      use mprotein

      implicit none
      include 'mpif.h' ! MPI libraries
      include 'MPI.h'
 
      integer k,vx(4),vy(4)

      integer total,ix(3)
      integer *1 displ_temp
      integer i,il,u1,u2,iii,ii,ll, jj
      integer j,ncha

      real*8 Rb, Zb, tip2bias, accept, rateBias

      real*8 indax, inday, indaz
      
      real*8 altx,alty, altz
      real*8 rij,theta,theta1, rn1, rn2

      real*8 rpos, zpos
           
      integer total1,iglobal

      integer spp ! sitios por procesador

      integer cuantas_p(2)

      real*8 q_tosend
      integer*1 displ_one(0:7)

      integer*1 displacement(maxlong+1,2)
      integer*1 binary(int(maxlong/2))
     
      real*8 xend(3,maxlong)
      integer chains(ncha_max, maxlong)
      real*8 endtoendtemp(10000)
      real*8 x(maxlong)
      real*8 y(maxlong)
      real*8 z(maxlong)

      real*8 in1(maxlong, 3)
      common /endtoend/ endtoendtemp

      seed = readseed
 
      rpos = posgraft(1)
      zpos = posgraft(2)

      ii = rank+1
    
      newcuantas(ii) = 0

      il=0
      iglobal=1
 
      if(calq.eq.1)ncha=-1

      do while (il.lt.cuantas)
 
         if(cadenastype.eq.1)
     & call cadenas72mr(chains,ncha, rpos, zpos, long(ii))
         if(cadenastype.eq.2)
     & call cadenas_mk(chains,ncha, rpos, zpos, long(ii))
      
         do i=1,ncha

!         il = il + 1
 
         if(il.ge.cuantas) goto 100

! Check collision with protein...

         do j = 1, long(ii)
         if(proteinC(chains(i,j)).ne.0.0)goto 200 ! collides with protein, don't use this chain
         enddo

         Rb = (dble(indexa(chains(i,long(ii)),1))-0.5)*delta
         Zb = (dble(indexa(chains(i,long(ii)),2))-0.5)*delta
         tip2bias = sqrt((Zb-zpos-Zbias)**2)
!         tip2bias = sqrt((Rb-Rbias)**2+(Zb-zpos-Zbias)**2)
         accept = dexp(-kBias*tip2bias)
         rateBias=rands(seed)
         if(rateBias .gt. accept) goto 200

         newcuantas(ii) = newcuantas(ii)+1
         il = il + 1

         endtoend(newcuantas(ii)) = endtoendtemp(i)

          fs(newcuantas(ii))=chains(i,1)
          pbias(ii,newcuantas(ii))=1.0
          pbias(ii,newcuantas(ii))=1.0/accept

          if(fs(newcuantas(ii)).eq.0) then 
          print*, 'Error', il, i, rank, ncha
          stop
          endif


          displacement(1,1) = 0
          displacement(1,2) = 0

          do j=2,long(ii)        

          displacement(j,1)= 
     &  (indexa(chains(i,j),1)-indexa(chains(i,j-1),1)) ! displacement in R
          
          displacement(j, 2) =
     &  (indexa(chains(i,j),2)-indexa(chains(i,j-1),2)) ! displacement in Z

          enddo ! 

          call encode(displacement,binary, long(ii))
          displ(newcuantas(ii),:) = binary(:)

 200      enddo !ncha
          enddo ! il 

      print*, 'Processor ', rank+1, 'has',newcuantas(ii),'conformations'
      if(newcuantas(ii).eq.0)stop
          
 100  return
      end

