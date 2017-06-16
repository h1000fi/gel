!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Subrutina que llama a kinsol
!

      subroutine call_kinsol(x1, xg1, ier)

      use mncells
      use mparameters
      use mparameters_chain
      use mparameters_monomer

      implicit none

      real*8 x1((2+N_poorsol)*ncells), xg1((2+N_poorsol)*ncells)

      integer*8 iout(15) ! Kinsol additional output information
      real*8 rout(2) ! Kinsol additional out information

      integer*8 msbpre
      real*8 fnormtol, scsteptol
      real*8 scale((2+N_poorsol)*ncells)
      real*8 constr((2+N_poorsol)*ncells)

      integer*4  globalstrat, maxl, maxlrst, max_niter

      integer *4 ier ! Kinsol error flag
      integer *8 neq ! Kinsol number of equations

      integer i, ierr

      common /psize/ neq ! Kinsol


! INICIA KINSOL


      neq = (2+N_poorsol)*ncells

     

      msbpre  = 20 !10 ! maximum number of iterations without prec. setup (?)
      fnormtol = tol !1.0d-9 ! Function-norm stopping tolerance
      scsteptol = tol !1.0d-9 ! Function-norm stopping tolerance

      maxl = 500 ! maximum Krylov subspace dimesion (?!?!?!) ! Esto se usa para el preconditioner
      maxlrst = 100 !20 ! maximum number of restarts

      max_niter = 2000

      globalstrat = 0

      call fnvinits(3, neq, ier) ! fnvinits inits NVECTOR module
      if (ier .ne. 0) then       ! 3 for Kinsol, neq ecuantion number, ier error flag (0 is OK)
      print*, 'SUNDIALS_ERROR: FNVINITS returned IER = ', ier
      call MPI_FINALIZE(ierr) ! finaliza MPI
         stop

      endif

      call fkinmalloc(iout, rout, ier)    ! Allocates memory and output additional information
      if (ier .ne. 0) then
         print*, 'SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
         call MPI_FINALIZE(ierr) ! finaliza MPI
         stop
      endif


      call fkinsetiin('MAX_SETUPS', msbpre, ier)  ! Additional input information
      call fkinsetrin('FNORM_TOL', fnormtol, ier)
      call fkinsetrin('SSTEP_TOL', scsteptol, ier)
      call fkinsetiin('MAX_NITER', max_niter, ier)
                                                      
      do i = 1, ncells  !constraint vector
         constr(i) = 2.0
      enddo
      do i = ncells+1, 2*ncells  !constraint vector
         constr(i) = 0.0
      enddo
      do i = 2*ncells+1, neq
         constr(i) = 0.0
      enddo 

      call fkinsetvin('CONSTR_VEC', constr, ier) ! constraint vector

c       CALL FKINSPTFQMR (MAXL, IER)
      call fkinspgmr(maxl, maxlrst, ier) !  Scale Preconditioned GMRES solution of linear system (???)
c       call fkinspbcg(maxl, ier) !  Scale Preconditioned BCG


      if (ier .ne. 0) then
      print*, 'SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
         call fkinfree ! libera memoria
         call MPI_FINALIZE(ierr) ! finaliza MPI
         stop
      endif
      call fkinspilssetprec(1, ier) ! preconditiones

      do i = 1, neq ! scaling vector
       scale(i) = 1.0
      enddo

      do i = 1, neq ! Initial guess
      x1(i) = xg1(i)
      enddo


      call fkinsol(x1, globalstrat, scale, scale, ier)         ! Llama a kinsol
      if (ier .lt. 0) then
           print*, 'SUNDIALS_ERROR: FKINSOL returned IER = ', ier
           print*, 'Linear Solver returned IER = ', iout(9)
           call fkinfree
      endif
      call fkinfree

      return
      end
