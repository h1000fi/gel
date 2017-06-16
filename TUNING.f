! ================================================================== MK SPECIAL REQUIRES s and p, finds q
        subroutine TUNING(N,rpos,zpos)

        use mlookup
        use mparameters
        use mporesystem
        use mrands
        use mparameters_chain
        use mparameters_monomer

        implicit none

        integer N,ncha
        real*8 tic,q_best
        integer cuantas_done,check_deaths,check_rounds
        real*8 cputotal,check_DAR,check_ok
        integer check_CPT
        real*8 rpos,zpos
        integer j
        integer chains(ncha_max,maxlong)
        integer*4 deathtotal
        common /crit/ deathtotal
        integer spacer,branch_p
        common /com1/ spacer,branch_p
        common /qprob/ qprob

        real*8 qprob,deltaq
        logical converged
        print '(A,I6,A)','+++ TUNING wanted ',wantedCPUsecs,' secs'
        ! print *,'+++ mk-tuning cpu ',wantedCPUsecs/dble(1000)
c        print 101,'N','s','p','q','DAR','CPT','cpu','cpu-total','ok%'
        qprob  = 1.0            ! DO NOT EDIT
        deltaq = 0.1            ! DO NOT EDIT
        converged = .false.
        ncha   = 1
        do while (.not.converged)
         tic = secnds(0.)
         cuantas_done = 0
         check_deaths = 0
         check_rounds = 0
         do while ((secnds(0.).lt.tic+wantedCPUsecs/1000).      ! DO NOT EDIT
     >             or.(check_rounds.lt.40))                     ! DO NOT EDIT
          call cadenas_mk(chains,ncha,rpos,zpos,N)
          check_rounds = check_rounds + 1
          cuantas_done = cuantas_done + ncha
          check_deaths = check_deaths + deathtotal
         enddo
         tic = secnds(0.)-tic
         cputotal = tic*cuantas/dble(1e-6+cuantas_done)
         tic = 10**6*tic/dble(max(1,cuantas_done))
         check_DAR = check_deaths/dble(max(1,cuantas_done))
         check_CPT = cuantas_done/dble(check_rounds)
         check_ok  = 100*abs(cputotal-wantedCPUsecs)/dble(wantedCPUsecs)
c         print 100,N,spacer,branch_p,qprob,check_DAR,
c     >         check_CPT,tic,cputotal,check_ok
         if (int(cputotal).gt.wantedCPUsecs) goto 1
         if (cuantas_done.eq.0) goto 1
         qprob = qprob - deltaq
         goto 2
1        continue
         qprob = min(1.D0,qprob + deltaq)
         deltaq = deltaq/1.5
         if (deltaq.lt.1e-4) goto 3
2        continue
        enddo
3       continue
        print *,'now operating at N, q ',N,qprob
        

c100     format("+++ ",3(I6,1x),F8.6,1x,F8.3,1x,I8,1x,F8.1,2(1x,F8.1))
c101     format("+++ ",3(A6,1x),A8,1x,A8,1x,A8,1x,A8,1x,A8,1x,A8)

        return
        end

   
C****************************************************************
C **********************************************************************
        double precision FUNCTION RANDS (SEED)
C **********************************************************************

C-----  this is a special function for random number generation
C        on 32-bit machines that do not support long integer
C        multiplication and truncation.  the technique used is to do
C        the multiplication and addition in parts, by splitting all
C       integers in a 'high' and a 'low' part.  the algorithm is
C-----        exact, and should give machine-independent results.

C-----        the algorithm implemented is (following d.e. knuth):
C        seed = seed*1592653589 + 453816691
C        if (seed.lt.0) seed = seed + 1 + 2147483647
C-----        note that 1592653589 = 48603*2**15 + 30485

C 32768 = 2^15, 65536 = 2^16, 4.65661287308E-10 = 2^(-31)

        INTEGER SEED, I1, I2

        I2 = MOD (SEED, 32768) * 30485
        I1 = MOD (SEED / 32768 * 30485, 65536) + MOD (MOD (SEED, 32768)
     X    * 48603, 65536) + I2 / 32768 + MOD (SEED / 32768, 2) * 32768 +
     X     13849
        I2 = MOD (I2, 32768) + 12659
        I1 = I1 + I2 / 32768
        I2 = MOD (I2, 32768)
        I1 = MOD (I1, 65536)
        SEED = I1 * 32768 + I2
        RANDS = REAL(I1 * 32768 + I2) * 4.65661287308E-10

        RETURN
        END

