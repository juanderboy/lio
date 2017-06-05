!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE TD()
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
c REAL TIME-TDDFT
c 
c Dario Estrin, 1992
c Nano, Dario, Uriel, Damian 2012
c
c  This subrutine takes the converged density matrix from an SCF calculation
c  and evolves it in time. In the input file the total number of propagation
c  steps is specified (nstep) as well as the time of each evolution step 
c  (tdstep). 
c  This implementation has two alternatives to evolve the density in time. The 
c  first one (propagator=1) is the Verlet algorithm that uses a convination of 
c  Liouville von Newmann expresion for the time derivative of the density matrix 
c  and a first order Taylor expansion of the density matrix. The second one 
c  (propagator=2) is the Magnus propagation scheme that uses Backer Campbell
c  Hausdorff (BCH) formula. For this reason when Magnus is used the number of 
c  total conmutators in the BCH espansion has to be specified (NBCH, default=10). 
c  A narrow gaussian type electric field can be introduced during the time 
c  evolution in order to excite all electronic frequencies with the same intensity.
c  Once this perturbation is turned on (Field=t, exter=t) each component of the
c  external electric field has to be specified in the input file (Fx,Fy,Fz).
c  In each step of the propagation the cartesian components of the sistems dipole
c  are stored in files x.dip, y.dip, z.dip.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
c       USE latom
       USE garcha_mod
       use mathsubs
#ifdef CUBLAS
       use cublasmath
#endif
       use general_module
       IMPLICIT REAL*8 (a-h,o-z)
       INTEGER :: istep
       REAL*8 :: t,E2
       REAL*8,ALLOCATABLE,DIMENSION(:,:) :: 
     >   xnano2,xmm,xtrans,ytrans,Y,fock,
     >   F1a,F1b,overlap,rhoscratch
       real*8, dimension (:,:), ALLOCATABLE :: elmu
#ifdef TD_SIMPLE
       COMPLEX*8 :: Im,Ix
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:) ::
     >   rho,rhonew,rhold,xnano,rho1
#else
       COMPLEX*16 :: Im,Ix
       COMPLEX*16,ALLOCATABLE,DIMENSION(:,:) ::
     >   rho,rhonew,rhold,xnano,rho1
#endif
       DIMENSION q(natom)
       REAL*8,dimension(:),ALLOCATABLE :: factorial
       INTEGER            :: LWORK,ii,jj
       REAL*8,ALLOCATABLE :: WORK(:)
!!------------------------------------!!
!! FFR ADD
       INTEGER ::
     >   pert_steps,lpfrg_steps,chkpntF1a,chkpntF1b
       REAL*8 ::
     >   dt_magnus,dt_lpfrg
        logical :: just_int3n,ematalloct,lpop
!! CUBLAS
#ifdef CUBLAS
      integer sizeof_real
      parameter(sizeof_real=8)
      integer sizeof_complex
#ifdef TD_SIMPLE
      parameter(sizeof_complex=8)
#else
      parameter(sizeof_complex=16)
#endif
      integer stat
      integer*8 devPtrX, devPtrY,devPtrXc,devPtrS
      external CUBLAS_INIT, CUBLAS_SET_MATRIX,CUBLAS_FREE
      external CUBLAS_SHUTDOWN, CUBLAS_ALLOC,CUBLAS_GET_MATRIX
      integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX,CUBLAS_GET_MATRIX
#endif
!!   GROUP OF CHARGES
       LOGICAL             :: groupcharge
       INTEGER             :: ngroup
       INTEGER,ALLOCATABLE :: group(:)
       REAL*8,ALLOCATABLE  :: qgr(:)
!! TRANSPORT - ELECTROSTAT -
       COMPLEX*8,ALLOCATABLE,DIMENSION(:,:) :: rho_et
       LOGICAL :: ET
       LOGICAL :: GATEFIELD
       REAL*8 :: scratchgamma
       REAL*8 :: fxx, fyy, fzz
!       LOGICAL :: TRANSPORT_CALC
#ifdef TD_SIMPLE
       COMPLEX*8, allocatable :: rhofirst(:,:)
#else
       COMPLEX*16, allocatable :: rhofirst(:,:)
#endif
       INTEGER,ALLOCATABLE,DIMENSION(:,:) :: mapmat
!FFR!
       logical             :: dovv
       real*8              :: weight
       integer,allocatable :: atom_group(:)
       integer,allocatable :: orb_group(:)
       integer,allocatable :: orb_selection(:)

       real*8,dimension(:,:),allocatable :: fockbias
       real*8,dimension(:,:),allocatable :: sqsm
       real*8,dimension(:,:),allocatable :: Vmat
       real*8,dimension(:),allocatable   :: Dvec
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       call g2g_timer_start('TD')
       call g2g_timer_start('inicio')
       just_int3n = .false.
       ALLOCATE(factorial(NBCH))
!!------------------------------------!!
! Mulliken
       ipop=1
! Group of charges
       groupcharge=.true.
!!------------------------------------!!
#ifdef CUBLAS
       write(*,*) 'USING CUBLAS'
       stat=CUBLAS_INIT()
       if (stat.NE.0) then
           write(*,*) "initialization failed -TD"
           call CUBLAS_SHUTDOWN
           stop
       endif
#endif
#ifdef TD_SIMPLE
        write(*,*) 'simple presition complex'
#else
        write(*,*) 'double presition complex'
#endif
       if(propagator.eq.2) then
          dt_magnus=tdstep
          dt_lpfrg=tdstep*0.10D0
          factorial(1)=1.0D0
#ifdef CUBLAS
          DO ii=1,NBCH
             factorial(ii)=1.0D0/ii
          ENDDO
#else     
       DO ii=2,NBCH
         factorial(ii)=factorial(ii-1)/ii
       ENDDO
#endif
       endif
       if(propagator.eq.1) then
          dt_lpfrg=tdstep
       endif
!!------------------------------------!!
!! FFR ADD:
       pert_steps=100
       lpfrg_steps=200
       chkpntF1a=185
       chkpntF1b=195
! CASO RHOFIRST
       ALLOCATE(rhofirst(M,M))
! LOWDIN POPULATION
       lpop=.false.
! TRANSPORT
       GATEFIELD=.true.
       IF (TRANSPORT_CALC) THEN
       ALLOCATE(mapmat(M,M))
       ENDIF
       if (groupcharge) then
           inquire(file='atomgroup',exist=exists)
           if (.not.exists) then
               write(*,*) 'ERROR CANNOT FIND atomgroup file'
               stop
           endif
           open(unit=678,file='atomgroup')
           allocate(group(natom))
           ngroup=0
           do n=1,natom
               read(678,*) kk
               group(n)=kk
               if (kk.gt.ngroup) ngroup=kk
           enddo
           allocate(qgr(ngroup))
           if(ngroup.gt.3) write(*,*) 'if the number of group
     > is greater than 3 then group 1 should be the donor and 2 the 
     > acceptor'
           close(unit=678)
       endif
       IF(TRANSPORT_CALC) THEN
          call mat_map(group,mapmat)
       ENDIF
! Pointers -
       Ndens=1
       E=0.0D0
       E1=0.0D0
       En=0.0D0
       E2=0.0D0
       idip=1
       ngeo=ngeo+1
       Im=(0.0D0,2.0D0)
       sq2=sqrt(2.D0)
       MM=M*(M+1)/2 
       MM2=M**2
       MMd=Md*(Md+1)/2
       Md2=2*Md
       M2=2*M
!
       ALLOCATE(xnano(M,M),xnano2(M,M),fock(M,M),rhonew(M,M),
     >   rhold(M,M),rho(M,M),xmm(M,M),Y(M,M),ytrans(M,M),xtrans(M,M),
     >   rho1(M,M))
!
      if(propagator.eq.2) allocate(F1a(M,M),F1b(M,M))
      if (tdrestart) then 
         inquire(file='rho.restart',exist=exists)
         if (.not.exists) then
             write(*,*) 'ERROR CANNOT FIND rho.restart'
             write(*,*) '(if you are not restarting a previous 
     > run set tdrestart= false)'
             stop
         endif
         open(unit=1544,file='rho.restart',status='old')
         do j=1,M
            do k=1,M
               read(1544,*) rho1(j,k)
            enddo
         enddo
         do j=1,M
              do k=j,M
                  if(j.eq.k) then
                    RMM(k+(M2-j)*(j-1)/2)=REAL(rho1(j,k))
                  else
                    RMM(k+(M2-j)*(j-1)/2)=(REAL(rho1(j,k)))*2
                  endif
              enddo
         enddo
         rho=rho1
         if (propagator .eq. 2) then
            inquire(file='F1a.restart',exist=exists)
            if (.not.exists) then
               write(*,*) 'ERROR CANNOT FIND F1a.restart'
               write(*,*) '(if you are not restarting a 
     > previous run set tdrestart= false)'
               stop
            endif
            inquire(file='F1b.restart',exist=exists)
            if (.not.exists) then
               write(*,*) 'ERROR CANNOT FIND F1b.restart'
               write(*,*) '(if you are not restarting a
     > previous run set tdrestart= false)'
               stop
            endif
            open(unit=7777,file='F1a.restart',status='old')
            do i=1,M
               do j=1,M
                  read(7777,*) F1a(i,j)
               enddo
            enddo
            open(unit=7399,file='F1b.restart',status='old')
            do i=1,M
               do j=1,M
                  read(7399,*) F1b(i,j)
               enddo
            enddo
         endif 
!--------------------------------------------------------------------!
! We read the density matrix stored in RMM(1,2,3,...,MM) and it is copied in rho matrix.
         else
            call spunpack_rtc('L',M,RMM,rho)
         endif
!----------------------------------------------------------------------!
c first i
            M1=1
c now Fold
            M3=M1+MM
c now S, F also uses the same position after S was used
            M5=M3+MM
c now G
            M7=M5+MM
c now Gm
            M9=M7+MMd
c now H
            M11=M9+MMd
c W ( eigenvalues ), also this space is used in least squares
            M13=M11+MM
c aux ( vector for ESSl)
            M15=M13+M
c Least squares
            M17=M15+MM
c vectors of MO
            M18=M17+MMd
c weights (in case of using option )
            M19=M18+M*NCO
c RAM storage of two-electron integrals (if MEMO=T)
            M20 = M19 + natom*50*Nang   
c
            Nel=2*NCO+Nunp
c Initializations/Defaults
       write(*,*) ' TD CALCULATION  '
!--------------------------------------!
           niter=0
           D1=1.D0
           D2=1.D0
!--------------------------------------!
           Qc=0.0D0
           do i=1,natom
             Qc=Qc+Iz(i)
           enddo
           Qc=Qc-Nel
           Qc2=Qc**2
! FFR: Variable Allocation
!--------------------------------------------------------------------!
       allocate(Vmat(M,M),Dvec(M))
       allocate(sqsm(M,M))
       allocate(fockbias(M,M))

       dovv=.true.
       if (dovv.eq..true.) then
        if (.not.allocated(orb_group)) then
          allocate(orb_group(M))
          call atmorb(group,nuc,orb_group)
        endif
        if (.not.allocated(orb_selection)) then
          allocate(orb_selection(M))
        endif
       endif
!------------------------------------------------------------------------------!
! Two electron integral with neighbor list.
!
            do i=1,natom
              natomc(i)=0
!
              do j=1,natom
                d(i,j)=(r(i,1)-r(j,1))**2+(r(i,2)-r(j,2))**2+
     >            (r(i,3)-r(j,3))**2
                zij=atmin(i)+atmin(j)
                ti=atmin(i)/zij
                tj=atmin(j)/zij
                alf=atmin(i)*tj
                rexp=alf*d(i,j)
                if (rexp.lt.rmax) then
                  natomc(i)=natomc(i)+1
                  jatc(natomc(i),i)=j
                endif 
              enddo
            enddo
            do ii=nshell(0),1,-1
              nnps(nuc(ii))=ii
            enddo
            do ii=nshell(0)+nshell(1),nshell(0)+1,-1
              nnpp(nuc(ii))=ii
            enddo
            do ii=M,nshell(0)+nshell(1)+1,-1
              nnpd(nuc(ii))=ii
            enddo
!------------------------------------------------------------------!
c
c Create integration grid for XC here
c Assign points to groups (spheres/cubes)
c Assign significant functions to groups
c -Calculate point weights
c
      call g2g_timer_sum_start('Exchange-correlation grid setup')
      call g2g_reload_atom_positions(igrid2)
      call g2g_timer_sum_stop('Exchange-correlation grid setup')

      call aint_query_gpu_level(igpu)
      if (igpu.gt.1) call aint_new_step()

      if (predcoef.and.npas.gt.3) then

c        if (.not.OPEN) then
c          if(verbose) write(*,*) 'prediciendo densidad'
c          do i=1,MM
c            RMM(i)=(3*old1(i))-(3*old2(i))+(old3(i))
c          enddo
c         endif
       endif
!------------------------------------------------------------------------------!
! H H core, 1 electron matrix elements
      call g2g_timer_sum_start('1-e Fock')
      call g2g_timer_sum_start('Nuclear attraction')
      call int1(En)
      call g2g_timer_sum_stop('Nuclear attraction')
      if(nsol.gt.0.or.igpu.ge.4) then
          call g2g_timer_sum_start('QM/MM')
        if (igpu.le.1) then
          call g2g_timer_start('intsol')
          call intsol(E1s,Ens,.true.)
          call g2g_timer_stop('intsol')
        else
          call aint_qmmm_init(nsol,r,pc)

          call g2g_timer_start('aint_qmmm_fock')
          call aint_qmmm_fock(E1s,Ens)
          call g2g_timer_stop('aint_qmmm_fock')
        endif
          call g2g_timer_sum_stop('QM/MM')
      endif
!--------------------------------------!
            E1=0.D0
            do k=1,MM
              E1=E1+RMM(k)*RMM(M11+k-1)
            enddo
            call g2g_timer_sum_stop('1-e Fock')
!            if(ipop.eq.1)then
              allocate(overlap(M,M),rhoscratch(M,M))
              call spunpack('L',M,RMM(M5),overlap)
!            endif
! LED-DEVICE ORTHOGONALIZATION TEST
!            OPEN(unit=232323232,file='overlap')
!            OPEN(unit=131313131, file='mapmat')
!            DO i=1,M
!            DO j=1,M
!            write(232323232,*) overlap(i,j)
!            write(131313131,*) mapmat(i,j)
!            ENDDO
!            ENDDO
!            stop 'Ya escribi overlap y mapmat'
!
#ifdef CUBLAS
! Copy overlap matrix to device for Mulliken population analysis
            stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrS)
            if (stat.NE.0) then
               write(*,*) "S memory allocation failed"
               call CUBLAS_SHUTDOWN
               stop
            endif
            stat=CUBLAS_SET_MATRIX(M,M,sizeof_real,overlap,M,devPtrS,M)
            if (stat.NE.0) then
               write(*,*) "S matrix setting failed"
               call CUBLAS_SHUTDOWN
               stop
            endif
#endif
! FFR: Canonical Diagonalization of Overlap
!--------------------------------------------------------------------!
! I am keeping Y,Ytrans and Xtrans but they should be replaced
! by the much nicer Ymat,Ytrp,Xtrp (and X by Xmat). Also, copy
! into RMM.
!
         call sdiag_canonical(overlap,Dvec,Vmat,Xmm,Xtrans,Y,Ytrans)
         sqsm=matmul(Vmat,Ytrans)
         X=Xmm
         do kk=1,M
           RMM(M13+kk-1)=Dvec(kk)
         enddo
!------------------------------------------------------------------------------!
! Here rho1 is used as an auxiliar matrix. Then we need to restore its original value
#ifdef CUBLAS
            DO i=1,M
               DO j=1,M
                  rho1(i,j)=dcmplx(X(i,j),0.0D0)
               ENDDO
            ENDDO
            stat = CUBLAS_ALLOC(M*M, sizeof_real, devPtrX)
            stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrXc)
            stat = CUBLAS_ALLOC(M*M, sizeof_complex, devPtrY)
            if (stat.NE.0) then
            write(*,*) "X and/or Y memory allocation failed"
            call CUBLAS_SHUTDOWN
            stop
            endif
            stat=CUBLAS_SET_MATRIX(M,M,sizeof_complex,rho1,M,devPtrXc,M)
            stat=CUBLAS_SET_MATRIX(M,M,sizeof_real,x,M,devPtrX,M)
            DO i=1,M
               DO j=1,M
                  rho1(i,j)=dcmplx(Y(i,j),0.0D0)
               ENDDO
            ENDDO
            stat=CUBLAS_SET_MATRIX(M,M,sizeof_complex,rho1,M,devPtrY,M)
            if (stat.NE.0) then
            write(*,*) "X and/or Y setting failed"
            call CUBLAS_SHUTDOWN
            stop
            endif
#endif
! we have just used rho1 as an auxiliar matrix. Now, we restore it to the previous value:
            rho1=rho
!------------------------------------------------------------------------------!
! the transformation matrices is copied in xmm
!
            do i=1,M
               do j=1,M
                  xmm(i,j)=X(i,j)
               enddo
            enddo
!            do i=1,M
!               do j=1,M
!                 xtrans(j,i)=X(i,j)
!                 ytrans(j,i)=Y(i,j)
!               enddo
!            enddo
!------------------------------------------------------------------------------!
! External Electric Field components
!
!       write(*,*) 'fx =', fx
!       write(*,*) 'fy =', fy
!       write(*,*) 'fz =', fz
!------------------------------------------------------------------------------!
! Rho is transformed to the orthonormal basis
! TRANSPORT
             IF(TRANSPORT_CALC) THEN
                open(unit=100000,file='rhofirst')
                IF(generate_rho0) then
                   DO i=1,M
                      DO j=1,M
                         write(100000,*) rho(i,j)
                      ENDDO
                   ENDDO
                   rhofirst=rho
                   write(*,*) 'HE ESCRITO RHOFIRST'
                else            
                   DO i=1,M
                      DO j=1,M
                         read(100000,*) rhofirst(i,j)
                      ENDDO
                   ENDDO
                   write(*,*) ' HE LEIDO RHOFIRST '
                ENDIF
             ENDIF
! with matmul:
#ifdef CUBLAS
       call g2g_timer_start('complex_rho_ao_to_on-cu')
       rho=basechange_cublas(M,rho,devPtrY,'dir')
!       rho=rho1
       call g2g_timer_stop('complex_rho_ao_to_on-cu')
#else
       rho1=matmul(ytrans,rho)
       rho=matmul(rho1,y)
#endif
! TRANSPORT
       IF (TRANSPORT_CALC) THEN
         DO i=1,M
            traza0=traza0+rho(i,i)
         ENDDO
         write(*,*) 'traza0 =', traza0
       ENDIF
!
c Precalculate three-index (two in MO basis, one in density basis) matrix
c used in density fitting / Coulomb F element calculation here
c (t_i in Dunlap)
c
      call aint_query_gpu_level(igpu)
      if (igpu.gt.2) then
        call aint_coulomb_init()
      endif
      if (igpu.eq.5) MEMO = .false.
      !MEMO=.true.
      if (MEMO) then
         call g2g_timer_start('int3mem')
         call g2g_timer_sum_start('Coulomb precalc')
c Large elements of t_i put into double-precision cool here
c Size criteria based on size of pre-factor in Gaussian Product Theorem
c (applied to MO basis indices)
         call int3mem()
c Small elements of t_i put into single-precision cools here
c         call int3mems()
         call g2g_timer_stop('int3mem')
         call g2g_timer_sum_stop('Coulomb precalc')
      endif
#ifdef CUBLAS
            if(.not.TRANSPORT_CALC) call CUBLAS_FREE(devPtrY)
#endif
!TRANSPORT
          if(TRANSPORT_CALC) then
             traza10=0
             traza100=0
             alfa=1.0D0
!            GammaMagnus=0.036285
!            GammaMagnus=0.0048380
!             GammaMagnus=0.000012095
!            GammaMagnus=0.00184 ! GammaMagnus es la Gamma que se usa en los pasos de Magnus
!            GammaMagnus=0.0048380
             GammaMagnus=driving_rate
!            GammaMagnus=0.012
!            GammaMagnus=0.05
             GammaVerlet=GammaMagnus*0.1 ! GammaVerlet se usa en los primeros 200 pasos en los que se hace Verlet
             write(*,*) ' Driving Rate =', GammaMagnus
           endif
!
            call g2g_timer_stop('inicio')
!##############################################################################!
! HERE STARTS THE TIME EVOLUTION
!##############################################################################!
            write(*,*) 'PROPAGATION'
            do 999 istep=1, ntdstep
!--------------------------------------!
              call g2g_timer_start('TD step')
              if ((propagator.eq.2).and.(istep.lt.lpfrg_steps)
     >      .and. (.not.tdrestart)) then
                 t=(istep-1)*tdstep*0.1
              elseif(.not.tdrestart) then
                 t=20*tdstep
                 t=t+(istep-200)*tdstep
              else
                 t=tdstep*(istep-1)
              endif
              if (propagator.eq.1) then
                 t=(istep-1)*tdstep
              endif
              t=t*0.02419
              write(*,*) 'evolution time (fs)  =', t
!----------------------------------------------------!
            if ((propagator.eq.1).or.
     >      (((propagator.eq.2).and.(istep.lt.lpfrg_steps))
     >      .and. (.not.tdrestart))) then
            call int3lu(E2)
            call g2g_solve_groups(0,Ex,0)
            endif

         write(*,*) '! step & energy', istep,E
         E1=0.0D0
c ELECTRIC FIELD CASE - Type=gaussian (ON)
         fxx=0.0D0
         fyy=0.0D0
         fzz=0.0D0
         if(.not.TRANSPORT_CALC) then
            if(istep.lt.pert_steps) then
                if (field) then
                  write(*,*) 'FIELD ADD'
                  call dip(ux,uy,uz)
                  if (exter) then
                     g=1.0D0
                     factor=2.54D0
                     fxx=fx*exp(-0.2*(real(istep-50))**2)
                     fyy=fy*exp(-0.2*(real(istep-50))**2)
                     fzz=fz*exp(-0.2*(real(istep-50))**2)
                     write(*,*) fxx,fyy,fzz
!
                  else
                  g=2.0D0*(epsilon-1.0D0)/((2.0D0*epsilon+1.0D0)*a0**3)
                     Fx=ux/2.54D0
                     Fy=uy/2.54D0
                     Fz=uz/2.54D0
                     factor=(2.54D0*2.00D0)
!
                  endif
                     call intfld(g,Fxx,Fyy,Fzz)
                     E1=-1.00D0*g*(Fx*ux+Fy*uy+Fz*uz)/factor -
     >            0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0
                endif
            endif
         endif
         if(GATEFIELD) then 
              if((istep.le.15000).and.(istep.gt.10000))  then
                 call dip(ux,uy,uz)
                 g=1.0D0
                 factor=2.54D0
                 fxx=fx*exp(-0.00001*(dble(istep-15000))**2)
                 fyy=fy*exp(-0.00001*(dble(istep-15000))**2)
                 fzz=fz*exp(-0.00001*(dble(istep-15000))**2)
              elseif(istep.gt.15000) then
                 call dip(ux,uy,uz)
                 g=1.0D0
                 factor=2.54D0
                 fxx=fx
                 fyy=fy
                 fzz=fz
              endif  
              if(abs(fxx).lt.1E-16) fxx=0.0D0
              if(abs(fyy).lt.1E-16) fyy=0.0D0
              if(abs(fzz).lt.1E-16) fzz=0.0D0
              if((fxx.eq.0.0D0).and.(fyy.eq.0.0D0)
     >        .and.(fzz.eq.0.0D0)) then 
              else
                 call intfld(g,Fxx,Fyy,Fzz)
              endif
!              E1=-1.00D0*g*(Fxx*ux+Fyy*uy+Fzz*uz)/factor -
!     >        0.50D0*(1.0D0-1.0D0/epsilon)*Qc2/a0
              write(*,*) 'a0 =', a0
              write(*,*) 'Qc2 =', Qc2
              write(*,*) 'ux,uy,uz =', ux,uy,uz
              write(*,*) 'GATEFIELD =', fxx,fyy,fzz
         endif   
!------------------------------------------------------------------------------!
! E1 includes solvent 1 electron contributions
         do k=1,MM
           E1=E1+RMM(k)*RMM(M11+k-1)
         enddo
!        write(*,*) '1 electron contribution',E1
!------------------------------------------------------------------------------!
! Here we obtain the fock matrix in the molecular orbital (MO) basis.
! where U matrix with eigenvectors of S , and s is vector with
! eigenvalues
            if ((propagator.eq.1).or.
     >      (((propagator.eq.2).and.(istep.lt.lpfrg_steps))
     >      .and. (.not.tdrestart))) then
             call g2g_timer_start('fock')
             call spunpack('L',M,RMM(M5),fock)
#ifdef CUBLAS
!            call cumxtf(fock,devPtrX,fock,M)
!            call cumfx(fock,DevPtrX,fock,M)
!             call fock_ao_to_on(fock,devPtrX,fock,M)
             xnano2 = basechange_cublas(M,fock,devPtrX,'dir')
             fock=xnano2
#else
! This should be replaced by dgemm in order to avoid having xmm matrix which is a dummy copy of the first M*M elements of X
            xnano2=matmul(xtrans,fock)
            fock=matmul(xnano2,xmm)
!             call fock_ao_to_on(fock,x,fock,M)
#endif
             call sprepack('L',M,RMM(M5),fock)
             call g2g_timer_stop('fock')
             endif
c Now fock is stored in molecular orbital basis.
c
!  stores F1a and F1b for magnus propagation
            if((propagator.eq.2) .and. (.not.tdrestart)) then
               if(istep.eq.chkpntF1a) then
                  F1a=fock         
               endif
               if(istep.eq.chkpntF1b) then
                  F1b=fock         
               endif
            endif
            if(writedens .and. propagator.eq.1) then
               kk=istep+5
               ii=istep+15
               if(mod (kk,500) == 0) then
                 open(unit=7624,file='F1b.restart')
                 rewind 7624
                 do i=1,M
                    do j=1,M
                       write(7624,*) fock(i,j)
                    enddo
                 enddo
               endif 
               if(mod (ii,500) == 0) then
                 open(unit=7625,file='F1a.restart')
                 rewind 7625
                 do i=1,M
                    do j=1,M
                       write(7625,*) fock(i,j)
                    enddo
                 enddo
               endif
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            E=E1+E2+En
            if (sol) then
                E=E+Es
            endif
!--------------------------------------------------------------------!
            if ((propagator.eq.1).or.
     >      (((propagator.eq.2).and.(istep.lt.lpfrg_steps))
     >      .and. (.not.tdrestart))) then
           write(*,*) 'Verlet E-Propagator'
c In the first step of the propagation we extrapolate rho back in time
c using Verlet algorithm to calculate rhold.
c using matmul 
!           if(istep.eq.1) then
!             rhold=rho+(dt_lpfrg*Im*(matmul(fock,rho)))
!             rhold=rhold-(dt_lpfrg*Im*(matmul(rho,fock)))
!           endif
c using commutator
              if(istep.eq.1) then
#ifdef CUBLAS
                call g2g_timer_start('cuconmut')
                rhold=commutator_cublas(fock,rho)
                rhold=rho+dt_lpfrg*(Im*rhold)
                call g2g_timer_stop('cuconmut')
#else
                call g2g_timer_start('conmutc')
                 rhold=commutator(fock,rho)
                 rhold=rho+dt_lpfrg*(Im*rhold)
                call g2g_timer_stop('conmutc')
#endif
              endif
!####################################################################!
! DENSITY MATRIX PROPAGATION USING VERLET ALGORITHM
! using matmul:
!           rhonew=rhold-(dt_lpfrg*Im*(matmul(fock,rho)))
!           rhonew=rhonew+(dt_lpfrg*Im*(matmul(rho,fock)))
c--------------------------------------c
! TRANSPORT
             if(TRANSPORT_CALC) then
               call g2g_timer_start('TRANSPORT - b Verlet -')
               if(istep.eq.1) then 
                   open(unit=55555,file='DriveMul')
!                   open(unit=51515,file='DriveMulAtom')
               endif
               if(istep.ge.3) then
! compute the driving term for transport properties
             scratchgamma=GammaVerlet*exp(-0.0001*(dble(istep-1000))**2)    
             call ELECTROSTAT(rho1,mapmat,overlap,rhofirst,gammascratch)
                 re_traza=0.0D0
! Mulliken population analysis for the driving term
                 if((ipop.eq.1).and.
     >             (mod(istep-1,save_charge_freq*10)==0)) then
                    rhoscratch=REAL(rho1)
                    call g2g_timer_start('Mulliken charges - cu -')
                    call cumsp_r(rhoscratch,devPtrS,rhoscratch,M)
                    call g2g_timer_stop('Mulliken charges - cu -')
                    do n=1,natom
                       q(n)=0.0D0
                    enddo
                    do i=1,M
                       q(Nuc(i))=q(Nuc(i))-rhoscratch(i,i)
                    enddo
!                   do i=1,natom
!                       write(51515,*) i,i,q(i)
!                    enddo
                    if(groupcharge) then
                       qgr=0.0d0
                       do n=1,natom
                          qgr(group(n))=qgr(group(n))+q(n)
                       enddo
                       do n=1,ngroup
                          write(55555,*) n,n,qgr(n)
                          re_traza=re_traza+qgr(n)
                       enddo
                       write(55555,*) 'tot=',re_traza
                       re_traza=0
                       write(55555,*) '-------------------------'
                    endif
                 endif
               endif
! Lowdin Population
               if(lpop) then
                  if(istep.eq.1) then
                     open(unit=525252,file='DriveLowd')
!                     open(unit=535353,file='DriveLowdAtom')
                  endif
                  if((istep.ge.3).and.
     >               (mod(istep-1,save_charge_freq*10)==0)) then
                    rhoscratch=REAL(rho1)
                    do n=1,natom
                        q(n)=0.0D0
                    enddo
                    call lowdinpop(M,natom,rhoscratch,sqsm,Nuc,q)
!                 do i=1,natom
!                    write(535353,*) i,i,q(i)
!                 enddo
                       if(groupcharge) then
                      qgr=0.0d0
                      do n=1,natom
                         qgr(group(n))=qgr(group(n))+q(n)
                      enddo
                      do n=1,ngroup
                         write(525252,*) n,n,qgr(n)
                     enddo
                     write(525252,*) '-------------------------' 
                   endif
                endif    
             endif
!            rho1=matmul(ytrans,rho1)
!            rho1=matmul(rho1,y)
#ifdef CUBLAS
            call g2g_timer_start('complex_rho_ao_to_on-cu')
            rho1=basechange_cublas(M,rho1,devPtrY,'dir')
            call g2g_timer_stop('complex_rho_ao_to_on-cu')
#endif
            call g2g_timer_stop('TRANSPORT - b Verlet -')
          endif
          call g2g_timer_start('Verlet')
          call g2g_timer_start('commutator')
#ifdef CUBLAS
          rhonew=commutator_cublas(fock,rho)
          rhonew=rhold-dt_lpfrg*(Im*rhonew)
#else
          rhonew=commutator(fock,rho)
          rhonew=rhold-dt_lpfrg*(Im*rhonew)
#endif
          call  g2g_timer_stop('commutator')
          if((istep.ge.3).and.(TRANSPORT_CALC)) then
! Add the driving term to the propagation
               write(*,*) 'adding driving term to the density'
               rhonew=rhonew-rho1
          endif
c Density update (rhold-->rho, rho-->rhonew)
                    rhold=rho
                    rho=rhonew
            call g2g_timer_stop('Verlet')
! END OF VERLET PROPAGATOR
!####################################################################!
          else
          write(*,*) 'Magnus E-Propagator'
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
! DENSITY MATRIX PROPAGATION USING MAGNUS ALGORITHM
! compute the driving term for transport properties
! TRANSPORT
          if(TRANSPORT_CALC) then
              call g2g_timer_start('TRANSPORT - b magnus -')
              if(istep.eq.1) then
                  open(unit=55555,file='DriveMul')
!                  open(unit=51515,file='DriveMulAtom')
              endif
              if(istep.le.1000) then
            scratchgamma=GammaMagnus*exp(-0.0001*(dble(istep-1000))**2)
              else
            scratchgamma=GammaMagnus
              endif
            call ELECTROSTAT(rho1,mapmat,overlap,rhofirst,scratchgamma)
            re_traza=0.0D0
! Mulliken population analysis for the driving term
              if((ipop.eq.1).and.
     >          (mod(istep-1,save_charge_freq)==0)) then
                rhoscratch=REAL(rho1)
                call g2g_timer_start('Mulliken charges - cu -')
                call cumsp_r(rhoscratch,devPtrS,rhoscratch,M)
                call g2g_timer_stop('Mulliken charges - cu -')
                do n=1,natom
                   q(n)=0.0D0
                enddo
                do i=1,M
                   q(Nuc(i))=q(Nuc(i))-rhoscratch(i,i)
                enddo
!               do i=1,natom
!                 write(51515,*) i,i,q(i)
!               enddo
                if(groupcharge) then
                   qgr=0.0d0
                   do n=1,natom
                      qgr(group(n))=qgr(group(n))+q(n)
                   enddo
                   do n=1,ngroup
                       write(55555,*) n,n,qgr(n)
                       re_traza=re_traza+qgr(n)
                   enddo
                   write(55555,*) 'tot=',re_traza
                   re_traza=0
                   write(55555,*) '-------------------------'
                 endif
               endif
! Lowdin Population
              if(lpop) then
                  if(istep.eq.1) then
                     open(unit=525252,file='DriveLowd')
!                    open(unit=535353,file='DriveLowdAtom')
                  endif
                  if((istep.ge.3).and.
     >              (mod(istep-1,save_charge_freq)==0)) then
                    rhoscratch=REAL(rho1)
                    do n=1,natom
                       q(n)=0.0D0
                    enddo
                    call lowdinpop(M,natom,rhoscratch,sqsm,Nuc,q)
!              do i=1,natom
!                 write(535353,*) i,i,q(i)
!              enddo
                    if(groupcharge) then
                       qgr=0.0d0
                       do n=1,natom
                         qgr(group(n))=qgr(group(n))+q(n)
                       enddo
                       do n=1,ngroup
                          write(525252,*) n,n,qgr(n)
                       enddo
                       write(525252,*) '-------------------------'
                    endif
                  endif
              endif      
!              rho1=matmul(ytrans,rho1)
!              rho1=matmul(rho1,y)
#ifdef CUBLAS
              call g2g_timer_start('complex_rho_ao_to_on-cu')
              rho1=basechange_cublas(M,rho1,devPtrY,'dir')
              call g2g_timer_stop('complex_rho_ao_to_on-cu')
#endif
              call g2g_timer_stop('TRANSPORT - b magnus -')
          endif
#ifdef CUBLAS
                call g2g_timer_start('cupredictor')
                call cupredictor(F1a,F1b,fock,rho,devPtrX,factorial,
     > scratchgamma,fyy,fzz,g,devPtrXc) 
                call g2g_timer_stop('cupredictor')
                call g2g_timer_start('cumagnus')
                call cumagnusfac(fock,rho,rhonew,M,NBCH,dt_magnus,
     >factorial)
                call g2g_timer_stop('cumagnus')
#else
                call g2g_timer_start('predictor')
                call predictor(F1a,F1b,fock,rho,xtrans,factorial)
                call g2g_timer_stop('predictor')
                call g2g_timer_start('magnus')
               call magnus(fock,rho,rhonew,M,NBCH,dt_magnus,factorial)
                call g2g_timer_stop('magnus')
#endif
            if(TRANSPORT_CALC == .true.) then
! Add the driving term to the propagation
               write(*,*) 'adding driving term to the density'
               rhonew=rhonew-rho1
            endif
! density update and fock storage
               F1a=F1b
               F1b=fock
               rho=rhonew
! END OF MAGNUS PROPAGATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
              endif
!####################################################################!
c Here we transform the density to the atomic orbital basis and take the real part of it. The imaginary part of the density 
c can be descarted since for a basis set of purely real functions the fock matrix is real and symetric and depends only on 
c the real part of the complex density matrix. (This wont be true in the case of hybrid functionals)
c with matmul:
#ifdef CUBLAS
             call g2g_timer_start('complex_rho_on_to_ao-cu')
!             call cumxp(rho,devPtrX,rho1,M)
!             call cumpxt(rho1,devPtrX,rho1,M)
             rho1=basechange_cublas(M,rho,devPtrXc,'inv')
             call g2g_timer_stop('complex_rho_on_to_ao-cu')
#else
             call g2g_timer_start('complex_rho_on_to_ao')
!             rho1=matmul(x,rho)
!             rho1=matmul(rho1,xtrans)
              rho1=basechange(M,X,rho,Xtrans)
             call g2g_timer_stop('complex_rho_on_to_ao')
#endif
!       rho1=REAL(rho1)
c with matmulnanoc:
c (NO LONGER AVAILABLE; USE BASECHANGE INSTEAD)
c          call matmulnanoc(rho,xtrans,rho1,M)
c          rho1=basechange(M,X,rho,Xtrans)
c          rho1 = REAL(rho1)
c The real part of the density matrix in the atomic orbital basis is copied in RMM(1,2,3,...,MM) to compute the corresponding fock matrix.
              do j=1,M
                  do k=j,M
                      if(j.eq.k) then
                        RMM(k+(M2-j)*(j-1)/2)=REAL(rho1(j,k))
                      else
                        RMM(k+(M2-j)*(j-1)/2)=(REAL(rho1(j,k)))*2
                      endif
                  enddo
              enddo
! Stores the density matrix each 500 steps to be able to restart the dynamics
              if(writedens) then
                 if(mod (istep,500) == 0) then
                     open(unit=5374,file='rho.restart')
                     rewind 5374  
                     do j=1,M
                        do k=1,M
                           write(5374,*) rho1(j,k)   
                        enddo
                     enddo
                     open(unit=7624,file='F1b.restart')
                     open(unit=7625,file='F1a.restart')
                     rewind 7624
                     rewind 7625
                     do i=1,M
                        do j=1,M
                           write(7624,*) F1b(i,j)
                           write(7625,*) F1a(i,j)
                        enddo
                     enddo
                 endif
! In the last step density matrix is stored
                  if (istep.eq.ntdstep) then
                    open(unit=44,file='rholast')
                    do j=1,M
                       do k=1,M
                          write(44,*) rho1(j,k)
                       enddo
                    enddo
                  endif
              endif
!Compute the trace of the density matrix for population analysis
              if(TRANSPORT_CALC) then
              traza=dcmplx(0.0D0,0.0D0)
              DO i=1,M
                 traza=traza+rho(i,i)
              ENDDO
              write(*,*) 'TRAZA =', traza
              endif
!###################################################################!
!# DIPOLE MOMENT CALCULATION
              if(ipop.eq.1) then
              call g2g_timer_start('DIPOLE')
              if(istep.eq.1) then
                open(unit=134,file='x.dip')
                open(unit=135,file='y.dip')
                open(unit=136,file='z.dip')
                open(unit=13600,file='abs.dip')
        write(134,*) '#Time (fs) vs DIPOLE MOMENT, X COMPONENT (DEBYES)'
        write(135,*) '#Time (fs) vs DIPOLE MOMENT, Y COMPONENT (DEBYES)'
        write(136,*) '#Time (fs) vs DIPOLE MOMENT, Z COMPONENT (DEBYES)'
        write(13600,*) '#Time (fs) vs DIPOLE MOMENT (DEBYES)'
              endif
              if ((propagator.eq.2).and.(istep.lt.lpfrg_steps)
     >      .and. (.not.tdrestart)) then
                  if(mod ((istep-1),10) == 0) then
                     call dip(ux,uy,uz)
                     write(134,901) t,ux
                     write(135,901) t,uy
                     write(136,901) t,uz
                  endif
              else
                  call dip(ux,uy,uz)
                  write(134,901) t,ux
                  write(135,901) t,uy
                  write(136,901) t,uz
              endif
              call g2g_timer_stop('DIPOLE')
              endif
c u in Debyes
!# END OF DIPOLE MOMENT CALCULATION
c-------------------------MULLIKEN CHARGES-----------------------------------------------!
             if(ipop.eq.1) then  
                call g2g_timer_start('Mulliken Population')   
! open files to store Mulliken Population Analysis each step of the dynamics
                if(istep.eq.1) then
                   open(unit=1111111,file='Mulliken')
                   if (groupcharge) then
                      open(unit=678,file='MullikenGroup')
                   endif
                endif
                if ((propagator.eq.2).and.(istep.lt.lpfrg_steps)
     >          .and. (.not.tdrestart)) then
                      if(mod ((istep-1),save_charge_freq*10) == 0) then 
                         rhoscratch=REAL(rho1)
#ifdef CUBLAS
                         call cumsp_r(rhoscratch,devPtrS,rhoscratch,M)
#else
                         rhoscratch=matmul(overlap,rhoscratch)
#endif
                         do n=1,natom
                            q(n)=Iz(n)
                         enddo
                         do i=1,M
                            q(Nuc(i))=q(Nuc(i))-rhoscratch(i,i)
                         enddo
                         if(groupcharge) qgr=0.0d0
                         do n=1,natom
                            write(1111111,*) n,Iz(n),q(n)
                            if(groupcharge) then
                               qgr(group(n))=qgr(group(n))+q(n)
                            endif
                         enddo
                         if(groupcharge) then
                             do n=1,ngroup
                                write(678,*) n,n,qgr(n)
                             enddo
                             write(678,*) '-------------------------'
                         endif
                      endif
                   elseif(mod ((istep-1),save_charge_freq) == 0) then
                   rhoscratch=REAL(rho1)
#ifdef CUBLAS
                   call cumsp_r(rhoscratch,devPtrS,rhoscratch,M)
#else
                   rhoscratch=matmul(overlap,rhoscratch)
#endif
                   do n=1,natom
                      q(n)=Iz(n)
                   enddo
                   do i=1,M
                      q(Nuc(i))=q(Nuc(i))-rhoscratch(i,i)
                   enddo
                   if(groupcharge) qgr=0.0d0
                   do n=1,natom
                      write(1111111,*) n,Iz(n),q(n)
                      if(groupcharge) then
                         qgr(group(n))=qgr(group(n))+q(n)
                      endif
                   enddo
                   if(groupcharge) then
                      do n=1,ngroup
                         write(678,*) n,n,qgr(n)
                      enddo
                     write(678,*) '------------------------------------'
                   endif
               endif
               call g2g_timer_stop('Mulliken Population')
             endif
c-------------------END OF MULLIKEN CHARGES-----------------------------------------------
!-------------------LOWDIN POPULATIONS----------------------------------------------------
! open files to store Mulliken Population Analysis each step of the dynamics
                if(lpop) then
                   call g2g_timer_start('Lowdin Population')
                   if(istep.eq.1) then
!                      open(unit=12121212,file='Lowdin')
                       if (groupcharge) then
                          open(unit=13131313,file='LowdinGroup')
                       endif
                    endif
                    if ((propagator.eq.2).and.(istep.lt.lpfrg_steps)
     >                 .and. (.not.tdrestart)) then
                       if(mod ((istep-1),save_charge_freq*10) == 0) then
                          rhoscratch=REAL(rho1)
                          do n=1,natom
                            q(n)=Iz(n)
                          enddo
                          call lowdinpop(M,natom,rhoscratch,sqsm,Nuc,q)
                          if(groupcharge) qgr=0.0d0
                          do n=1,natom
!                             write(12121212,*) n,Iz(n),q(n)
                             if(groupcharge) then
                                qgr(group(n))=qgr(group(n))+q(n)
                             endif
                          enddo
                          if(groupcharge) then
                             do n=1,ngroup
                                write(13131313,*) n,n,qgr(n)
                             enddo
                             write(13131313,*) '----------------------'
                         endif
                      endif
                   elseif(mod ((istep-1),save_charge_freq) == 0) then
                   rhoscratch=REAL(rho1)
                   do n=1,natom
                         q(n)=Iz(n)
                   enddo
                   call lowdinpop(M,natom,rhoscratch,sqsm,Nuc,q)
                   if(groupcharge) qgr=0.0d0
                        do n=1,natom
!                            write(12121212,*) n,Iz(n),q(n)
                            if(groupcharge) then
                               qgr(group(n))=qgr(group(n))+q(n)
                            endif
                         enddo
                         if(groupcharge) then
                             do n=1,ngroup
                                write(13131313,*) n,n,qgr(n)
                             enddo
                             write(13131313,*) '----------------------'
                         endif
                    endif
                    call g2g_timer_stop('Lowdin Population')
               endif
!!-----------------------------------------------------------------------------------------
               write(*,*) 'Coulomb E',E2-Ex,Ex
               call g2g_timer_stop('TD step')
               write(*,*)
               if((istep.ge.10).and.(generate_rho0)) 
     >         stop 'RHO0 GENERATED'
 999           continue
!
!##############################################################################!
! HERE FINISHES THE PROPAGATION
!##############################################################################!

 995   continue
c
c
#ifdef CUBLAS
      call CUBLAS_FREE(devPtrX)
      call CUBLAS_FREE(devPtrS)
      call CUBLAS_FREE(devPtrXc)
      call CUBLAS_FREE(devPtrY)
#endif
         if (memo) then
            deallocate (kkind,kkinds)
            deallocate(cool,cools)
         endif
         if(propagator.eq.2) then
           deallocate (F1a,F1b)
         endif
         if (GRAD) then
            if(nopt.eq.0) then
              write(*,*)
              write(*,600)
              write(*,610)
              write(*,620) E1,E2-Ex,En
              if (sol) then
                 write(*,615)
                 write(*,625) Es
              endif
            endif
            write(*,*)
            write(*,450) E
         else
            E=E-Ex
         endif
c calculation of energy weighted density matrix
c
          kk=0
          do 307 j=1,M
             do 307 i=j,M
                kk=kk+1
                RMM(M15+kk-1)=0.D0
                if(i.eq.j) then
                    ff=2.D0
                else
                    ff=4.D0
                endif
                do 309 k=1,NCO
                   RMM(M15+kk-1)=RMM(M15+kk-1)-RMM(M13+k-1)
     >  *ff*X(i,M2+k)*X(j,M2+k)
 309  continue
 307   continue
c
          if (nopt.eq.0) then
c calculates Mulliken poputations
             if (ipop.eq.1) then
                call int1(En)
                do n=1,natom
                   q(n)=Iz(n)
                enddo
                do i=1,M
                   do j=1,i-1
                      kk=i+(M2-j)*(j-1)/2
                      t0=RMM(kk)*RMM(M5+kk-1)/2.D0
                      q(Nuc(i))=q(Nuc(i))-t0
                   enddo
                   kk=i+(M2-i)*(i-1)/2
                   t0=RMM(kk)*RMM(M5+kk-1)
                   q(Nuc(i))=q(Nuc(i))-t0
                   do j=i+1,M
                      kk=j+(M2-i)*(i-1)/2
                      t0=RMM(kk)*RMM(M5+kk-1)/2.D0
                      q(Nuc(i))=q(Nuc(i))-t0
                   enddo
                 enddo
                 write(*,*) 'MULLIKEN POPULATION ANALYSIS'
                 write(*,770)
                 do n=1,natom
                    write(*,760) n,Iz(n),q(n)
                 enddo
                 write(*,*)
             endif
c ELECTRICAL POTENTIAL AND POINT CHARGES EVALUATION
c
c        if (icharge.eq.1) then
c          Q1=-(2*NCO+Nunp)
c         do n=1,natom
c          Q1=Q1+Iz(n)
c         enddo
c         call charge(NORM,natom,r,Nuc,Iz,M,Md,ncont,nshell,
c     >            c,a,RMM,map,Q1)
c        endif
c
c outputs final  MO ---------------------
!      do l=1,M
c      do n=1,NCO+3
!      do n=1,M
!        X(indexii(l),M+n)=X(l,M2+n)
!      enddo
!      enddo
!
!      do l=1,M
!        write(2,400) (X(l,M+n),n=1,NCO)
!      enddo
!--------------------------------------!
! Writes down MO coefficients and orbital energies
!       write(29,*) 'ORBITAL COEFFICIENTS AND ENERGIES, CLOSED SHELL'
!       do n=1,NCO
!         write(29,850) n,RMM(M13+n-1)
!         write(29,400) (X(l,M+n),l=1,M)
!       enddo
!       do n=NCO+1,M
!         write(29,851) n,RMM(M13+n-1)
!         write(29,400) (X(l,M+n),l=1,M)
!       enddo
!       close(29)
!--------------------------------------!
      endif
c
c
c---- DEBUGGINGS
c      write(*,*) 'Exc, integrated and calculated',Exc,Ex
c      write(*,*) 'Coulomb energy',E2-Ex
c
       call g2g_timer_stop('TD')
      DEALLOCATE(xnano,xnano2,fock,rhonew,
     >   rhold,rho,xmm,xtrans,Y,ytrans,
     >   rho1)
       DEALLOCATE(factorial)
!------------------------------------------------------------------------------!
 500  format('SCF TIME ',I6,' sec')
 450  format ('FINAL ENERGY = ',F19.12)
 400  format(4(E14.7E2,2x))
 300  format(I3,E14.6,2x,F14.7)
 600  format('  ENERGY CONTRIBUTIONS IN A.U.')
 610  format(2x,'ONE ELECTRON',9x,'COULOMB',11x,'NUCLEAR')
 615  format(2x,'SOLVENT')
 620  format(F14.7,4x,F14.7,4x,F14.7)
 625  format(F14.7)
 760  format(I3,9x,I3,6x,F10.4)
 761  format(I3,9x,I3,6x,F10.4)
 770  format('ATOM #',4x,'ATOM TYPE',4x,'POPULATION')
 850  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7)
 851  format('MOLECULAR ORBITAL #',2x,I3,3x,'ORBITAL ENERGY ',F14.7,
     >    '(NON OCC.)')
 900  format(F15.9,2x,3(F15.9,2x),2x,F15.9)
 901  format(F15.9,2x,F15.9)
 777  format(4(F8.4,2x))
 776  format (3(F8.4,2x))

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
