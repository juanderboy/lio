! ---------------------------------------------------
!   Solve  Ax^2 + Bx + C = 0 
! ---------------------------------------------------
       SUBROUTINE  Quad_Solv(a,b,c,root1,root2)
       IMPLICIT  NONE
#ifdef TD_SIMPLE
       COMPLEX*8,intent(in)  :: a, b, c
       COMPLEX*8  :: d
       COMPLEX*8,intent(out)  :: root1, root2
#else
       COMPLEX*16,intent(in)  :: a, b, c
       COMPLEX*16  :: d
       COMPLEX*16,intent(out)  :: root1, root2
#endif
       d = b*b
       d= d- 4.0D0*a*c
       d     = SQRT(d)
       root1 = (-b + d)/(2.0D0*a)     ! first root
       root2 = (-b - d)/(2.0D0*a)     ! second root
       RETURN;END
! ---------------------------------------------------
!   ROOT SELECTION SUBROUTINE 
!   This subroutine selects one of the two roots given by quad solv to obtain the coeficient of density matrix rescaling in TE calculations.
!   The coeficient should be non-negative. If both roots are positive, we select the smaller one.
! ---------------------------------------------------
       SUBROUTINE  ROOT_SELECT(root1,root2,coef)
       IMPLICIT  NONE
       REAL*8 :: x1,x2
#ifdef TD_SIMPLE
       COMPLEX*8,intent(in)  :: root1,root2
       COMPLEX*8,intent(out)  :: coef
#else
       COMPLEX*16,intent(in)  :: root1,root2
       COMPLEX*16,intent(out)  :: coef
#endif
       x1=abs((root1)**2)
       x1=abs(1-x1)
       x2=abs((root2)**2)
       x2=abs(1-x2)        
       IF(x2.gt.x1) THEN
          coef=root1
       ELSE
          coef=root2
       ENDIF
       RETURN;END
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       subroutine delta_edens(group,density,edens_d,edens_a,
     > delta_edens_a, delta_edens_d,mapmat,overlap,istep,edens_d_old,
     > edens_a_old)
! calcula delta_edens_A & delta_edens_D
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       use garcha_mod
       integer :: i,j,k,kk
       integer,intent(in) :: mapmat(M,M)
       integer,intent(in) :: group(natom)
       integer,intent(in) :: istep
#ifdef TD_SIMPLE
       COMPLEX*8,intent(inout) :: edens_d,edens_a
       COMPLEX*8 :: t0,edens_prev_a,edens_prev_d
       COMPLEX*8,intent(out) :: delta_edens_a,delta_edens_d
       COMPLEX*8,intent(inout) :: edens_d_old,edens_a_old
       COMPLEX*8,intent(in) :: density(M,M)
       COMPLEX*8,allocatable :: scratch(:,:)
       COMPLEX*8 :: TOTAL
#else
       COMPLEX*16,intent(inout) :: edens_d,edens_a
       COMPLEX*16 :: t0,edens_prev_a,edens_prev_d
       COMPLEX*16,intent(out) :: delta_edens_a,delta_edens_d
       COMPLEX*16,intent(inout) :: edens_d_old,edens_a_old
       COMPLEX*16,intent(in) :: density(M,M)
       COMPLEX*16,allocatable :: scratch(:,:)
       COMPLEX*16 :: TOTAL
#endif
!#ifdef CUBLAS
!       integer*8,intent(in) :: devPtrS
!#else
       REAL*8,intent(in) :: overlap(M,M)
!#endif
       allocate(scratch(M,M))
       edens_prev_a=edens_a
       edens_prev_d=edens_d
       edens_a=dcmplx(0.0D0,0.0D0)
       edens_d=dcmplx(0.0D0,0.0D0)
!#ifdef CUBLAS
!                 call g2g_timer_start('Electrostatic density - cu -')
!                 call cumxp(density,devPtrS,scratch,M)
!                 call g2g_timer_stop('Electrostatic density - cu -')
!#else
!                 call g2g_timer_start('Electrostatic density')
!                 scratch=matmul(density,overlap)
!                 call g2g_timer_stop('Electrostatic density')
!#endif
!       do i=1,M
!             if(group(nuc(i)).eq.1) then
!                 edens_d=edens_d+scratch(i,i)
!             endif
!             if(group(nuc(i)).eq.2) then
!                 edens_a=edens_a+scratch(i,i)
!             endif
!        enddo
           DO i=1,M
             DO k=1,M
                t0=density(i,k)*overlap(k,i)
                IF(mapmat(i,k).eq.1) THEN
                   edens_d=edens_d+t0
!                   write(999,*) mapmat(i,k), density(i,k), density(k,i)
                ELSE IF(mapmat(i,k).eq.2) THEN
!                   edens_d=edens_d+t0
!                   write(999,*) mapmat(i,k), density(i,k), density(k,i)
                ELSE IF(mapmat(i,k).eq.3) THEN
!                   edens_d=edens_d+t0
!                   write(999,*) mapmat(i,k), density(i,k), density(k,i)
                ELSE IF(mapmat(i,k).eq.4) THEN
!                   edens_a=edens_a+t0
!                   write(999,*) mapmat(i,k), density(i,k), density(k,i)
                ELSE IF(mapmat(i,k).eq.5) THEN
                   edens_a=edens_a+t0
!                   write(999,*) mapmat(i,k), density(i,k), density(k,i)
                ELSE IF(mapmat(i,k).eq.6) THEN
!                   edens_a=edens_a+t0
!                   write(999,*) mapmat(i,k), density(i,k), density(k,i)
                ENDIF
             ENDDO
           ENDDO
!         call HERMITICITY(density,MAPMAT,M)
!        TOTAL=cmplx(0.0D0,0.0D0)
!       DO i=1,M
!           if(mapmat(i,i).eq.1) then
!              edens_d=edens_d+scratch(i,i)
!           endif
!           if(mapmat(i,i).eq.5) then
!              edens_a=edens_a+scratch(i,i)
!           endif
!           TOTAL=TOTAL+scratch(i,i)
!        ENDDO
!        write(*,*) 'TOTAL =', TOTAL
        if(istep.eq.1) then
           edens_a_old=edens_a
           edens_d_old=edens_d
        endif
        write(*,*) 'edens_a =', edens_a
        write(*,*) 'edens_d =', edens_d
        delta_edens_a=edens_a-edens_a_old
        delta_edens_d=edens_d-edens_d_old
        write(*,*) 'edens_d_old =', edens_d_old
        write(*,*) 'edens_a_old =', edens_a_old
        write(*,*) 'delta_edens_a =', delta_edens_a
        write(*,*) 'delta_edens_d =', delta_edens_d
        write(*,*) 'current_a =', edens_a-edens_prev_a
        write(*,*) 'current_d =', edens_d-edens_prev_d
        deallocate(scratch)
        return;end
! ---------------------------------------------------
!   MAT_MAP  SUBROUTINE 
! ---------------------------------------------------
       SUBROUTINE  MAT_MAP(group,mapmat)
       USE garcha_mod , ONLY: Nuc,M,natom
       integer,intent(out) :: mapmat(M,M)
       integer,intent(in) :: group(natom)
       integer :: nn,n,k
!====================================================!
       write(*,*) 'M =', M
       write(*,*) 'natoms=', natom
       mapmat=0
       n=0
       nn=0
       k=0
       DO i=1,M
        DO j=1,M
         IF((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.1)) mapmat(i,j)=1
         IF((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.2)) mapmat(i,j)=2
         IF((group(nuc(i)).eq.1).and.(group(nuc(j)).eq.3)) mapmat(i,j)=3
         IF((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.1)) mapmat(i,j)=4
         IF((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.2)) mapmat(i,j)=5
         IF((group(nuc(i)).eq.2).and.(group(nuc(j)).eq.3)) mapmat(i,j)=6
         IF((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.1)) mapmat(i,j)=7
         IF((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.2)) mapmat(i,j)=8
         IF((group(nuc(i)).eq.3).and.(group(nuc(j)).eq.3)) mapmat(i,j)=9
        ENDDO
       ENDDO
!TEST
       DO i=1,M
             if(mapmat(i,i).eq.1) k=k+1
             if(mapmat(i,i).eq.5) nn=nn+1
             if(mapmat(i,i).eq.9) n=n+1
       ENDDO
       write(*,*) 'Basis from group 1 =', k
       write(*,*) 'Basis from group 2 =', nn
       write(*,*) 'Basis from group 3 =', n
       RETURN;END
!---------------------------------------------------!
!  HERMITICITY TEST SUBROUTINE
!---------------------------------------------------!
       SUBROUTINE HERMITICITY(MATRIX,MAPMAT,M)
       COMPLEX*16,intent(in) :: MATRIX(M,M)
       integer,intent(in) :: mapmat(M,M)
       COMPLEX*16 :: scratch
       DO i=1,M
          DO j=i,M
             IF (MAPMAT(i,j).eq.1) THEN
                write(999,*) MATRIX(i,j)+MATRIX(j,i)
             ENDIF
             IF (MAPMAT(i,j).eq.5) THEN
                write(1000,*) MATRIX(i,j)+MATRIX(j,i)
             ENDIF
          ENDDO
       ENDDO
       write(999,*) '<------------------------------>'
       write(1000,*) '<----------------------------->'
       RETURN;END
!----------------------------------------------------!
! ELECTROSTAT SUBROUTINE
! computes the driven term needed for transport calculations
! rho1 is the variable where the density is stored and where the driven term is retrieved
!----------------------------------------------------!
!       SUBROUTINE ELECTROSTAT(rho1,mapmat,overlap,rhofirst,Gamma0)
!       use garcha_mod
!       integer, intent(in) :: mapmat(M,M)
!       REAL*8,intent(in)  :: overlap(M,M)
!       REAL*8 :: deltaa, deltad,dif
!       logical :: TRACE_INVARIANT
!       integer :: iter
!       REAL*8 :: c1,c2,re_traza
!       REAL*8,intent(in)  :: Gamma0
!       real*8 :: GammaIny, GammaAbs
!#ifdef TD_SIMPLE
!       COMPLEX*8  :: delta_a,delta_d
!       COMPLEX*8,intent(inout) :: rho1(M,M)
!       COMPLEX*8  :: AA,BB,CC,DD,G,H,Q,RR,c01,c02,root1,root2
!       COMPLEX*8,ALLOCATABLE :: rho_scratch(:,:,:)
!       COMPLEX*8 :: t0, AAA, BBB, CCC,traza,suma
!       COMPLEX*8 :: edens_d,edens_a
!       COMPLEX*8,intent(in) :: rhofirst(M,M)
!#else
!       COMPLEX*16  :: delta_a,delta_d !,trazain,traza0
!       COMPLEX*16, intent(inout) :: rho1(M,M)
!       COMPLEX*16  :: AA,BB,CC,DD,G,H,Q,RR,c01,c02,root1,root2
!       COMPLEX*16,ALLOCATABLE :: rho_scratch(:,:,:)
!       COMPLEX*16 :: t0, AAA, BBB, CCC,traza,suma
!       COMPLEX*16 :: edens_d,edens_a
!       COMPLEX*16,intent(in) :: rhofirst(M,M)
!#endif
!!====================================================!
!       call g2g_timer_start('electrostat')
!         ALLOCATE(rho_scratch(M,M,2))
!         rho_scratch=0
!         DO i=1,M
!            DO j=1,M
!               IF(mapmat(i,j).eq.0) THEN
!                 rho_scratch(i,j,1)=dcmplx(0.0D0,0.0D0)
!                 rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0)
!               ENDIF
!               IF(mapmat(i,j).eq.9) THEN
!                 rho_scratch(i,j,1)=dcmplx(0.0D0,0.0D0)
!                 rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0)
!               ENDIF
!               IF(mapmat(i,j).eq.1) THEN
!                  rho_scratch(i,j,1)=(rho1(i,j))
!                  rho_scratch(i,j,2)=(rhofirst(i,j))
!               ENDIF
!               IF((mapmat(i,j).eq.3).or.(mapmat(i,j).eq.7)) THEN
!                  rho_scratch(i,j,1)=(0.50D0*rho1(i,j))
!                  rho_scratch(i,j,2)=(0.50D0*rhofirst(i,j))
!!                  rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0) !mini-prueba
!               ENDIF
!               IF((mapmat(i,j).eq.2).or.(mapmat(i,j).eq.4)) THEN
!                  rho_scratch(i,j,1)=rho1(i,j)
!                  rho_scratch(i,j,2)=(rhofirst(i,j))
!!                  rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0) !mini-prueba
!               ENDIF
!               IF(mapmat(i,j).eq.5) THEN
!                  rho_scratch(i,j,1)=rho1(i,j)
!                  rho_scratch(i,j,2)=rhofirst(i,j)
!               ENDIF
!               IF((mapmat(i,j).eq.6).or.(mapmat(i,j).eq.8)) THEN
!                  rho_scratch(i,j,1)=(0.50D0*rho1(i,j))
!                  rho_scratch(i,j,2)=(0.50D0*rhofirst(i,j))
!!                  rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0) !mini-prueba 
!               ENDIF
!            ENDDO
!         ENDDO
!! sacamos lo que sigue para que quede como el metodo de Hodded Hod
!         edens_a=dcmplx(0.0D0,0.0D0)
!         edens_d=dcmplx(0.0D0,0.0D0)
!         DO i=1,M
!            DO k=1,M
!               t0=rho_scratch(i,k,1)*overlap(k,i)
!               edens_d=edens_d+t0
!               t0=rho_scratch(i,k,2)*overlap(k,i)
!               edens_a=edens_a+t0
!            ENDDO
!         ENDDO
!         write(*,*) 'Absorption =', edens_d 
!         write(*,*) 'Inyection =', edens_a
!         GammaAbs=Gamma0*edens_a
!         GammaAbs=GammaAbs/(edens_a+edens_d)
!         GammaIny=Gamma0-GammaAbs
!         write(*,*) 'GammaAbs,GammaIny =',GammaAbs,GammaIny
!! Usamos GammaAbs=GammaIny=Gamma0
!!        Gammainy=Gamma0
!!         GammaAbs=Gamma0
!!
!         DO i=1,M
!            DO j=1,M
!               rho1(i,j)= (GammaAbs*rho_scratch(i,j,1))-
!     >         (GammaIny*rho_scratch(i,j,2))
!            ENDDO
!         ENDDO
!!----PRUEBA - NOBORRAR -------!
!!         t0=cmplx(0.0D0,0.0D0)
!!         edens_a=cmplx(0.0d0,0.0d0)
!!         edens_d=cmplx(0.0d0,0.0d0)
!!         DO i=1,M
!!            DO k=1,M
!!               t0=rho(i,k)*overlap(k,i)
!!               IF(mapmat(i,k).eq.1) THEN
!!                  edens_d=edens_d+t0
!!               ELSE IF(mapmat(i,k).eq.5) THEN
!!                  edens_a=edens_a+t0
!!               ENDIF
!!            ENDDO
!!         ENDDO
!!         write(*,*) 'EDENS_A - DESPUES - =', edens_a
!!         write(*,*) 'EDENS_D - DESPUES - =', edens_d
!!---- Paramos si encontramos NaN -----------!
!!         DO i=1,M
!!            DO j=1,M
!!               write(999999,*) mapmat(i,j), rho1(i,j)
!!               IF(rho1(i,j).ne.rho1(i,j)) THEN
!!               stop 'Huston, we have a problem'
!!               ENDIF
!!            ENDDO
!!         ENDDO
!         DEALLOCATE(rho_scratch)
!        call g2g_timer_stop('electrostat')
!        RETURN;END
!!============================================================================================!
       SUBROUTINE ELECTROSTAT(rho1,mapmat,overlap,rhofirst,Gamma0)
       use garcha_mod
       integer, intent(in) :: mapmat(M,M)
       REAL*8,intent(in)  :: overlap(M,M)
       REAL*8 :: deltaa, deltad,dif
       logical :: TRACE_INVARIANT
       integer :: iter
       REAL*8 :: c1,c2,re_traza
       REAL*8,intent(in)  :: Gamma0
       real*8 :: GammaIny, GammaAbs
#ifdef TD_SIMPLE
       COMPLEX*8  :: delta_a,delta_d
       COMPLEX*8,intent(inout) :: rho1(M,M)
       COMPLEX*8  :: AA,BB,CC,DD,G,H,Q,RR,c01,c02,root1,root2
       COMPLEX*8,ALLOCATABLE :: rho_scratch(:,:,:)
       COMPLEX*8 :: t0, AAA, BBB, CCC,traza,suma
       COMPLEX*8 :: edens_d,edens_a
       COMPLEX*8,intent(in) :: rhofirst(M,M)
#else
       COMPLEX*16  :: delta_a,delta_d !,trazain,traza0
       COMPLEX*16, intent(inout) :: rho1(M,M)
       COMPLEX*16  :: AA,BB,CC,DD,G,H,Q,RR,c01,c02,root1,root2
       COMPLEX*16,ALLOCATABLE :: rho_scratch(:,:,:)
       COMPLEX*16 :: t0, AAA, BBB, CCC,traza,suma
       COMPLEX*16 :: edens_d,edens_a
       COMPLEX*16,intent(in) :: rhofirst(M,M)
#endif
!====================================================!
       call g2g_timer_start('electrostat')
         ALLOCATE(rho_scratch(M,M,2))
         rho_scratch=0
         DO i=1,M
            DO j=1,M
               IF(mapmat(i,j).eq.0) THEN
                 rho_scratch(i,j,1)=dcmplx(0.0D0,0.0D0)
                 rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0)
               ENDIF
               IF(mapmat(i,j).eq.9) THEN
                 rho_scratch(i,j,1)=dcmplx(0.0D0,0.0D0)
                 rho_scratch(i,j,2)=dcmplx(0.0D0,0.0D0)
               ENDIF
               IF(mapmat(i,j).eq.1) THEN
                  rho_scratch(i,j,1)=(rho1(i,j))
                  rho_scratch(i,j,2)=(rhofirst(i,j))
               ENDIF
               IF((mapmat(i,j).eq.3).or.(mapmat(i,j).eq.7)) THEN
                  rho_scratch(i,j,1)=(0.50D0*rho1(i,j))
                  rho_scratch(i,j,2)=(0.50D0*rhofirst(i,j))
               ENDIF
               IF((mapmat(i,j).eq.2).or.(mapmat(i,j).eq.4)) THEN
                  rho_scratch(i,j,1)=rho1(i,j)
                  rho_scratch(i,j,2)=(rhofirst(i,j))
               ENDIF
               IF(mapmat(i,j).eq.5) THEN
                  rho_scratch(i,j,1)=rho1(i,j)
                  rho_scratch(i,j,2)=rhofirst(i,j)
               ENDIF
               IF((mapmat(i,j).eq.6).or.(mapmat(i,j).eq.8)) THEN
                  rho_scratch(i,j,1)=(0.50D0*rho1(i,j))
                  rho_scratch(i,j,2)=(0.50D0*rhofirst(i,j))
               ENDIF
            ENDDO
         ENDDO
!         edens_a=dcmplx(0.0D0,0.0D0)
!         edens_d=dcmplx(0.0D0,0.0D0)
!         DO i=1,M
!            DO k=1,M
!               t0=rho_scratch(i,k,1)*overlap(k,i)
!               edens_d=edens_d+t0
!               t0=rho_scratch(i,k,2)*overlap(k,i)
!               edens_a=edens_a+t0
!            ENDDO
!         ENDDO
!         write(*,*) 'Absorption =', edens_d
!         write(*,*) 'Inyection =', edens_a
!         GammaAbs=Gamma0*edens_a
!         GammaAbs=GammaAbs/(edens_a+edens_d)
!         GammaIny=Gamma0-GammaAbs
         GammaIny=Gamma0*0.5D0
         GammaAbs=GammaIny
         write(*,*) 'GammaAbs,GammaIny =',GammaAbs,GammaIny
         DO i=1,M
            DO j=1,M
               rho1(i,j)= (GammaAbs*rho_scratch(i,j,1))-
     >         (GammaIny*rho_scratch(i,j,2))
            ENDDO
         ENDDO
!----PRUEBA - NOBORRAR -------!
!         t0=cmplx(0.0D0,0.0D0)
!         edens_a=cmplx(0.0d0,0.0d0)
!         edens_d=cmplx(0.0d0,0.0d0)
!         DO i=1,M
!            DO k=1,M
!               t0=rho(i,k)*overlap(k,i)
!               IF(mapmat(i,k).eq.1) THEN
!                  edens_d=edens_d+t0
!               ELSE IF(mapmat(i,k).eq.5) THEN
!                  edens_a=edens_a+t0
!               ENDIF
!            ENDDO
!         ENDDO
!         write(*,*) 'EDENS_A - DESPUES - =', edens_a
!         write(*,*) 'EDENS_D - DESPUES - =', edens_d
!---- Paramos si encontramos NaN -----------!
         DO i=1,M
            DO j=1,M
!               write(999999,*) mapmat(i,j), rho1(i,j)
               IF(rho1(i,j).ne.rho1(i,j)) THEN
               stop 'Huston, we have a problem'
               ENDIF
            ENDDO
         ENDDO
         DEALLOCATE(rho_scratch)
        call g2g_timer_stop('electrostat')
        RETURN;END
!============================================================================================!
                                            
