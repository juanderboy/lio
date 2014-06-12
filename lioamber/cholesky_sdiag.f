!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       SUBROUTINE cholesky_sdiag(NM,SMvec,L,Li,U,Ui)
!------------------------------------------------------------------------------!
!
! This subroutine receives a simetric matrix in packed
! format inside 'SMvec', and returns the matrixes of
! the basis change that diagonalizes it via Cholesky
! LU decomposition, such that:
!
! S = L * U = Y * Id * Yt
!
! Id = Xt * S * X
!
!
! For more information, read:
!   (*) JChemPhys 130 224106 (2009)
!   (*) Modern Quantum Chemistry, Attila Szabo
!
! 04/2014 || F.F.R
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!
       IMPLICIT NONE
       INTEGER,INTENT(IN)                  :: NM
       REAL*8,INTENT(IN)                   :: SMvec(NM*(NM+1)/2)
       REAL*8,INTENT(OUT),DIMENSION(NM,NM) :: L,Li,U,Ui
       REAL*8,ALLOCATABLE                  :: SMcpy(:)
       INTEGER                             :: NV,ErrID,ii,jj
!
!------------------------------------------------------------------------------!
!
       NV=NM*(NM+1)/2
       ALLOCATE(SMcpy(NV))
       SMcpy = SMvec
       L     = 0.0d0
       Li    = 0.0d0
       U     = 0.0d0
       Ui    = 0.0d0

       CALL DPPTRF   ('L',NM,SMcpy,ErrID)
       CALL spunpack ('L',NM,SMcpy,L)
       CALL DTPTRI   ('L','N',NM,SMcpy,ErrID)
       CALL spunpack ('L',NM,SMcpy,Li)

       DO ii=1,NM
       DO jj=1,NM
         IF (ii.LT.jj) THEN
           L(ii,jj)  = 0.0d0
           Li(ii,jj) = 0.0d0
         ENDIF
         U(jj,ii)  =  L(ii,jj)
         Ui(jj,ii) = Li(ii,jj)
       ENDDO
       ENDDO

       DEALLOCATE(SMcpy)
       RETURN;END SUBROUTINE
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
