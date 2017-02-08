!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!  CONMUTATOR - batched version -
! 
!  CONMUTATOR(MA,MB)=MC=[MA*MB-MB*MA]
!====================================================================!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
       function commutator_batched_zd(MA,MB)
     > result(MC)
       use ISO_C_BINDING
       use cublas_f
       implicit none
       complex*16,intent(in)  :: MA(:,:)
       real*8,intent(in)      :: MB(:,:)
       complex*16,allocatable :: MC(:,:)
       complex*16 :: alpha,beta
       integer                :: nn
       integer sizeof_complex
       integer*8 devPtrA
       integer*8 devPtrB
       integer*8 devPtrC
       integer*8 devPtrScratch
       integer i,j,stat
       external CUBLAS_INIT, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       external CUBLAS_SHUTDOWN, CUBLAS_ALLOC
       integer CUBLAS_ALLOC, CUBLAS_SET_MATRIX, CUBLAS_GET_MATRIX
       integer CUBLAS_INIT
       COMPLEX*16, dimension (:,:), ALLOCATABLE :: scratch
       parameter(sizeof_complex=16)
!
       stat=CUBLAS_INIT()
       if (stat.NE.0) then
           write(*,*) "initialization failed -commutator_zd"
           call CUBLAS_SHUTDOWN
           stop
       endif
       nn=size(MA,1)
       allocate(MC(nn,nn))
       allocate(scratch(nn,nn))
       alpha=(1.0D0,0.0D0)
       beta=(0.0D0,0.0D0)
       MC=(0.0D0,0.0D0)
       do i=1,nn
          do j=1,nn
             scratch(i,j)=CMPLX(MB(i,j),0.0D0)
          enddo
      enddo
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrA)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_zd"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrB)
      if (stat.NE.0) then
         write(*,*) "allocation failed -commutator_zd"
         call CUBLAS_SHUTDOWN
         stop
      endif
      stat= CUBLAS_ALLOC(nn*nn, sizeof_complex, devPtrC)
      if (stat.NE.0) then
        write(*,*) "allocation failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      
      stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MA,nn,devPtrA,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,scratch,nn,devPtrB,nn)
      if (stat.NE.0) then
        write(*,*) "matrix setting failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat = CUBLAS_SET_MATRIX(nn,nn,sizeof_complex,MC,nn,devPtrC,nn)
      if (stat.NE.0) then
        call CUBLAS_FREE( devPtrA )
        write(*,*) "matrix setting failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      call CUBLAS_ZGEMM ('N','N',nn,nn,nn,alpha,devPtrB
     > ,nn ,devPtrA,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      beta=(-1.0D0,0.0D0)
      call CUBLAS_ZGEMM ('N','N',nn,nn,nn,alpha,
     > devPtrA,nn ,devPtrB,nn, beta, devPtrC,nn)
      if (stat.NE.0) then
        write(*,*) "ZGEMM failed -conmutator_zd"
        call CUBLAS_SHUTDOWN
        stop
      endif
      stat=CUBLAS_GET_MATRIX(nn, nn, sizeof_complex, devPtrC, nn, MC,nn)
      if (stat.NE.0) then
         write(*,*) "data upload failed"
         call CUBLAS_FREE ( devPtrA )
         call CUBLAS_FREE ( devPtrB )
         call CUBLAS_FREE ( devPtrC )
         call CUBLAS_SHUTDOWN
         stop
      endif
      call CUBLAS_FREE ( devPtrA )
      call CUBLAS_FREE ( devPtrB )
      call CUBLAS_FREE ( devPtrC )
      call CUBLAS_SHUTDOWN
      DEALLOCATE(scratch)
      return;end function
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
