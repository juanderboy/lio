!½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½!
       subroutine dft_get_qm_forces(dxyzqm)
!--------------------------------------------------------------------!
!
!
!
!½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½!
       use garcha_mod,only:natom
       implicit none
       real*8,intent(out) :: dxyzqm(3,natom)
       real*8,allocatable :: ff1G(:,:),ffSG(:,:),ff3G(:,:)
       real*8             :: factor
       integer            :: fileunit,kk,ii
       logical            :: print_forces

!--------------------------------------------------------------------!
       allocate(ff1G(natom,3),ffSG(natom,3),ff3G(natom,3))

       call g2g_timer_start('int1G')
       ff1G=0.0d0
       call int1G(ff1G)
       call g2g_timer_stop('int1G')

       call g2g_timer_start('intSG')
       ffSG=0.0d0
       call intSG(ffSG)
       call g2g_timer_stop('intSG')

       call g2g_timer_start('int3G')
       ff3G=0.0d0
       call int3G(ff3G,.true.)
       call g2g_timer_stop('int3G')

       factor=1.D0
c       factor=627.509391D0/0.5291772108D0
       do kk=1,natom
       do ii=1,3
         dxyzqm(ii,kk)=ff1G(kk,ii)+ffSG(kk,ii)+ff3G(kk,ii)
         dxyzqm(ii,kk)=dxyzqm(ii,kk)*factor
       enddo
       enddo

!--------------------------------------------------------------------!
       print_forces=.true.
       if (print_forces) then
         fileunit=3242
         open(unit=fileunit,file='Forces.log',access='APPEND')

         write(fileunit,'(A)')
     >   '------------------------------------------------------------'
         do kk=1,natom
           write(fileunit,200) 'TOTS',kk,
     >       ff1G(kk,1)+ffSG(kk,1)+ff3G(kk,1),
     >       ff1G(kk,2)+ffSG(kk,2)+ff3G(kk,2),
     >       ff1G(kk,3)+ffSG(kk,3)+ff3G(kk,3)
         enddo
         write(fileunit,'(A)')
     >   '------------------------------------------------------------'
         write(fileunit,*)

         do kk=1,natom
           write(fileunit,200) 
     >       'FF1G',kk,ff1G(kk,1),ff1G(kk,2),ff1G(kk,3)
           write(fileunit,200)
     >       'FFSG',kk,ffSG(kk,1),ffSG(kk,2),ffSG(kk,3)
           write(fileunit,200)
     >       'FF3G',kk,ff3G(kk,1),ff3G(kk,2),ff3G(kk,3)
           write(fileunit,*)
         enddo
         write(fileunit,*)

         close(fileunit)
       endif

!--------------------------------------------------------------------!
       deallocate(ff1G,ffSG,ff3G)
 200   format(1X,A4,1X,I4,3(2X,E14.7))
       return;end subroutine
!½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½½!
