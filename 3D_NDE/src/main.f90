program main
    use m_iter
    use m_subcritical
    
    implicit none
    integer :: PCG_mode = PCG_PRECON_NONE
    integer :: max_iter = 200000

    call iter_xyz_multi(PCG_mode, max_iter, fixed_source=.false.)

    ! call iter_xyz_multi(PCG_mode, max_iter, adjoint=.true., subcritical=.true.)

 end program main

!rm -rf build
!mkdir build && cd build
!cmake ..
!make -j
