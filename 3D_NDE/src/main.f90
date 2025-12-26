program main
    use m_iter
    use m_subcritical
    
    implicit none
    integer :: PCG_mode = PCG_PRECON_NONE
    integer :: max_iter = 200000
    real(8) :: splits(3)

    real(8) :: t_start, t_end
    
    ! call cpu_time(t_start)
    ! call iter_xyz_multi(PCG_mode, max_iter, adjoint=.false.)
    ! call cpu_time(t_end)

    ! splits(1) = t_end - t_start

    ! call cpu_time(t_start)
    ! call iter_xyz_multi(PCG_PRECON_NONE, max_iter, fixed_source=.false.)
    ! call cpu_time(t_end)

    ! splits(2) = t_end - t_start

    ! call cpu_time(t_start)
    ! call iter_xyz_multi(PCG_mode, max_iter, fixed_source=.false.)
    ! call cpu_time(t_end)

    ! splits(3) = t_end - t_start

    call cpu_time(t_start)
    call iter_xyz_multi(PCG_mode, max_iter, fixed_source=.true.)
    call cpu_time(t_end)

    call cpu_time(t_start)
    call iter_xyz_multi(PCG_mode, max_iter, adjoint=.true.)
    call cpu_time(t_end)

    call cpu_time(t_start)
    call iter_xyz_multi(PCG_mode, max_iter, adjoint=.false., subcritical=.true., kinetic=.false.)
    call cpu_time(t_end)

    ! call cpu_time(t_start)
    ! call iter_xyz_multi(PCG_mode, max_iter, adjoint=.false., subcritical=.true., kinetic=.true.)
    ! call cpu_time(t_end)

    ! print '(A17, F8.2, A)', 'PETsc Operations : ', splits(1), "s"
    ! print '(A17, F8.2, A)', 'PCG_PRECON_NONE  : ', splits(2), 's'
    ! print '(A18, F8.2, A)', 'PCG_PRECON_JACOBI: ', splits(3), 's'

 end program main

!rm -rf build
!mkdir build && cd build
!cmake ..
!make -j
