program main
    use m_iter
    use m_driver
    implicit none
    integer :: PCG_mode = PCG_PRECON_JACOBI
    integer :: max_iter = 1000
    call iter_xyz(PCG_mode, max_iter)
    
    !multigroup_diffusion_iter(4, 3)

    !call iter(PCG_mode, max_iter)
    !call PCG_solver !multigroup_diffusion_iter() !CSR_rebuild !sparse_matrix_multiplication !iterative_power_driver() !heterogeneous_driver_m2() !MMS_driver() !vacuum_BC_driver()
end program main

!rm -rf build
!mkdir build && cd build
!cmake ..
!make -j

!1.0d-10