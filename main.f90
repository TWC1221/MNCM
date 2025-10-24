program main
    use m_driver
    use m_sparse_mm
    implicit none
    call PCG_solver !CSR_rebuild !sparse_matrix_multiplication !multigroup_diffusion_iter() !iterative_power_driver() !heterogeneous_driver_m2() !MMS_driver() !vacuum_BC_driver()
end program main

!rm -rf build
!mkdir build && cd build
!cmake ..
!make -j
