program main
    use m_driver
    implicit none
    call multigroup_diffusion_iter() !iterative_power_driver() !heterogeneous_driver_m2() !MMS_driver() !vacuum_BC_driver()
end program main

!rm -rf build
!mkdir build && cd build
!cmake ..
!make -j
