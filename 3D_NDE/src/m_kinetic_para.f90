
!------------------------------------------------------------------------!
! Purpose:                                                              -!
!  Calculate and print neutron kinetic parameters                       -!
!    - read in data from VTK files output by m_iter                     -!
!    - calculate phi, phi_adj and normalise                             -!
!    - output: mean generation time, delayed fraction                   -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 12/11/25     T. Charlton    Implemented 2D matrix                     -!
! 15/11/25     T. Charlton    Implemented rth matrix                    -!
! 17/11/25     T. Charlton    Implemented rz matrix                     -!
! 27/11/25     T. Charlton    Implemented 3D matrix                     -!
! 07/12/25     T. Charlton    Implemented 3D multigroup                 -!
! 10/12/25     T. Charlton    Implemented 3D current calcs              -! 
!------------------------------------------------------------------------! 
module m_kinetic_para
    use mesh_types
    implicit none
    contains
    subroutine kinetic_para(mesh, chi, G)
        type(MeshGrid), intent(in) :: mesh
        real(8), intent(in), allocatable :: chi(:)
        integer, intent(in) :: G
        
        real(8), allocatable :: pla1(:), pla2(:), pla3(:)  
        real(8), allocatable :: phi_adj(:,:), phi_for(:,:), src_ext(:,:)
        real(8), allocatable :: nu_Sigma_f(:), nu_d_sigma_f(:), v(:)

        real(8) :: meanlt, delayed, NF, NS, adj_norm
        integer :: gg, gp

        allocate(pla1(mesh%N),pla2(mesh%N),pla3(mesh%N))
        allocate(phi_adj(3,mesh%N),phi_for(3,mesh%N),src_ext(3,mesh%N))
        allocate(nu_sigma_f(3), nu_d_sigma_f(3), v(3))

        nu_sigma_f = [0.02d0, 0.1d0, 0.35d0]
        nu_d_sigma_f(1:3) = nu_sigma_f * 0.0065
        v = [1.4d9, 5.0d6, 2.2d5]   ! cm/s

        call read_vtk_scalars('../critical_flux_adjoint.vtk', mesh%nx*mesh%ny*mesh%nz, pla1, pla2, pla3)
            phi_adj(1,:) = pla1
            phi_adj(2,:) = pla2
            phi_adj(3,:) = pla3
        call read_vtk_scalars('../critical_flux.vtk', mesh%nx*mesh%ny*mesh%nz, pla1, pla2, pla3)
            phi_for(1,:) = pla1
            phi_for(2,:) = pla2
            phi_for(3,:) = pla3

        adj_norm = 0.0d0
        do gg = 1, G
            adj_norm = adj_norm + sum( phi_adj(gg,:) * phi_for(gg,:) * nu_sigma_f(gg) * mesh%dV )
        end do

        phi_adj = phi_adj / adj_norm


        NS = 0.0d0
        NF = 0.0d0

        do gg = 1, G
            NS = NS + sum( phi_adj(gg,:) * phi_for(gg,:) / v(gg) * mesh%dV )
            NF = NF + sum( phi_adj(gg,:) * phi_for(gg,:) * chi(gg) * nu_sigma_f(gg) * mesh%dV )
        end do

        meanlt = NS / NF

        NS = 0.0d0
        NF = 0.0d0

        do gp = 1, G
            do gg = 1, G
                NS = NS + sum( phi_adj(gp,:) * nu_d_sigma_f(gg) * phi_for(gg,:) * mesh%dV )
                NF = NF + sum( phi_adj(gp,:) * nu_sigma_f(gg) * phi_for(gg,:) * mesh%dV )
            end do
        end do

        !print*,sum(matmul(transpose(phi_adj * phi_for * mesh%dV), nu_sigma_f))

        delayed = NS/NF

        print*, "KINETIC PARAMETERS"

        print*, "Mean Neutron Generation Time:", meanlt
        print*, "Delayed Fraction:", delayed

    end subroutine

    subroutine read_vtk_scalars(filename, ncell, phi_g1, phi_g2, phi_g3)
        implicit none
        
        ! ---- Arguments ----
        character(len=*), intent(in) :: filename
        integer, intent(in) :: ncell
        real(8), intent(out) :: phi_g1(ncell), phi_g2(ncell), phi_g3(ncell)

        ! ---- Local variables ----
        character(len=256) :: line
        integer :: i, unit
        logical :: found1, found2, found3

        found1 = .false.
        found2 = .false.
        found3 = .false.

        unit = 21
        open(unit, file=filename, status="old", action="read")

        do
            read(unit,'(A)', end=100) line

            ! ----------- find phi_g1 block -----------
            if (index(line, "SCALARS phi_g1") > 0) then
                read(unit,'(A)') line  ! skip LOOKUP_TABLE
                do i = 1, ncell
                    read(unit,*) phi_g1(i)
                end do
                found1 = .true.
            end if

            ! ----------- find phi_g2 block -----------
            if (index(line, "SCALARS phi_g2") > 0) then
                read(unit,'(A)') line
                do i = 1, ncell
                    read(unit,*) phi_g2(i)
                end do
                found2 = .true.
            end if

            ! ----------- find phi_g3 block -----------
            if (index(line, "SCALARS phi_g3") > 0) then
                read(unit,'(A)') line
                do i = 1, ncell
                    read(unit,*) phi_g3(i)
                end do
                found3 = .true.
            end if

            if (found1 .and. found2 .and. found3) exit
        end do

    100 continue
        close(unit)

    end subroutine read_vtk_scalars
end module
