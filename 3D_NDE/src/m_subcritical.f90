module m_subcritical
    use mesh_types
    use CSR_types
    implicit none
    contains
subroutine k_subcritical(mesh, MATs, K_eff)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    type(MatMatrix), intent(in) :: MATs
    type(MeshGrid),  intent(in) :: mesh
    real(dp),        intent(in) :: K_eff

    ! Data
    real(dp), allocatable :: phi_adj(:,:), phi_crit(:,:), phi_s(:,:), src_ext(:,:)
    real(dp), allocatable :: nu_Sigma_f(:)

    integer  :: G, ncell, gg

    ! Integrals for φ*
    real(dp) :: A, B, C, D           ! A=<φ0†, s>, B=<s>, C=<φ0†,Pφ0>, D=<Pφ0>
    real(dp) :: phi_star, ks_ads

    ! Normalisation
    real(dp) :: norm_adj

    ! Diagnostics
    real(dp) :: sum_phi0, sum_phi0_adj, sum_src, sum_Pphi0
    real(dp) :: max_phi0, max_phi0_adj, max_src
    real(dp) :: L2_phi0, L2_phi0_adj, L2_src

    real(dp), allocatable :: tmp(:)

    !------------------------------------------------------------------
    ! Dimensions and material data
    !------------------------------------------------------------------
    G          = MATs%G
    nu_Sigma_f = MATs%nu_Sigma_f

    ncell = mesh%nx * mesh%ny * mesh%nz
    
    allocate(phi_adj(G,ncell), phi_crit(G,ncell), phi_s(G,ncell), src_ext(G,ncell), tmp(ncell))

    !------------------------------------------------------------------
    ! Read data
    !   phi_crit : critical eigenflux φ0
    !   phi_adj  : adjoint eigenflux φ0†
    !   phi_s    : subcritical source-driven flux (not used in φ* now)
    !   src_ext  : real external source s_r
    !------------------------------------------------------------------
    call read_vtk_scalars('../critical_flux.vtk',         G, ncell, phi_crit)
    call read_vtk_scalars('../critical_flux_adjoint.vtk', G, ncell, phi_adj)
    call read_vtk_scalars('../fixed_flux.vtk',            G, ncell, phi_s)
    call read_vtk_scalars('../src_external.vtk',          G, ncell, src_ext)

    !==================================================================
    ! 0. Basic diagnostics BEFORE any scaling
    !==================================================================
    sum_phi0     = 0.0_dp
    sum_phi0_adj = 0.0_dp
    sum_src      = 0.0_dp
    sum_Pphi0    = 0.0_dp

    max_phi0     = 0.0_dp
    max_phi0_adj = 0.0_dp
    max_src      = 0.0_dp

    L2_phi0      = 0.0_dp
    L2_phi0_adj  = 0.0_dp
    L2_src       = 0.0_dp

    do gg = 1, G
        sum_phi0     = sum_phi0     + sum( phi_crit(gg,:) * mesh%dV )
        sum_phi0_adj = sum_phi0_adj + sum( phi_adj(gg,:)  * mesh%dV )
        sum_src      = sum_src      + sum( src_ext(gg,:)  * mesh%dV )

        sum_Pphi0    = sum_Pphi0    + sum( nu_Sigma_f(gg) * phi_crit(gg,:) * mesh%dV )

        max_phi0     = max( max_phi0,     maxval(phi_crit(gg,:)) )
        max_phi0_adj = max( max_phi0_adj, maxval(phi_adj(gg,:)) )
        max_src      = max( max_src,      maxval(src_ext(gg,:)) )

        L2_phi0      = L2_phi0      + sum( (phi_crit(gg,:)**2) * mesh%dV )
        L2_phi0_adj  = L2_phi0_adj  + sum( (phi_adj(gg,:)**2)  * mesh%dV )
        L2_src       = L2_src       + sum( (src_ext(gg,:)**2)  * mesh%dV )
    end do

    write(*,*) '========== RAW FIELD DIAGNOSTICS (BEFORE NORMALIZATION) =========='
    write(*,'(A,1PE16.8)') 'sum_phi0 (∫ φ0 dV)           = ', sum_phi0
    write(*,'(A,1PE16.8)') 'sum_phi0_adj (∫ φ0† dV)      = ', sum_phi0_adj
    write(*,'(A,1PE16.8)') 'sum_src (∫ s dV)             = ', sum_src
    write(*,'(A,1PE16.8)') 'sum_Pphi0 (∫ νΣf φ0 dV)      = ', sum_Pphi0

    write(*,'(A,1PE16.8)') 'max_phi0                     = ', max_phi0
    write(*,'(A,1PE16.8)') 'max_phi0_adj                 = ', max_phi0_adj
    write(*,'(A,1PE16.8)') 'max_src                      = ', max_src

    write(*,'(A,1PE16.8)') 'L2_phi0  (||φ0||^2)          = ', L2_phi0
    write(*,'(A,1PE16.8)') 'L2_phi0_adj (||φ0†||^2)       = ', L2_phi0_adj
    write(*,'(A,1PE16.8)') 'L2_src  (||s||^2)            = ', L2_src

    !==================================================================
    ! 1. Normalise adjoint so that <φ0†, Pφ0> = 1
    !==================================================================
    norm_adj = 0.0_dp
    do gg = 1, G
        norm_adj = norm_adj + sum( phi_adj(gg,:) * nu_Sigma_f(gg) * phi_crit(gg,:) * mesh%dV )
    end do

    write(*,'(A,1PE16.8)') '<φ0†, Pφ0> BEFORE norm        = ', norm_adj

    if (norm_adj /= 0.0_dp) then
        phi_adj = phi_adj / norm_adj
    else
        write(*,*) 'WARNING: <φ0†, Pφ0> = 0; cannot normalise adjoint.'
    end if

    ! Check <φ0†, Pφ0> after normalisation
    C = 0.0_dp
    do gg = 1, G
        C = C + sum( phi_adj(gg,:) * nu_Sigma_f(gg) * phi_crit(gg,:) * mesh%dV )
    end do
    write(*,'(A,1PE16.8)') '<φ0†, Pφ0> AFTER norm         = ', C

    !==================================================================
    ! 2. Compute A,B,C,D for φ*
    !==================================================================
    A = 0.0_dp   ! <φ0†, s>
    B = 0.0_dp   ! <s>
    C = 0.0_dp   ! <φ0†, Pφ0>
    D = 0.0_dp   ! <Pφ0>

    do gg = 1, G
        A = A + sum( phi_adj(gg,:) * src_ext(gg,:) * mesh%dV )
        B = B + sum( src_ext(gg,:) * mesh%dV )

        C = C + sum( phi_adj(gg,:) * nu_Sigma_f(gg) * phi_crit(gg,:) * mesh%dV )
        D = D + sum( nu_Sigma_f(gg) * phi_crit(gg,:) * mesh%dV )
    end do

    write(*,*) '================ φ* INTEGRALS ================'
    write(*,'(A,1PE16.8)') 'A = <φ0†, s>                = ', A
    write(*,'(A,1PE16.8)') 'B = <s>                     = ', B
    write(*,'(A,1PE16.8)') 'C = <φ0†, Pφ0>              = ', C
    write(*,'(A,1PE16.8)') 'D = <Pφ0>                   = ', D

    if (B /= 0.0_dp) then
        write(*,'(A,1PE16.8)') 'A/B                         = ', A/B
    else
        write(*,*) 'WARNING: B=<s>=0'
    end if

    if (D /= 0.0_dp) then
        write(*,'(A,1PE16.8)') 'C/D                         = ', C/D
    else
        write(*,*) 'WARNING: D=<Pφ0>=0'
    end if

    if (B == 0.0_dp.or. D == 0.0_dp) then
        phi_star = 0.0_dp
        write(*,*) 'ERROR: cannot compute φ*; denominator zero.'
    else
        phi_star = (A/B) / (C/D)
    end if

    write(*,'(A,1PE16.8)') 'phi_star (φ*)               = ', phi_star

    !==================================================================
    ! 3. ks from ADS relation:
    !    (1 - 1/Keff) = φ* (1 - 1/ks)
    !    => ks = 1 / [ 1 - (1/φ*) (1 - 1/Keff) ]
    !==================================================================
    if (K_eff == 0.0_dp) then
        write(*,*) 'ERROR: K_eff = 0 in ks calculation.'
        ks_ads = 0.0_dp
    else
        ks_ads = 1.0_dp / ( 1.0_dp - (1.0_dp/phi_star) * ( 1.0_dp - 1.0_dp/K_eff ) )
    end if

    write(*,*) '================ RESULT ================'
    write(*,'(A,1PE16.8)') 'Keff                         = ', K_eff
    write(*,'(A,1PE16.8)') 'ks (ADS)                     = ', ks_ads

deallocate(phi_adj, phi_crit, phi_s, src_ext, tmp)

end subroutine k_subcritical

    subroutine read_vtk_scalars(filename, G, ncell, phi)
        implicit none

        ! ---- Arguments ----
        character(len=*), intent(in) :: filename
        integer, intent(in) :: G, ncell
        real(8), intent(out) :: phi(G, ncell)

        ! ---- Local variables ----
        character(len=256) :: line, tag
        integer :: i, gg, unit
        logical :: found(G)

        found = .false.
        unit = 21

        open(unit, file=filename, status="old", action="read")

        do
            read(unit,'(A)', end=100) line

            do gg = 1, G
                write(tag, '(A,I0)') 'SCALARS phi_g', gg
                if (.not. found(gg) .and. index(line, trim(tag)) > 0) then
                    read(unit,'(A)') line   ! skip LOOKUP_TABLE
                    do i = 1, ncell
                        read(unit,*) phi(gg,i)
                    end do
                    found(gg) = .true.
                end if
            end do

            if (all(found)) exit
        end do

    100 continue
        close(unit)
    end subroutine read_vtk_scalars

end module
