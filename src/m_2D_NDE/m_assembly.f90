module m_diffusion_matrix
  use m_constants
  implicit none
  contains

  subroutine assemble_2D_diffusion_matrix(nx, ny, x_domain, y_domain, D, Sigma_a, ia, ja, a, dx, dy)
    implicit none
    !---------------------------------------------------------------
    ! Assembles CSR matrix for 2D diffusion operator:
    !   -D*(d2/dx2 + d2/dy2) + Sigma_a
    ! Dirichlet boundary conditions (phi=0) on all edges.
    !
    ! Input:
    !   nx, ny     - grid points (including boundaries)
    !   dx, dy     - grid spacing
    !   D          - diffusion coefficient
    !   Sigma_a    - absorption coefficient
    !
    ! Output:
    !   ia(:)      - kk pointer (CSR)
    !   ja(:)      - column indices (CSR)
    !   a(:)       - nonzero values (CSR)
    !
    ! Notes:
    !   Lexicographic ordering: k = ii + (jj-1)*nx
    !---------------------------------------------------------------
    integer, intent(in) :: nx, ny
    real(8), intent(in) :: x_domain, y_domain, D, Sigma_a
    real(8), intent(out) :: dx, dy
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(8), allocatable, intent(out) :: a(:)

    integer :: N, ii, jj, kk, ptr, est_nnz
    real(8) :: aL, aR, aB, aT, aC

    !-----------------------------------------------------------
    ! For box method: dx,dy = cell sizes
    !-----------------------------------------------------------
    dx = x_domain / real(nx)
    dy = y_domain / real(ny)

    N = nx * ny
    est_nnz = 5*N

    allocate(ia(N+1), ja(est_nnz), a(est_nnz))
    ptr = 1
    ia(1) = 1

    do jj = 1, ny
      do ii = 1, nx
        kk = ii + (jj-1)*nx

        ! Base interior coefficients
        aL = -D * dy / dx
        aR = -D * dy / dx
        aB = -D * dx / dy
        aT = -D * dx / dy

        !-----------------------------------------------------
        ! Apply Dirichlet BCs via zero-value ghost cells
        ! Remove outside neighbors and add to diagonal
        !-----------------------------------------------------

        ! Left boundary face
        if (ii == 1) then
            ! Φ_g = 0 → diagonal += -aL ; no left entry
            aC = aC - aL
            aL = 0.0
        else
            aC = 0.0
        end if

        ! Right boundary
        if (ii == nx) then
            aC = aC - aR
            aR = 0.0
        end if

        ! Bottom
        if (jj == 1) then
            aC = aC - aB
            aB = 0.0
        end if

        ! Top
        if (jj == ny) then
            aC = aC - aT
            aT = 0.0
        end if

        ! Complete diagonal
        aC = aC + -(aL + aR + aB + aT) + Sigma_a*dx*dy

        !-----------------------------------------------------
        ! Write row in CSR
        !-----------------------------------------------------

        ! Center
        ja(ptr) = kk
        a(ptr)  = aC
        ptr = ptr + 1

        ! Left neighbor
        if (aL /= 0.0) then
            ja(ptr) = kk - 1
            a(ptr)  = aL
            ptr = ptr + 1
        end if

        ! Right
        if (aR /= 0.0) then
            ja(ptr) = kk + 1
            a(ptr)  = aR
            ptr = ptr + 1
        end if

        ! Bottom
        if (aB /= 0.0) then
            ja(ptr) = kk - nx
            a(ptr)  = aB
            ptr = ptr + 1
        end if

        ! Top
        if (aT /= 0.0) then
            ja(ptr) = kk + nx
            a(ptr)  = aT
            ptr = ptr + 1
        end if

        ia(kk+1) = ptr

      end do
    end do

    ia(N+1) = ptr

    ! Trim storage
    a  = a(1:ptr-1)
    ja = ja(1:ptr-1)

  end subroutine assemble_2D_diffusion_matrix

  subroutine assemble_rz_diffusion_matrix(nr, nz, r_domain, z_domain, D, Sigma_a, ia, ja, a, dr, dz)
    implicit none
    !---------------------------------------------------------------
    ! Assembles CSR matrix for 2D diffusion operator:
    !   -D*(d2/dx2 + d2/dy2) + Sigma_a
    ! Dirichlet boundary conditions (phi=0) on all edges.
    !
    ! Input:
    !   nx, ny     - grid points (including boundaries)
    !   dx, dy     - grid spacing
    !   D          - diffusion coefficient
    !   Sigma_a    - absorption coefficient
    !
    ! Output:
    !   ia(:)      - kk pointer (CSR)
    !   ja(:)      - column indices (CSR)
    !   a(:)       - nonzero values (CSR)
    !
    ! Notes:
    !   Lexicographic ordering: k = ii + (jj-1)*nx
    !---------------------------------------------------------------
    integer, intent(in) :: nr, nz
    real(8), intent(in) :: r_domain, z_domain, D, Sigma_a
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(8), allocatable, intent(out) :: a(:)
    real(8), intent(out) ::  dr, dz

    integer :: N, ii, jj, kk, ptr, est_nnz
    real(8) :: aL, aB, aR, aT, aC
    real(8), allocatable :: ri(:)

    dr = r_domain/real(nr) ; dz = z_domain/real(nz) !Linearly spaced dx, dy
    
    N = nr * nz
    est_nnz = 5 * N
    allocate(ia(N+1), ja(est_nnz), a(est_nnz), ri(nr))

    ptr = 1
    ia(1) = 1

    do ii = 1, nr
      ri(ii) = (ii-0.5)*dr
    end do

    do jj = 1, nz
      do ii = 1, nr

        kk = ii + (jj-1)*nr

        aL = -D*(ri(ii)-dr/2)/(dr**2*ri(ii)) !-D / dr**2 + D/(2*(ri(ii)-0.5*ii*dr)*dr)
        aR = -D*(ri(ii)+dr/2)/(dr**2*ri(ii)) !-D / dr**2 - D/(2*(ri(ii)+0.5*ii*dr)*dr) 
        aB = -D / dz**2
        aT = -D / dz**2
        
        aC = Sigma_a - (aL + aR + aB + aT)

        !-------------------------
        ! Boundary conditions
        !-------------------------

        !------------------------------------------------------------
        ! Reflective at r = 0  (Neumann: ∂φ/∂r = 0)
        ! → Left flux = 0 → aL = 0, but DO NOT change aC
        !------------------------------------------------------------
        if (ii == 1) then
            aC = aC - aL 
            aL = 0
        end if

        !------------------------------------------------------------
        ! Dirichlet at r = r_max
        !------------------------------------------------------------
        if (ii == nr) then
            aC = aC - aR
            aR = 0.0d0
        end if

        !------------------------------------------------------------
        ! Dirichlet at z = 0
        !------------------------------------------------------------
        if (jj == 1) then
            aC = aC - aB
            aB = 0.0d0
        end if

        !------------------------------------------------------------
        ! Dirichlet at z = z_max
        !------------------------------------------------------------
        if (jj == nz) then
            aC = aC - aT
            aT = 0.0d0
        end if

        ! Center
        ja(ptr) = kk
        a(ptr) = aC
        ptr = ptr + 1

        ! Left
        if (ii > 1) then
            ja(ptr) = kk-1
            a(ptr)  = aL
            ptr = ptr + 1
        endif
        ! Right
        if (ii < nr) then
            ja(ptr) = kk+1
            a(ptr)  = aR
            ptr = ptr + 1
        endif
        ! Down
        if (jj > 1) then
            ja(ptr) = kk-nr
            a(ptr)  = aB
            ptr = ptr + 1
        endif
        ! Up
        if (jj < nz) then
            ja(ptr) = kk+nr
            a(ptr)  = aT
            ptr = ptr + 1
        endif
        ia(kk+1) = ptr
      end do
    end do

    ! Finalize CSR structure
    ia(N+1) = ptr
    if (ptr-1 < size(a)) then
      a = a(1:ptr-1)
      ja = ja(1:ptr-1)
    end if

  end subroutine assemble_rz_diffusion_matrix
  
  subroutine assemble_rth_diffusion_matrix(nr, nth, ri, dr, dth, D, Sigma_a, ia, ja, a)
    implicit none
    !---------------------------------------------------------------
    ! Assembles CSR matrix for 2D diffusion operator:
    !   -D*(d2/dr2 + 1/r d/dr + 1/r2 d2/dth2) + Sigma_a
    !
    ! Input:
    !   nr, nth     - grid points (including boundaries)
    !   dr, dth     - grid spacing
    !   D          - diffusion coefficient
    !   Sigma_a    - absorption coefficient
    !
    ! Output:
    !   ia(:)      - kk pointer (CSR)
    !   ja(:)      - column indices (CSR)
    !   a(:)       - nonzero values (CSR)
    !
    ! Notes:
    !   Lexicographic ordering: k = ii + (jj-1)*nx
    !---------------------------------------------------------------
    integer, intent(in) :: nr, nth
    real(8), intent(in) :: D, Sigma_a, ri(:)
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(8), allocatable, intent(out) :: a(:)
    real(8), intent(out) :: dr, dth

    integer :: ii, jj, kk, row, nnz, est_nnz, N
    real(8) :: aN, aE, aS, aW, a0

    N = nr * nth
    est_nnz = 5 * N
    nnz = 1

    allocate(ia(nr*nth+1),ja(est_nnz),a(est_nnz))
    ia(1) = 1

    do ii = 1, nr
      do jj = 1, nth
        row = (ii-1)*nth + jj

        aN = -D*(2*ri(ii)+dr)/(2*ri(ii)*dr**2)
        aS = -D*(2*ri(ii)-dr)/(2*ri(ii)*dr**2)
        aE = -D/(ri(ii)**2 * dth**2)
        aW = aE

        a0 = -(aN+aS)-2.0d0*aE + Sigma_a

        if (ii == 1) then
          a0 = a0 - aS
        else 
          kk = (ii-2)*nth + jj
          ja(nnz) = kk
          a(nnz) = aS
          nnz = nnz + 1
        end if

        if (ii == nr) then
          kk = (ii-1)*nth + jj
          ja(nnz) = kk
          a(nnz) = 1.0d0
          nnz = nnz + 1
        else
          kk = ii*Nth + jj
          ja(nnz) = kk
          a(nnz) = aN
          nnz = nnz + 1
        end if

        kk = (ii-1)*nth + merge(jj-1,nth,jj>1)
        ja(nnz) = kk
        a(nnz) = aW
        nnz = nnz + 1

        kk = (ii-1)*nth + merge(jj+1,1,jj<nth)
        ja(nnz) = kk
        a(nnz) = aE
        nnz = nnz + 1
        
        ja(nnz) = row
        a(nnz) = a0
        nnz = nnz + 1

        ia(row+1) = nnz
      end do 
    end do 
    
    ! Finalize CSR structure
    if (nnz-1 < size(a)) then
      a = a(1:nnz-1)
      ja = ja(1:nnz-1)
    end if

  end subroutine assemble_rth_diffusion_matrix     

end module m_diffusion_matrix
