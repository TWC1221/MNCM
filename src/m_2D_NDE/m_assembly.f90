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
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(8), allocatable, intent(out) :: a(:), dx, dy

    integer :: N, ii, jj, kk, ptr, est_nnz
    real(8) :: aL, aB, aR, aT, aC

    dx = x_domain/real(nx) ; dy = y_domain/real(ny) !Linearly spaced dx, dy
    
    N = nx * ny
    est_nnz = 5 * N
    allocate(ia(N+1), ja(est_nnz), a(est_nnz))

    aL = -D / dx**2
    aR = -D / dx**2
    aB = -D / dy**2
    aT = -D / dy**2

    ptr = 1
    ia(1) = 1

    do jj = 1, ny
      do ii = 1, nx
        kk = ii + (jj-1)*nx !Flattening index, assign each (ii,jj) to a kk
        
        aC = Sigma_a - (aL + aR + aB + aT)
          if (ii == 1)   aC = aC - aL   ! has left neighbor
          if (ii == nx)  aC = aC - aR   ! has right neighbor
          if (jj == 1)   aC = aC - aB   ! has bottom neighbor
          if (jj == ny)  aC = aC - aT   ! has top neighbor

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
        if (ii < nx) then
            ja(ptr) = kk+1
            a(ptr)  = aR
            ptr = ptr + 1
        endif
        ! Down
        if (jj > 1) then
            ja(ptr) = kk-nx
            a(ptr)  = aB
            ptr = ptr + 1
        endif
        ! Up
        if (jj < ny) then
            ja(ptr) = kk+nx
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
    real(8), allocatable, intent(out) :: a(:), dr, dz

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

        aL = -D / dr**2 + D/(2*(ri(ii)-0.5*ii)*dr)
        aR = -D / dr**2 - D/(2*(ri(ii)+0.5*ii)*dr) 
        aB = -D / dz**2
        aT = -D / dz**2

        kk = ii + (jj-1)*nr !Flattening index, assign each (ii,jj) to a kk
        
        aC = Sigma_a - (aL + aR + aB + aT)
          if (ii == 1)   aL = aC - aL   ! has left neighbor
          if (ii == nr)  aC = aC - aR  ! has right neighbor
          if (jj == 1)   aC = aC - aB  ! has bottom neighbor
          if (jj == nz)  aC = aC - aT  ! has top neighbor

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
