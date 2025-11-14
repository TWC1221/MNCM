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
    real(8), allocatable :: AA(:,:)

    dx = x_domain/real(nx) ; dy = y_domain/real(ny) !Linearly spaced dx, dy
    
    N = nx * ny
    est_nnz = 5 * N
    allocate(ia(N+1), ja(est_nnz), a(est_nnz), AA(N,N))

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

    AA(1:N,1:N) = 0.0d0

    do ii = 1, N
      do ptr = ia(ii), ia(ii+1)-1
        AA(ii, ja(ptr)) = a(ptr)
      end do
    end do

    !print '(25F5.2)', ((AA(ii,jj), jj=1,N), ii=1,N)

    ! Open output file
    open(unit=10, file="matrix_points.csv", status="replace", action="write")

    ! Write CSV header
    write(10,'(A)') 'X,Y,Z'

    ! Loop over sparse matrix structure
    do ii = 1, N
        do ptr = ia(ii), ia(ii+1)-1
          AA(ii, ja(ptr)) = a(ptr)
        end do
    end do

    close(10)

  end subroutine assemble_2D_diffusion_matrix
  
  subroutine assemble_rth_diffusion_matrix(nr, nth, R_domain, D, Sigma_a, ia, ja, a, dr, dth, ri)
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
    real(8), intent(in) :: R_domain, D, Sigma_a
    integer, allocatable, intent(out) :: ia(:), ja(:)
    real(8), allocatable, intent(out) :: a(:), ri(:)
    real(8), intent(out) :: dr, dth

    integer :: ii, jj, kk, row, nnz
    real(8) :: aN, aE, aS, aW, a0
    real(8), allocatable :: AA(:,:)

    allocate(AA(1:nth,1:nr), ri(nr))
    dr = R_domain/nr ; dth = 2*PI/nth

    nnz = 1
    ia(1) = 1

    do ii = 1, nr
      ri(ii) = ii*dr
      do jj = 1, nth
        row = (ii-1)*nth + jj

        aN = D*(2*ii*dr+dr)/(2*ii*dr*dr*dr)
        aS = D*(2*ii*dr-dr)/(2*ii*dr*dr*dr)
        aE = D/(ii*dr*ii*dr*dth*dth)
        aW = aE
        a0 = -(aN+aS)-2.0d0*aE + Sigma_a

        if (ii == 1) then
          a0 = a0 + aS
        else 
          kk = (ii-2)*nth + jj
          ja(nnz) = kk
          a(nnz) = aS
          nnz = nnz + 1
        end if

        if (ii == nr) then
          a0 = a0 + aN
        else
          kk = ii*Nth + jj
          ja(nnz) = kk
          a(nnz) = aN
          nnz = nnz + 1
        end if

        kk = (ii-1)*nth + merge(jj-1,Nth,jj>1)
        ja(nnz) = kk
        a(nnz) = aW

        kk = (ii-1)*nth + merge(jj-1,1,jj<nth)
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
  
    do ii = 1, nnz
      do kk = ia(ii), ia(ii+1)-1
        AA(ii, ja(kk)) = a(kk)
      end do
    end do

    print '(25F5.2)', ((AA(ii,jj), jj=1,nth), ii=1,nr)

  end subroutine assemble_rth_diffusion_matrix     

end module m_diffusion_matrix
