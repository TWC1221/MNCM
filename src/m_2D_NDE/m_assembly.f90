module m_diffusion_matrix
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
            ! Write each nonzero as (column, row, value)
            write(10,'(I0,",",I0,",",ES16.8)') ja(ptr), ii, a(ptr)
        end do
    end do

    close(10)

  end subroutine assemble_2D_diffusion_matrix

  subroutine assemble_2D_g_diffusion_matrix(nx, ny, x_domain, y_domain, D, Sigma_a, ia, ja, a, dx, dy)
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
            ! Write each nonzero as (column, row, value)
            write(10,'(I0,",",I0,",",ES16.8)') ja(ptr), ii, a(ptr)
        end do
    end do

    close(10)

  end subroutine assemble_2D_g_diffusion_matrix

end module m_diffusion_matrix
