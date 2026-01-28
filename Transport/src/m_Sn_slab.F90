program sn_slab
  use m_constants
  use m_quadrature
  implicit none

  integer :: Nx, N, it, max_iters, i, nn, f
  real(8) :: a, b, tol, omega, eps, res
  type(t_Quadrature) :: Quad

  real(8), allocatable :: dx(:), sig_t(:), sig_s(:)
  real(8), allocatable :: mu(:), w(:)
  real(8), allocatable :: S(:,:), q(:,:)
  real(8), allocatable :: psi(:,:), phi(:), phi_new(:), phi_old(:)
  real(8), allocatable :: psi_in_left(:), psi_in_right(:)
  real(8) :: up, down, tau, denom, psi_i, alpha
  
  call system('pkill gnuplot')

  ! ------------- Setup -------------
  Nx = 50
  a  = 0.0_dp; b = 10.0_dp
  allocate(dx(Nx))
  dx = (b - a) / real(Nx, 8)

  alpha = 1.0 ! step

  allocate(sig_t(Nx), sig_s(Nx))
  sig_t = 0.8_dp
  sig_s = 0.7_dp

  N = 4
  call GetLineQuad(Quad, N-1)     ! IntegOrder=7 -> N=8 (no mu=0)
  allocate(mu(N), w(N))
  mu = Quad%Xi
  w  = Quad%w

  allocate(S(Nx, N))
  S = 0.0_dp
  do i = 1, Nx
    S(i, :) = 1.0_dp            ! isotropic-in-angle source, example
  end do

  allocate(psi(Nx, N)); psi = 0.0_dp
  allocate(phi(Nx));     phi = 0.0_dp
  allocate(phi_new(Nx)); phi_new = 0.0_dp
  allocate(phi_old(Nx)); phi_old = 0.0_dp
  allocate(q(Nx, N));    q = 0.0_dp
  allocate(psi_in_left(N));  psi_in_left  = 0.0_dp   ! vacuum
  allocate(psi_in_right(N)); psi_in_right = 0.0_dp   ! vacuum

  max_iters = 5000
  tol  = 1.0e-6_dp
  omega = 0.9_dp
  eps  = 1.0e-14_dp

do alpha = 0.5_dp, 1.0_dp, 0.5_dp

  ! ------------- Source iteration -------------
  do it = 1, max_iters

    do i = 1, Nx
      do nn = 1, N
        q(i, nn) = 0.5_dp * sig_s(i) * phi(i) + S(i, nn)
      end do
    end do

    ! 2) Angle sweeps
    do nn = 1, N
      if (mu(nn) > 0.0_dp) then
        up = psi_in_left(nn)
        do i = 1, Nx
          tau   = dx(i) * alpha / abs(mu(nn))
          denom = 1.0_dp + tau * sig_t(i)
          psi_i = (up + tau * q(i, nn)) / denom
          down  = (psi_i - alpha * up)/(1-alpha)

          psi(i, nn) = psi_i
          up = down  !Required for weighted Diamond equations
          if (alpha == 1.0_dp) up = psi_i !Required for step
        end do

      else if (mu(nn) < 0.0_dp) then
        up = psi_in_right(nn)
        do i = Nx, 1, -1
          tau   = dx(i) * alpha / abs(mu(nn))
          denom = 1.0_dp + tau * sig_t(i)
          psi_i = (up + tau * q(i, nn)) / denom
          down  = (psi_i - alpha * up)/(1-alpha)

          psi(i, nn) = psi_i
          up = down  !Required for weighted Diamond equations
          if (alpha == 1.0_dp) up = psi_i !Required for step
        end do
      else
        print*, "mu(nn) = 0"
        stop
      end if
    end do

    ! 3) Scalar flux (unrelaxed)
    phi_new = 0.0_dp
    do i = 1, Nx
      do nn = 1, N
        phi_new(i) = phi_new(i) + w(nn) * psi(i, nn)
      end do
    end do

    ! 4) Residual on unrelaxed phi_new vs phi_old

    res = 0.0_dp
    do i = 1, Nx
      res = max(res, abs(phi_new(i) - phi_old(i)) / (abs(phi_new(i)) + eps))
    end do

    ! 5) Underrelax phi (store into phi)
    do i = 1, Nx
      phi(i) = omega * phi_new(i) + (1.0_dp - omega) * phi(i)
    end do


    ! 6) Convergence check
    if (res < tol) exit

    phi_old = phi_new
  end do


  write(*,*) "Converged? ", (res < tol), " in iters:", it, " residual:", res

  if (alpha == 1.0_dp) open (action='write', file='SN_phi_1.0.txt', newunit=f, status='replace')
  if (alpha == 0.5_dp) open (action='write', file='SN_phi_0.5.txt', newunit=f, status='replace')

  do i = 1,Nx
    write (f,'(ES16.8,1X,ES16.8)') (i-1)*dx(1)*Nx/(Nx-1), phi_new(i)
  end do 

  close(f)

end do

call execute_command_line( &
  "gnuplot -persist -e ""set grid; set xlabel 'x'; set ylabel 'phi'; " // &
  "plot 'SN_phi_0.5.txt' using 1:2 with lines lw 2 title 'diamond difference', " // &
  "'SN_phi_1.0.txt' using 1:2 with lines lw 2 title 'step'""")

! open (action='write', file='SN_phi.txt', newunit=f, status='replace')

! do i = 1, Nx
!   write (f,'(21(ES16.8,2X))') (i-1)*dx(1)*Nx/(Nx-1), psi(i,1:N)
! end do

! close(f)

!call execute_command_line( &
!  "gnuplot -persist -e ""set grid; set xlabel 'x'; set ylabel 'psi'; " // &
!  "plot for [i=2:5] 'SN_phi.txt' using 1:i with lines lw 2 title sprintf('psi(%d)',i-1)""" )

end program sn_slab