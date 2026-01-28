!------------------------------------------------------------------------!
! Purpose:                                                              -!
!  1D spherical S_N discrete ordinates transport solver (steady state)  -!
!  - Conservative shell balance with spherical face areas and volumes   -!
!  - Gauss–Legendre quadrature (N-1 to avoid μ=0)                       -!
!  - Angular differencing coefficients (EE Lewis CMNT formulation)      -!
!  - Isotropic external source and isotropic scattering                 -!
!  - Source iteration with inward/outward diamond-difference sweeps     -!
!  - Vacuum outer boundary; center regularity via angular differencing  -!
!  - Outputs phi.dat (r, φ, ψ_n) and a gnuplot script                   -!
!  - Colors ψ_n curves by their μ with colorbar scaled to [-1, 1]       -!
!------------------------------------------------------------------------!
!! Record of revisions:                                                 -!
!   Date       Programmer     Description of change                     -!
!   ====       ==========     =====================                     -!
! 01/14/26     T. Charlton    Initial 1D spherical S_N solver; file     -!
!                             output and μ-colored gnuplot visualization-!
!------------------------------------------------------------------------!

program sn_XY
  use m_constants
  use m_quadrature
  implicit none
  type(t_sn_quadrature) :: sn_Quad

  integer :: N, I, J, nn, ii, jj, Nang

  real(8), allocatable :: X, Y, x_edge(:), y_edge(:), dx(:), dy(:)
  
  real(8), allocatable :: mu(:), eta(:), w(:), scale(:)

  real(8), allocatable :: sigs(:,:), sigt(:,:)

  real(8), allocatable :: psi(:,:,:), psi_i(:,:,:), psi_j(:,:,:)

  real(8), allocatable :: S(:,:,:), q(:,:,:), phi(:,:), phi_prime(:,:)

  real(8), allocatable :: q_MMS(:,:,:), S_MMS(:,:,:), x_MMS(:), y_MMS(:)

  integer :: iter, max_iter

  integer :: udat, uscr
  real(dp), allocatable :: x_center(:), y_center(:)
  character(len=256) :: cmd

  call system('pkill gnuplot')

  N = 2
  I = 100
  J = 100
  X = 1.0
  Y = 1.0

  allocate(x_edge(I+1), Y_edge(J+1), dx(I), dy(J))
  x_edge = [((X/I) * ii, ii = 0,I)]
  y_edge = [((Y/J) * jj, jj = 0,J)]

  dx = [(x_edge(ii+1)-x_edge(ii),ii=1,I)]
  dy = [(y_edge(jj+1)-y_edge(jj),jj=1,J)]

  call Get2DAngleQuadrature(sn_Quad, N)
  Nang = sn_Quad%NoAngles 
  mu = sn_Quad%Angles(1:Nang,1)
  eta = sn_Quad%Angles(1:Nang,2)  
  w = sn_Quad%w(1:Nang)

  ! print'(40F5.2)',mu
  ! print'(40F5.2)',eta
  ! print'(40F5.2)',w

  allocate(sigs(I,J), sigt(I,J))
  sigs = 0.9
  sigt = 1.0
  
  ! allocate(x_MMS(I), y_MMS(J), q_MMS(Nang,I,J),S_MMS(Nang,I,J))
  ! do nn = 1,Nang
  !   do ii = 1,I
  !     do jj = 1,J
  !       x_MMS(ii) = (x_edge(ii) + x_edge(ii+1))/2.0
  !       y_MMS(jj) = (y_edge(jj) + y_edge(jj+1))/2.0
  !       q_MMS(nn,ii,jj) = mu(nn)*(y_MMS(jj)-sin(2.0*pi*x_MMS(ii)) - 2.0*pi*(x_MMS(ii)+cos(2.0*pi*y_MMS(jj)))*cos(2.0*pi*x_MMS(ii))) + eta(nn)*(x_MMS(ii) + cos(2.0*pi*y_MMS(jj)) - 2.0*pi*(y_MMS(jj)-sin(2.0*pi*x_MMS(ii)))*sin(2.0*pi*y_MMS(jj)))
  !       S_MMS(nn,ii,jj) = q_MMS(1,ii,jj) + (sigt(ii,jj)-sigs(ii,jj)) * (x_MMS(ii)+cos(2.0*pi*y_MMS(jj)))*(y_MMS(jj)-sin(2.0*pi*x_MMS(ii)))
  !     end do
  !   end do
  ! end do

  allocate(psi(Nang,I,J),psi_i(Nang,I+1,J),psi_j(Nang,I,J+1), q(Nang,I,J), S(Nang,I,J), phi(I,J), phi_prime(I,J))
  q = 0.0
  S(:,:,:) = 1.0 
  phi_prime = 0.0

  do iter = 1, 100

    q = 0.0
    do nn = 1, Nang
      do jj = 1, J
        do ii = 1, I
          q(nn,ii,jj) = (sigs(ii,jj) * phi_prime(ii,jj)) + S(nn,ii,jj)
        end do
      end do
    end do
    
    phi = 0.0

    do nn = 1,Nang
      if (mu(nn) >= 0 .and. eta(nn) >=0) then
        psi_i(nn,1,1:J) = 0.0
        psi_j(nn,1:I,1) = 0.0
        do jj = 1,J
          do ii = 1,I
            psi(nn,ii,jj) = (q(nn,ii,jj) + 2.0*abs(mu(nn))/dx(ii) * psi_i(nn,ii,jj) + 2.0*abs(eta(nn))/dy(jj)*psi_j(nn,ii,jj))/(2.0*abs(mu(nn))/dx(ii)+2.0*abs(eta(nn))/dy(jj) + sigt(ii,jj))
            psi_i(nn,ii+1,jj) = 2.0 * psi(nn,ii,jj) - psi_i(nn,ii,jj)
            psi_j(nn,ii,jj+1) = 2.0 * psi(nn,ii,jj) - psi_j(nn,ii,jj)
          end do
        end do
            
        do jj = 1,J
          do ii = 1,I
            phi(ii,jj) = phi(ii,jj) + 1.0/(4.0) * w(nn)*psi(nn,ii,jj)
          end do
        end do
      
      else if (mu(nn) >= 0 .and. eta(nn) <= 0) then
        psi_i(nn,1,1:J) = 0.0
        psi_j(nn,1:I,J+1) = 0.0
        do jj = J, 1, -1
          do ii = 1,I
            psi(nn,ii,jj) = (q(nn,ii,jj) + 2.0*abs(mu(nn))/dx(ii) * psi_i(nn,ii,jj) + 2.0*abs(eta(nn))/dy(jj)*psi_j(nn,ii,jj+1))/(2.0*abs(mu(nn))/dx(ii)+2.0*abs(eta(nn))/dy(jj) + sigt(ii,jj))
            psi_i(nn,ii+1,jj) = 2.0 * psi(nn,ii,jj) - psi_i(nn,ii,jj)
            psi_j(nn,ii,jj) = 2.0 * psi(nn,ii,jj) - psi_j(nn,ii,jj+1)
          end do
        end do

        do jj = 1,J
          do ii = 1,I
            phi(ii,jj) = phi(ii,jj) + 1.0/(4.0) * w(nn)*psi(nn,ii,jj)
          end do
        end do
      
      else if (mu(nn) <= 0 .and. eta(nn) <= 0) then
        psi_i(nn,I+1,1:J) = 0.0
        psi_j(nn,1:I,J+1) = 0.0
        do jj = J, 1, -1
          do ii = I, 1, -1
            psi(nn,ii,jj) = (q(nn,ii,jj) + 2.0*abs(mu(nn))/dx(ii) * psi_i(nn,ii+1,jj) + 2.0*abs(eta(nn))/dy(jj)*psi_j(nn,ii,jj+1))/(2.0*abs(mu(nn))/dx(ii)+2.0*abs(eta(nn))/dy(jj) + sigt(ii,jj))
            psi_i(nn,ii,jj) = 2.0*psi(nn,ii,jj) - psi_i(nn,ii+1,jj)
            psi_j(nn,ii,jj) = 2.0 * psi(nn,ii,jj) - psi_j(nn,ii,jj+1)
          end do
        end do
            
        do jj = 1,J
          do ii = 1,I
            phi(ii,jj) = phi(ii,jj) + 1.0/(4.0) * w(nn)*psi(nn,ii,jj)
          end do
        end do

      else if (mu(nn) <= 0 .and. eta(nn) >= 0) then
        psi_i(nn,I+1,1:J) = 0.0
        psi_j(nn,1:I,1) = 0.0
        do jj = 1,J
          do ii = I, 1, -1
            psi(nn,ii,jj) = (q(nn,ii,jj) + 2.0*abs(mu(nn))/dx(ii) * psi_i(nn,ii+1,jj) + 2.0*abs(eta(nn))/dy(jj)*psi_j(nn,ii,jj))/(2.0*abs(mu(nn))/dx(ii)+2.0*abs(eta(nn))/dy(jj) + sigt(ii,jj))
            psi_i(nn,ii,jj) = 2.0 * psi(nn,ii,jj) - psi_i(nn,ii+1,jj)
            psi_j(nn,ii,jj+1) = 2.0 * psi(nn,ii,jj) - psi_j(nn,ii,jj)
          end do
        end do
        
        do jj = 1,J
          do ii = 1,I
            phi(ii,jj) = phi(ii,jj) + 1.0/(4.0) * w(nn)*psi(nn,ii,jj)
          end do
        end do
      end if

    end do
    
    if (sum(phi)-sum(phi_prime) <= 10E-8) print*,iter
    if ( maxval(abs(phi - phi_prime)) < 1d-8 ) exit
    phi_prime = phi

  end do 

print*,phi(1,1),phi(I/2,J/2)

! -----------------------------
! Write phi(x,y) as grid for pm3d
! -----------------------------

allocate(x_center(I), y_center(J))
do ii = 1, I
  x_center(ii) = 0.5_dp * ( x_edge(ii) + x_edge(ii+1) )
end do
do jj = 1, J
  y_center(jj) = 0.5_dp * ( y_edge(jj) + y_edge(jj+1) )
end do

udat = 10
open(unit=udat, file='phi_xy.dat', status='replace', action='write')
! pm3d expects rows separated by blank lines; we write one row per j
! Columns: x  y  phi
do jj = 1, J
  do ii = 1, I
    write(udat,'(F14.7,1X,F14.7,1X,ES16.8)') x_center(ii), y_center(jj), phi(ii,jj)
  end do
  write(udat,*)  ! blank line to separate rows
end do
close(udat)

! -----------------------------
! Build gnuplot script
! -----------------------------
uscr = 20
open(unit=uscr, file='plot_phi.gp', status='replace', action='write')

write(uscr,'(A)') "set term qt"
write(uscr,'(A)') "set grid"
write(uscr,'(A)') "set xlabel 'x'"
write(uscr,'(A)') "set ylabel 'y'"
write(uscr,'(A)') "set cblabel 'phi'"
write(uscr,'(A)') "set colorbox vertical"
write(uscr,'(A)') "set xrange ["//trim(adjustl(to_str(x_edge(1))))//":"//trim(adjustl(to_str(x_edge(I+1))))//"]"
write(uscr,'(A)') "set yrange ["//trim(adjustl(to_str(y_edge(1))))//":"//trim(adjustl(to_str(y_edge(J+1))))//"]"
write(uscr,'(A)') "set size ratio -1"            ! equal scaling x:y

! Heatmap with pm3d
write(uscr,'(A)') "set view map"
write(uscr,'(A)') "set pm3d at b"
write(uscr,'(A)') "set palette defined (0 '#000004', 0.25 '#3b0f70', 0.5 '#8c2981', 0.75 '#de4968', 1 '#fcfdbf')"
write(uscr,'(A)') "splot 'phi_xy.dat' using 1:2:3 with pm3d notitle"

! Optional: contours overlay (uncomment next lines if desired)
write(uscr,'(A)') "# set contour base"
write(uscr,'(A)') "# set cntrparam levels 10"
write(uscr,'(A)') "# replot 'phi_xy.dat' using 1:2:3 with lines lc rgb 'black' notitle"

! Optional: 3D surface view (uncomment next lines if desired)
write(uscr,'(A)') "# unset view; set view 60,30"
write(uscr,'(A)') "# unset pm3d"
write(uscr,'(A)') "# splot 'phi_xy.dat' using 1:2:3 with lines lw 1 notitle"

close(uscr)

call execute_command_line("gnuplot -persist plot_phi.gp")

contains
  function to_str(x) result(s)
    real(dp), intent(in) :: x
    character(len=32) :: s
    write(s,'(F14.7)') x
  end function to_str

end program sn_XY