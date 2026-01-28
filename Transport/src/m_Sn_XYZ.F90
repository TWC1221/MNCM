program sn_XYZ
  use m_constants
  use m_quadrature
  use m_outpVTK
  implicit none
  type(t_sn_quadrature) :: sn_Quad

  integer :: N, I, J, K, nn, ii, jj, kk, Nang

  real(8), allocatable :: X, Y, Z, x_edge(:), y_edge(:), z_edge(:), dx(:), dy(:), dz(:)
  
  real(8), allocatable :: mu(:), eta(:), zeta(:), w(:)

  real(8), allocatable :: sigs(:,:,:), sigt(:,:,:), nu_sigf(:,:,:)

  real(8), allocatable :: psi(:,:,:,:), psi_i(:,:,:,:), psi_j(:,:,:,:), psi_k(:,:,:,:)

  real(8), allocatable :: S_ext(:,:,:,:), q(:,:,:,:), phi(:,:,:), phi_prime(:,:,:), denom
  
  integer :: istart, iend, istep, jstart, jend, jstep, kstart, kend, kstep, ifor, jfor, kfor, irev, jrev, krev

  real(8) :: k_eff, k_eff_prime

  N = 16
  I = 40
  J = 40
  K = 40
  X = 1.0
  Y = 1.0
  Z = 1.0

  allocate(x_edge(I+1), Y_edge(J+1), Z_edge(K+1), dx(I), dy(J), dz(K))
  x_edge = [((X/I) * ii, ii = 0,I)]
  y_edge = [((Y/J) * jj, jj = 0,J)]
  z_edge = [((Z/K) * kk, kk = 0,K)]

  dx = [(x_edge(ii+1)-x_edge(ii),ii=1,I)]
  dy = [(y_edge(jj+1)-y_edge(jj),jj=1,J)]
  dz = [(z_edge(kk+1)-z_edge(kk),kk=1,K)]

  call Get3DAngleQuadrature(sn_Quad, N, flag_adjoint = .false.)
  Nang = sn_Quad%NoAngles 
  mu = sn_Quad%Angles(1:Nang,1)
  eta = sn_Quad%Angles(1:Nang,2)
  zeta = sn_Quad%Angles(1:Nang,3)
  w = sn_Quad%w(1:Nang)

  allocate(sigs(I,J,K), sigt(I,J,K), nu_sigf(I,J,K))
  sigs = 0.9
  sigt = 1.0
  nu_sigf = 0.05

  allocate(psi(Nang,I,J,K),psi_i(Nang,I+1,J,K),psi_j(Nang,I,J+1,K), psi_k(Nang,I,J,K+1), q(Nang,I,J,K), S_ext(Nang,I,J,K), phi(I,J,K), phi_prime(I,J,K))
  q = 0.0
  S_ext = 0.0 
  phi = 0.0
  phi_prime = 1.0
  k_eff = 0.0
  k_eff_prime = 1.0

  do
    do nn = 1,Nang
      do kk = 1,K
        do jj = 1,J
          do ii = 1,I
            q(nn,ii,jj,kk) = (sigs(ii,jj,kk) + nu_sigf(ii,jj,kk)/k_eff_prime) * phi_prime(ii,jj,kk) + S_ext(nn,ii,jj,kk)
          end do
        end do
      end do
    end do

    do nn = 1,Nang

      istep  = merge(-1, 1, mu(nn) <= 0)
      istart = merge(I, 1, mu(nn)  <= 0)
      iend   = merge(1, I, mu(nn)  <= 0)

      jstep  = merge(-1, 1, eta(nn) <= 0)
      jstart = merge(J, 1, eta(nn)  <= 0)
      jend   = merge(1, J, eta(nn)  <= 0)

      kstep  = merge(-1, 1, zeta(nn) <= 0)
      kstart = merge(K, 1, zeta(nn)  <= 0)
      kend   = merge(1, K, zeta(nn)  <= 0)

      ifor = merge(1, 0, mu(nn) >=0)
      jfor = merge(1, 0, eta(nn) >=0)
      kfor = merge(1, 0, zeta(nn) >=0)
      irev = merge(1, 0, mu(nn) <=0)
      jrev = merge(1, 0, eta(nn) <=0)
      krev = merge(1, 0, zeta(nn) <=0)

      psi_i(nn, merge(1, I+1, mu(nn)  >= 0), :, :) = 0.0
      psi_j(nn, :, merge(1, J+1, eta(nn) >= 0), :) = 0.0
      psi_k(nn, :, :, merge(1, K+1, zeta(nn) >= 0)) = 0.0

      do kk = kstart,kend,kstep
        do jj = jstart,jend,jstep
          do ii = istart,iend,istep
            denom = sigt(ii,jj,kk) + (2.0*abs(mu(nn))/dx(ii)) + (2.0*abs(eta(nn))/dy(jj)) + (2.0*abs(zeta(nn))/dz(kk))
            psi(nn,ii,jj,kk) = ((2.0*abs(mu(nn))/dx(ii))*psi_i(nn,ii+irev,jj,kk) + (2.0*abs(eta(nn))/dy(jj))*psi_j(nn,ii,jj+jrev,kk) + (2.0*abs(zeta(nn))/dz(kk))*psi_k(nn,ii,jj,kk+krev) + q(nn,ii,jj,kk))/denom

            psi_i(nn,ii+ifor,jj,kk) = 2.0 * psi(nn,ii,jj,kk) - psi_i(nn,ii+irev,jj,kk)
            psi_j(nn,ii,jj+jfor,kk) = 2.0 * psi(nn,ii,jj,kk) - psi_j(nn,ii,jj+jrev,kk) 
            psi_k(nn,ii,jj,kk+kfor) = 2.0 * psi(nn,ii,jj,kk) - psi_k(nn,ii,jj,kk+krev) 

            phi(ii,jj,kk) = phi(ii,jj,kk) + 0.125*w(nn)*psi(nn,ii,jj,kk)
          end do
        end do           
      end do
    end do
    
    !if ( maxval(abs(phi - phi_prime)) < 1d-8 ) exit !FIXED SOURCE
    if ( maxval(abs(phi - phi_prime)) < 1d-8 .and. abs(k_eff - k_eff_prime) < 1d-10 ) exit

    k_eff = k_eff_prime * sum(phi * nu_sigf)/sum(phi_prime * nu_sigf)
    k_eff_prime = k_eff
    print*,k_eff


    !phi_prime = phi !FIXED SOURCE
    phi_prime = phi / sum(phi)
    phi = 0.0

  end do

  call outpVTK_xyz_transport('SN_phi', phi, I, J, K, dx, dy, dz)

  call outpVTK_xyz_vector('SN_psi_field', psi, I, J, K, dx, dy, dz, Nang, mu, eta, zeta)

  call outpVTK_xyz_vector_allangles('SN_psi_decomposition', psi, I, J, K, dx, dy, dz, Nang, mu, eta, zeta)

  print*,maxval(phi),minval(phi), k_eff

end program sn_XYZ