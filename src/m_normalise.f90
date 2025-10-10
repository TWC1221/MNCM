module m_Normalise

    use m_constants
    use m_Material
    use m_Element1D
    use m_Elements
    use m_Element3D
    use m_SphElement1D
    use m_CylElement1D
    use m_ElementsCyl

    implicit none

    contains

        subroutine Normalise1D(ScalarFlux, FissionCX, keff)
            real(dp), dimension(:,:)  ::  ScalarFlux
            real(dp), dimension(:,:)  ::  FissionCX
            real(dp)                  ::  keff, Norm

            norm = sum(ScalarFlux * FissionCX) / keff
            ScalarFlux = ScalarFlux / Norm

        end subroutine Normalise1D


        subroutine Normalise2DXY(ScalarFlux, FissionCX, dx, dy, keff)

            real(dp), dimension(:,:,:)  ::  ScalarFlux
            real(dp), dimension(:,:,:)  ::  FissionCX

            real(dp)                    ::  dx, dy
            real(dp)                    ::  keff, NormalisationFactor
            integer                     ::  gindex

            NormalisationFactor = 0.0_dp
            do gindex = 1, size(ScalarFlux, dim=1)
                NormalisationFactor = NormalisationFactor + sum(ScalarFlux(gindex,:,:) * FissionCX(gindex,:,:))
            end do

            NormalisationFactor = NormalisationFactor * dx * dy / keff
            ScalarFlux = ScalarFlux / NormalisationFactor

        end subroutine Normalise2DXY

        subroutine NormaliseFEM1DTransport(NoGroups,Material, Elements, keff, Adjoint)

            type(t_material), dimension(:)   ::  Material
            type(t_Element1D), dimension(:)     ::  Elements
            real(dp)                            ::  keff    
            real(dp)                            ::  ScalarFluxAvg
            real(dp)                            ::  NormalisationFactor
            logical                             ::  Adjoint
            integer                             ::  gindx, i, NoGroups

            NormalisationFactor = 0.0_dp
            if (Adjoint) then
                do i = 1, size(Elements)
                    if (maxval(Elements(i)%ScalarFluxIter) > NormalisationFactor) then
                        NormalisationFactor = maxval(Elements(i)%ScalarFluxIter)
                    end if
                end do

                do i = 1, size(Elements)
                    Elements(i)%ScalarFluxIter = Elements(i)%ScalarFluxIter / NormalisationFactor
                end do
            
            else
                NormalisationFactor = 0.0_dp
                do i = 1, size(Elements)
                    do gindx = 1, NoGroups
                        ! Average the scalar flux over the element
                        ScalarFluxAvg = sum(Elements(i)%ScalarFluxIter(gindx,:)) / real(Elements(i)%PolyOrder+1, dp)
                        NormalisationFactor = NormalisationFactor + ScalarFluxAvg * Material(Elements(i)%MatID)%nusigmaF(gindx) * Elements(i)%Length
                    end do
                end do

                NormalisationFactor = NormalisationFactor / keff
                do i = 1, size(Elements)
                    Elements(i)%ScalarFluxIter = Elements(i)%ScalarFluxIter / NormalisationFactor
                end do
            end if

        end subroutine NormaliseFEM1DTransport

        subroutine NormaliseFEM2DTransport(NoGroups,Material, Elements, keff, Adjoint)

            type(t_material), dimension(:)   ::  Material
            type(t_Element2D), dimension(:)     ::  Elements
            real(dp)                            ::  keff    
            real(dp)                            ::  ScalarFluxAvg
            real(dp)                            ::  NormalisationFactor
            logical                             ::  Adjoint
            integer                             ::  gindx, i, NoGroups

            NormalisationFactor = 0.0_dp
            if (Adjoint) then
                do i = 1, size(Elements)
                    if (maxval(Elements(i)%ScalarFluxIter) > NormalisationFactor) then
                        NormalisationFactor = maxval(Elements(i)%ScalarFluxIter)
                    end if
                end do

                do i = 1, size(Elements)
                    Elements(i)%ScalarFluxIter = Elements(i)%ScalarFluxIter / NormalisationFactor
                end do
            
            else
                NormalisationFactor = 0.0_dp
                do i = 1, size(Elements)
                    do gindx = 1, NoGroups
                        ! Average the scalar flux over the element
                        ScalarFluxAvg = sum(Elements(i)%ScalarFluxIter(gindx,:)) / real(Elements(i)%NoNodes, dp)
                        NormalisationFactor = NormalisationFactor + ScalarFluxAvg * Material(Elements(i)%MatID)%nusigmaF(gindx) * Elements(i)%Area
                    end do
                end do

                NormalisationFactor = NormalisationFactor / keff
                do i = 1, size(Elements)
                    Elements(i)%ScalarFluxIter = Elements(i)%ScalarFluxIter / NormalisationFactor
                end do
            end if

        end subroutine NormaliseFEM2DTransport

        subroutine NormaliseFEMQuad(Material, Elements, ScalarFlux, keff, Adjoint)

            type(t_material), dimension(:)   ::  Material
            type(t_Element2D), dimension(:)     ::  Elements
            real(dp), dimension(:,:)            ::  ScalarFlux
            real(dp)                            ::  keff    
            real(dp)                            ::  ScalarFluxAvg
            real(dp)                            ::  NormalisationFactor
            logical                             ::  Adjoint
            integer                             ::  gindx, eindx

            if (Adjoint) then
                ScalarFlux = ScalarFlux / maxval(ScalarFlux)
            else
                NormalisationFactor = 0.0_dp
                do eindx = 1, size(Elements)
                    do gindx = 1, size(ScalarFlux,dim=1)
                        ! Average the scalar flux over the element
                        ScalarFluxAvg = sum(ScalarFlux(gindx,Elements(eindx)%EN_List+1)) / real(Elements(eindx)%NoNodes, dp)
                        NormalisationFactor = NormalisationFactor + ScalarFluxAvg * Material(Elements(eindx)%MatID)%NuSigmaF(gindx) * Elements(eindx)%Area
                    
                    end do
                end do

                NormalisationFactor = NormalisationFactor / keff
                write(*,*) 'Normalisation Factor = ', NormalisationFactor
                ScalarFlux = ScalarFlux / NormalisationFactor
            end if

        end subroutine NormaliseFEMQuad


        subroutine NormaliseFEM3DTransport(NoGroups,Material, Elements, keff, Adjoint)

            type(t_material), dimension(:)   ::  Material
            type(t_Element3D), dimension(:)     ::  Elements
            real(dp)                            ::  keff    
            real(dp)                            ::  ScalarFluxAvg
            real(dp)                            ::  NormalisationFactor
            logical                             ::  Adjoint
            integer                             ::  gindx, NoGroups, i


            if (Adjoint) then
                
            else
                NormalisationFactor = 0.0_dp
                do i = 1, size(Elements)
                    do gindx = 1, NoGroups
                        ! Average the scalar flux over the element
                        ScalarFluxAvg = sum(Elements(i)%ScalarFluxIter(gindx,:)) / real(Elements(i)%NoNodes, dp)
                        NormalisationFactor = NormalisationFactor + ScalarFluxAvg * Material(Elements(i)%MatID)%nusigmaF(gindx) * Elements(i)%Volume
                    end do
                end do

                NormalisationFactor = NormalisationFactor / keff
                do i = 1, size(Elements)
                    Elements(i)%ScalarFluxIter = Elements(i)%ScalarFluxIter / NormalisationFactor
                end do
            end if

        end subroutine NormaliseFEM3DTransport

        subroutine NormaliseFEMTri(Material, Elements, ScalarFlux, keff)

            type(t_material), dimension(:)   ::  Material
            type(t_Element2D), dimension(:)     ::  Elements
            real(dp), dimension(:,:)            ::  ScalarFlux
            real(dp)                            ::  keff    
            real(dp)                            ::  ScalarFluxAvg
            real(dp)                            ::  NormalisationFactor
            integer                             ::  gindx, eindx

            NormalisationFactor = 0.0_dp
            do eindx = 1, size(Elements)
                do gindx = 1, size(ScalarFlux,dim=1)
                    ! Average the scalar flux over the element
                    ScalarFluxAvg = sum(ScalarFlux(gindx,Elements(eindx)%EN_List+1)) / real(Elements(eindx)%NoNodes, dp)
                    NormalisationFactor = NormalisationFactor + ScalarFluxAvg * Material(Elements(eindx)%MatID)%NuSigmaF(gindx) * Elements(eindx)%Area
                end do
            end do

            NormalisationFactor = NormalisationFactor / keff
            ScalarFlux = ScalarFlux / NormalisationFactor

        end subroutine NormaliseFEMTri


        subroutine Normalise1DSphFE(Material, Elements, keff, Adjoint)

            type(t_material), dimension(:)   ::  Material
            type(t_SphElement1D), dimension(:)  ::  Elements
            real(dp)                            ::  keff
            real(dp)                            ::  ScalarFluxAvg
            real(dp)                            ::  NormalisationFactor
            logical                             ::  Adjoint
            integer                             ::  gindx, i

            if (Adjoint) then

            else
                NormalisationFactor = 0.0_dp
                do i = 1, size(Elements)
                    do gindx = 1, size(Elements(i)%ScalarFluxIter, dim=1)
                        ! Average the scalar flux over the element
                        ScalarFluxAvg = sum(Elements(i)%ScalarFluxIter(gindx,:)) / real(Elements(i)%NoNodes, dp)
                        NormalisationFactor = NormalisationFactor + ScalarFluxAvg * Material(Elements(i)%MatID)%nusigmaF(gindx) * Elements(i)%Delta_r * 4.0_dp * pi * Elements(i)%r_bar**2 / keff
                    end do
                end do

                
                write(*,*) 'Normalisation Factor = ', NormalisationFactor
                do i = 1, size(Elements)
                    Elements(i)%ScalarFluxIter = Elements(i)%ScalarFluxIter / NormalisationFactor
                end do
            end if
        end subroutine Normalise1DSphFE

        subroutine Normalise1DCylFE(Material, Elements, keff, Adjoint)

            type(t_material), dimension(:)   ::  Material
            type(t_CylElement1D), dimension(:)  ::  Elements
            real(dp)                            ::  keff
            real(dp)                            ::  ScalarFluxAvg
            real(dp)                            ::  NormalisationFactor
            logical                             ::  Adjoint
            integer                             ::  gindx, i

            if (Adjoint) then

            else
                NormalisationFactor = 0.0_dp
                do i = 1, size(Elements)
                    do gindx = 1, size(Elements(i)%ScalarFluxIter, dim=1)
                        ! Average the scalar flux over the element
                        ScalarFluxAvg = sum(Elements(i)%ScalarFluxIter(gindx,:)) / real(Elements(i)%NoNodes, dp)
                        NormalisationFactor = NormalisationFactor + ScalarFluxAvg * Material(Elements(i)%MatID)%nusigmaF(gindx) * Elements(i)%Delta_rho * 2.0_dp * pi * Elements(i)%rho_bar / keff
                    end do
                end do

                
                write(*,*) 'Normalisation Factor = ', NormalisationFactor
                do i = 1, size(Elements)
                    Elements(i)%ScalarFluxIter = Elements(i)%ScalarFluxIter / NormalisationFactor
                end do
            end if
        end subroutine Normalise1DCylFE

        subroutine Normalise2DCylFE(Material, Elements, keff, Adjoint)

            type(t_material), dimension(:)   ::  Material
            type(t_Element2DCyl), dimension(:)  ::  Elements
            real(dp)                            ::  keff
            real(dp)                            ::  ScalarFluxAvg
            real(dp)                            ::  NormalisationFactor
            logical                             ::  Adjoint
            integer                             ::  gindx, i

            if (Adjoint) then

            else
                NormalisationFactor = 0.0_dp
                do i = 1, size(Elements)
                    do gindx = 1, size(Elements(i)%ScalarFluxIter, dim=1)
                        ! Average the scalar flux over the element
                        ScalarFluxAvg = sum(Elements(i)%ScalarFluxIter(gindx,:)) / real(Elements(i)%NoNodes, dp)
                        NormalisationFactor = NormalisationFactor + ScalarFluxAvg * Material(Elements(i)%MatID)%nusigmaF(gindx) * Elements(i)%Area * 2.0_dp * pi * Elements(i)%r / keff
                    end do
                end do

                
                write(*,*) 'Normalisation Factor = ', NormalisationFactor
                do i = 1, size(Elements)
                    Elements(i)%ScalarFluxIter = Elements(i)%ScalarFluxIter / NormalisationFactor
                end do
            end if
        end subroutine Normalise2DCylFE


end module m_Normalise
