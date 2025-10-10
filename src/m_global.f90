module m_global

    implicit none
    private 
    public :: t_global

    type  :: t_global
        integer                             :: NoDimension, Coordsys, NoGroups, NoMaterials
        integer                             :: Preconditioner, SN, Anisotropy, MaxIter, NoCores, FE_SE_Order, NURB_Order, IntegOrder, OutIntegOrder, Refinement
        character(20)                       :: Method, Source, Discretisation
        character(5)                        :: problem_disc
        logical                             :: Adjoint
        logical :: OutputVTK
        logical :: CalcCondition
        logical :: TerminalOutput
        logical :: OutputData
        logical :: OutMaterial
        logical :: CalcDiffCoeff
        logical :: OutFluxInt
        logical :: Outkeff
        logical :: RUN_ANALYSIS
        real(8)                             :: Tolerance, keff, flux_ratio
        real(8), dimension(:), allocatable  :: cond_number, lambda_min, lambda_max
        character(128)                      :: MeshPath, MatPath, OutVTKPath, OutDataPath
        
        contains
            procedure :: Print => print_globals
            
    end type t_global

    contains

        subroutine print_globals(this)
            class(t_global)         :: this

            write(*,'(A)') '|------------------------------------------|'
            if (this%Coordsys == 0) then
                write(*,'(A)') '| Coordinate system:  Cartesian            |'
            elseif (this%Coordsys == 1) then
                write(*,'(A)') '| Coordinate system:  Cylindrical          |'
            elseif (this%Coordsys == 2) then
                write(*,'(A)') '| Coordinate system:  Spherical            |'
            end if
            
            write(*,'(A,I1,A21)') '| Dimenension:        ', this%NoDimension, '|'
            write(*,'(A,A5,A17)') '| Problem Type:       ', this%problem_disc, '|'
            write(*,'(A,A15,A7)') '| Method:             ', this%Method, '|'

            if (this%method == 'Diffusion') then
                write(*,'(A,I2,A20)') '| Preconditioner:     ', this%Preconditioner, '|'
            end if

            write(*,'(A,A10,A12)') '| Source:             ', this%Source, '|'
            if (this%Adjoint) then
                write(*,'(A)') '| Adjoint:            True                 |'
            else
                write(*,'(A)') '| Adjoint:            False                |'
            end if
            write(*,'(A,I2,A20)') '| Cores:              ', this%NoCores, '|'
            write(*,'(A,I1,A21)') '| Groups:             ', this%NoGroups, '|'
            write(*,'(A,I1,A21)') '| Materials:          ', this%NoMaterials, '|'
            write(*,'(A,I1,A21)') '| FE/SE order:        ', this%FE_SE_Order, '|'

            if (this%IntegOrder >= 10) then
                write(*,'(A,I2,A21)') '| Integration order:  ', this%IntegOrder, '|'
            else
                write(*,'(A,I1,A21)') '| Integration order:  ', this%IntegOrder, '|'
            end if
            if (this%Method == 'Transport' .or. this%Method == 'FEM_Transport') then
                write(*,'(A,I2,A21)') '| SN quadrature:     ', this%SN, '|'
                write(*,'(A,I2,A21)') '| Anisotropy order:  ', this%Anisotropy, '|'
            end if
            write(*,'(A)') '|------------------------------------------|'

        end subroutine print_globals

end module m_global