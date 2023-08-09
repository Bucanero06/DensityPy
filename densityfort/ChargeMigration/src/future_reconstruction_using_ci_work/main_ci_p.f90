! {{{ Detailed description

!> \mainpage Program <ProgramName> <Insert here what the program does>
!!
!! Synopsis:
!! ---------
!!
!!     <Program Name> <mandatory run-time parameters (RTP)> [<optional RTP>]
!!
!! ___
!! Description:
!! ------------
!!
!! Input parameters:      {#Input_Parameters}
!! =================
!! [...Input](@ref ...Input) as specified in the command line.
!!
!> \file
!!
!!
! }}}
program ChargeMigration

    use, intrinsic :: ISO_FORTRAN_ENV
    use ModuleErrorHandling
    use ModuleSystemUtils
    use ModuleString
    use ModuleIO
    use ModuleConstants
    use ModuleDiagonalize
    use ModulePulses_3D
    !
    !.. Local modules
    use ModuleRTP_p
    use Module_CD_IO

    implicit none

    external zgemm, zgemv, system !Explicit declaration of the EXTERNAL attribute is required.
    !.. Run-time parameters
    !..
    character(len = :), allocatable :: input_directory
    character(len = :), allocatable :: output_directory
    character(len = :), allocatable :: molecular_geometry_file
    integer :: n_times
    real(kind(1d0)) :: t_min
    real(kind(1d0)) :: t_max
    character(len = :), allocatable :: Ext_field_file
    character(len = :), allocatable :: Weight_File
    logical :: Verbous
    logical :: read_precomputed_weights_flag
    logical :: save_charge_migration_flag
    integer, allocatable :: ivorb(:)
    integer :: counted_number_of_orbitals
    real(kind(1d0)) :: volume
    real(kind(1d0)) :: dephasing_factor
    real(kind(1d0)) :: relaxation_factor
    real(kind(1d0)) :: bath_temperature
    integer, parameter :: GS_IDX = 1

    !.. Local parameters
    integer :: nStates, number_of_orbitals, i
    real(kind(1d0)), allocatable :: Evec(:)
    !.. Dmat(i,j,alpha) = $ \langle \varphi_i | \hat{\mu}_\alpha$ | \varphi_j \rangle $
    real(kind(1d0)), allocatable :: Dmat(:, :, :)
    !.. TDM(j,i,B,A) = $ \langle A | \hat{a}_i^\dagger \hat{a}_j | B \rangle $
    real(kind(1d0)), allocatable :: TDM(:, :, :, :)
    !.. Expectation Value of the Dipole Moment (Mu)
    complex(kind(1d0)), allocatable :: zMuEV(:, :), zDmat_t(:, :, :)

    integer :: npts, nAtoms, nxpoints
    real(kind(1d0)), allocatable :: gridv(:, :), AtCoord(:, :) ! 3 x npts
    real(kind(1d0)), allocatable :: OrbTab(:, :) ! 3 x npts
    !    real(kind(1d0)) :: volume

    !.. Statistical Density Matrix
    complex(kind(1d0)), allocatable :: zStatRho(:, :), zStatRho0(:, :)
    real   (kind(1d0)) :: t, dt
    real   (kind(1d0)), allocatable :: ChDen(:), WEIGHTV(:, :)
    real   (kind(1d0)), allocatable :: AtomicChargeVec(:), AtomicChargeEvolution(:, :)

    character(len = 30) :: istrn
    character(len = 16), allocatable :: atom_names(:)

    integer :: iPts, it, iPol, iAtom
    integer :: uid_AtomicCharge, uid_dipole

    complex(kind(1d0)), allocatable :: Liouvillian0(:, :), L0_LEvec(:, :), L0_REvec(:, :), L0_Eval(:)

    !.. Pulse parameters
    integer :: iSim, N_Simulations, uid
    character(len = 64), pointer :: Simulation_Tagv(:)
    character(len = 1000) :: strn
    type(pulse_train), pointer :: train(:)

    !..
    !    B_{ij}^{\alpha} = \int d^3r \phi_i(\vec{r})\phi_j(\vec{r}) w_\alpha(\vec{r})
    !                    = BeckeMatrix(i,j,alpha)
    !..
    real(kind(1d0)), allocatable :: BeckeMatrix   (:, :, :)
    real(kind(1d0)), allocatable :: OrbitalDensity(:, :)
    real(kind(1d0)), allocatable :: Radius_BS(:)
    complex(kind(1d0)), external :: zdotu

    !..Test
    real(kind(1d0)), allocatable :: ChargeTotal1(:), ChargeTotal2(:), Becke_new(:, :, :, :), Dipole_new(:, :), Dipole_new_(:), Becke_new1(:, :, :, :), data(:, :, :, :)
    real(kind(1d0)) :: data1, data2, data3, Orbital_overlap_self, Orbital_overlap_other, Computed_volume, MASK
    real   (kind(1d0)), allocatable :: AtomicChargeVec_new(:, :), AtomicChargeEvolution_new(:, :, :)
    integer :: iOrb, jOrb, number_of_orbitalss, i_excitation, i_epsilon
    real(kind(1d0)), allocatable :: QchargeVec_new    (:, :), R_el(:, :)
    character(len = 64) :: excitation_string

    call GetRunTimeParameters(input_directory, output_directory, molecular_geometry_file, n_times, t_min, t_max, Ext_field_file, Verbous, Weight_File, &
            read_precomputed_weights_flag, save_charge_migration_flag, ivorb, counted_number_of_orbitals, volume, dephasing_factor, relaxation_factor, &
            bath_temperature, i_excitation, i_epsilon)

    call Set_CD_IO_Verbous(Verbous)

    open(newunit = uid, &
            file = Ext_field_file, &
            form = "formatted", &
            status = "old")
    call Parse_Simulation_File(uid, OUTPUT_UNIT, N_Simulations, Simulation_Tagv, train)
    close(uid)
    write(*, *)"N_Simulations=", N_Simulations

    !    !.. Print the pulse
    !    !..
    dt = (t_max - t_min) / dble(n_times - 1)
    if (i_excitation ==-1 .and. i_epsilon==-1) then
        call system("mkdir -p " // trim(output_directory))
        call system("mkdir -p " // output_directory // "/Dipole")
        call system("mkdir -p " // output_directory // "/AtomicCharge")
        call system("mkdir -p " // output_directory // "/ChargeDensity")
        call system("mkdir -p " // output_directory // "/Pulses")
        call system("cp " // Ext_field_file // " " // output_directory // "/")
        do i = 1, N_Simulations
            write(*, *) "Writing pulse " // trim(Simulation_Tagv(i))
            strn = trim(output_directory) // "/Pulses/pulse" // trim(Simulation_Tagv(i))
            call train(i)%Write(strn, t_min, t_max, dt)
            strn = trim(output_directory) // "/Pulses/FTpulse" // trim(Simulation_Tagv(i))
            call train(i)%WriteFTA(strn)
        enddo
    end if
!    stop!stop here to look at pulses
    !
    !
    call LoadEnergies(input_directory // "/ROOT_ENERGIES", nStates, Evec)

    !.. At the moment, we assume that the TDM are defined as
    !   $\rho^{JI}_{nm} = \langle J | a_n^\dagger a_m | I \rangle $
    !..
    call LoadTDMs    (input_directory // "/DENSITY_MATRIX", input_directory // "/TRANSITION_DENSITY_MATRIX", &
            nStates, number_of_orbitals, TDM)

    !.. Load Dipole matrix elements $\mu_{IJ}$
    call LoadDipoleME(Dmat, input_directory, nStates)
    allocate(zDmat_t(nStates, nStates, 3))
    do i = 1, 3
        zDmat_t(:, :, i) = Z1 * transpose(Dmat(:, :, i))
    enddo

    if (number_of_orbitals .ne. counted_number_of_orbitals) then
        write(*, *) "Number of orbitals given by '-iOrb' do not match number of orbitals in active the active space"
        stop
    endif

    !.. Load Geometry,volume, Grid and Orbitals
    call LoadGeometry(nAtoms, AtCoord, molecular_geometry_file, atom_names)
    call LoadGrid(input_directory // "/gridcoord", npts, gridv)! $G=\{\vec{r}_i, i=1,\ldots,N_{G}\}$
    call Computevolume(nPts, Computed_volume, gridv)
    call LoadOrbitals(input_directory, number_of_orbitals, ivOrb, npts, OrbTab)! $\varphi_n(\vec{r}_i)$
    !.. Compute Liouvillian and Diagonalizes it
    call ComputeLiouvillian01(Evec, Dmat, Liouvillian0, bath_temperature, dephasing_factor, relaxation_factor)
    call DiagonalizeLindblad0(Liouvillian0, L0_Eval, L0_LEvec, L0_REvec)

    !..Computes Real and Imaginary Components of the Dipole for each Molecular Excitation and the Ground State with no Excitations
    if (i_excitation ==0 .and. i_epsilon==0) then !Non_Parallel
        Single_Excitation_Loop : do i_excitation = 1, nStates
            write(excitation_string, '(a,i0.2)') 'mkdir -p ' // output_directory // '/Dipole/Dipole_State_', i_excitation
            call system(excitation_string)
            write(excitation_string, '(a,i0.2)') output_directory // '/Dipole/Dipole_State_', i_excitation
            call EMULATE(excitation_string, N_Simulations, zMuEV, zDmat_t, L0_Eval, L0_LEvec, L0_REvec, GS_IDX, i_excitation, n_times, t_min, dt, nStates)
        end do Single_Excitation_Loop
    else !Run Parallel
        write(excitation_string, '(a,i0.2)') 'mkdir -p ' // output_directory // '/Dipole/Dipole_State_', i_excitation
        call system(excitation_string)
        write(excitation_string, '(a,i0.2)') output_directory // '/Dipole/Dipole_State_', i_excitation
        call EMULATE_Parallel(excitation_string, N_Simulations, zMuEV, zDmat_t, L0_Eval, L0_LEvec, L0_REvec, GS_IDX, i_excitation, n_times, t_min, dt, nStates, i_epsilon)
    end if
    stop
    !
    !    if(.not. save_charge_migration_flag)then
    !        call Write_Summary(output_directory // "/Simulation_Summary", nPts, nAtoms, volume, Computed_volume, n_times, t_min, t_max, atom_names, Radius_BS, number_of_orbitals, OrbTab)
    !    end if

    stop


contains

    !###############################>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine EMULATE(FileName, N_Simulations, zMuEV, zDmat_t, L0_Eval, L0_LEvec, L0_REvec, GS_IDX, STATE_TO_EXCITE, n_times, t_min, dt, nStates)
        character(len = *), intent(in) :: FileName
        real   (kind(1d0)), intent(in) :: t_min, dt
        complex(kind(1d0)), intent(in) :: zDmat_t(:, :, :), L0_LEvec(:, :), L0_REvec(:, :), L0_Eval(:)
        !
        integer, intent(in) :: N_Simulations, GS_IDX, STATE_TO_EXCITE, n_times, nStates
        !
        complex(kind(1d0)), allocatable :: zStatRho(:, :), zMuEV(:, :)
        complex   (kind(1d0)) :: EPSILON
        character(len = 64) :: PART, excitation_string, file
        integer :: i_epsilon
        !
        real   (kind(1d0)) :: Mean, Variance, StdDev   ! results
        integer :: n, i                   ! actual array size
        !
        !        !!!!!!!!!!!!!!!!!!!!!!!
        !        real   (kind(1d0)) :: Mean, Variance, StdDev   ! results
        !        Mean = 0.0                           ! compute mean
        !        do i = 1, N_Simulations

        !            Mean = Mean + i     !Simulation_tagv(i)
        !        end do
        !        Mean = Mean / N_Simulations
        !        Variance = 0.0                       ! compute variance
        !        do i = 1, N_Simulations
        !            Variance = Variance + (i - Mean)**2
        !        end do
        !        Variance = Variance / (N_Simulations - 1)
        !        StdDev = SQRT(Variance)            ! compute standard deviation
        !
        !
        !                    MASK = (iSim - 1) / StdDev
        !                    zMuEV(iPol, it) = zMuEV(iPol, it) * MASK
        !        !!!!!!!!!!!!!!!!!!!!!!!!
        !
        EPSILON = 10e-5
        !
        EPSILON_loop : do i_epsilon = 1, 2
            !
            if (i_epsilon == 1) then
                PART = "Real"
            elseif (i_epsilon == 2) then
                PART = "Imaginary"
                EPSILON = EPSILON * Zi
            end if
            !
            if ((STATE_TO_EXCITE == GS_IDX) .and. (i_epsilon == 2)) cycle
            !
            if(allocated(zStatRho))deallocate(zStatRho)
            allocate(zStatRho(nStates, nStates))
            if(allocated(zMuEV))deallocate(zMuEV)
            allocate(zMuEV(3, n_times))
            zMuEV = Z0
            !
            write(*, *)
            write(*, '(a,i0.2,a)') "Starting Sim Loop for Molecular Excitation ", STATE_TO_EXCITE, " " // trim(PART)
            !
            Sim_loop : do iSim = 1, N_Simulations
                !
                zStatRho = Z0
                zStatRho(GS_IDX, GS_IDX) = 1.d0
                !
                write(*, *) "Time_Delay = ", iSim, "/", N_Simulations!!!print to screen
                !.. Time cycle
                time_loop : do it = 1, n_times
                    !
                    t = t_min + dt * dble(it - 1)
                    !
                    !..Excites State at t = 0
                    if (t==0.d0) then
                        if (STATE_TO_EXCITE .ne. GS_IDX) then
                            zStatRho(STATE_TO_EXCITE, STATE_TO_EXCITE) = EPSILON
                        end if
                    end if
                    !
                    !.. Computes the excitation density wrt ground state
                    zStatRho(GS_IDX, GS_IDX) = zStatRho(GS_IDX, GS_IDX) - 1.d0
                    !
                    !.. Evaluates the expectation value of the dipole as a function of time
                    do iPol = 1, 3
                        zMuEV(iPol, it) = zdotu(nStates * nStates, zStatRho, 1, zDmat_t(1, 1, iPol), 1)
                    enddo
                    !
                    zStatRho(GS_IDX, GS_IDX) = zStatRho(GS_IDX, GS_IDX) + 1.d0
                    !
                    if(it == n_times) exit time_loop
                    !
                    call LiouvillianPropagator(L0_Eval, L0_LEvec, L0_REvec, zStatRho)
                    !
                enddo time_loop
                !
                !.. Save Dipole to file
                excitation_string = trim(FileName) // "/Dipole" // trim(Simulation_tagv(iSim)) // "_" // trim(PART)
                call Write_Dipole(excitation_string, zMuEV, n_times, t_min, dt)
                !
            end do Sim_loop
            !
        end do EPSILON_loop
        !
    end subroutine EMULATE

    subroutine EMULATE_Parallel(FileName, N_Simulations, zMuEV, zDmat_t, L0_Eval, L0_LEvec, L0_REvec, GS_IDX, STATE_TO_EXCITE, n_times, t_min, dt, nStates, i_epsilon)
        character(len = *), intent(in) :: FileName
        real   (kind(1d0)), intent(in) :: t_min, dt
        complex(kind(1d0)), intent(in) :: zDmat_t(:, :, :), L0_LEvec(:, :), L0_REvec(:, :), L0_Eval(:)
        !
        integer, intent(in) :: N_Simulations, GS_IDX, STATE_TO_EXCITE, n_times, nStates
        !
        complex(kind(1d0)), allocatable :: zStatRho(:, :), zMuEV(:, :)
        complex   (kind(1d0)) :: EPSILON
        character(len = 64) :: PART, excitation_string, file
        integer, intent(in) :: i_epsilon
        !
        real   (kind(1d0)) :: Mean, Variance, StdDev   ! results
        integer :: n, i                   ! actual array size
        !
        EPSILON = 10e-5
        !
        !
        if (i_epsilon == 1) then
            PART = "Real"
        elseif (i_epsilon == 2) then
            PART = "Imaginary"
            EPSILON = EPSILON * Zi
        end if
        !
        if(allocated(zStatRho))deallocate(zStatRho)
        allocate(zStatRho(nStates, nStates))
        if(allocated(zMuEV))deallocate(zMuEV)
        allocate(zMuEV(3, n_times))
        zMuEV = Z0
        !
        write(*, *)
        write(*, '(a,i0.2,a)') "Starting Sim Loop for Molecular Excitation ", STATE_TO_EXCITE, " " // trim(PART)
        !
        Sim_loop : do iSim = 1, N_Simulations
            !
            zStatRho = Z0
            zStatRho(GS_IDX, GS_IDX) = 1.d0
            !
            write(*, *) "Time_Delay = ", iSim, "/", N_Simulations!!!print to screen
            !.. Time cycle
            time_loop : do it = 1, n_times
                !
                t = t_min + dt * dble(it - 1)
                !
                !..Excites State at t = 0
                if (t==0.d0) then
                    if (STATE_TO_EXCITE .ne. GS_IDX) then
                        zStatRho(STATE_TO_EXCITE, STATE_TO_EXCITE) = EPSILON
                    end if
                end if
                !
                !.. Computes the excitation density wrt ground state
                zStatRho(GS_IDX, GS_IDX) = zStatRho(GS_IDX, GS_IDX) - 1.d0
                !
                !.. Evaluates the expectation value of the dipole as a function of time
                do iPol = 1, 3
                    zMuEV(iPol, it) = zdotu(nStates * nStates, zStatRho, 1, zDmat_t(1, 1, iPol), 1)
                enddo
                !
                zStatRho(GS_IDX, GS_IDX) = zStatRho(GS_IDX, GS_IDX) + 1.d0
                !
                if(it == n_times) exit time_loop
                !
                call LiouvillianPropagator(L0_Eval, L0_LEvec, L0_REvec, zStatRho)
                !
            enddo time_loop
            !
            !.. Save Dipole to file
            excitation_string = trim(FileName) // "/Dipole" // trim(Simulation_tagv(iSim)) // "_" // trim(PART)
            call Write_Dipole(excitation_string, zMuEV, n_times, t_min, dt)
            !
        end do Sim_loop
        !
        !
    end subroutine EMULATE_Parallel

    subroutine Compute_R_el(gridv, WeightV, OrbTab, R_el)
        real(kind(1d0)), allocatable, intent(in) :: OrbTab (:, :), gridv(:, :), WeightV(:, :)
        real(kind(1d0)), allocatable, intent(out) :: R_el(:, :)

        real(kind(1d0)) :: sum, m, r, sum1
        integer :: nAtoms, nPts, number_of_orbitalss, iOrb, jOrb, uid
        nAtoms = size(WeightV, 2)
        nPts = size(WeightV, 1)
        number_of_orbitalss = size(OrbTab, 2)

        !
        write(*, "(a)") "Computing Barycenters of the Atomic Charges"
        !
        allocate(R_el(3, nAtoms))

        do iPol = 1, 3
            do iAtom = 1, nAtoms
                sum = 0.d0
                sum1 = 0.d0
                do iOrb = 1, number_of_orbitals
                    do jOrb = 1, number_of_orbitals
                        do iPts = 1, nPts
                            m = WeightV(iPts, iAtom) * OrbTab(iPts, iOrb) * OrbTab(iPts, iOrb)
                            r = gridv(iPol, iPts)
                            sum = sum + ((m * r))
                            sum1 = sum1 + m
                        end do
                    end do
                end do
                !
                !
                R_el(iPol, iAtom) = sum / sum1
            end do
        end do

        open(newunit = uid, &
                file = "R_el_bc", &
                form = "formatted", &
                status = "unknown")
        do iAtom = 1, nAtoms
            write(uid, "(*(x,e24.14e3))") (R_el(iPol, iAtom), iPol = 1, 3)
        end do
        close(uid)
    end subroutine Compute_R_el
    !
    subroutine ComputeNewBeckeMatrix1(WeightV, OrbTab, Becke_new, Bary_center)
        real(kind(1d0)), intent(in) :: WeightV(:, :)
        real(kind(1d0)), intent(in) :: OrbTab (:, :)
        real(kind(1d0)), allocatable, intent(out) :: Becke_new(:, :, :, :)
        real(kind(1d0)), intent(in) :: Bary_center(:, :)

        integer :: nPts, number_of_orbitalss, nAtoms
        integer :: iPts, iOrb1, iOrb2, iAtom, iPol
        real(kind(1d0)) :: dsum

        !
        write(*, "(a)") "Computing Becke's Matrix"
        !
        nAtoms = size(WeightV, 2)
        nPts = size(WeightV, 1)
        number_of_orbitalss = size(OrbTab, 2)

        !        if(allocated(Becke_new))deallocate(Becke_new)
        allocate(Becke_new(3, number_of_orbitalss, number_of_orbitalss, nAtoms))

        Becke_new = 0.d0
        do iPol = 1, 3
            do iAtom = 1, nAtoms
                do iOrb2 = 1, number_of_orbitalss
                    do iOrb1 = 1, number_of_orbitalss
                        do iPts = 1, nPts
                            Becke_new(iPol, iOrb1, iOrb2, iAtom) = Becke_new(iPol, iOrb1, iOrb2, iAtom) + &
                                    WeightV(iPts, iAtom) * OrbTab(iPts, iOrb1) * OrbTab(iPts, iOrb2) * &
                                            Bary_center(iPol, iAtom)

                        enddo
                    enddo
                enddo
            enddo
        end do
    end subroutine ComputeNewBeckeMatrix1
    !
    subroutine ComputeNewBeckeMatrix2(WeightV, OrbTab, Becke_new, gridv)
        real(kind(1d0)), intent(in) :: WeightV(:, :)
        real(kind(1d0)), intent(in) :: OrbTab (:, :)
        real(kind(1d0)), allocatable, intent(out) :: Becke_new(:, :, :, :)
        real(kind(1d0)), intent(in) :: gridv(:, :)

        integer :: nPts, number_of_orbitalss, nAtoms
        integer :: iPts, iOrb1, iOrb2, iAtom, iPol
        real(kind(1d0)) :: dsum

        !
        write(*, "(a)") "Computing Becke's Matrix"
        !
        nAtoms = size(WeightV, 2)
        nPts = size(WeightV, 1)
        number_of_orbitalss = size(OrbTab, 2)

        if(allocated(Becke_new))deallocate(Becke_new)
        allocate(Becke_new(3, number_of_orbitalss, number_of_orbitalss, nAtoms))

        Becke_new = 0.d0
        do iPol = 1, 3
            do iAtom = 1, nAtoms
                do iOrb2 = 1, number_of_orbitalss
                    do iOrb1 = 1, number_of_orbitalss
                        do iPts = 1, nPts
                            Becke_new(iPol, iOrb1, iOrb2, iAtom) = Becke_new(iPol, iOrb1, iOrb2, iAtom) + &
                                    WeightV(iPts, iAtom) * OrbTab(iPts, iOrb1) * OrbTab(iPts, iOrb2) * &
                                            gridv(iPol, iPts)
                        enddo
                    enddo
                enddo
            enddo
        end do
    end subroutine ComputeNewBeckeMatrix2
    !
    subroutine ComputeNewAtomicCharges(OrbitalDensity, Becke_new, QchargeVec_new)
        real(kind(1d0)), intent(in) :: OrbitalDensity(:, :)
        real(kind(1d0)), intent(in) :: Becke_new(:, :, :, :)
        real(kind(1d0)), allocatable, intent(out) :: QchargeVec_new    (:, :)

        integer :: iAtom, nAtoms, number_of_orbitalss, jOrb, iOrb, nPts, iPol
        real(kind(1d0)), external :: DDOT
        real(kind(1d0)) :: dsum
        nPts = size(Becke_new, 1)
        number_of_orbitalss = size(Becke_new, 2)
        nAtoms = size(Becke_new, 4)
        allocate(QchargeVec_new(3, nAtoms))

        QchargeVec_new = 0.d0
        do iPol = 1, 3
            do iAtom = 1, nAtoms
                do jOrb = 1, number_of_orbitalss
                    do iOrb = 1, number_of_orbitalss
                        QchargeVec_new(iPol, iAtom) = QchargeVec_new(iPol, iAtom) + OrbitalDensity(iOrb, jOrb) &
                                * Becke_new(iPol, iOrb, jOrb, iAtom)
                        !                QchargeVec_new(iPol, iAtom) = ddot(number_of_orbitalss * number_of_orbitalss, OrbitalDensity, 1, Becke_new(iPol, 1, 1, iAtom), 1)
                    end do
                end do
            enddo
        end do
    end subroutine ComputeNewAtomicCharges


    subroutine ComputeOrbitalDensity(zStatRho, TDM, Amat)
        complex(kind(1d0)), intent(in) :: zStatRho(:, :)
        real   (kind(1d0)), intent(in) :: TDM(:, :, :, :)
        real   (kind(1d0)), allocatable, intent(out) :: Amat  (:, :)
        integer :: number_of_orbitals, ii_orb, ij_orb, ii_state, ij_state

        number_of_orbitals = size(TDM, 1)

        if(allocated(Amat))then
            if(abs(size(Amat, 1) - number_of_orbitals) + abs(size(Amat, 2) - number_of_orbitals)>0)then
                deallocate(Amat)
                allocate(Amat(number_of_orbitals, number_of_orbitals))
            endif
        else
            allocate(Amat(number_of_orbitals, number_of_orbitals))
        endif
        !
        !.. Build the expansion Matrix
        !   $A_{nm}(t) = \sum_{IJ} P_{IJ}(t) \rho^{JI}_{nm}$
        !   where $P_{IJ}(t)$ is the solution of a suitable Master equation, starting from $P_{IJ}(0)$
        !   Here, we will neglect relaxation and decoherence and hence $P_{IJ}(t) = P_{IJ}(0) e^{-i\omega_{IJ}t}$
        !   where $\omega_{IJ}=E_I-E_J$
        Amat = 0.d0
        do ii_orb = 1, number_of_orbitals
            do ij_orb = 1, number_of_orbitals
                !
                do ii_state = 1, nStates
                    do ij_state = 1, nStates
                        Amat(ii_orb, ij_orb) = Amat(ii_orb, ij_orb) + &
                                dble(zStatRho(ii_state, ij_state) * TDM(ii_orb, ij_orb, ij_state, ii_state))
                    enddo
                enddo
                !
            enddo
        enddo

    end subroutine ComputeOrbitalDensity


    subroutine TabulateChargeDensity(Amat, OrbTab, ChDen)
        real(kind(1d0)), intent(in) :: Amat(:, :)
        real(kind(1d0)), intent(in) :: OrbTab(:, :)
        real(kind(1d0)), intent(out) :: ChDen(:)

        integer :: number_of_orbitals, ii_orb, ij_orb

        number_of_orbitals = size(OrbTab, 2)
        !.. Tabulate Charge Density
        ChDen = 0.d0
        do iPts = 1, nPts
            do ii_orb = 1, number_of_orbitals
                do ij_orb = 1, number_of_orbitals
                    ChDen(iPts) = ChDen(iPts) + OrbTab(iPts, ii_orb) * Amat(ii_orb, ij_orb) * OrbTab(iPts, ij_orb)
                enddo
            enddo
        enddo
    end subroutine TabulateChargeDensity
    !!!t
    subroutine TabulateChargeDensity1(Amat, OrbTab, ChDen)
        real(kind(1d0)), intent(in) :: Amat(:, :)
        real(kind(1d0)), intent(in) :: OrbTab(:, :)
        real(kind(1d0)), intent(out) :: ChDen(:)

        integer :: number_of_orbitals, ii_orb, ij_orb

        number_of_orbitals = size(OrbTab, 2)
        !.. Tabulate Charge Density
        ChDen = 0.d0
        do iPts = 1, nPts
            do ii_orb = 1, number_of_orbitals
                do ij_orb = 1, number_of_orbitals
                    ChDen(iPts) = ChDen(iPts) + OrbTab(iPts, ii_orb) * OrbTab(iPts, ij_orb)
                enddo
            enddo
        enddo
    end subroutine TabulateChargeDensity1


    complex(kind(1d0)) function zTraceFunction(zA) result(zTrace)
        integer :: i
        complex(kind(1d0)), intent(in) :: zA(:, :)
        zTrace = 0.d0
        do i = 1, size(zA, 1)
            zTrace = zTrace + zA(i, i)
        enddo
    end function zTraceFunction


    subroutine LiouvillianPropagator(L0_Eval, L0_LEvec, L0_REvec, zStatRho)

        complex(kind(1d0)), intent(in) :: L0_LEvec(:, :), L0_REvec(:, :), L0_Eval(:)
        complex(kind(1d0)), intent(inout) :: zStatRho(:, :)

        logical, save :: FIRST_CALL = .TRUE.
        integer :: nLiou
        integer :: i_state, j
        real   (kind(1d0)) :: EFieldX, EFieldY, EFieldZ
        complex(kind(1d0)), allocatable, save :: zmat1(:, :), zmat2(:, :)
        complex(kind(1d0)), allocatable, save :: zPureMat(:, :)
        real   (kind(1d0)), allocatable, save :: Epure(:)
        complex(kind(1d0)), allocatable, save :: zUMAT_DX(:, :), zUMAT_DY(:, :), zUMAT_DZ(:, :)
        complex(kind(1d0)), allocatable, save :: zUMAT_DYs_DX(:, :), zUMAT_DZs_DY(:, :)
        complex(kind(1d0)), allocatable, save :: RhoVec2(:), RhoVec(:)
        real   (kind(1d0)), allocatable, save :: EVEC_DX(:), EVEC_DY(:), EVEC_DZ(:)

        nLiou = size(L0_Eval, 1)

        if(FIRST_CALL)then
            allocate(RhoVec (nLiou))
            allocate(RhoVec2(nLiou))
            RhoVec = Z0
            RhoVec2 = Z0
        endif
        call HilbertToLiouvilleMatrix(zStatRho, RhoVec)

        !.. First half of the field-free propagation in Liouville Space across an interval dt/2
        call ZGEMV("C", nLiou, nLiou, Z1, L0_LEvec, nLiou, RhoVec, 1, Z0, RhoVec2, 1)
        RhoVec2 = RhoVec2 * exp(-Zi * L0_Eval * dt / 2.d0)
        call ZGEMV("N", nLiou, nLiou, Z1, L0_REvec, nLiou, RhoVec2, 1, Z0, RhoVec, 1)
        call LiouvilleToHilbertMatrix(RhoVec, zStatRho)

        nStates = size(zStatRho, 1)
        if(FIRST_CALL)then
            allocate(zmat1(nStates, nStates), zmat2(nStates, nStates), Epure(nStates), zPureMat(nStates, nStates))
            zmat1 = 0.d0
            zmat2 = 0.d0
            Epure = 0.d0
            zPureMat = 0.d0

            !.. Diagonalizes the dipole matrices
            call DiagonalizeDipole(Dmat, &
                    EVEC_DX, zUMAT_DX, &
                    EVEC_DY, zUMAT_DY, &
                    EVEC_DZ, zUMAT_DZ)
            !
            !.. Computes combined matrices
            allocate(zUMAT_DYs_DX(nStates, nStates))
            zUMAT_DYs_DX = Z0
            call ZGEMM("C", "N", nStates, nStates, nStates, Z1, zUMAT_DY, nStates, &
                    zUMAT_DX, nStates, Z0, zUMAT_DYs_DX, nStates)
            allocate(zUMAT_DZs_DY(nStates, nStates))
            zUMAT_DZs_DY = Z0
            call ZGEMM("C", "N", nStates, nStates, nStates, Z1, zUMAT_DZ, nStates, &
                    zUMAT_DY, nStates, Z0, zUMAT_DZs_DY, nStates)

            FIRST_CALL = .FALSE.
        endif

        !.. Diagonalizes the density matrix to determine the pure states it is composed of
        zPureMat = zStatRho
        call Short_ZHEEV(nStates, zPureMat, EPure)
        !
        !.. Propagates the individual states in the presence of the field

        EFieldX = ExternalElectricFieldCart(t + dt / 2.d0, train(iSim), 1)
        EFieldY = ExternalElectricFieldCart(t + dt / 2.d0, train(iSim), 2)
        EFieldZ = ExternalElectricFieldCart(t + dt / 2.d0, train(iSim), 3)
        zmat1 = Z0
        zmat2 = Z0
        call ZGEMM("C", "N", nStates, nStates, nStates, Z1, zUMAT_DX, nStates, zPureMat, nStates, Z0, zmat1, nStates)
        do i = 1, nStates
            zmat1(i, :) = exp(-Zi * EfieldX * EVEC_DX(i) * dt / 2.d0) * zmat1(i, :)
        enddo
        call ZGEMM("N", "N", nStates, nStates, nStates, Z1, zUMAT_DYs_DX, nStates, zmat1, nStates, Z0, zmat2, nStates)
        do i = 1, nStates
            zmat2(i, :) = exp(-Zi * EfieldY * EVEC_DY(i) * dt / 2.d0) * zmat2(i, :)
        enddo
        call ZGEMM("N", "N", nStates, nStates, nStates, Z1, zUMAT_DZs_DY, nStates, zmat2, nStates, Z0, zmat1, nStates)
        do i = 1, nStates
            zmat1(i, :) = exp(-Zi * EfieldZ * EVEC_DZ(i) * dt) * zmat1(i, :)
        enddo
        call ZGEMM("C", "N", nStates, nStates, nStates, Z1, zUMAT_DZs_DY, nStates, zmat1, nStates, Z0, zmat2, nStates)
        do i = 1, nStates
            zmat2(i, :) = exp(-Zi * EfieldY * EVEC_DY(i) * dt / 2.d0) * zmat2(i, :)
        enddo
        call ZGEMM("C", "N", nStates, nStates, nStates, Z1, zUMAT_DYs_DX, nStates, zmat2, nStates, Z0, zmat1, nStates)
        do i = 1, nStates
            zmat1(i, :) = exp(-Zi * EfieldX * EVEC_DX(i) * dt / 2.d0) * zmat1(i, :)
        enddo
        call ZGEMM("N", "N", nStates, nStates, nStates, Z1, zUMAT_DX, nStates, zmat1, nStates, Z0, zPureMat, nStates)
        !LiouvillianPropagator
        !.. Form the density matrix again
        zStatRho = Z0
        do i_state = 1, nStates
            do j = 1, nStates
                zStatRho(:, j) = zStatRho(:, j) + EPure(i_state) * zPureMat(:, i_state) * conjg(zPureMat(j, i_state))
            enddo
        enddo
        call HilbertToLiouvilleMatrix(zStatRho, RhoVec)

        !.. Second half of the field-free propagation in Liouville Space across an interval dt/2
        call ZGEMV("C", nLiou, nLiou, Z1, L0_LEvec, nLiou, RhoVec, 1, Z0, RhoVec2, 1)
        RhoVec2 = RhoVec2 * exp(-Zi * L0_Eval * dt / 2.d0)
        call ZGEMV("N", nLiou, nLiou, Z1, L0_REvec, nLiou, RhoVec2, 1, Z0, RhoVec, 1)
        call LiouvilleToHilbertMatrix(RhoVec, zStatRho)
    end subroutine LiouvillianPropagator


    subroutine HilbertToLiouvilleMatrix(RhoMat, RhoVec)
        complex(kind(1d0)), intent(in) :: RhoMat(:, :)
        complex(kind(1d0)), intent(out) :: RhoVec(:)
        integer :: i, j, n, iPair
        n = size(RhoMat, 1)
        do i = 1, n
            do j = 1, n
                call HilbertToLiouvilleIndexes(i, j, iPair)
                RhoVec(iPair) = RhoMat(i, j)
            enddo
        enddo
    end subroutine HilbertToLiouvilleMatrix

    subroutine LiouvilleToHilbertMatrix(RhoVec, RhoMat)
        complex(kind(1d0)), intent(in) :: RhoVec(:)
        complex(kind(1d0)), intent(out) :: RhoMat(:, :)
        integer :: i, j, n, iPair
        n = size(RhoMat, 1)
        RhoMat = Z0
        do j = 1, n
            do i = 1, n
                call HilbertToLiouvilleIndexes(i, j, iPair)
                RhoMat(i, j) = RhoVec(iPair)
            enddo
        enddo
    end subroutine LiouvilleToHilbertMatrix

    subroutine HilbertToLiouvilleIndexes(i, j, iPair)
        integer, intent(in) :: i, j
        integer, intent(out) :: iPair
        if(i==j)then
            iPair = i**2
        elseif(i<j)then
            iPair = (j - 1)**2 + i
        else
            iPair = (i - 1) * i + j
        endif
    end subroutine HilbertToLiouvilleIndexes

    subroutine LiouvilleToHilbertIndexes(iPair, i, j)
        integer, intent(in) :: iPair
        integer, intent(out) :: i, j
        integer :: n, n2
        n = int(sqrt(1.d0 * iPair))
        n2 = n**2
        if(n2==iPair)then
            i = n
            j = n
        elseif(iPair - n2<=n)then
            i = iPair - n2
            j = n + 1
        else
            i = n + 1
            j = iPair - n2 - n
        endif
    end subroutine LiouvilleToHilbertIndexes

    !    !####################
    !!
    !    !..FUNCTIONS
    !    !
    !    !.. Eq. 3 $w_a \vec(r)$  (WEIGHTS)
    !    !Result =P_a(\vec{r})/\sum_bP_b(\vec{r}) $
    real(kind(1d0)) function wkfun1(rvec, iAtom, AtCoord, &
            nAtoms, k) result(res)
        !
        integer, intent(in) :: iAtom, k, nAtoms
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms)
        !
        res = Pkfuna1(rvec, iAtom, AtCoord, nAtoms, k) / PkfunTot1(rvec, AtCoord, nAtoms, k)
    end function wkfun1
    !
    !    !..Eq. 3 \sum_bP_b(\vec{r}) in the denominator
    !    !
    real(kind(1d0)) function PkFunTot1(rvec, AtCoord, nAtoms, k) result(res)
        integer, intent(in) :: nAtoms
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms)
        integer, intent(in) :: k
        integer :: iAtom
        res = 0.d0
        do iAtom = 1, nAtoms
            res = res + Pkfuna1(rvec, iAtom, AtCoord, nAtoms, k)
        enddo
    end function PkFunTot1
    !
    !    !.. Eq. 2 P_a(\vec{r}) in and nominator of Eq. 3
    !    !
    real(kind(1d0)) function Pkfuna1(rvec, iAtom, AtCoord, nAtoms, k) result(res)
        integer, intent(in) :: nAtoms, iAtom
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms)
        integer, intent(in) :: k
        integer :: jAtom
        res = 1.d0
        do jAtom = 1, nAtoms
            if(jAtom == iAtom) cycle
            res = res * skfunab1(rvec, AtCoord(:, iAtom), AtCoord(:, jAtom), k)
        enddo
    end function Pkfuna1
    !
    !    !..
    !    !
    real(kind(1d0)) function skfunab1(rvec, avec, bvec, k) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        integer, intent(in) :: k
        real(kind(1d0)) :: mu
        mu = EllipticalCoord1(rvec, avec, bvec)
        res = skfun1(mu, k)
    end function skfunab1
    !
    !    !.. $s(\mu): [-1,1]  \rightarrow [0,1]$  (Step Function)
    !    !
    real(kind(1d0)) function skfun1(mu, k) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu
        integer, intent(in) :: k
        res = 0.5d0 * (1.d0 - fkfun1(mu, k))
    end function skfun1
    !
    !    !.. Eq. 4 Recursive $f_k(\mu)$
    !    !
    real(kind(1d0)) recursive function fkfun1(mu, k) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu
        integer, intent(in) :: k
        res = 0.5d0 * mu * (3.d0 - mu**2)
        if(k==1)return
        res = fkfun1(res, k - 1)
    end function fkfun1
    !
    !    !.. Eq. 1 $\mu_{ab}(\vec{r})$
    !    !
    real(kind(1d0)) function EllipticalCoord1(rvec, avec, bvec) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        real(kind(1d0)) :: r_a, r_b, R_ab ! $r_a$, $r_b$, $R_ab$
        r_a = EuclDist1(rvec, avec)
        r_b = EuclDist1(rvec, bvec)
        R_ab = EuclDist1(avec, bvec)
        res = (r_a - r_b) / R_ab
    end function EllipticalCoord1
    !
    !    !.. Distance
    !    !
    real(kind(1d0)) function EuclDist1(avec, bvec) result(res)
        real(kind(1d0)), intent(in) :: avec(3), bvec(3)
        real(kind(1d0)) :: cvec(3)
        cvec = (avec - bvec)
        res = sqrt(sum(cvec * cvec))
    end function EuclDist1
    !    !
    !    !END OF FUNCTIONS
    !
    subroutine ComputeAtomicWeights1(nPts, gridv, nAtoms, AtCoord, WeightV)

        integer, intent(in) :: nPts, nAtoms
        real(kind(1d0)), intent(in) :: gridv(:, :), AtCoord(:, :)
        real(kind(1d0)), allocatable, intent(out) :: WeightV(:, :)

        integer, parameter :: k = 4
        integer :: iPts, iAtom
        real(kind(1d0)) :: rvec(3)

        !
        write(*, "(a)") "Computing Becke's Weights"
        !

        allocate(WEIGHTV(nPts, nAtoms))
        do iPts = 1, nPts
            rvec = gridv(:, iPts)
            do iAtom = 1, nAtoms
                WeightV(iPts, iAtom) = wkfun1(rvec, iAtom, AtCoord, nAtoms, k)
            enddo
        enddo
    end subroutine ComputeAtomicWeights1


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!>>>>>>>>>>>>>>
    subroutine Computevolume(nPts, volume, gridv)
        integer, intent(in) :: nPts
        real(kind(1d0)), intent(out) :: volume
        real(kind(1d0)), allocatable, intent(in) :: gridv(:, :)
        real(kind(1d0)) :: coord1, coord2, delta
        real(kind(1d0)), parameter :: threshold = 1d-10

        integer :: iPts, iCoord
        volume = 1.d0

        do iCoord = 1, 3
            coord1 = gridv(iCoord, 1)
            delta = 1d10
            do iPts = 2, nPts
                coord2 = gridv(iCoord, iPts)
                if (abs(coord1 - coord2)>threshold) then
                    delta = min(delta, abs(coord1 - coord2))
                endif
            enddo
            volume = volume * delta
        enddo

    end subroutine Computevolume

    subroutine AtomicRadius_Bragg_Slater_Becke(atom_names, nAtom, Radius_BS)
        real(kind(1d0)), allocatable, intent(out) :: Radius_BS(:)
        integer, intent(in) :: nAtom
        character(len = 16), intent(in) :: atom_names(:)
        integer :: iAtom
        allocate(Radius_BS(nAtom))

        do iAtom = 1, nAtom
            Radius_BS(iAtom) = Radius_Table(atom_names(iAtom))
        end do
    end subroutine AtomicRadius_Bragg_Slater_Becke


    real(kind(1d0)) function Radius_Table(Atomic_Name) result(Atomic_Radius)
        character(len = 16), intent(in) :: Atomic_Name

        if (Atomic_Name == "H") then
            Atomic_Radius = 0.35

        elseif (Atomic_Name == "C") then
            Atomic_Radius = 0.70

        elseif (Atomic_Name == "N") then
            Atomic_Radius = 0.65

        elseif (Atomic_Name == "O") then
            Atomic_Radius = 0.60

        elseif (Atomic_Name == "Mg") then
            Atomic_Radius = 1.50
        end if

    end function Radius_Table

    subroutine ComputeAtomicCharges(OrbitalDensity, BeckeMatrix, QchargeVec)
        real(kind(1d0)), intent(in) :: OrbitalDensity(:, :)
        real(kind(1d0)), intent(in) :: BeckeMatrix   (:, :, :)
        real(kind(1d0)), allocatable, intent(out) :: QchargeVec    (:)
        integer :: iAtom, nAtoms, number_of_orbitalss, jOrb, iOrb
        real(kind(1d0)), external :: DDOT
        number_of_orbitalss = size(BeckeMatrix, 1)
        nAtoms = size(BeckeMatrix, 3)
        if(.not.allocated(QchargeVec))allocate(QchargeVec(nAtoms))
        QchargeVec = 0.d0
        do iAtom = 1, nAtoms
            !            do jOrb = 1, number_of_orbitalss
            !                do iOrb = 1, number_of_orbitalss
            !                    QchargeVec(iAtom) = QchargeVec(iAtom) + OrbitalDensity(iOrb, jOrb) * BeckeMatrix(iOrb, jOrb, iAtom)
            !                end do
            !            end do
            QchargeVec(iAtom) = ddot(number_of_orbitalss * number_of_orbitalss, OrbitalDensity, 1, BeckeMatrix(1, 1, iAtom), 1)
        enddo
    end subroutine ComputeAtomicCharges

    subroutine ComputeBeckeMatrix(WeightV, OrbTab, Becke)
        real(kind(1d0)), allocatable, intent(in) :: WeightV(:, :)
        real(kind(1d0)), allocatable, intent(in) :: OrbTab (:, :)
        real(kind(1d0)), allocatable, intent(out) :: Becke(:, :, :)

        integer :: nPts, number_of_orbitalss, nAtoms
        integer :: iPts, iOrb1, iOrb2, iAtom
        real(kind(1d0)) :: dsum

        nAtoms = size(WeightV, 2)
        nPts = size(WeightV, 1)
        number_of_orbitalss = size(OrbTab, 2)

        if(allocated(Becke))deallocate(Becke)
        allocate(Becke(number_of_orbitalss, number_of_orbitalss, nAtoms))

        Becke = 0.d0
        do iAtom = 1, nAtoms
            do iOrb2 = 1, number_of_orbitalss
                do iOrb1 = 1, number_of_orbitalss
                    dsum = 0.d0
                    do iPts = 1, nPts
                        dsum = dsum + WeightV(iPts, iAtom) * OrbTab(iPts, iOrb1) * OrbTab(iPts, iOrb2)
                    enddo
                    Becke(iOrb1, iOrb2, iAtom) = dsum
                enddo
            enddo
        enddo
    end subroutine ComputeBeckeMatrix


    subroutine ComputeAtomicWeights(nPts, gridv, nAtoms, AtCoord, WeightV, Radius_BS)
        integer, intent(in) :: nPts, nAtoms
        real(kind(1d0)), intent(in) :: gridv(:, :), AtCoord(:, :), Radius_BS(:)
        real(kind(1d0)), allocatable, intent(out) :: WeightV(:, :)

        integer, parameter :: k = 4
        integer :: iPts, iAtom
        real(kind(1d0)) :: rvec(3)

        allocate(WEIGHTV(nPts, nAtoms))
        do iPts = 1, nPts
            rvec = gridv(:, iPts)
            do iAtom = 1, nAtoms
                WeightV(iPts, iAtom) = wkfun(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS)
            enddo
        enddo
    end subroutine ComputeAtomicWeights
    !..FUNCTIONS
    !
    !.. Eq. 3 $w_a \vec(r)$  (WEIGHTS)
    !Result =P_a(\vec{r})/\sum_bP_b(\vec{r}) $
    real(kind(1d0)) function wkfun(rvec, iAtom, AtCoord, &
            nAtoms, k, Radius_BS) result(res)
        !
        integer, intent(in) :: iAtom, k, nAtoms
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms), Radius_BS(:)
        !
        res = Pkfuna(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS) / PkfunTot(rvec, AtCoord, nAtoms, k, Radius_BS)
    end function wkfun

    !..Eq. 3 \sum_bP_b(\vec{r}) in the denominator
    !
    real(kind(1d0)) function PkFunTot(rvec, AtCoord, nAtoms, k, Radius_BS) result(res)
        integer, intent(in) :: nAtoms
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms), Radius_BS(:)
        integer, intent(in) :: k
        integer :: iAtom
        res = 0.d0
        do iAtom = 1, nAtoms
            res = res + Pkfuna(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS)
        enddo
    end function PkFunTot

    !.. Eq. 2 P_a(\vec{r}) in and nominator of Eq. 3
    !
    real(kind(1d0)) function Pkfuna(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS) result(res)
        integer, intent(in) :: nAtoms, iAtom
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms), Radius_BS(:)
        integer, intent(in) :: k
        integer :: jAtom
        real(kind(1d0)) :: a_ij, R_i, R_j
        res = 1.d0
        do jAtom = 1, nAtoms
            if(jAtom == iAtom) cycle
            R_i = Radius_BS(iAtom)
            R_j = Radius_BS(jAtom)
            a_ij = get_param_a(R_i, R_j)
            res = res * skfunab(rvec, AtCoord(:, iAtom), AtCoord(:, jAtom), k, a_ij)
        enddo
    end function Pkfuna


    !..
    !
    real(kind(1d0)) function skfunab(rvec, avec, bvec, k, a_ij) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        integer, intent(in) :: k
        real(kind(1d0)) :: mu_old, mu_new, a_ij
        mu_old = EllipticalCoord(rvec, avec, bvec)
        mu_new = new_mu_transformation(mu_old, a_ij)
        res = skfun(mu_new, k)
    end function skfunab

    !.. $s(\mu): [-1,1]  \rightarrow [0,1]$  (Step Function)
    !
    real(kind(1d0)) function skfun(mu, k) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu
        integer, intent(in) :: k
        real(kind(1d0)) :: r_m

        res = 0.5d0 * (1.d0 - fkfun(mu, k))
    end function skfun

    !.. Eq. 4 Recursive $f_k(\mu)$
    !
    real(kind(1d0)) recursive function fkfun(mu, k) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu
        integer, intent(in) :: k
        res = 0.5d0 * mu * (3.d0 - mu**2)
        if(k==1)return
        res = fkfun(res, k - 1)
    end function fkfun

    !.. Eq. 1 $\mu_{ab}(\vec{r})$
    !
    real(kind(1d0)) function EllipticalCoord(rvec, avec, bvec) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        real(kind(1d0)) :: r_a, r_b, R_ab ! $r_a$, $r_b$, $R_ab$
        r_a = EuclDist(rvec, avec)
        r_b = EuclDist(rvec, bvec)
        R_ab = EuclDist(avec, bvec)
        res = (r_a - r_b) / R_ab
    end function EllipticalCoord

    !.. Distance
    !
    real(kind(1d0)) function EuclDist(avec, bvec) result(res)
        real(kind(1d0)), intent(in) :: avec(3), bvec(3)
        real(kind(1d0)) :: cvec(3)
        cvec = (avec - bvec)
        res = sqrt(sum(cvec * cvec))
    end function EuclDist
    !
    !END OF FUNCTIONS


    !..get_param_aij
    !..
    real(kind(1d0)) function get_param_a(R_i, R_j) result(a_ij)
        real(kind(1d0)), intent(in) :: R_i, R_j
        real(kind(1d0)) :: Chi
        Chi = R_i / R_j
        a_ij = (Chi - 1) / (Chi + 1)
        if (abs(a_ij)>(0.5d0)) then
            a_ij = 0.5d0
        endif
    end function get_param_a

    real(kind(1d0)) function new_mu_transformation(mu, a_ij) result(res)
        real(kind(1d0)), intent(in) :: mu, a_ij
        res = mu + a_ij * (1 - mu**2)
    end function new_mu_transformation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!>>>>>>>>>>>>>>



    subroutine DiagonalizeDipole(Dmat, &
            EVEC_DX, zUMAT_DX, &
            EVEC_DY, zUMAT_DY, &
            EVEC_DZ, zUMAT_DZ)
        real   (kind(1d0)), intent(in) :: Dmat(:, :, :)
        complex(kind(1d0)), allocatable, intent(out) :: zUMAT_DX(:, :), zUMAT_DY(:, :), zUMAT_DZ(:, :)
        real   (kind(1d0)), allocatable, intent(out) :: EVEC_DX(:), EVEC_DY(:), EVEC_DZ(:)
        !
        real   (kind(1d0)), allocatable :: rmat(:, :)
        integer :: nStates

        nStates = size(Dmat, 1)
        allocate(rmat(nStates, nStates))
        !
        !.. X dipole
        if(allocated(zUMAT_DX)) deallocate(zUMAT_DX)
        if(allocated(EVEC_DX)) deallocate(EVEC_DX)
        allocate(zUMAT_DX(nStates, nStates), EVEC_DX(nStates))
        rmat = Dmat(:, :, 1)
        EVEC_DX = 0.d0
        call Short_Diag(nStates, rmat, EVEC_DX)
        zUMAT_DX = Z1 * rmat
        !
        !.. Y dipole
        if(allocated(zUMAT_DY)) deallocate(zUMAT_DY)
        if(allocated(EVEC_DY)) deallocate(EVEC_DY)
        allocate(zUMAT_DY(nStates, nStates), EVEC_DY(nStates))
        rmat = Dmat(:, :, 2)
        EVEC_DY = 0.d0
        call Short_Diag(nStates, rmat, EVEC_DY)
        zUMAT_DY = Z1 * rmat
        !
        !.. Z dipole
        if(allocated(zUMAT_DZ)) deallocate(zUMAT_DZ)
        if(allocated(EVEC_DZ)) deallocate(EVEC_DZ)
        allocate(zUMAT_DZ(nStates, nStates), EVEC_DZ(nStates))
        rmat = Dmat(:, :, 3)
        EVEC_DZ = 0.d0
        call Short_Diag(nStates, rmat, EVEC_DZ)
        zUMAT_DZ = Z1 * rmat
        !
        deallocate(rmat)

    end subroutine DiagonalizeDipole

    !.. Build the time-independent Lindblad superoperator
    !..
    !.. Linblad equation (Schroedinger representation)
    !.. d/dt \rho(t) = -i [H,\rho(t)]+
    !                     \sum_{mn} [
    !                                 L_{mn} \rho(t) L_{mn}^\dagger -
    !                    -\frac{1}{2} L_{mn}^\dagger L_{mn} \rho(t)
    !                    -\frac{1}{2} \rho(t) L_{mn}^\dagger L_{mn}
    !
    !   H = H_0 + (\hat{\epsilon}\cdot\vec{\mu}) E(t)
    !
    !.. For pure dephasing without relaxation:
    !
    !   L_{nm} = \delta_{mn} \sqrt{\gamma_m/2} | m \rangle \langle m |
    !
    !.. For relaxation without pure dephasing:
    !
    !   L_{mn} = (1-\delta_{nm}) \sqrt{\gamma_{mn}} |m\rangle \langle n|
    !
    !   subject to the constraint \gamma_{mn} = e^{-\beta \omega_{mn}} \gamma_{nm}
    !
    !   where \beta = 1/(k_B T)
    !
    !################################################
    subroutine ComputeLiouvillian0(Evec, Dmat, Liouvillian0)
        real   (kind(1d0)), intent(in) :: Evec(:)
        real   (kind(1d0)), intent(in) :: Dmat(:, :, :)
        complex(kind(1d0)), allocatable, intent(out) :: Liouvillian0(:, :)

        real(kind(1d0)), parameter :: WATER_MELTING_POINT = 273.15
        real(kind(1d0)), parameter :: bath_temperature = WATER_MELTING_POINT + 3000
        real(kind(1d0)), parameter :: BETA = 1.d0 / (BOLTZMANN_CONSTANT_auK * bath_temperature)
        real(kind(1d0)), parameter :: dephasing_factor = 1.d-3
        real(kind(1d0)), parameter :: relaxation_factor = 1.d-3

        integer :: nLiou, nStates, i, j, iLiou, i1, i2, iLiou1, iLiou2
        !.. Dephasing and Relaxation Constants
        real   (kind(1d0)), allocatable :: PairGamma(:, :), TotalGamma(:)
        real   (kind(1d0)) :: AverageDipole, dBuf

        nStates = size(Evec, 1)

        nLiou = nStates ** 2
        allocate(Liouvillian0(nLiou, nLiou))
        Liouvillian0 = Z0
        !
        !.. Unitary component
        do iLiou = 1, nLiou
            call LiouvilleToHilbertIndexes(iLiou, i, j)
            Liouvillian0(iLiou, iLiou) = Evec(i) - Evec(j)
        end do
        !
        !.. Determine the factors \gamma_{mn}
        allocate(PairGamma(nStates, nStates))
        PairGamma = 0.d0
        do i = 1, nStates
            !
            dBuf = 0.d0
            do j = 1, nStates
                if(j==i)cycle
                !
                AverageDipole = sqrt(sum(Dmat(i, j, :)**2) / 3.d0)
                dBuf = dBuf + AverageDipole
                !
                !.. Relaxation Factor
                PairGamma(i, j) = relaxation_factor * AverageDipole
                if(Evec(i)>Evec(j)) PairGamma(i, j) = PairGamma(i, j) * exp(-BETA * (Evec(i) - Evec(j)))
                !
            enddo
            !
            !.. Dephasing Factor
            PairGamma(i, i) = dephasing_factor * sqrt(dBuf)
            !
        enddo
        !
        !.. Determines \Gamma_i = \sum_m \gamma_{mi}
        allocate(TotalGamma(nStates))
        do i = 1, nStates
            TotalGamma(i) = sum(PairGamma(:, i))
        enddo
        !
        !.. Dissipative component
        do iLiou = 1, nLiou
            call LiouvilleToHilbertIndexes(iLiou, i, j)
            Liouvillian0(iLiou, iLiou) = Liouvillian0(iLiou, iLiou) - Zi * (TotalGamma(i) + TotalGamma(j)) / 2.d0
        enddo
        do i1 = 1, nStates
            call HilbertToLiouvilleIndexes(i1, i1, iLiou1)
            do i2 = 1, nStates
                call HilbertToLiouvilleIndexes(i2, i2, iLiou2)
                Liouvillian0(iLiou1, iLiou2) = Liouvillian0(iLiou1, iLiou2) + Zi * PairGamma(i1, i2)
            enddo
        enddo

    end subroutine ComputeLiouvillian0
    !###########################################

    subroutine ComputeLiouvillian01(Evec, Dmat, Liouvillian0, bath_temperature, dephasing_factor, relaxation_factor)
        real   (kind(1d0)), intent(in) :: Evec(:)
        real   (kind(1d0)), intent(in) :: Dmat(:, :, :)
        complex(kind(1d0)), allocatable, intent(out) :: Liouvillian0(:, :)

        real(kind(1d0)), intent(in) :: dephasing_factor
        real(kind(1d0)), intent(in) :: relaxation_factor
        real(kind(1d0)), intent(in) :: bath_temperature

        integer :: nLiou, nStates, i, j, iLiou, i1, i2, iLiou1, iLiou2
        !.. Dephasing and Relaxation Constants
        real   (kind(1d0)), allocatable :: PairGamma(:, :), TotalGamma(:)
        real   (kind(1d0)) :: AverageDipole, dBuf
        real(kind(1d0)) :: BETA

        BETA = 1.d0 / (BOLTZMANN_CONSTANT_auK * bath_temperature)

        nStates = size(Evec, 1)

        nLiou = nStates ** 2
        allocate(Liouvillian0(nLiou, nLiou))
        Liouvillian0 = Z0
        !
        !.. Unitary component
        do iLiou = 1, nLiou
            call LiouvilleToHilbertIndexes(iLiou, i, j)
            Liouvillian0(iLiou, iLiou) = Evec(i) - Evec(j)
        end do
        !
        !.. Determine the factors \gamma_{mn}
        allocate(PairGamma(nStates, nStates))
        PairGamma = 0.d0
        do i = 1, nStates
            !
            dBuf = 0.d0
            do j = 1, nStates
                if(j==i)cycle
                !
                AverageDipole = sqrt(sum(Dmat(i, j, :)**2) / 3.d0)
                dBuf = dBuf + AverageDipole
                !
                !.. Relaxation Factor
                PairGamma(i, j) = relaxation_factor * AverageDipole
                if(Evec(i)>Evec(j)) PairGamma(i, j) = PairGamma(i, j) * exp(-BETA * (Evec(i) - Evec(j)))
                !
            enddo
            !
            !.. Dephasing Factor
            PairGamma(i, i) = dephasing_factor * sqrt(dBuf)
            !
        enddo
        !
        !.. Determines \Gamma_i = \sum_m \gamma_{mi}
        allocate(TotalGamma(nStates))
        do i = 1, nStates
            TotalGamma(i) = sum(PairGamma(:, i))
        enddo
        !
        !.. Dissipative component
        do iLiou = 1, nLiou
            call LiouvilleToHilbertIndexes(iLiou, i, j)
            Liouvillian0(iLiou, iLiou) = Liouvillian0(iLiou, iLiou) - Zi * (TotalGamma(i) + TotalGamma(j)) / 2.d0
        enddo
        do i1 = 1, nStates
            call HilbertToLiouvilleIndexes(i1, i1, iLiou1)
            do i2 = 1, nStates
                call HilbertToLiouvilleIndexes(i2, i2, iLiou2)
                Liouvillian0(iLiou1, iLiou2) = Liouvillian0(iLiou1, iLiou2) + Zi * PairGamma(i1, i2)
            enddo
        enddo

    end subroutine ComputeLiouvillian01


    !.. Diagonalize the time-independent Lindblad superoperator
    subroutine DiagonalizeLindblad0(Liouvillian0, L0_Eval, L0_LEvec, L0_REvec)
        complex(kind(1d0)), intent(inout) :: Liouvillian0(:, :)
        complex(kind(1d0)), allocatable, intent(out) :: L0_LEvec(:, :), L0_REvec(:, :), L0_Eval(:)
        complex(kind(1d0)), allocatable :: zmat(:, :)
        integer :: nLiou, iLiou, n, i, j

        nLiou = size(Liouvillian0, 1)

        allocate(L0_LEvec(nLiou, nLiou))
        allocate(L0_REvec(nLiou, nLiou))
        allocate(L0_Eval (nLiou))
        L0_LEvec = Z0
        L0_REvec = Z0
        L0_Eval = Z0
        call Short_Diag(nLiou, Liouvillian0, L0_Eval, L0_LEvec, L0_REvec)
        do iLiou = 1, nLiou
            write(*, *) iLiou, L0_Eval(iLiou)
        enddo
        !.. Test that U_R U_L^\dagger is indeed the identity
        n = nLiou
        allocate(zmat(n, n))
        zmat = Z0
        call ZGEMM("N", "T", n, n, n, Z1, L0_REvec, n, L0_LEvec, n, Z0, zmat, n)
        do i = 1, n
            zmat(i, i) = zmat(i, i) - Z1
        enddo
        write(*, *)"|U_R U_L^+ -1|_1 = ", sum(abs(zmat))

    end subroutine DiagonalizeLindblad0

    !>>
    subroutine Write_Dipole(FileName, Dipole, n_times, t_min, dt)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), intent(in) :: Dipole(:, :)
        real   (kind(1d0)), intent(in) :: t_min, dt
        integer, intent(in) :: n_times

        real   (kind(1d0)) :: t
        integer :: uid_dipole, iPol, it

        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do it = 1, n_times
            t = t_min + dt * dble(it - 1)
            write(uid_dipole, "(i4,*(x,E24.16))") it, t, ((dble(Dipole(iPol, it)), aimag(Dipole(iPol, it))), iPol = 1, 3)
        enddo
        close(uid_dipole)
    end subroutine Write_Dipole
    !
    !
    subroutine Write_Dipole1(FileName, Dipole, n_times, t_min, dt)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(in) :: Dipole(:, :)
        !        complex(kind(1d0)), intent(in) :: Dipole(:, :)
        real   (kind(1d0)), intent(in) :: t_min, dt
        integer, intent(in) :: n_times

        real   (kind(1d0)) :: t
        integer :: uid_dipole, iPol, it

        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do it = 1, n_times
            t = t_min + dt * dble(it - 1)
            !            write(uid_dipole, "(i4,*(x,E24.16))") it, t, ((dble(Dipole(iPol, it)), aimag(Dipole(iPol, it))), iPol = 1, 3)
            write(uid_dipole, "(i4,*(x,E24.16))") it, t, ((Dipole(iPol, it), 0.d0), iPol = 1, 3)
        enddo
        close(uid_dipole)
    end subroutine Write_Dipole1

    subroutine Write_Q_Charge(FileName, Charge, n_times, t_min, dt, nAtoms)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(in) :: Charge(:, :)
        real   (kind(1d0)), intent(in) :: t_min, dt
        integer, intent(in) :: n_times, nAtoms

        real   (kind(1d0)) :: t
        integer :: uid_AtomicCharge, iPol, it, iAtom

        !.. Save Q_Charge
        open(newunit = uid_AtomicCharge, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        !..Regular
        do it = 1, n_times
            t = t_min + dt * dble(it - 1)
            write(uid_AtomicCharge, "(*(x,e24.16))") t, sum(Charge(:, it)), (Charge(iAtom, it), iAtom = 1, nAtoms)
        enddo
        !
        !
        !..New
        !        do it = 1, n_times
        !            t = t_min + dt * dble(it - 1)
        !            write(uid_AtomicCharge, "(*(x,e24.16))") t, sum(Charge(:, :, it)), ((Charge(iPol, iAtom, it), iPol = 1, 3), iAtom = 1, nAtoms)
        !        enddo
        !        !
        close(uid_AtomicCharge)
    end subroutine Write_Q_Charge
    !
    !
    subroutine Write_Q_Charge1(FileName, Charge, n_times, t_min, dt, nAtoms)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(in) :: Charge(:, :, :)
        real   (kind(1d0)), intent(in) :: t_min, dt
        integer, intent(in) :: n_times, nAtoms

        real   (kind(1d0)) :: t
        integer :: uid_AtomicCharge, iPol, it, iAtom

        !.. Save Q_Charge
        open(newunit = uid_AtomicCharge, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")

        !..New
        do it = 1, n_times
            t = t_min + dt * dble(it - 1)
            write(uid_AtomicCharge, "(*(x,e24.16))") t, sum(Charge(:, :, it)), ((Charge(iPol, iAtom, it), iPol = 1, 3), iAtom = 1, nAtoms)
        enddo
        !
        close(uid_AtomicCharge)
    end subroutine Write_Q_Charge1

    !    !*** APPARENTLY, IT IS NOT WORKING YET!
    !    subroutine ComputeDipoleFTplus(&
    !            L0_Eval, L0_LEvec, L0_REvec, &
    !            OmegaVec, zStatRho0, DipoleFTplus, &
    !            StepTime, StepWidth)
    !
    !        complex(kind(1d0)), intent(in) :: L0_LEvec(:, :), L0_REvec(:, :), L0_Eval(:), zStatRho0(:, :)
    !        real   (kind(1d0)), intent(in) :: OmegaVec(:)
    !        complex(kind(1d0)), intent(out) :: DipoleFTplus(:, :)
    !        real   (kind(1d0)), intent(in) :: StepTime, StepWidth
    !
    !        real(kind(1d0)), parameter :: DIPOLE_THRESHOLD = 1.d-49
    !
    !        logical, save :: FIRST_CALL = .TRUE.
    !        integer, save :: nliou, nStates, nOmegas
    !        complex(kind(1d0)), allocatable, save :: zmat1(:, :), zmat2(:, :), zvec1(:), zvec2(:), zvec3(:)
    !        integer :: iPol, iOmega, iLiou, i_state
    !        real(kind(1d0)) :: w, sigma
    !
    !        if(FIRST_CALL)then
    !
    !            nliou = size(L0_Eval, 1)
    !            nStates = size(zStatRho0, 1)
    !            nOmegas = size(OmegaVec)
    !            allocate(zvec1(nLiou))
    !            allocate(zvec2(nLiou))
    !            allocate(zvec3(nLiou))
    !            allocate(zmat1(nStates, nStates))
    !            allocate(zmat2(nStates, nStates))
    !            zmat1 = Z0
    !            zmat2 = Z0
    !            zvec1 = Z0
    !            zvec2 = Z0
    !            zvec3 = Z0
    !
    !            FIRST_CALL = .FALSE.
    !
    !        endif
    !
    !        sigma = StepWidth / sqrt(2.d0)
    !
    !        do iPol = 1, 3
    !            call ZGEMM("N", "N", nStates, nStates, nStates, Z1, zStatRho0, nStates, &
    !                    Dmat(1, 1, iPol), nStates, Z0, zmat1, nStates)
    !            call HilbertToLiouvilleMatrix(zmat1, zvec1)
    !            call ZGEMV("C", nLiou, nLiou, Z1, L0_LEvec, nLiou, zvec1, 1, Z0, zvec2, 1)
    !            do iOmega = 1, nOmegas
    !                w = OmegaVec(iOmega)
    !                do iLiou = 1, nLiou
    !                    !*** CHECK CONSISTENCY OF SIGMA AND STEPWIDTH
    !                    zvec1(iLiou) = exp(-sigma**2 * (w - L0_Eval(iLiou))**2) / (w - L0_Eval(iLiou)) * zvec2(iLiou)
    !                    !***
    !                enddo
    !                call ZGEMV("N", nLiou, nLiou, Z1, L0_REvec, nLiou, zvec1, 1, Z0, zvec3, 1)
    !                call LiouvilleToHilbertMatrix(zvec3, zmat2)
    !                DipoleFTplus(iPol, iOmega) = 0.d0
    !                do i_state = 1, nStates
    !                    DipoleFTplus(iPol, iOmega) = DipoleFTplus(iPol, iOmega) + zmat2(i_state, i_state)
    !                enddo
    !                DipoleFTplus(iPol, iOmega) = DipoleFTplus(iPol, iOmega) * exp(Zi * w * StepTime) / (2.d0 * PI)
    !
    !                if(abs(DipoleFTplus(iPol, iOmega)) < DIPOLE_THRESHOLD) DipoleFTplus(iPol, iOmega) = 0.d0
    !
    !            enddo
    !        enddo
    !
    !    end subroutine ComputeDipoleFTplus


end program ChargeMigration


!*** IDEALLY, SHOULD COMPUTE THE DIPOLE FROM THE AO - AO DIPOLES.
!!$  !.. Load Dipole matrix element between active MOs
!!$  call LoadDipoleMO( input_directory, number_of_orbitals, ivOrb, MuOrb )
!!$
!!$  !.. Check Dipole Matrix Elements between States
!!$  do ii_state = 1, nStates
!!$     do ij_state = 1, nStates
!!$        TME = 0.d0
!!$        do i = 1, number_of_orbitals
!!$           do j = 1, number_of_orbitals
!!$              TME = TME + TDM(i,j,ii_state,ij_state) * MuOrb(j,i)
!!$           enddo
!!$        enddo
!!$        write(*,*) ii_state, ij_state, TME
!!$     enddo
!!$  enddo


!        if (R_i < 0.36) then
!            r_m = R_i
!        else
!            r_m = R_i * 0.5d0
!        end if