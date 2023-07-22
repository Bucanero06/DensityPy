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
    use ModuleRTP
    use Module_CD_IO

    implicit none

    external zgemm, zgemv, system !Explicit declaration of the EXTERNAL attribute is required.
    !.. Run-time parameters
    !..
    character(len = :), allocatable :: InpDir
    character(len = :), allocatable :: OutDir
    character(len = :), allocatable :: FileGeometry
    integer :: nTimes
    real(kind(1d0)) :: Tmin
    real(kind(1d0)) :: Tmax
    character(len = :), allocatable :: Ext_field_file
    character(len = :), allocatable :: Weight_File
    logical :: Verbous
    logical :: Weight_File_flag
    logical :: SaveDensity
    integer, allocatable :: ivorb(:)
    integer :: OrbNumber
    real(kind(1d0)) :: Volume
    real(kind(1d0)) :: dephasing_factor
    real(kind(1d0)) :: RELAXATION_FACTOR
    real(kind(1d0)) :: BATH_TEMPERATURE
    integer, parameter :: GS_IDX = 1

    !.. Local parameters
    integer :: nStates, nOrb, i
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
    !    real(kind(1d0)) :: Volume

    !.. Statistical Density Matrix
    complex(kind(1d0)), allocatable :: zStatRho(:, :), zStatRho0(:, :)
    real   (kind(1d0)) :: t, dt
    real   (kind(1d0)), allocatable :: ChDen(:), WEIGHTV(:, :)
    real   (kind(1d0)), allocatable :: AtomicChargeVec(:), AtomicChargeEvolution(:, :)

    character(len = 30) :: istrn
    character(len = 16), allocatable :: AtomName(:)

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
    real(kind(1d0)) :: data1, data2, data3, Orbital_overlap_self, Orbital_overlap_other, Computed_Volume
    real   (kind(1d0)), allocatable :: AtomicChargeVec_new(:, :), AtomicChargeEvolution_new(:, :, :)
    integer :: iOrb, jOrb, nOrbs
    real(kind(1d0)), allocatable :: QchargeVec_new    (:, :), R_el(:, :)

    call GetRunTimeParameters(InpDir, OutDir, FileGeometry, &
            nTimes, Tmin, Tmax, Ext_field_file, Verbous, Weight_File, Weight_File_flag, &
            SaveDensity, ivorb, OrbNumber, Volume, dephasing_factor, RELAXATION_FACTOR, BATH_TEMPERATURE)

    call system("mkdir -p " // trim(OutDir))
    call system("mkdir -p " // OutDir // "/Dipole")
    call system("mkdir -p " // OutDir // "/AtomicCharge")
    call system("mkdir -p " // OutDir // "/ChargeDensity")
    call system("mkdir -p " // OutDir // "/Pulses")
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
    dt = (tmax - tmin) / dble(nTimes - 1)
    do i = 1, N_Simulations
        write(*, *) "Writing pulse " // trim(Simulation_Tagv(i))
        strn = trim(OUTDir) // "/Pulses/pulse" // trim(Simulation_Tagv(i))
        call train(i)%Write(strn, Tmin, Tmax, dt)  !>>pulse_trainWrite
        strn = trim(OUTDir) // "/Pulses/FTpulse" // trim(Simulation_Tagv(i))
        call train(i)%WriteFTA(strn) !>>pulse_trainWriteFTA
    enddo
    !    stop!stop here to look at pulses
    !
    !
    call LoadEnergies(InpDir // "/ROOT_ENERGIES", nStates, Evec)

    !.. At the moment, we assume that the TDM are defined as
    !   $\rho^{JI}_{nm} = \langle J | a_n^\dagger a_m | I \rangle $
    call LoadTDMs    (InpDir // "/DENSITY_MATRIX", InpDir // "/TRANSITION_DENSITY_MATRIX", nStates, nOrb, TDM)

    !.. Load Dipole matrix elements $\mu_{IJ}$
    call LoadDipoleME(Dmat, InpDir, nStates)
    allocate(zDmat_t(nStates, nStates, 3))
    do i = 1, 3
        zDmat_t(:, :, i) = Z1 * transpose(Dmat(:, :, i))
    enddo

    if (nOrb .ne. OrbNumber) then
        write(*, *) "Number of orbitals given by '-iOrb' do not match number of orbitals in active the active space"
        stop
    endif

    !.. Load Geometry,Volume, Grid and Orbitals
    call LoadGeometry(nAtoms, AtCoord, FileGeometry, AtomName)
    call LoadGrid(InpDir // "/gridcoord", npts, gridv)! $G=\{\vec{r}_i, i=1,\ldots,N_{G}\}$
    !        nxpoints = int(exp(log(nPts + 0.1d0) / 3.d0))
    call ComputeVolume(nPts, Computed_Volume, gridv) ! Should probably check that this computed volume is the same as the one given in the input file

    !.. Load Orbitals
    call LoadOrbitals(InpDir, nOrb, ivOrb, npts, OrbTab)! $\varphi_n(\vec{r}_i)$

    !.. Compute and Write | or Read Becke's Weights
    if(.not. SaveDensity)then
        !
        call AtomicRadius_Bragg_Slater_Becke(AtomName, nAtoms, Radius_BS)
        if (Weight_File_flag == .True.) then
            call Read_Weights(Weight_File, WEIGHTV, nAtoms, nPts)
        else
            call ComputeAtomicWeights(nPts, gridv, nAtoms, AtCoord, WeightV, Radius_BS)
            call Write_Weights(OUTDir // "/" // Weight_File // "_" // OUTDir, WEIGHTV, gridv, nAtoms, nPts)
        endif
        !.. Compute Berycenter of Atmoic Charges
        call Compute_R_el(gridv, WeightV, OrbTab, R_el)
        call Write_R_el_bc(OutDir, nAtoms, R_el)

        !.. Compute Becke's Matrix
        call ComputeNewBeckeMatrix(WeightV, OrbTab, Becke_new, R_el) !!using electronic barycenter or AtCoord for the nuclear barycenter
        !    BeckeMatrix = BeckeMatrix * Computed_Volume
    end if

    allocate(AtomicChargeVec(nAtoms))
    allocate(AtomicChargeEvolution(nAtoms, nTimes))
    allocate(AtomicChargeEvolution_new(3, nAtoms, nTimes))


    !.. Compute Liouvillian and Diagonalizes it
    call ComputeLiouvillian_0(Evec, Dmat, Liouvillian0, BATH_TEMPERATURE, dephasing_factor, RELAXATION_FACTOR)
    call DiagonalizeLindblad_0(Liouvillian0, L0_Eval, L0_LEvec, L0_REvec)

    allocate(zStatRho (nStates, nStates))
    allocate(zStatRho0(nStates, nStates))

    allocate(ChDen(nPts))

    allocate(zMuEV(3, nTimes))
    zMuEV = Z0

    allocate(Dipole_new(3, nTimes))

    write(*, *) "Starting Sim Loop"
    Sim_loop : do iSim = 1, N_Simulations

        zStatRho = Z0
        zStatRho(GS_IDX, GS_IDX) = 1.d0
        ! Lets keep a percentage write of the simulation only 2 digits after the decimal point
        write(*, '(A, F0.2, A, I0, A)') "Simulation progress: ", 100.d0 * iSim / (N_Simulations)

        !.. Time cycle
        time_loop : do it = 1, nTimes
            !
            t = tmin + dt * dble(it - 1)
            !.. Computes the excitation density wrt ground state
            zStatRho(GS_IDX, GS_IDX) = zStatRho(GS_IDX, GS_IDX) - 1.d0

            !.. Evaluates the expectation value of the dipole as a function of time
            !            write(*, *) it, iSim
            if(.not. SaveDensity)then
                do iPol = 1, 3
                    zMuEV(iPol, it) = zdotu(nStates * nStates, zStatRho, 1, zDmat_t(1, 1, iPol), 1)
                enddo
            end if
            !
            call ComputeOrbitalDensity(zStatRho, TDM, OrbitalDensity)

            if(.not. SaveDensity)then
                !
                call ComputeNewAtomicCharges(OrbitalDensity, Becke_new, AtomicChargeVec_new)
                AtomicChargeVec_new = AtomicChargeVec_new * Computed_Volume
                AtomicChargeEvolution_new(:, :, it) = AtomicChargeVec_new
            end if

            !
            if(SaveDensity)then
                call TabulateChargeDensity(OrbitalDensity, OrbTab, ChDen)
                write(istrn, "(f12.4)")t
                !..Makes new directory inside of ChargeDensity directory for each simulation
                call system("mkdir -p " // OutDir // "/ChargeDensity/ChDenSim" // trim(Simulation_tagv(iSim)))
                call Write_Charge_Density(OutDir // "/ChargeDensity/ChDenSim" // trim(Simulation_tagv(iSim)) // "/ChDen" // trim(adjustl(istrn)), &
                        nPts, gridv, ChDen, Weightv, nAtoms)
            endif
            !
            zStatRho(GS_IDX, GS_IDX) = zStatRho(GS_IDX, GS_IDX) + 1.d0
            !
            if(it == nTimes) exit time_loop
            !
            call LiouvillianPropagator(L0_Eval, L0_LEvec, L0_REvec, zStatRho)
            !
        enddo time_loop
        !.. Save Dipole to filerc
        call Write_Dipole(OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(iSim)), zMuEV, nTimes, tmin, dt)
        !        call Write_Dipole1(OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(iSim)), Dipole_new, nTimes, tmin, dt)

        !.. Save Q_Charge
        !        call Write_Q_Charge(OutDir // "/AtomicCharge/AtomicCharge" // trim(Simulation_tagv(iSim)), AtomicChargeEvolution, nTimes, tmin, dt, nAtoms)
        call Write_Q_Charge1(OutDir // "/AtomicCharge/AtomicCharge" // trim(Simulation_tagv(iSim)), AtomicChargeEvolution_new, nTimes, tmin, dt, nAtoms)

    end do Sim_loop
    !
    if(.not. SaveDensity)then
        call Write_Summary(OutDir // "/Simulation_Summary", nPts, nAtoms, Volume, Computed_Volume, nTimes, tmin, tmax, AtomName, Radius_BS, nOrb, OrbTab)
    end if
    stop


contains

    ! ################################################
    !... Lindblad Equation - LiouvillianPropagator
    ! ----------------------------------------------------
    !.. Build the time-independent Lindblad superoperator
    !..
    !.. Linblad equation (Schroedinger representation)
    !.. d/dt \rho(t) = -i [H,\rho(t)]+
    !                     \sum_{mn} [ L_{mn} \rho(t) L_{mn}^\dagger -
    !                    -\frac{1}{2} L_{mn}^\dagger L_{mn} \rho(t)
    !                    -\frac{1}{2} \rho(t) L_{mn}^\dagger L_{mn}
    !   H = H_0 + (\hat{\epsilon}\cdot\vec{\mu}) E(t)
    !
    !.. For pure dephasing without relaxation:
    !   L_{nm} = \delta_{mn} \sqrt{\gamma_m/2} | m \rangle \langle m |
    !
    !.. For relaxation without pure dephasing:
    !   L_{mn} = (1-\delta_{nm}) \sqrt{\gamma_{mn}} |m\rangle \langle n|
    !   subject to the constraint \gamma_{mn} = e^{-\beta \omega_{mn}} \gamma_{nm}
    !   where \beta = 1/(k_B T)
    !
    !.. LiouvillianPropagator: The central routine in simulating the dynamics of a quantum system subject to an external field. This routine is named after the concept of the Liouville operator, which in quantum mechanics, describes the time evolution of a quantum system in Liouville space.
    !
    ! Input parameters:
    !
    ! L0_Eval: a complex array holding the eigenvalues of the Liouville operator.
    ! L0_LEvec: a complex array holding the left eigenvectors of the Liouville operator.
    ! L0_REvec: a complex array holding the right eigenvectors of the Liouville operator.
    ! zStatRho: a complex array holding the initial state density matrix of the quantum system. The array is updated within the subroutine, and hence is an input/output parameter.
    ! Important variables used within the function:
    !
    ! nLiou: the size of the L0_Eval array, which represents the dimension of the Liouville space.
    ! EFieldX, EFieldY, EFieldZ: real numbers representing the components of the external electric field at time t + dt/2.
    ! nStates: the size of the zStatRho array, which represents the number of quantum states.
    ! Various complex arrays (zmat1, zmat2, zPureMat, RhoVec, RhoVec2, etc.) and real arrays (Epure, EVEC_DX, EVEC_DY, EVEC_DZ, etc.) for temporary storage and computation.
    ! Key steps in the subroutine:
    !
    !   Transform the initial density matrix zStatRho into Liouville space.
    !   Propagate the transformed matrix in the absence of external fields for a time interval dt/2.
    !   Convert the propagated matrix back into Hilbert space.
    !   Consider the system's interaction with an external electric field by applying rotations in the Hilbert space.
    !   Diagonalize the density matrix to determine the pure states it comprises.
    !   Propagate each individual state under the influence of the external field.
    !   Re-form the density matrix from the propagated pure states.
    !   Repeat steps 1-3 for the re-formed density matrix to complete the propagation process over the time interval dt.
    !   This subroutine heavily relies on various matrix operations and makes use of other subroutines such as HilbertToLiouvilleMatrix, ZGEMV, LiouvilleToHilbertMatrix, DiagonalizeDipole, ZGEMM, Short_ZHEEV, and ExternalElectricFieldCart, which are not defined within this subroutine and are expected to be defined elsewhere.
    !
    ! Important to note is that certain variables used within this subroutine, such as dt, t, iSim, train, Z1, Zi, Z0, nStates, and Dmat, are not defined or passed as arguments, implying that they are either global variables or passed from a higher-level routine that calls this function.
    !
    ! The FIRST_CALL logical variable is used to execute certain lines of code only during the first invocation of the subroutine, primarily for the allocation and initialization of arrays.
    !
    !This subroutine effectively encapsulates a single iteration or time step of the simulation. For a complete simulation over a given time period, this subroutine would typically be called repeatedly in a loop.
    subroutine LiouvillianPropagator(L0_Eval, L0_LEvec, L0_REvec, zStatRho)
        complex(kind(1d0)), intent(in) :: L0_LEvec(:, :), L0_REvec(:, :), L0_Eval(:)
        complex(kind(1d0)), intent(inout) :: zStatRho(:, :)

        logical, save :: FIRST_CALL = .TRUE.
        integer :: nLiou
        integer :: iState, j
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
        do iState = 1, nStates
            do j = 1, nStates
                zStatRho(:, j) = zStatRho(:, j) + EPure(iState) * zPureMat(:, iState) * conjg(zPureMat(j, iState))
            enddo
        enddo
        call HilbertToLiouvilleMatrix(zStatRho, RhoVec)

        !.. Second half of the field-free propagation in Liouville Space across an interval dt/2
        call ZGEMV("C", nLiou, nLiou, Z1, L0_LEvec, nLiou, RhoVec, 1, Z0, RhoVec2, 1)
        RhoVec2 = RhoVec2 * exp(-Zi * L0_Eval * dt / 2.d0)
        call ZGEMV("N", nLiou, nLiou, Z1, L0_REvec, nLiou, RhoVec2, 1, Z0, RhoVec, 1)
        call LiouvilleToHilbertMatrix(RhoVec, zStatRho)
    end subroutine LiouvillianPropagator
    !.. ComputeLiouvillian_0:  Calculates the time-independent Lindblad superoperator (Liouvillian), a crucial element in the Lindblad master equation, used for describing the time evolution of the density matrix of a quantum system in an open quantum system framework.
    !
    ! Here's a breakdown of the subroutine:
    !
    ! Parameters:
    ! Evec (input, real array): The array of energy eigenvalues of the system's Hamiltonian.
    ! Dmat (input, real array): The dipole moment matrix of the system.
    ! Liouvillian0 (output, complex array): The computed Liouvillian superoperator. It's an output parameter that's updated in the subroutine.
    ! BATH_TEMPERATURE (input, real): The temperature of the bath or environment interacting with the quantum system.
    ! dephasing_factor (input, real): The factor representing dephasing effects in the system.
    ! RELAXATION_FACTOR (input, real): The factor representing relaxation effects in the system.
    ! Important Variables:
    ! nLiou and nStates: The size of the Liouvillian and the number of quantum states in the system, respectively.
    ! PairGamma, TotalGamma: The arrays containing the pair relaxation rates and the total relaxation rates, respectively.
    ! BETA: The inverse temperature factor (1/kT) with k being the Boltzmann constant and T being the temperature.
    ! Detailed Steps:
    ! First, the subroutine initializes parameters and allocates space for the Liouvillian.
    ! The unitary part of the Liouvillian is computed, based on the system Hamiltonian eigenvalues.
    ! The subroutine calculates the pair relaxation rates and the dephasing rates for all pairs of states, taking into account the energy levels, dipole moment matrix, and input parameters related to dephasing and relaxation.
    ! The total relaxation rate for each state is computed as the sum of the pair relaxation rates involving the state.
    ! Finally, the dissipative part of the Liouvillian, capturing the effects of the interaction of the system with its environment, is computed based on these relaxation rates.
    ! Please note that this subroutine computes the Lindblad superoperator for a specific case where the relaxation and dephasing rates are assumed to be constants and the Lindblad operators are assumed to be proportional to the dipole moment operator. In more complex or specific cases, modifications would be required.
    subroutine ComputeLiouvillian_0(Evec, Dmat, Liouvillian0, BATH_TEMPERATURE, dephasing_factor, RELAXATION_FACTOR)
        real   (kind(1d0)), intent(in) :: Evec(:)
        real   (kind(1d0)), intent(in) :: Dmat(:, :, :)
        complex(kind(1d0)), allocatable, intent(out) :: Liouvillian0(:, :)

        real(kind(1d0)), intent(in) :: dephasing_factor
        real(kind(1d0)), intent(in) :: RELAXATION_FACTOR
        real(kind(1d0)), intent(in) :: BATH_TEMPERATURE

        integer :: nLiou, nStates, i, j, iLiou, i1, i2, iLiou1, iLiou2
        !.. Dephasing and Relaxation Constants
        real   (kind(1d0)), allocatable :: PairGamma(:, :), TotalGamma(:)
        real   (kind(1d0)) :: AverageDipole, dBuf
        real(kind(1d0)) :: BETA

        BETA = 1.d0 / (BOLTZMANN_CONSTANT_auK * BATH_TEMPERATURE)

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
                PairGamma(i, j) = RELAXATION_FACTOR * AverageDipole
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

    end subroutine ComputeLiouvillian_0
    !
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Hilbert and Liouville space transformations subroutines
    ! ----------------------------------------------------
    !
    !*************************************************************************************************
    !.. DiagonalizeLindblad_0: This subroutine performs the diagonalization of the time-independent
    ! Lindblad superoperator (Liouvillian0). The diagonalization is performed using the
    ! Short_Diag subroutine, and the results are printed to the standard output.
    !
    ! This routine also checks that the product of right and conjugate transpose of left eigenvectors
    ! is an identity matrix (a necessary condition for a valid diagonalization).
    !
    ! Input:
    ! Liouvillian0 - A complex square matrix representing the Lindblad superoperator to be diagonalized.
    !
    ! Output:
    ! L0_Eval - A complex vector that will hold the eigenvalues of Liouvillian0.
    ! L0_LEvec - A complex square matrix that will hold the left eigenvectors of Liouvillian0.
    ! L0_REvec - A complex square matrix that will hold the right eigenvectors of Liouvillian0.
    !
    ! Note: The input matrix Liouvillian0 is also modified during the execution of this subroutine.
    !*************************************************************************************************
    subroutine DiagonalizeLindblad_0(Liouvillian0, L0_Eval, L0_LEvec, L0_REvec)
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
    end subroutine DiagonalizeLindblad_0
    !
    !*************************************************************************************************
    !..  HilbertToLiouvilleMatrix: This subroutine transforms a matrix from Hilbert space into a vector
    ! in Liouville space. The transformation is performed by looping through each element of the input
    ! matrix and mapping it to the corresponding position in the output vector using the subroutine
    ! HilbertToLiouvilleIndexes. This transformation is required for the application of the Liouville operator.
    !
    ! Input:
    ! RhoMat - A complex matrix representing the state in Hilbert space.
    !
    ! Output:
    ! RhoVec - A complex vector representing the transformed state in Liouville space.
    !*************************************************************************************************
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
    !
    !*************************************************************************************************
    ! LiouvilleToHilbertMatrix: This subroutine performs the inverse transformation of
    ! HilbertToLiouvilleMatrix, transforming a vector in Liouville space back into a matrix in
    ! Hilbert space. The transformation is performed by looping through each element of the output matrix
    ! and assigning it the corresponding value from the input vector using the subroutine
    ! HilbertToLiouvilleIndexes.
    !
    ! Input:
    ! RhoVec - A complex vector representing the state in Liouville space.
    !
    ! Output:
    ! RhoMat - A complex matrix representing the transformed state in Hilbert space.
    !*************************************************************************************************
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
    !
    !*************************************************************************************************
    ! HilbertToLiouvilleIndexes: This subroutine transforms matrix indices in Hilbert space into a
    ! single index in Liouville space. The transformation is dependent on the relative values of the
    ! input indices. The resulting single index is required for the transformation from Hilbert to
    ! Liouville space.
    !
    ! Input:
    ! i, j - Integers representing the indices in Hilbert space.
    !
    ! Output:
    ! iPair - An integer representing the corresponding single index in Liouville space.
    !*************************************************************************************************
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
    !
    !*************************************************************************************************
    ! LiouvilleToHilbertIndexes: This subroutine takes in a single integer index in Liouville space
    ! and converts it into a pair of matrix indices in Hilbert space. The indices are calculated by
    ! decomposing the given Liouville index into a square plus a remainder, which are then used to
    ! calculate the corresponding Hilbert indices.
    !
    ! Input:
    ! iPair - An integer representing a single index in Liouville space.
    !
    ! Output:
    ! i, j - Integers representing the corresponding indices in Hilbert space.
    !*************************************************************************************************
    subroutine LiouvilleToHilbertIndexes(iPair, i, j)
        ! Input integer (Liouville space index)
        integer, intent(in) :: iPair

        ! Output integers (Hilbert space indices)
        integer, intent(out) :: i, j

        ! Intermediate variables to hold the square root and square of iPair
        integer :: n, n2

        ! Calculate the integer square root of iPair
        n = int(sqrt(1.d0 * iPair))

        ! Calculate the square of the integer square root
        n2 = n**2

        ! Check if iPair is a perfect square
        if(n2==iPair)then
            i = n
            j = n

            ! Check if the difference between iPair and n2 is less than or equal to n
        elseif(iPair - n2<=n)then
            i = iPair - n2
            j = n + 1

            ! Otherwise, the difference is greater than n
        else
            i = n + 1
            j = iPair - n2 - n
        endif
    end subroutine LiouvilleToHilbertIndexes
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Called ... rn not going to organize these yet
    ! ----------------------------------------------------
    subroutine DiagonalizeDipole(Dmat, EVEC_DX, zUMAT_DX, EVEC_DY, zUMAT_DY, EVEC_DZ, zUMAT_DZ)
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
    subroutine ComputeVolume(nPts, Volume, gridv)
        integer, intent(in) :: nPts
        real(kind(1d0)), intent(out) :: Volume
        real(kind(1d0)), allocatable, intent(in) :: gridv(:, :)
        real(kind(1d0)) :: coord1, coord2, delta
        real(kind(1d0)), parameter :: threshold = 1d-10

        integer :: iPts, iCoord
        Volume = 1.d0

        do iCoord = 1, 3
            coord1 = gridv(iCoord, 1)
            delta = 1d10
            do iPts = 2, nPts
                coord2 = gridv(iCoord, iPts)
                if (abs(coord1 - coord2)>threshold) then
                    delta = min(delta, abs(coord1 - coord2))
                endif
            enddo
            Volume = Volume * delta
        enddo

    end subroutine ComputeVolume
    subroutine Compute_R_el(gridv, WeightV, OrbTab, R_el) !!!***** refactor writing to Module RTP, hard coded... *****
        real(kind(1d0)), allocatable, intent(in) :: OrbTab (:, :), gridv(:, :), WeightV(:, :)
        real(kind(1d0)), allocatable, intent(out) :: R_el(:, :)

        real(kind(1d0)) :: sum, m, r, sum1
        integer :: nAtoms, nPts, nOrbs, iOrb, jOrb, uid
        nAtoms = size(WeightV, 2)
        nPts = size(WeightV, 1)
        nOrbs = size(OrbTab, 2)

        !
        write(*, "(a)") "Computing Barycenters of the Atomic Charges"
        !
        allocate(R_el(3, nAtoms))

        do iPol = 1, 3
            do iAtom = 1, nAtoms
                sum = 0.d0
                sum1 = 0.d0
                do iOrb = 1, nOrb
                    do jOrb = 1, nOrb
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

        !        open(newunit = uid, &
        !                file = "R_el_bc", & !!!!!!! hard coded
        !                form = "formatted", &
        !                status = "unknown")
        !        do iAtom = 1, nAtoms
        !            write(uid, "(*(x,e24.14e3))") (R_el(iPol, iAtom), iPol = 1, 3)
        !        end do
        !        close(uid)
    end subroutine Compute_R_el
    subroutine AtomicRadius_Bragg_Slater_Becke(AtomName, nAtom, Radius_BS)
        real(kind(1d0)), allocatable, intent(out) :: Radius_BS(:)
        integer, intent(in) :: nAtom
        character(len = 16), intent(in) :: AtomName(:)
        integer :: iAtom
        allocate(Radius_BS(nAtom))

        do iAtom = 1, nAtom
            Radius_BS(iAtom) = Radius_Table(AtomName(iAtom))
        end do
    end subroutine AtomicRadius_Bragg_Slater_Becke
    subroutine ComputeNewAtomicCharges(OrbitalDensity, Becke_new, QchargeVec_new)
        real(kind(1d0)), intent(in) :: OrbitalDensity(:, :)
        real(kind(1d0)), intent(in) :: Becke_new(:, :, :, :)
        real(kind(1d0)), allocatable, intent(out) :: QchargeVec_new    (:, :)

        integer :: iAtom, nAtoms, nOrbs, jOrb, iOrb, nPts, iPol
        real(kind(1d0)), external :: DDOT
        real(kind(1d0)) :: dsum
        nPts = size(Becke_new, 1)
        nOrbs = size(Becke_new, 2)
        nAtoms = size(Becke_new, 4)
        allocate(QchargeVec_new(3, nAtoms))

        QchargeVec_new = 0.d0
        do iPol = 1, 3
            do iAtom = 1, nAtoms
                do jOrb = 1, nOrbs
                    do iOrb = 1, nOrbs
                        QchargeVec_new(iPol, iAtom) = QchargeVec_new(iPol, iAtom) + OrbitalDensity(iOrb, jOrb) &
                                * Becke_new(iPol, iOrb, jOrb, iAtom)
                        !                QchargeVec_new(iPol, iAtom) = ddot(nOrbs * nOrbs, OrbitalDensity, 1, Becke_new(iPol, 1, 1, iAtom), 1)
                    end do
                end do
            enddo
        end do
    end subroutine ComputeNewAtomicCharges
    subroutine ComputeOrbitalDensity(zStatRho, TDM, Amat)
        complex(kind(1d0)), intent(in) :: zStatRho(:, :)
        real   (kind(1d0)), intent(in) :: TDM(:, :, :, :)
        real   (kind(1d0)), allocatable, intent(out) :: Amat  (:, :)
        integer :: nOrb, iOrbi, iOrbj, iStatei, iStatej

        nOrb = size(TDM, 1)

        if(allocated(Amat))then
            if(abs(size(Amat, 1) - nOrb) + abs(size(Amat, 2) - nOrb)>0)then
                deallocate(Amat)
                allocate(Amat(nOrb, nOrb))
            endif
        else
            allocate(Amat(nOrb, nOrb))
        endif
        !
        !.. Build the expansion Matrix
        !   $A_{nm}(t) = \sum_{IJ} P_{IJ}(t) \rho^{JI}_{nm}$
        !   where $P_{IJ}(t)$ is the solution of a suitable Master equation, starting from $P_{IJ}(0)$
        !   Here, we will neglect relaxation and decoherence and hence $P_{IJ}(t) = P_{IJ}(0) e^{-i\omega_{IJ}t}$
        !   where $\omega_{IJ}=E_I-E_J$
        Amat = 0.d0
        do iOrbi = 1, nOrb
            do iOrbj = 1, nOrb
                !
                do iStatei = 1, nStates
                    do iStatej = 1, nStates
                        Amat(iOrbi, iOrbj) = Amat(iOrbi, iOrbj) + &
                                dble(zStatRho(iStatei, iStatej) * TDM(iOrbi, iOrbj, iStatej, iStatei))
                    enddo
                enddo
                !
            enddo
        enddo

    end subroutine ComputeOrbitalDensity
    subroutine ComputeAtomicWeights(nPts, gridv, nAtoms, AtCoord, WeightV, Radius_BS)
        integer, intent(in) :: nPts, nAtoms
        real(kind(1d0)), intent(in) :: gridv(:, :), AtCoord(:, :), Radius_BS(:)
        real(kind(1d0)), allocatable, intent(out) :: WeightV(:, :)

        integer, parameter :: k = 4
        integer :: iPts, iAtom
        real(kind(1d0)) :: rvec(3)

        write(*, "(a)") "Computing Weights"

        allocate(WEIGHTV(nPts, nAtoms))
        do iPts = 1, nPts
            rvec = gridv(:, iPts)
            do iAtom = 1, nAtoms
                WeightV(iPts, iAtom) = wkfun(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS)
            enddo
        enddo
    end subroutine ComputeAtomicWeights
    subroutine TabulateChargeDensity(Amat, OrbTab, ChDen)
        real(kind(1d0)), intent(in) :: Amat(:, :)
        real(kind(1d0)), intent(in) :: OrbTab(:, :)
        real(kind(1d0)), intent(out) :: ChDen(:)

        integer :: nOrb, iOrbi, iOrbj

        nOrb = size(OrbTab, 2)
        !.. Tabulate Charge Density
        ChDen = 0.d0
        do iPts = 1, nPts
            do iOrbi = 1, nOrb
                do iOrbj = 1, nOrb
                    ChDen(iPts) = ChDen(iPts) + OrbTab(iPts, iOrbi) * Amat(iOrbi, iOrbj) * OrbTab(iPts, iOrbj)
                enddo
            enddo
        enddo
    end subroutine TabulateChargeDensity
    subroutine ComputeNewBeckeMatrix(WeightV, OrbTab, Becke_new, Bary_center)
        real(kind(1d0)), intent(in) :: WeightV(:, :)
        real(kind(1d0)), intent(in) :: OrbTab (:, :)
        real(kind(1d0)), allocatable, intent(out) :: Becke_new(:, :, :, :)
        real(kind(1d0)), intent(in) :: Bary_center(:, :)

        integer :: nPts, nOrbs, nAtoms
        integer :: iPts, iOrb1, iOrb2, iAtom, iPol
        real(kind(1d0)) :: dsum

        !
        write(*, "(a)") "Computing Becke's Matrix"
        !
        nAtoms = size(WeightV, 2)
        nPts = size(WeightV, 1)
        nOrbs = size(OrbTab, 2)

        !        if(allocated(Becke_new))deallocate(Becke_new)
        allocate(Becke_new(3, nOrbs, nOrbs, nAtoms))

        Becke_new = 0.d0
        do iPol = 1, 3
            do iAtom = 1, nAtoms
                do iOrb2 = 1, nOrbs
                    do iOrb1 = 1, nOrbs
                        do iPts = 1, nPts
                            Becke_new(iPol, iOrb1, iOrb2, iAtom) = Becke_new(iPol, iOrb1, iOrb2, iAtom) + &
                                    WeightV(iPts, iAtom) * OrbTab(iPts, iOrb1) * OrbTab(iPts, iOrb2) * &
                                            Bary_center(iPol, iAtom)

                        enddo
                    enddo
                enddo
            enddo
        end do
    end subroutine ComputeNewBeckeMatrix
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Not used in the code ...?
    ! ----------------------------------------------------
    subroutine ComputeAtomicCharges(OrbitalDensity, BeckeMatrix, QchargeVec)
        real(kind(1d0)), intent(in) :: OrbitalDensity(:, :)
        real(kind(1d0)), intent(in) :: BeckeMatrix   (:, :, :)
        real(kind(1d0)), allocatable, intent(out) :: QchargeVec    (:)
        integer :: iAtom, nAtoms, nOrbs, jOrb, iOrb
        real(kind(1d0)), external :: DDOT
        nOrbs = size(BeckeMatrix, 1)
        nAtoms = size(BeckeMatrix, 3)
        if(.not.allocated(QchargeVec))allocate(QchargeVec(nAtoms))
        QchargeVec = 0.d0
        do iAtom = 1, nAtoms
            !            do jOrb = 1, nOrbs
            !                do iOrb = 1, nOrbs
            !                    QchargeVec(iAtom) = QchargeVec(iAtom) + OrbitalDensity(iOrb, jOrb) * BeckeMatrix(iOrb, jOrb, iAtom)
            !                end do
            !            end do
            QchargeVec(iAtom) = ddot(nOrbs * nOrbs, OrbitalDensity, 1, BeckeMatrix(1, 1, iAtom), 1)
        enddo
    end subroutine ComputeAtomicCharges
    subroutine ComputeDipole_new(Amat, OrbTab, Dipole_new_, gridv)
        real(kind(1d0)), intent(in) :: Amat(:, :)
        real(kind(1d0)), intent(in) :: OrbTab (:, :)
        real(kind(1d0)), intent(in) :: gridv(:, :)
        real(kind(1d0)), allocatable, intent(out) :: Dipole_new_(:)
        integer :: nPts, nOrbs
        integer :: iPts, iOrb, jOr

        nPts = size(OrbTab, 1)
        nOrbs = size(OrbTab, 2)

        if (.not.allocated(Dipole_new_))allocate(Dipole_new_(3))
        Dipole_new_ = 1
        do iPol = 1, 3
            do iOrb = 1, nOrb
                do jOrb = 1, nOrb
                    do iPts = 1, nPts
                        Dipole_new_(iPol) = Dipole_new_(iPol) + OrbTab(iPts, iOrb) * OrbTab(iPts, jOrb) * Amat(iOrb, jOrb) * gridv(iPol, iPts)
                    enddo
                enddo
            enddo
        end do

    end subroutine ComputeDipole_new
    complex(kind(1d0)) function zTraceFunction(zA) result(zTrace)
        integer :: i
        complex(kind(1d0)), intent(in) :: zA(:, :)
        zTrace = 0.d0
        do i = 1, size(zA, 1)
            zTrace = zTrace + zA(i, i)
        enddo
    end function zTraceFunction
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! FUNCTIONS ... **notes can be found right below, while  documentation is present right above each function.**
    !Radius_Table: Given the name of an atom, returns its atomic radius. Only supports H, C, N, O, and Mg.
    !get_param_a: Given two radii R_i and R_j, returns a transformed parameter a_ij. If the absolute value of a_ij is greater than 0.5, it's set to 0.5.
    !new_mu_transformation: Returns a new mu parameter, given the original mu and a_ij.
    !wkfun: Calculates a weighted function for atom i.
    !PkFunTot: Sums up the function Pkfuna over all atoms.
    !Pkfuna: Calculates a product of skfunab over all atoms excluding atom i.
    !skfunab: Transforms an elliptical coordinate to a new mu, then returns the result of skfun with this new mu.
    !skfun: Transforms a mu parameter into a step function.
    !fkfun: Recursively calculates a function of mu.
    !EllipticalCoord: Given coordinates rvec, avec, and bvec, calculates an elliptical coordinate.
    !EuclDist: Given two points in space, calculates their Euclidean distance.
    !The main considerations:
    !This code does not handle the case where the input atomic name does not match one of the provided options in Radius_Table. This could lead to an error or unexpected behavior.
    !The recursion in fkfun does not seem to have a base case when k is less than 1, which could potentially lead to a stack overflow error.
    !There are no checks for zero or negative distances in EuclDist which might cause issues if not handled properly upstream.
    !The input values are not validated for possible errors. Adding error handling and input validation could make the code more robust.
    !The functions are tightly coupled, which might make modifications and debugging more difficult. Consider making them more modular.
    ! ----------------------------------------------------
    !Result =P_a(\vec{r})/\sum_bP_b(\vec{r}) $
    ! ----------------------------------------------------
    !.. Euclidean distance
    ! This function calculates the Euclidean distance between two points in 3D space. The inputs are two 3-component vectors representing the coordinates of the points, and the output is a scalar representing the distance between them.
    real(kind(1d0)) function EuclDist(avec, bvec) result(res)
        real(kind(1d0)), intent(in) :: avec(3), bvec(3)
        real(kind(1d0)) :: cvec(3)
        cvec = (avec - bvec)
        res = sqrt(sum(cvec * cvec))
    end function EuclDist
    ! ----------------------------------------------------
    !.. Eq. 1 $\mu_{ab}(\vec{r})$
    ! This function calculates the elliptical coordinate ($\mu_{ab}$) of a point relative to two other points (defined as atoms 'a' and 'b'). This coordinate is a measure of how much closer the point is to atom 'a' compared to atom 'b'.
    real(kind(1d0)) function EllipticalCoord(rvec, avec, bvec) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        real(kind(1d0)) :: r_a, r_b, R_ab ! $r_a$, $r_b$, $R_ab$
        r_a = EuclDist(rvec, avec)
        r_b = EuclDist(rvec, bvec)
        R_ab = EuclDist(avec, bvec)
        res = (r_a - r_b) / R_ab
    end function EllipticalCoord
    ! ----------------------------------------------------
    !..Radius_Table for atomic radius based on atomic name
    ! This function provides the atomic radius based on the atomic name. The input is a string representing the atomic name and the output is the corresponding atomic radius. This is useful for defining the size of atoms in molecular modelling.
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
    ! ----------------------------------------------------
    !..get_param_aij
    ! This function calculates the parameter 'a' (a_ij) which depends on the radii of two atoms. This parameter is used in the transformation of the elliptical coordinate.
    real(kind(1d0)) function get_param_a(R_i, R_j) result(a_ij)
        real(kind(1d0)), intent(in) :: R_i, R_j
        real(kind(1d0)) :: Chi
        Chi = R_i / R_j
        a_ij = (Chi - 1) / (Chi + 1)
        if (abs(a_ij)>(0.5d0)) then
            a_ij = 0.5d0
        endif
    end function get_param_a
    ! ----------------------------------------------------
    !..new_mu_transformation
    ! This function transforms the elliptical coordinate using the parameter 'a' (a_ij). This transformation is used to handle cases where the point is much closer to one atom than the other.
    real(kind(1d0)) function new_mu_transformation(mu, a_ij) result(res)
        real(kind(1d0)), intent(in) :: mu, a_ij
        res = mu + a_ij * (1 - mu**2)
    end function new_mu_transformation
    ! ----------------------------------------------------
    !.. Eq. 4 Recursive $f_k(\mu)$
    ! This function calculates a recursive function of the transformed elliptical coordinate. The output of this function will be used in the calculation of the step function skfun.
    real(kind(1d0)) recursive function fkfun(mu, k) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu
        integer, intent(in) :: k
        res = 0.5d0 * mu * (3.d0 - mu**2)
        if(k==1)return
        res = fkfun(res, k - 1)
    end function fkfun
    ! ----------------------------------------------------
    !.. $s(\mu): [-1,1]  \rightarrow [0,1]$  (Step Function)
    ! This function calculates a step function of the transformed elliptical coordinate. This function maps the transformed elliptical coordinate from the range [-1,1] to [0,1].
    real(kind(1d0)) function skfun(mu, k) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu
        integer, intent(in) :: k
        real(kind(1d0)) :: r_m

        res = 0.5d0 * (1.d0 - fkfun(mu, k))
    end function skfun
    ! ----------------------------------------------------
    !.. This function combines the calculations of the transformed elliptical coordinate and the step function. The output of this function will be used in the calculation of the partial partition function Pkfuna.
    real(kind(1d0)) function skfunab(rvec, avec, bvec, k, a_ij) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        integer, intent(in) :: k
        real(kind(1d0)) :: mu_old, mu_new, a_ij
        mu_old = EllipticalCoord(rvec, avec, bvec)
        mu_new = new_mu_transformation(mu_old, a_ij)
        res = skfun(mu_new, k)
    end function skfunab
    ! ----------------------------------------------------
    !.. Eq. 2 P_a(\vec{r}) in and nominator of Eq. 3
    ! This function calculates the partial partition function for a specific atom. The output of this function represents the contribution of a specific atom to the total partition function.
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
    ! ----------------------------------------------------
    !..Eq. 3 \sum_bP_b(\vec{r}) in the denominator
    ! This function calculates the total partition function. This is done by summing the partial partition functions for all atoms.
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
    ! ----------------------------------------------------
    !.. Eq. 3 $w_a \vec(r)$  (WEIGHTS)
    ! This function calculates the weight of an atom. The weight is defined as the ratio of the partial partition function of the atom to the total partition function. This weight represents the relative importance of the atom.
    real(kind(1d0)) function wkfun(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS) result(res)
        !
        integer, intent(in) :: iAtom, k, nAtoms
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms), Radius_BS(:)
        !
        res = Pkfuna(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS) / PkfunTot(rvec, AtCoord, nAtoms, k, Radius_BS)
    end function wkfun
    ! ----------------------------------------------------
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


end program ChargeMigration


