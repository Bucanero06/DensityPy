! {{{ Detailed description

!> \mainpage Program <ProgramName> <Insert here what the program does>
!!
!! Synopsis: <Program Name> <mandatory run-time parameters (RTP)> [<optional RTP>]
!! ---------
!!
!! ChargeMigration simulates charge migration in molecular systems under external fields.
!!
!! ___
!! Description:
!! ------------
!!
!!
!!
!! Key Subroutines:
!!
!! LiouvillianPropagator: This central routine simulates the dynamics of a quantum system subject to an external field. It updates the state density matrix zStatRho using the eigenvalues L0_Eval and eigenvectors L0_LEvec, L0_REvec of the Liouvillian operator.
!!
!! ComputeLiouvillian_0: Calculates the time-independent Lindblad superoperator (Liouvillian), crucial in describing the time evolution of the density matrix in an open quantum system. It uses system energy eigenvalues Evec, dipole moment matrix Dmat, and parameters like bath_temperature, dephasing_factor, and relaxation_factor.
!!
!! DiagonalizeLindblad_0: Performs diagonalization of the Lindblad superoperator Liouvillian0, determining its eigenvalues L0_Eval and eigenvectors L0_LEvec and L0_REvec.
!!
!! HilbertToLiouvilleMatrix and LiouvilleToHilbertMatrix: These subroutines transform between Hilbert space matrices and Liouville space vectors. They are crucial for applying the Liouville operator in the simulations.
!!
!! ComputeOrbitalDensity: Builds the orbital density matrix Amat from the state density matrix zStatRho and transition density matrices TDM.
!!
!! ComputeAtomicCharges: Calculates new atomic charges QchargeVec from the orbital density OrbitalDensity and Becke matrix BeckeMatrix.
!!
!! TabulateChargeDensity: Generates charge density ChDen from the orbital density matrix Amat and orbital table OrbTab.
!!
!! ComputeBeckeMatrix: Constructs the Becke matrix BeckeMatrix using weights WeightV, orbital table OrbTab, and barycenter coordinates Bary_center.
!!
!! ComputeAtomicWeights: Calculates atomic weights WeightV for a given set of points gridv, atomic coordinates AtCoord, and radii Radius_BS.
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
    use Module_Becke
    use omp_lib

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

    integer :: npts, nAtoms
    real(kind(1d0)), allocatable :: gridv(:, :), AtCoord(:, :) ! 3 x npts
    real(kind(1d0)), allocatable :: OrbTab(:, :) ! 3 x npts
    !    real(kind(1d0)) :: volume

    !.. Statistical Density Matrix
    complex(kind(1d0)), allocatable :: zStatRho(:, :)
    real   (kind(1d0)) :: t, dt
    real   (kind(1d0)), allocatable :: ChDen(:), WEIGHTV(:, :)
    real   (kind(1d0)), allocatable :: AtomicChargeVec(:, :), AtomicChargeEvolution(:, :, :)

    character(len = 30) :: istrn
    character(len = 16), allocatable :: atom_names(:)

    integer :: iPts, it, iPol, iAtom

    complex(kind(1d0)), allocatable :: Liouvillian0(:, :), L0_LEvec(:, :), L0_REvec(:, :), L0_Eval(:)

    !.. Pulse parameters
    integer :: iSim, N_Simulations, uid
    character(len = 64), pointer :: Simulation_Tagv(:)
    character(len = 1000) :: strn
    type(pulse_train), pointer :: train(:)
    !..
    !  B_{ij}^{\alpha} = \int d^3r \phi_i(\vec{r})\phi_j(\vec{r}) w_\alpha(\vec{r}) = BeckeMatrix(i,j,alpha)
    real(kind(1d0)), allocatable :: BeckeMatrix(:, :, :, :)
    real(kind(1d0)), allocatable :: OrbitalDensity(:, :)
    real(kind(1d0)), allocatable :: Radius_BS(:)
    complex(kind(1d0)), external :: zdotu

    !..Test
    real(kind(1d0)) :: Computed_volume
    integer :: iOrb, jOrb
    real(kind(1d0)), allocatable :: QchargeVec(:, :), R_el(:, :)

    call GetRunTimeParameters(input_directory, output_directory, molecular_geometry_file, &
            n_times, t_min, t_max, Ext_field_file, Verbous, Weight_File, read_precomputed_weights_flag, &
            save_charge_migration_flag, ivorb, counted_number_of_orbitals, volume, dephasing_factor, relaxation_factor, bath_temperature)

    call system("mkdir -p " // trim(output_directory))
    call system("mkdir -p " // output_directory // "/Dipole")
    call system("mkdir -p " // output_directory // "/AtomicCharge")
    call system("mkdir -p " // output_directory // "/ChargeDensity")
    call system("mkdir -p " // output_directory // "/Pulses")
    call Set_CD_IO_Verbous(Verbous)

    open(newunit = uid, &
            file = Ext_field_file, &
            form = "formatted", &
            status = "old")
    call Parse_Simulation_File(uid, OUTPUT_UNIT, N_Simulations, Simulation_Tagv, train)
    close(uid)
    write(*, *) "N_Simulations=", N_Simulations

    !.. Write and Print the pulse
    !..
    dt = (t_max - t_min) / dble(n_times)
    write(*, *) "fortran t_min=", t_min
    write(*, *) "fortran t_min=", t_min
    write(*, *) "fortran t_max=", t_max
    write(*, *) "fortran n_times=", n_times
    write(*, *) "fortran dble(n_times - 1)=", dble(n_times)
    write(*, *) "fortran dt=", dt

    call omp_set_num_threads(30) ! fixme do not hard code, use number of passed or allowed threads

    !$OMP PARALLEL DO PRIVATE(strn, i)
    do i = 1, N_Simulations
        write(*, *) "Writing pulse " // trim(Simulation_Tagv(i))
        strn = trim(output_directory) // "/Pulses/pulse" // trim(Simulation_Tagv(i))
        call train(i)%Write(strn, t_min, t_max, dt)  !>>pulse_trainWrite
        strn = trim(output_directory) // "/Pulses/FTpulse" // trim(Simulation_Tagv(i))
        call train(i)%WriteFTA(strn) !>>pulse_trainWriteFTA
    end do
    !$OMP END PARALLEL DO
    !

    call LoadEnergies(input_directory // "/ROOT_ENERGIES", nStates, Evec)

    !.. At the moment, we assume that the TDM are defined as
    !   $\rho^{JI}_{nm} = \langle J | a_n^\dagger a_m | I \rangle $
    call LoadTDMs    (input_directory // "/DENSITY_MATRIX", input_directory // "/TRANSITION_DENSITY_MATRIX", nStates, number_of_orbitals, TDM)

    !.. Load Dipole matrix elements $\mu_{IJ}$
    call LoadDipoleME(Dmat, input_directory, nStates)

    allocate(zDmat_t(nStates, nStates, 3))
    do i = 1, 3
        zDmat_t(:, :, i) = Z1 * transpose(Dmat(:, :, i))
    enddo

    if (number_of_orbitals .ne. counted_number_of_orbitals) then
        ! TODO change to assert method of ASTRA
        write(*, *) "Number of orbitals given by '-iOrb' do not match number of orbitals in active the active space"
        stop
    endif

    !.. Load Geometry,volume, Grid and Orbitals
    call LoadGeometry(nAtoms, AtCoord, molecular_geometry_file, atom_names)
    call LoadGrid(input_directory // "/gridcoord.csv", npts, gridv)! $G=\{\vec{r}_i, i=1,\ldots,N_{G}\}$
    call Computevolume(nPts, Computed_volume, gridv) ! Should probably check that this computed volume is the same as the one given in the input file

    !.. Load Orbitals
    call LoadOrbitals(input_directory, number_of_orbitals, npts, OrbTab)! $\varphi_n(\vec{r}_i)$

    !.. Compute and Write | or Read Becke's Weights
    !
    call AtomicRadius_Bragg_Slater_Becke(atom_names, nAtoms, Radius_BS)
    if (read_precomputed_weights_flag == .True.) then
        call Read_Weights(Weight_File, WEIGHTV, nAtoms, nPts)
    else
        call ComputeAtomicWeights(nPts, gridv, nAtoms, AtCoord, WeightV, Radius_BS)
        call Write_Weights(output_directory // "/" // Weight_File // "_" // output_directory // ".csv", &
                WEIGHTV, gridv, nAtoms, nPts, atom_names)
    endif
    !.. Compute Berycenter of Atmoic Charges
    call Compute_R_el(gridv, WeightV, OrbTab, R_el)
    call Write_R_el_bc(output_directory, atom_names, nAtoms, R_el)

    !.. Compute Becke's Matrix
    call ComputeNewBeckeMatrix(WeightV, OrbTab, Becke_new, R_el) !!using electronic barycenter or AtCoord for the nuclear barycenter
    !    BeckeMatrix = BeckeMatrix * Computed_volume
    write(*, *) "Becke_new", sum(Becke_new(1, :, :, 1))
    write(*, *) "Becke_new", sum(Becke_new(2, :, :, 1))
    write(*, *) "Becke_new", sum(Becke_new(3, :, :, 1))

    write(*, *) "allocating Atomic Charge matrixes"
    allocate(AtomicChargeEvolution(3, nAtoms, n_times))

    !.. Compute Liouvillian and Diagonalizes it
    write(*, *) "Computing Liouvillian"
    call ComputeLiouvillian_0(Evec, Dmat, Liouvillian0, bath_temperature, dephasing_factor, relaxation_factor)
    call DiagonalizeLindblad_0(Liouvillian0, L0_Eval, L0_LEvec, L0_REvec)

    allocate(zStatRho (nStates, nStates))
    allocate(ChDen(nPts))
    allocate(zMuEV(3, n_times))
    zMuEV = Z0

    write(*, *) "Starting Sim Loop"
    Sim_loop : do iSim = 1, N_Simulations
        ! Lets keep a percentage write of the simulation only 2 digits after the decimal point
        write(*, '(A, F0.2, A, I0, A)') "Simulation progress: ", 100.d0 * iSim / (N_Simulations)
        if(save_charge_migration_flag)then
            write(*, *) "Computing Charge Migration for simulation ", iSim, trim(Simulation_tagv(iSim))
        end if

        !.. Load Initial State
        zStatRho = Z0
        zStatRho(GS_IDX, GS_IDX) = 1.d0

        !.. Time cycle
        time_loop : do it = 1, n_times
            !
            t = t_min + dt * dble(it - 1)

            !.. Computes the excitation density wrt ground state
            zStatRho(GS_IDX, GS_IDX) = zStatRho(GS_IDX, GS_IDX) - 1.d0

            !.. Evaluates the expectation value of the dipole as a function of time
            !            write(*, *) it, iSim, 100.d0 * iSim / (N_Simulations)
            do iPol = 1, 3
                zMuEV(iPol, it) = zdotu(nStates * nStates, zStatRho, 1, zDmat_t(1, 1, iPol), 1)
            enddo
            !
            call ComputeOrbitalDensity(zStatRho, TDM, OrbitalDensity)

            !
            call ComputeAtomicCharges(OrbitalDensity, BeckeMatrix, AtomicChargeVec)
            AtomicChargeVec = AtomicChargeVec * Computed_volume
            AtomicChargeEvolution(:, :, it) = AtomicChargeVec


            !$> this needs to be be moved outside of this module or loop since is costly and can be computed when needed
            if(save_charge_migration_flag)then
                call TabulateChargeDensity(OrbitalDensity, OrbTab, ChDen)
                write(istrn, "(f12.4)")t
                !..Makes new directory inside of ChargeDensity directory for each simulation
                call system("mkdir -p " // output_directory // "/ChargeDensity/ChDenSim" // trim(Simulation_tagv(iSim)))
                call Write_Charge_Density(&
                        output_directory // "/ChargeDensity/ChDenSim" // trim(Simulation_tagv(iSim)) &
                                // "/ChDen" // trim(adjustl(istrn)) // ".csv", &
                        nPts, gridv, ChDen, Weightv, nAtoms, atom_names)
            endif
            !
            zStatRho(GS_IDX, GS_IDX) = zStatRho(GS_IDX, GS_IDX) + 1.d0
            !
            if(it == n_times) exit time_loop
            !
            call LiouvillianPropagator(L0_Eval, L0_LEvec, L0_REvec, zStatRho)
            !
        enddo time_loop

        !.. Save Dipole to filerc
        call Write_Dipole(output_directory // "/Dipole/Dipole" // trim(Simulation_tagv(iSim)) // ".csv", zMuEV, n_times, t_min, dt)

        !.. Save Q_Charge
        call Write_Q_Charge(output_directory // "/AtomicCharge/AtomicCharge" // trim(Simulation_tagv(iSim)) // ".csv", &
                AtomicChargeEvolution, n_times, t_min, dt, nAtoms, atom_names)
    end do Sim_loop
    !
    call Write_Summary(output_directory // "/Simulation_Summary", nPts, nAtoms, volume, Computed_volume, n_times, t_min, t_max, atom_names, Radius_BS, number_of_orbitals, OrbTab)
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
    !*************************************************************************************************
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
    !
    !*************************************************************************************************
    !.. ComputeLiouvillian_0:  Calculates the time-independent Lindblad superoperator (Liouvillian), a crucial element in the Lindblad master equation, used for describing the time evolution of the density matrix of a quantum system in an open quantum system framework.
    !
    ! Here's a breakdown of the subroutine:
    !
    ! Parameters:
    ! Evec (input, real array): The array of energy eigenvalues of the system's Hamiltonian.
    ! Dmat (input, real array): The dipole moment matrix of the system.
    ! Liouvillian0 (output, complex array): The computed Liouvillian superoperator. It's an output parameter that's updated in the subroutine.
    ! bath_temperature (input, real): The temperature of the bath or environment interacting with the quantum system.
    ! dephasing_factor (input, real): The factor representing dephasing effects in the system.
    ! relaxation_factor (input, real): The factor representing relaxation effects in the system.
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
    !*************************************************************************************************
    subroutine ComputeLiouvillian_0(Evec, Dmat, Liouvillian0, bath_temperature, dephasing_factor, relaxation_factor)
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

    end subroutine ComputeLiouvillian_0
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
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Hilbert and Liouville space transformations subroutines
    ! ----------------------------------------------------
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
    ! Called ... rn not going to organize these yet ... figure out what module they go to
    ! ----------------------------------------------------
    !
    !*************************************************************************************************
    ! DiagonalizeDipole: Purpose: Diagonalizes the dipole matrices in the X, Y, and Z directions.
    !
    ! Inputs:
    ! Dmat(:, :, :): Real 3D array representing the dipole moment matrix.
    !
    ! Outputs:
    ! zUMAT_DX(:, :), zUMAT_DY(:, :), zUMAT_DZ(:, :): Complex 2D arrays for the diagonalized matrices in X, Y, and Z directions, respectively.
    ! EVEC_DX(:), EVEC_DY(:), EVEC_DZ(:): Real 1D arrays for the eigenvalues of the dipoles in X, Y, and Z directions.
    !
    ! Internal Operations: Allocation of matrices, diagonalization using Short_Diag, and assignment of the resultant matrices and eigenvalues.
    !*************************************************************************************************
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
    subroutine Compute_R_el(gridv, WeightV, OrbTab, R_el)
        real(kind(1d0)), allocatable, intent(in) :: OrbTab (:, :), gridv(:, :), WeightV(:, :)
        real(kind(1d0)), allocatable, intent(out) :: R_el(:, :)

        real(kind(1d0)) :: sum, m, r, sum1
        integer :: nAtoms, nPts, number_of_orbitals, iOrb, jOrb, uid
        nAtoms = size(WeightV, 2)
        nPts = size(WeightV, 1)
        number_of_orbitals = size(OrbTab, 2)

        !
        write(*, "(a)") "Computing Barycenters of the Atomic Charges"
        !
        allocate(R_el(3, nAtoms))

        do iAtom = 1, nAtoms
            do iPol = 1, 3
                sum = 0.d0
                sum1 = 0.d0
                do iOrb = 1, number_of_orbitals
                    do jOrb = 1, number_of_orbitals
                        do iPts = 1, nPts
                            m = WeightV(iPts, iAtom) * OrbTab(iPts, iOrb) * OrbTab(iPts, jOrb)
                            r = gridv(iPol, iPts)
                            sum = sum + ((m * r))
                            sum1 = sum1 + m
                        end do
                    end do
                end do
                !
                R_el(iPol, iAtom) = sum / sum1
            end do
        end do
    end subroutine Compute_R_el
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
    subroutine ComputeAtomicCharges(OrbitalDensity, BeckeMatrix, QchargeVec)
        real(kind(1d0)), intent(in) :: OrbitalDensity(:, :)
        real(kind(1d0)), intent(in) :: BeckeMatrix(:, :, :, :)
        real(kind(1d0)), allocatable, intent(out) :: QchargeVec    (:, :)

        integer :: iAtom, nAtoms, number_of_orbitals, jOrb, iOrb, iPol
        real(kind(1d0)), external :: DDOT
        real(kind(1d0)) :: dsum
        number_of_orbitals = size(BeckeMatrix, 2)
        nAtoms = size(BeckeMatrix, 4)
        allocate(QchargeVec(3, nAtoms))

        QchargeVec = 0.d0
        do iPol = 1, 3
            do iAtom = 1, nAtoms
                do jOrb = 1, number_of_orbitals
                    do iOrb = 1, number_of_orbitals
                        QchargeVec(iPol, iAtom) = QchargeVec(iPol, iAtom) + OrbitalDensity(iOrb, jOrb) &
                                * BeckeMatrix(iPol, iOrb, jOrb, iAtom)
                    end do
                end do
            enddo
        end do
    end subroutine ComputeAtomicCharges
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
    subroutine ComputeAtomicWeights(nPts, gridv, nAtoms, AtCoord, WeightV, Radius_BS)
        integer, intent(in) :: nPts, nAtoms
        real(kind(1d0)), intent(in) :: gridv(:, :), AtCoord(:, :), Radius_BS(:)
        real(kind(1d0)), allocatable, intent(out) :: WeightV(:, :)

        integer, parameter :: k = 4
        integer :: iPts, iAtom
        real(kind(1d0)) :: rvec(3)

        write(*, "(a)") "Computing Weights"

        allocate(WEIGHTV(nPts, nAtoms))
        !$OMP PARALLEL DO PRIVATE(iPts, rvec, iAtom)
        do iPts = 1, nPts
            rvec = gridv(:, iPts)
            do iAtom = 1, nAtoms
                WeightV(iPts, iAtom) = wkfun(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS)
            enddo
        enddo
        !$OMP END PARALLEL DO
    end subroutine ComputeAtomicWeights
    subroutine TabulateChargeDensity(Amat, OrbTab, ChDen)
        real(kind(1d0)), intent(in) :: Amat(:, :)
        real(kind(1d0)), intent(in) :: OrbTab(:, :)
        real(kind(1d0)), intent(out) :: ChDen(:)

        integer :: number_of_orbitals, ii_orb, ij_orb

        number_of_orbitals = size(OrbTab, 2)
        !.. Tabulate Charge Density
        ChDen = 0.d0

        !$OMP PARALLEL DO PRIVATE(iPts, ii_orb, ij_orb)
        do iPts = 1, nPts
            do ii_orb = 1, number_of_orbitals
                do ij_orb = 1, number_of_orbitals
                    ChDen(iPts) = ChDen(iPts) + OrbTab(iPts, ii_orb) * Amat(ii_orb, ij_orb) * OrbTab(iPts, ij_orb)
                enddo
            enddo
        enddo
        !$OMP END PARALLEL DO
    end subroutine TabulateChargeDensity
    subroutine ComputeBeckeMatrix(WeightV, OrbTab, BeckeMatrix, Bary_center)
        real(kind(1d0)), intent(in) :: WeightV(:, :)
        real(kind(1d0)), intent(in) :: OrbTab (:, :)
        real(kind(1d0)), allocatable, intent(out) :: BeckeMatrix(:, :, :, :)
        real(kind(1d0)), intent(in) :: Bary_center(:, :)

        integer :: nPts, number_of_orbitals, nAtoms
        integer :: iPts, iOrb1, iOrb2, iAtom, iPol
        real(kind(1d0)) :: dsum

        !
        write(*, "(a)") "Computing Becke's Matrix"
        !
        nAtoms = size(WeightV, 2)
        nPts = size(WeightV, 1)
        number_of_orbitals = size(OrbTab, 2)

        allocate(BeckeMatrix(3, number_of_orbitals, number_of_orbitals, nAtoms))
        BeckeMatrix = 0.d0
        !$omp parallel do private(iOrb2, iOrb1, iPts) shared(BeckeMatrix, WeightV, OrbTab, Bary_center)
        do iPol = 1, 3
            do iAtom = 1, nAtoms
                do iOrb2 = 1, number_of_orbitals
                    do iOrb1 = 1, number_of_orbitals
                        do iPts = 1, nPts
                            BeckeMatrix(iPol, iOrb1, iOrb2, iAtom) = BeckeMatrix(iPol, iOrb1, iOrb2, iAtom) + &
                                    WeightV(iPts, iAtom) * OrbTab(iPts, iOrb1) * OrbTab(iPts, iOrb2) * &
                                            Bary_center(iPol, iAtom) ! fixme????
                        end do
                    end do
                end do
            end do
        end do
        !$omp end parallel do

    end subroutine ComputeBeckeMatrix
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! Not used in the code ...?
    ! ----------------------------------------------------
    complex(kind(1d0)) function zTraceFunction(zA) result(zTrace)
        integer :: i
        complex(kind(1d0)), intent(in) :: zA(:, :)
        zTrace = 0.d0
        do i = 1, size(zA, 1)
            zTrace = zTrace + zA(i, i)
        enddo
    end function zTraceFunction
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


end program ChargeMigration


