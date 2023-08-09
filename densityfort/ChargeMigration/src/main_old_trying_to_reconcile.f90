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

    real(kind(1d0)), parameter :: WATER_MELTING_POINT = 273.15
    !    real(kind(1d0)), parameter :: bath_temperature = WATER_MELTING_POINT + 3000
    real(kind(1d0)) :: bath_temperature

    !    real(kind(1d0)), parameter :: BETA = 1.d0 / (BOLTZMANN_CONSTANT_auK * bath_temperature)

    real(kind(1d0)) :: dephasing_factor
    real(kind(1d0)) :: relaxation_factor


    !.. Run-time parameters
    !..
    character(len = :), allocatable :: input_directory
    character(len = :), allocatable :: output_directory
    character(len = :), allocatable :: molecular_geometry_file
    integer :: n_times
    real(kind(1d0)) :: t_min
    real(kind(1d0)) :: t_max
    real(kind(1d0)) :: PolVec(3)
    character(len = :), allocatable :: Ext_field_file
    logical :: Verbous
    integer, allocatable :: ivorb(:)
    character(len = :), allocatable :: Weight_File
    logical :: read_precomputed_weights_flag
    logical :: save_charge_migration_flag
    integer :: counted_number_of_orbitals

    real(kind(1d0)), parameter :: dephasing_constant = 0.d0!1.d-2
    integer, parameter :: GROUND_STATE_INDEX = 1

    !.. Local parameters
    integer :: nStates, number_of_orbitals, ii_state, i_state, ij_state, i, j
    real(kind(1d0)), allocatable :: Evec(:)
    !.. Dmat(i,j,alpha) = $ \langle \varphi_i | \hat{\mu}_\alpha$ | \varphi_j \rangle $
    real(kind(1d0)), allocatable :: Dmat(:, :, :)
    !.. TDM(j,i,B,A) = $ \langle A | \hat{a}_i^\dagger \hat{a}_j | B \rangle $
    real(kind(1d0)), allocatable :: TDM(:, :, :, :)
    !.. Expectation Value of the Dipole Moment (Mu)
    real(kind(1d0)) :: MuEV(3)

    integer :: npts, nAtoms, nxpoints
    real(kind(1d0)), allocatable :: gridv(:, :), AtCoord(:, :) ! 3 x npts
    real(kind(1d0)), allocatable :: OrbTab(:, :) ! 3 x npts
    real(kind(1d0)) :: Computed_volume

    real   (kind(1d0)), allocatable :: StatRho0(:, :) !.. Initial Statistical Density Matrix (real ...)
    complex(kind(1d0)), allocatable :: zStatRho(:, :) !.. Initial Statistical Density Matrix (real ...)
    real   (kind(1d0)) :: t, dt
    real   (kind(1d0)), allocatable :: Amat(:, :), ChDen(:), WEIGHTV(:, :), Q_Charge(:)
    complex(kind(1d0)) :: zTrace

    character(len = 30) :: istrn
    character(len = 16), allocatable :: atom_names(:)

    integer :: ii_orb, ij_orb, iPts, it, iPol, i1, j1, i2, j2, iLiou, iLiou1, iLiou2
    integer :: uid_ch

    integer :: nLiou, n
    complex(kind(1d0)), allocatable :: Liouvillian0(:, :), L0_LEvec(:, :), L0_REvec(:, :), L0_Eval(:)
    complex(kind(1d0)), allocatable :: RhoVec(:), RhoVec2(:), zmat(:, :), zmat1(:, :), zmat2(:, :)
    real   (kind(1d0)) :: AverageDipole, dBuf
    complex(kind(1d0)), allocatable :: zPureMat(:, :), zUMAT_DX(:, :), zUMAT_DY(:, :), zUMAT_DZ(:, :)
    complex(kind(1d0)), allocatable :: zUMAT_DYs_DX(:, :), zUMAT_DZs_DY(:, :)
    real   (kind(1d0)), allocatable :: EPure(:), rmat(:, :), EVEC_DX(:), EVEC_DY(:), EVEC_DZ(:)

    !.. Dephasing and Relaxation Constants
    real   (kind(1d0)), allocatable :: PairGamma(:, :), TotalGamma(:)

    !.. Pulse parameters
    integer :: iSim, N_Simulations, uid
    character(len = 64), pointer :: Simulation_Tagv(:)
    character(len = 1000) :: strn
    type(pulse_train), pointer :: train(:)
    real(kind(1d0)) :: EFieldX, EFieldY, EFieldZ
    real(kind(1d0)) :: BETA

    integer :: iOrb, jOrb, number_of_orbitalss
    real(kind(1d0)) :: Orbital_overlap_self, Orbital_overlap_other


    call GetRunTimeParameters(input_directory, output_directory, molecular_geometry_file, &
            n_times, t_min, t_max, Ext_field_file, Verbous, Weight_File, read_precomputed_weights_flag, &
            save_charge_migration_flag, ivorb, &
            counted_number_of_orbitals, &
            Computed_volume, &
            dephasing_factor, &
            relaxation_factor, &
            bath_temperature)

    BETA = 1.d0 / (BOLTZMANN_CONSTANT_auK * bath_temperature)

    call system("mkdir -p " // trim(output_directory))
    call Set_CD_IO_Verbous(Verbous)
    dt = (t_max - t_min) / dble(n_times - 1)

    open(newunit = uid, &
            file = Ext_field_file, &
            form = "formatted", &
            status = "old")
    call Parse_Simulation_File(uid, OUTPUT_UNIT, N_Simulations, Simulation_Tagv, train)
    close(uid)
    write(*, "(a,i3)")"N_Simulations=", N_Simulations

    !.. Print the pulse
    !..
    write(*, *) "=== WRITE PULSE TO FILE ==="
    do i = 1, N_Simulations
        write(*, *) "Writing pulse " // trim(Simulation_Tagv(i))
        strn = trim(output_directory) // "/pulse" // trim(Simulation_Tagv(i))
        call train(i)%Write(strn, t_min, t_max, dt)
        strn = trim(output_directory) // "/FTpulse" // trim(Simulation_Tagv(i))
        call train(i)%WriteFTA(strn)
    enddo

    call LoadEnergies(input_directory // "/ROOT_ENERGIES", nStates, Evec)   ! $E_I$

    do i = 1, nStates - 1
        do j = i + 1, nStates
            write(*, *) i, j, Evec(j) - Evec(i)
        enddo
    enddo

    !.. At the moment, we assume that the TDM are defined as
    !   $\rho^{JI}_{nm} = \langle J | a_n^\dagger a_m | I \rangle $
    !..
    call LoadTDMs    (input_directory // "/DENSITY_MATRIX", input_directory // "/TRANSITION_DENSITY_MATRIX", nStates, number_of_orbitals, TDM)

    !.. Load Dipole matrix elements $\mu_{IJ}$
    call LoadDipoleME(Dmat, input_directory, nStates)


    !.. Diagonalizes the dipole matrices
    allocate(rmat(nStates, nStates))
    !
    !.. X dipole
    allocate(zUMAT_DX(nStates, nStates), EVEC_DX(nStates))
    rmat = Dmat(:, :, 1)
    EVEC_DX = 0.d0
    call Short_Diag(nStates, rmat, EVEC_DX)
    zUMAT_DX = Z1 * rmat
    !
    !.. Y dipole
    allocate(zUMAT_DY(nStates, nStates), EVEC_DY(nStates))
    rmat = Dmat(:, :, 2)
    EVEC_DY = 0.d0
    call Short_Diag(nStates, rmat, EVEC_DY)
    zUMAT_DY = Z1 * rmat
    !
    !.. Z dipole
    allocate(zUMAT_DZ(nStates, nStates), EVEC_DZ(nStates))
    rmat = Dmat(:, :, 3)
    EVEC_DZ = 0.d0
    call Short_Diag(nStates, rmat, EVEC_DZ)
    zUMAT_DZ = Z1 * rmat
    !
    deallocate(rmat)
    !
    !.. Computes combined matrices
    allocate(zUMAT_DYs_DX(nStates, nStates))
    zUMAT_DYs_DX = Z0
    call ZGEMM("C", "N", nStates, nStates, nStates, Z1, zUMAT_DY, nStates, zUMAT_DX, nStates, Z0, zUMAT_DYs_DX, nStates)
    allocate(zUMAT_DZs_DY(nStates, nStates))
    zUMAT_DZs_DY = Z0
    call ZGEMM("C", "N", nStates, nStates, nStates, Z1, zUMAT_DZ, nStates, zUMAT_DY, nStates, Z0, zUMAT_DZs_DY, nStates)



    !..Load Geometry,Computed_volume, Grid and Orbitals
    call LoadGeometry(nAtoms, AtCoord, molecular_geometry_file, atom_names)
    call LoadGrid(input_directory // "/gridcoord", npts, gridv)! $G=\{\vec{r}_i, i=1,\ldots,N_{G}\}$
    nxpoints = int(exp(log(nPts + 0.1d0) / 3.d0))
    call Computevolume(nPts, Computed_volume, gridv)
    call LoadOrbitals(input_directory, number_of_orbitals, ivOrb, npts, OrbTab)! $\varphi_n(\vec{r}_i)$

    call ComputeAtomicWeights(nPts, gridv, nAtoms, AtCoord, WeightV)




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


    !.. Diagonalize the time-independent Lindblad superoperator
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

    call system("mkdir -p " // output_directory)

    !!$  call SetStatDenMat_ImpulsiveApp( nStates, GROUND_STATE_INDEX, Dmat, PolVec, StatRho0, Verbous )
    !!$  allocate(zStatRho(nStates,nStates))
    !!$  zStatRho=Z1*StatRho0

    allocate(zStatRho(nStates, nStates))
    zStatRho = Z0
    zStatRho(GROUND_STATE_INDEX, GROUND_STATE_INDEX) = 1.d0

    allocate(RhoVec (nLiou))
    allocate(RhoVec2(nLiou))
    RhoVec = Z0
    RhoVec2 = Z0
    call HilbertToLiouvilleMatrix(zStatRho, RhoVec)

    allocate(Amat(number_of_orbitals, number_of_orbitals))
    allocate(ChDen(nPts))

    allocate(zPureMat, source = zStatRho)
    allocate(EPure(nStates))
    allocate(zmat1(nStates, nStates))
    allocate(zmat2(nStates, nStates))

    Sim_loop : do iSim = 1, N_Simulations
        if (save_charge_migration_flag) then

            !.. Open Q_Charge File
            open(newunit = uid_ch, &
                    file = "Q_Charge" // trim(Simulation_tagv(iSim)), &
                    form = "formatted", &
                    status = "unknown", &
                    action = "write")
        endif
        !.. Time cycle
        time_loop : do it = 1, n_times
            !
            t = t_min + dt * dble(it - 1)

            !====== DIAGNOSTICS BEGIN =============================================

            !.. Computes the excitation density wrt ground state
            zStatRho(GROUND_STATE_INDEX, GROUND_STATE_INDEX) = zStatRho(GROUND_STATE_INDEX, GROUND_STATE_INDEX) - 1.d0
            !
            zTrace = 0.d0
            do i = 1, nStates
                zTrace = zTrace + zStatRho(i, i)
            enddo
            write(*, *) "trace", t, zTrace

            !.. Print the expectation value of the dipole as a function of time
            MuEV = 0.d0
            do iPol = 1, 3
                do i = 1, nStates
                    do j = 1, nStates
                        MuEV(iPol) = MuEV(iPol) + dble(zStatRho(i, j) * Dmat(j, i, iPol))
                    enddo
                enddo
            enddo
            !!$     write(*,"(i4,*(x,d14.6))") it, t, (MuEV(iPol),iPol=1,3)

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
            !!$     write(*,"(i4,*(x,d14.6))") it, t, (Amat(1,ij_orb),ij_orb=1,4)

            if (save_charge_migration_flag) then
                !.. Tabulate Charge Density
                ChDen = 0.d0
                do iPts = 1, nPts
                    do ii_orb = 1, number_of_orbitals
                        do ij_orb = 1, number_of_orbitals
                            ChDen(iPts) = ChDen(iPts) + OrbTab(iPts, ii_orb) * Amat(ii_orb, ij_orb) * OrbTab(iPts, ij_orb)
                        enddo
                    enddo
                enddo

                !.. Save Charge Density
                write(istrn, "(f12.4)")t
                call Write_Charge_Density(output_directory // "/ChDen" // trim(adjustl(istrn)), nPts, gridv, ChDen, Weightv, nAtoms)

                !..
                call Compute_Q_Charge(nAtoms, nPts, ChDen, Weightv, Computed_volume, Q_Charge)
                write(uid_ch, "(*(x,e24.16))") t, Q_Charge
            endif
            !======= DIAGNOSTIC ENDS =================================================

            if(it == n_times) exit time_loop

            !======= TIME-STEP PROPAGATION ===========================================

            !!$     !.. Compute StatRho(t)
            !!$     !
            !!$     do i = 1, nStates
            !!$        do j = 1, nStates
            !!$           zStatRho(i,j) = StatRho0(i,j) * exp( - Zi* ( Evec(i) - Evec(j) ) * ( t + dt - t_min ) )
            !!$           if (i .ne. j) then
            !!$              zStatRho(i,j) = zStatRho(i,j) * exp( - dephasing_constant * ( t + dt - t_min ) )
            !!$           endif
            !!$        enddo
            !!$     enddo


            !.. First half of the field-free propagation in Liouville Space across an interval dt/2
            call ZGEMV("C", nLiou, nLiou, Z1, L0_LEvec, nLiou, RhoVec, 1, Z0, RhoVec2, 1)
            RhoVec2 = RhoVec2 * exp(-Zi * L0_Eval * dt / 2.d0)
            call ZGEMV("N", nLiou, nLiou, Z1, L0_REvec, nLiou, RhoVec2, 1, Z0, RhoVec, 1)
            call LiouvilleToHilbertMatrix(RhoVec, zStatRho)


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
            !
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

            !====== END TIME-STEP PROPAGATION ========================================
        enddo time_loop

        if (save_charge_migration_flag) then
            close(uid_ch) !Close Q_Charge File
        end if
    end do Sim_loop

    write(*, *)
    write(*, *)
    write(*, *) "SUMMARY"
    write(*, *) "# of Points=", nPts
    write(*, *) "# of Atoms=", nAtoms
    write(*, *) "Computed_volume=", Computed_volume
    write(*, *) "n_time=", n_times
    write(*, *) "t_min=", t_min
    write(*, *) "t_max=", t_max

    do iOrb = 1, number_of_orbitals
        !Compute orbital norm
        Orbital_overlap_self = 0.d0
        do iPts = 1, nPts
            Orbital_overlap_self = Orbital_overlap_self + OrbTab(iPts, iOrb) * OrbTab(iPts, iOrb)
        end do
        Orbital_overlap_self = Orbital_overlap_self * Computed_volume

        !Compute Overlap with other Orbitals
        Orbital_overlap_other = 0.d0
        do jOrb = 1, number_of_orbitals
            if (jOrb==iOrb) cycle
            do iPts = 1, nPts
                Orbital_overlap_other = Orbital_overlap_other + OrbTab(iPts, iOrb) * OrbTab(iPts, jOrb)
            end do
        end do
        Orbital_overlap_other = Orbital_overlap_other * Computed_volume

        !Write Results to Screen
        write(*, "(a,i0,a,e14.6,a,e14.6)") "    Norm(", iOrb, ") =", Orbital_overlap_self, "    Overlap other =", Orbital_overlap_other
    end do

contains

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


    !> Computes Electron Density at Each Atom
    subroutine Compute_Q_Charge(nAtoms, nPts, ChDen, Weightv, Computed_volume, Q_Charge)
        !
        implicit none
        !
        integer, intent(in) :: nAtoms, nPts
        real(kind(1d0)), intent(in) :: Computed_volume
        real(kind(1d0)), allocatable, intent(in) :: Weightv(:, :)
        real(kind(1d0)), allocatable, intent(in) :: ChDen(:)
        real(kind(1d0)), allocatable, intent(out) :: Q_Charge(:)
        !z
        integer :: iAtom

        allocate(Q_Charge(nAtoms))

        do iAtom = 1, nAtoms
            Q_Charge(iAtom) = sum(Chden(1:nPts) * Weightv(1:nPts, iAtom)) * Computed_volume
        enddo
        !
    end subroutine Compute_Q_Charge

    !>Determines 1 Computed_volume Unit of the Grid
    !
    subroutine Computevolume(nPts, Computed_volume, gridv)
        integer, intent(in) :: nPts
        real(kind(1d0)), intent(out) :: Computed_volume
        real(kind(1d0)), allocatable, intent(in) :: gridv(:, :)
        real(kind(1d0)) :: coord1, coord2
        integer :: iPts, iCoord
        Computed_volume = 1.d0
        do iCoord = 1, 3
            coord1 = gridv(iCoord, 1)
            do iPts = 2, nPts
                coord2 = gridv(iCoord, iPts)
                if (coord1 .ne. coord2) then
                    Computed_volume = Computed_volume * abs(coord2 - coord1)
                    exit
                endif
            enddo
        enddo
    end subroutine Computevolume


    !..FUNCTIONS
    !
    !.. Eq. 3 $w_a \vec(r)$  (WEIGHTS)
    !Result =P_a(\vec{r})/\sum_bP_b(\vec{r}) $
    real(kind(1d0)) function wkfun(rvec, iAtom, AtCoord, &
            nAtoms, k) result(res)
        !
        integer, intent(in) :: iAtom, k, nAtoms
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms)
        !
        res = Pkfuna(rvec, iAtom, AtCoord, nAtoms, k) / PkfunTot(rvec, AtCoord, nAtoms, k)
    end function wkfun

    !..Eq. 3 \sum_bP_b(\vec{r}) in the denominator
    !
    real(kind(1d0)) function PkFunTot(rvec, AtCoord, nAtoms, k) result(res)
        integer, intent(in) :: nAtoms
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms)
        integer, intent(in) :: k
        integer :: iAtom
        res = 0.d0
        do iAtom = 1, nAtoms
            res = res + Pkfuna(rvec, iAtom, AtCoord, nAtoms, k)
        enddo
    end function PkFunTot

    !.. Eq. 2 P_a(\vec{r}) in and nominator of Eq. 3
    !
    real(kind(1d0)) function Pkfuna(rvec, iAtom, AtCoord, nAtoms, k) result(res)
        integer, intent(in) :: nAtoms, iAtom
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms)
        integer, intent(in) :: k
        integer :: jAtom
        res = 1.d0
        do jAtom = 1, nAtoms
            if(jAtom == iAtom) cycle
            res = res * skfunab(rvec, AtCoord(:, iAtom), AtCoord(:, jAtom), k)
        enddo
    end function Pkfuna

    !..
    !
    real(kind(1d0)) function skfunab(rvec, avec, bvec, k) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        integer, intent(in) :: k
        real(kind(1d0)) :: mu
        mu = EllipticalCoord(rvec, avec, bvec)
        res = skfun(mu, k)
    end function skfunab

    !.. $s(\mu): [-1,1]  \rightarrow [0,1]$  (Step Function)
    !
    real(kind(1d0)) function skfun(mu, k) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu
        integer, intent(in) :: k
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
        cvec = avec - bvec
        res = sqrt(sum(cvec * cvec))
    end function EuclDist
    !
    !END OF FUNCTIONS


    subroutine SetStatDenMat_ImpulsiveApp(nStates, InitialState, Dmat, PolVec, StatRho0, Verbous)
        implicit none
        integer, intent(in) :: nStates
        integer, intent(in) :: InitialState
        real(kind(1d0)), allocatable, intent(in) :: Dmat(:, :, :)
        real(kind(1d0)), intent(in) :: PolVec(3)
        real(kind(1d0)), allocatable, intent(out) :: StatRho0(:, :)
        logical, intent(in) :: Verbous
        !
        real(kind(1d0)), parameter :: EXCITATION_FACTOR = 1.d0
        !
        real(kind(1d0)), allocatable :: Amplitudev(:)
        real(kind(1d0)) :: NormFact, Amplitudei, Amplitudej
        integer :: ii_state, ij_state

        !.. Computes the expansion coefficients of the state immediately after the dipole excitation
        !   assuming the excitation to be instantaneous
        !..
        allocate(AmplitudeV(nStates))
        AmplitudeV = 1.d0
        do ii_state = 1, nStates
            if(ii_state==InitialState)cycle
            AmplitudeV(ii_state) = EXCITATION_FACTOR * dot_product(Dmat(ii_state, InitialState, :), PolVec)
        enddo

        !.. Compute the Statistical Density Matrix in the impulsive approximation
        !   $P_{IJ}(0) = \frac{ \mu_{Ig} \mu_{gJ}}{\sum_{I} |\mu_{Ig}|^2 }$, where $g$ refers to the ground state.
        allocate(StatRho0(nStates, nStates))
        StatRho0 = 0.d0
        NormFact = 0.d0
        do ii_state = 1, nStates
            do ij_state = 1, nStates
                StatRho0(ii_state, ij_state) = AmplitudeV(ii_state) * AmplitudeV(ij_state)
            enddo
            NormFact = NormFact + StatRho0(ii_state, ii_state)
        enddo
        StatRho0 = StatRho0 / NormFact
        if(Verbous)then
            write(*, *)
            write(*, *) "Statistical Density matrix"
            do ii_state = 1, nStates
                write(*, "(*(x,e14.6))")(StatRho0(ii_state, ij_state), ij_state = 1, nStates)
            enddo
        endif
    end subroutine SetStatDenMat_ImpulsiveApp

    subroutine ComputeAtomicWeights(nPts, gridv, nAtoms, AtCoord, WeightV)

        integer, intent(in) :: nPts, nAtoms
        real(kind(1d0)), intent(in) :: gridv(:, :), AtCoord(:, :)
        real(kind(1d0)), allocatable, intent(out) :: WeightV(:, :)

        integer, parameter :: k = 4
        integer :: iPts, iAtom
        real(kind(1d0)) :: rvec(3)

        allocate(WEIGHTV(nPts, nAtoms))
        do iPts = 1, nPts
            rvec = gridv(:, iPts)
            do iAtom = 1, nAtoms
                WEIGHTV(iPts, iAtom) = wkfun(rvec, iAtom, AtCoord, nAtoms, k)
            enddo
        enddo

    end subroutine ComputeAtomicWeights


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

