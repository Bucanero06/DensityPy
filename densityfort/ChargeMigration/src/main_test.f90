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
    logical :: Verbous
    logical :: SaveDensity
    integer, allocatable :: ivorb(:)
    integer :: OrbNumber
    real(kind(1d0)) :: Volume
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

    integer :: nPts, nAtoms, nxpoints
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

    complex(kind(1d0)), allocatable :: Liouvillian0(:, :),  L0_LEvec(:, :), L0_REvec(:, :), L0_Eval(:)

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
    real(kind(1d0)), allocatable :: OrbitalDensity(:, :), normalized_OrbTab(:, :), normalized_OrbTab1(:), normalized_OrbTab2(:)
    real(kind(1d0)), allocatable :: Radius_BS(:), ChargeTotal1(:), ChargeTotal2(:, :), timed_stat_rho(:), data(:, :), data1(:), data2(:), data4(:)
    real(kind(1d0)) :: data3, data5, trace, Progress, Computed_Volume, summation
    integer :: iOrb, jOrb, index

    real(kind(1d0)) :: Orbital_overlap_self, Orbital_overlap_other

    complex(kind(1d0)), external :: zdotu

    call GetRunTimeParameters(InpDir, OutDir, FileGeometry, &
            nTimes, Tmin, Tmax, Ext_field_file, Verbous, SaveDensity, ivorb, OrbNumber, Volume)

    call system("mkdir -p " // trim(OutDir))
    call system("mkdir -p " // OutDir // "/Dipole")
    call system("mkdir -p " // OutDir // "/AtomicCharge")
    call system("mkdir -p " // OutDir // "/ChargeDensity")
    call system("mkdir -p " // OutDir // "/Pulses")
    call system("cp " // Ext_field_file // " " // OutDir // "/")
    call system("cp chargemigration.ini " // OutDir // "/")
    call Set_CD_IO_Verbous(Verbous)

    open(newunit = uid, &
            file = Ext_field_file, &
            form = "formatted", &
            status = "old")
    call Parse_Simulation_File(uid, OUTPUT_UNIT, N_Simulations, Simulation_Tagv, train)
    close(uid)
    write(*, *)"N_Simulations=", N_Simulations
    !

    !    !.. Print the pulse
    !    !..
    write(*, *) "=== WRITE PULSE TO FILE ==="
    dt = (tmax - tmin) / dble(nTimes - 1)
    do i = 1, N_Simulations
        write(*, *) "Writing pulse " // trim(Simulation_Tagv(i))
        strn = trim(OUTDir) // "/Pulses/pulse" // trim(Simulation_Tagv(i))
        call train(i)%Write(strn, Tmin, Tmax, dt)
        strn = trim(OUTDir) // "/Pulses/FTpulse" // trim(Simulation_Tagv(i))
        call train(i)%WriteFTA(strn)
    enddo

    call LoadEnergies(InpDir // "/ROOT_ENERGIES", nStates, Evec)

    !.. At the moment, we assume that the TDM are defined as
    !   $\rho^{JI}_{nm} = \langle J | a_n^\dagger a_m | I \rangle $
    !..
    call LoadTDMs    (InpDir // "/DENSITY_MATRIX", InpDir // "/TRANSITION_DENSITY_MATRIX", &
            nStates, nOrb, TDM)

    !.. Load Dipole matrix elements $\mu_{IJ}$
    call LoadDipoleME(Dmat, InpDir, nStates)
    allocate(zDmat_t(nStates, nStates, 3))

    do i = 1, 3
        zDmat_t(:, :, i) = Z1 * transpose(Dmat(:, :, i))
    enddo
    write(*, *) nStates, nOrb

    !    !>.._(TEST) "Trace" of TDM, $\rho^{JI}_{nm} = \langle J | a_n^\dagger a_m | I \rangle $__
    !    !..                         trace($\rho^{JI}$)
    !    trace = 0.0d0
    !    allocate(data(nStates))
    !    do i = 1, nStates
    !        do iPts = 1, nOrb
    !            trace = trace + TDM(iPts, iPts, i, i)
    !        end do
    !        write(*, *) trace
    !        data(i) = trace
    !        trace = 0.0d0
    !    end do
    !    write(*, *)

    !    trace = 0.0d0
    !    do i = 1, nOrb
    !        trace = trace + TDM(:, :, 1, 1)
    !    end do
    !    write(*, *)trace
    !_________________________________________________________________________________________

    !    !>.._(TEST) "Trace" of ground, $\rho^{GG}__
    !    trace = 0.0d0
    !    do i = 1, nOrb
    !        trace = trace + TDM(i, i, 1, 1)
    !    end do
    !    data1 = trace
    !    write(*, *) "Trace(\rho^{GG}) = # electrons = ", trace
    !    write(*, *)
    !    write(*, *)
    !__________________________________________________________

    !    !>.. \int rho * dr
    !    data2 = 0.d0
    !    data3 = 0.d0
    !    do iPts = 1, nPts
    !        do iOrb = 1, nOrb
    !            do jOrb = 1, nOrb
    !                data2(iPts) = data2(iPts) + OrbTab(iPts, iOrb) * OrbTab(iPts, jOrb)
    !            end do
    !        end do
    !    end do
    !    data3 = sum(data2) * Volume


    !..Load Geometry,Volume, Grid and Orbitals
    call LoadGeometry(nAtoms, AtCoord, FileGeometry, AtomName)

    call LoadGrid(InpDir // "/gridcoord", npts, gridv)! $G=\{\vec{r}_i, i=1,\ldots,N_{G}\}$
    !    nxpoints = int(exp(log(nPts + 0.1d0) / 3.d0))
    call ComputeVolume(nPts, Computed_Volume, gridv)

    call LoadOrbitals(InpDir, nOrb, ivOrb, npts, OrbTab)! $\varphi_n(\vec{r}_i)$

    !    call ComputeAtomicWeights1(nPts, gridv, nAtoms, AtCoord, WeightV)

    !    call ComputeBeckeMatrix(WeightV, OrbTab, BeckeMatrix)
    !    BeckeMatrix = BeckeMatrix * Volume

    call ComputeLiouvillian0(Evec, Dmat, Liouvillian0)

    call DiagonalizeLindblad0(Liouvillian0, L0_Eval, L0_LEvec, L0_REvec)

    allocate(zStatRho (nStates, nStates))
    allocate(zStatRho0(nStates, nStates))

    !    allocate(AtomicChargeVec(nAtoms))
    !    allocate(AtomicChargeEvolution(nAtoms, nTimes))

    allocate(ChDen(nPts))


    !    allocate(normalized_OrbTab(nPts, nOrb), normalized_OrbTab1(nPts))
    !    call Normalize(OrbTab(:, 1), normalized_OrbTab1)
    !    do iOrb = 1, nOrb
    !        call Normalize(OrbTab(:, iO/home/rubenrb), normalized_OrbTab1)
    !        normalized_OrbTab(:, iOrb) = normalized_OrbTab1
    !    end do
    !
    !    do iPts = 1, nPts
    !        write(*, *) (normalized_OrbTab(iPts, iOrb), OrbTab(iPts, iOrb), iOrb = 1, nOrb)
    !    end do
    !>.. \int (|\varphi_n(\vec{r}_i)|^2) * dr

    !    do iPts = 1, nPts
    !        do iOrb = 1, nOrb
    !            data1(iOrb) = data1(iOrb) + OrbTab(iPts, iOrb) * OrbTab(iPts, iOrb)
    !        end do
    !    end do
    !    data1 = data1 * Volume
    !_________________________________________

    !    !>.. \int rho * dr
    !    data1 = 0.d0
    !    data3 = 0.d0
    !    index = 0
    !    do iOrb = 1, nOrb
    !        index = index + 1
    !        do jOrb = index, nOrb
    !            do iPts = 1, nPts
    !                if (iOrb .ne. jOrb) then
    !                    data3 = data3 + Modulus(normalized_OrbTab(iPts, iOrb)) * Modulus(normalized_OrbTab(iPts, jOrb))
    !                elseif (iOrb == jOrb) then
    !                    data1(iOrb) = data1(iOrb) + Modulus(normalized_OrbTab(iPts, iOrb)) * Modulus(normalized_OrbTab(iPts, jOrb))
    !                end if
    !            end do
    !        end do
    !    end do
    !    data1 = data1 * Volume
    !    data3 = data3 * Volume
    !    write(*, *) "Overlap \varphi_n(\vec{r}) w/ \varphi_n(\vec{r})"
    !    write(*, *) (data1(iOrb), iOrb = 1, nOrb)
    !    write(*, *)
    !    write(*, *)"Overlap \varphi_n(\vec{r}) w/ \varphi_m(\vec{r})"
    !    write(*, *) data3

    !>.. \int rho * dr
    allocate(data(nPts, nOrb))
    write(*, *)
    summation = sum(BeckeMatrix(:, :, :))
    write(*, "(a,e14.6)") "sum beckes", summation
    data = OrbTab * OrbTab
    summation = sum(data(:, :)) * Volume
    write(*, "(a,e14.6)") "sum Orbtab", summation
    write(*, *)

    do iOrb = 1, nOrb
        !Compute orbital norm
        Orbital_overlap_self = 0.d0
        do iPts = 1, nPts
            Orbital_overlap_self = Orbital_overlap_self + OrbTab(iPts, iOrb) * OrbTab(iPts, iOrb)
        end do
        Orbital_overlap_self = Orbital_overlap_self * Computed_Volume

        !Compute Overlap with other Orbitals
        Orbital_overlap_other = 0.d0
        do jOrb = 1, nOrb
            if (jOrb==iOrb) cycle
            do iPts = 1, nPts
                Orbital_overlap_other = Orbital_overlap_other + OrbTab(iPts, iOrb) * OrbTab(iPts, jOrb)
            end do
        end do
        Orbital_overlap_other = Orbital_overlap_other * Computed_Volume

        !Write Results to Screen
        write(*, "(a,i0,a,e14.6,a,e14.6)") "Norm(", iOrb, ") =", Orbital_overlap_self, " Overlap other = ", Orbital_overlap_other

    end do
    write(*, "(a,e14.6,a,e14.6,a,i0)") "Volume = ", Volume, " Compute Volume = ", Computed_Volume, " Number of Points = ", nPts
    !_________________________________________



    !    Sim_loop : do iSim = 1, N_Simulations
    !
    !        zStatRho = Z0
    !        zStatRho(GS_IDX, GS_IDX) = 1.d0
    !
    !        !.. Time cycle
    !        time_loop : do it = 1, nTimes
    !            !
    !            t = tmin + dt * dble(it - 1)
    !
    !            timed_stat_rho(it) = zTraceFunction(zStatRho)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !            !.. Computes the excitation density wrt ground state
    !            zStatRho(GS_IDX, GS_IDX) = zStatRho(GS_IDX, GS_IDX) - 1.d0
    !            !            write(*, *) "trace", t, zTraceFunction(zStatRho)
    !
    !            !>..\sum_{ij} \rho_{ij}^M(t) \rho_{nm}^{ji}-\rho_{nm}^{gg}
    !            call ComputeOrbitalDensity(zStatRho, TDM, OrbitalDensity)
    !
    !            !                do iPts = 1, nPts
    !            !                    do iOrb = 1, nOrb
    !            !                        data(iPts, iOrb) = data(iPts, iOrb) + OrbTab(iPts, iOrb) * OrbitalDensity(iOrb, iOrb) * OrbTab(iPts, iOrb)
    !            !                        do jOrb = 1, nOrb
    !            !                            ChDen(iPts) = ChDen(iPts) + OrbTab(iPts, iOrb) * OrbitalDensity(iOrb, jOrb) * OrbTab(iPts, jOrb)
    !            !                        enddo
    !            !                    enddo
    !            !                enddo
    !            !
    !            !                ChargeTotal1(it) = sum(ChDen) * Volume
    !            !                do iOrb = 1, nOrb
    !            !                    ChargeTotal2(it, iOrb) = sum(data(:, iOrb)) * Volume
    !            !                end do
    !
    !            !            !>.. \int rho * dr
    !            !            data1 = 0.d0
    !            !            data3 = 0.d0
    !            !            index = 0
    !            !            do iOrb = 1, nOrb
    !            !                index = index + 1
    !            !                do jOrb = 1, nOrb
    !            !                    do iPts = index, nPts
    !            !                        if (iOrb .ne. jOrb) then
    !            !                            data3 = data3 + OrbTab(iPts, iOrb) * OrbitalDensity(iOrb, jOrb) * OrbTab(iPts, jOrb)
    !            !                        elseif (iOrb == jOrb) then
    !            !                            data1(iOrb) = data1(iOrb) + OrbTab(iPts, iOrb) * OrbitalDensity(iOrb, jOrb) * OrbTab(iPts, jOrb)
    !            !                        end if
    !            !                    end do
    !            !                end do
    !            !            end do
    !            !            data1 = data1 * Volume
    !            !            data3 = data3 * Volume
    !            !            write(*, *) "Overlap \varphi_n(\vec{r}) w/ \varphi_n(\vec{r})"
    !            !            write(*, *) (data1(iOrb), iOrb = 1, nOrb)
    !            !            write(*, *)
    !            !            write(*, *)"Overlap \varphi_n(\vec{r}) w/ \varphi_m(\vec{r})"
    !            !            write(*, *) data3
    !            !            !_________________________________________
    !
    !
    !
    !
    !
    !
    !
    !
    !
    !            !            call ComputeAtomicCharges(OrbitalDensity, BeckeMatrix, AtomicChargeVec)
    !            !            AtomicChargeEvolution(:, it) = AtomicChargeVec
    !
    !            zStatRho(GS_IDX, GS_IDX) = zStatRho(GS_IDX, GS_IDX) + 1.d0
    !
    !            if(it == nTimes) exit time_loop
    !
    !            call LiouvillianPropagator(L0_Eval, L0_LEvec, L0_REvec, zStatRho)
    !
    !            write(*, *) it, nTimes
    !        enddo time_loop
    !        write(*, *)
    !
    !        !        !.. Save Q_Charge
    !        !        open(newunit = uid_AtomicCharge, &
    !        !                file = OutDir // "/AtomicCharge/AtomicCharge" // trim(Simulation_tagv(iSim)), &
    !        !                form = "formatted", &
    !        !                status = "unknown", &
    !        !                action = "write")
    !        !        do it = 1, nTimes
    !        !            t = tmin + dt * dble(it - 1)
    !        !            write(uid_AtomicCharge, "(*(x,e24.16))") t, sum(AtomicChargeEvolution(:, it)), (AtomicChargeEvolution(iAtom, it), iAtom = 1, nAtoms)!, ChargeTotal1(it), ChargeTotal2(it)
    !        !        enddo
    !        !        close(uid_AtomicCharge)
    !        !
    !        !        !.. Save file
    !        !        open(newunit = uid, &
    !        !                file = OutDir // "/test_file", &
    !        !                form = "formatted", &
    !        !                status = "unknown", &
    !        !                action = "write")
    !        !        !        write(uid, *) "Trace(\rho^{GG}) = #_electrons = ", data1
    !        !        !        write(uid, *)
    !        !        !        write(uid, *)
    !        !        !        write(uid, *) "Trace($\rho^{JI}$) = #_electrons * \delta_{ij} = ", (data(i), i = 1, nStates)
    !        !        !        write(uid, *)
    !        !        !        write(uid, *)
    !        !        !        write(uid, *) "   time  ", "            (Trace(\rho_{ij}^M(t))) ", "  (N_e[\sum_{i}\rho_{ii}^M(t)-1])"
    !        !        !        do it = 1, nTimes
    !        !        !            t = tmin + dt * dble(it - 1)
    !        !        !            write(uid, *) t, dble(timed_stat_rho(it)), ((data1) * ((dble(timed_stat_rho(it)) - 1)))
    !        !        !        enddo
    !        !
    !        !        write(uid, *) " time  ", " sum(ChargeDensity) ", "                sum(ChargeDensity(iOrb))"
    !        !        write(uid, *)
    !        !        do it = 1, nTimes
    !        !            t = tmin + dt * dble(it - 1)
    !        !            write(uid, "(i4,*(x,E24.16))") t, ChargeTotal1(it), (ChargeTotal2(it, iOrb), iOrb = 1, nOrb)
    !        !        end do
    !        !
    !        !        close(uid)
    !
    !    end do Sim_loop

    stop


contains

    real(kind(1d0)) function Modulus(avec) result(mod)
        real(kind(1d0)), intent(in) :: avec
        mod = sqrt(avec * avec)
    end function Modulus

    subroutine Normalize(amat, normalized_OrbTab)
        real(kind(1d0)), intent(in) :: amat(:)
        real(kind(1d0)), intent(out) :: normalized_OrbTab(:)
        normalized_OrbTab = amat / (sqrt(amat * amat))
    end subroutine Normalize


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
        Amat = 0.d0

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
                                dble(zStatRho(iStatei, iStatej) * TDM(iOrbi, iOrbj, iStatej, iStatei))!TDM(iOrbi, iOrbj, iStatei, iStatej))!TDM(iOrbi, iOrbj, iStatej, iStatei))!!!flip indexes?? (check)
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
    !    subroutine ComputeVolume(nPts, Volume, gridv)
    !        integer, intent(in) :: nPts
    !        real(kind(1d0)), intent(out) :: Volume
    !        real(kind(1d0)), allocatable, intent(in) :: gridv(:, :)
    !        real(kind(1d0)) :: coord1, coord2
    !        integer :: iPts, iCoord
    !        Volume = 1.d0
    !        do iCoord = 1, 3
    !            coord1 = gridv(iCoord, 1)
    !            do iPts = 2, nPts
    !                coord2 = gridv(iCoord, iPts)
    !                if (coord1 .ne. coord2) then
    !                    Volume = Volume * abs(coord2 - coord1)
    !                    exit
    !                endif
    !            enddo
    !        enddo
    !    end subroutine ComputeVolume
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
    !    subroutine ComputeBeckeMatrix(WeightV, OrbTab, Becke)
    !
    !        real(kind(1d0)), allocatable, intent(in) :: WeightV(:, :)
    !        real(kind(1d0)), allocatable, intent(in) :: OrbTab (:, :)
    !        real(kind(1d0)), allocatable, intent(out) :: Becke(:, :, :)
    !
    !        integer :: nPts, nOrbs, nAtoms
    !        integer :: iPts, iOrb1, iOrb2, iAtom
    !        real(kind(1d0)) :: dsum
    !
    !        nAtoms = size(WeightV, 2)
    !        nPts = size(WeightV, 1)
    !        nOrbs = size(OrbTab, 2)
    !
    !        if(allocated(Becke))deallocate(Becke)
    !        allocate(Becke(nOrbs, nOrbs, nAtoms))
    !
    !        Becke = 0.d0
    !        do iAtom = 1, nAtoms
    !            do iOrb2 = 1, nOrbs
    !                do iOrb1 = 1, nOrbs
    !                    dsum = 0.d0
    !                    do iPts = 1, nPts
    !                        dsum = dsum + WeightV(iPts, iAtom) * OrbTab(iPts, iOrb1) * OrbTab(iPts, iOrb2)
    !                    enddo
    !                    Becke(iOrb1, iOrb2, iAtom) = dsum
    !                enddo
    !            enddo
    !        enddo
    !
    !    end subroutine ComputeBeckeMatrix
    !
    !
    subroutine ComputeAtomicWeights1(nPts, gridv, nAtoms, AtCoord, WeightV)

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
                WeightV(iPts, iAtom) = wkfun1(rvec, iAtom, AtCoord, nAtoms, k)
            enddo
        enddo
    end subroutine ComputeAtomicWeights1
    !
    !    subroutine ComputeAtomicCharges(OrbitalDensity, BeckeMatrix, QchargeVec)
    !        real(kind(1d0)), intent(in) :: OrbitalDensity(:, :)
    !        real(kind(1d0)), intent(in) :: BeckeMatrix   (:, :, :)
    !        real(kind(1d0)), allocatable, intent(out) :: QchargeVec    (:)
    !        integer :: iAtom, nAtoms, nOrbs
    !        real(kind(1d0)), external :: DDOT
    !        nOrbs = size(BeckeMatrix, 1)
    !        nAtoms = size(BeckeMatrix, 3)
    !        if(.not.allocated(QchargeVec))allocate(QchargeVec(nAtoms))
    !        QchargeVec = 0.d0
    !        do iAtom = 1, nAtoms
    !            QchargeVec(iAtom) = ddot(nOrbs * nOrbs, OrbitalDensity, 1, BeckeMatrix(1, 1, iAtom), 1)
    !        enddo
    !    end subroutine ComputeAtomicCharges
    !    !####################


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!>>>>>>>>>>>>>>
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
        end if

    end function Radius_Table

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

    subroutine ComputeBeckeMatrix(WeightV, OrbTab, Becke)
        real(kind(1d0)), allocatable, intent(in) :: WeightV(:, :)
        real(kind(1d0)), allocatable, intent(in) :: OrbTab (:, :)
        real(kind(1d0)), allocatable, intent(out) :: Becke(:, :, :)

        integer :: nPts, nOrbs, nAtoms
        integer :: iPts, iOrb1, iOrb2, iAtom
        real(kind(1d0)) :: dsum

        nAtoms = size(WeightV, 2)
        nPts = size(WeightV, 1)
        nOrbs = size(OrbTab, 2)

        if(allocated(Becke))deallocate(Becke)
        allocate(Becke(nOrbs, nOrbs, nAtoms))

        Becke = 0.d0
        do iAtom = 1, nAtoms
            do iOrb2 = 1, nOrbs
                do iOrb1 = 1, nOrbs
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
    subroutine ComputeLiouvillian0(Evec, Dmat, Liouvillian0)
        real   (kind(1d0)), intent(in) :: Evec(:)
        real   (kind(1d0)), intent(in) :: Dmat(:, :, :)
        complex(kind(1d0)), allocatable, intent(out) :: Liouvillian0(:, :)

        real(kind(1d0)), parameter :: WATER_MELTING_POINT = 273.15
        real(kind(1d0)), parameter :: BATH_TEMPERATURE = WATER_MELTING_POINT + 3000
        real(kind(1d0)), parameter :: BETA = 1.d0 / (BOLTZMANN_CONSTANT_auK * BATH_TEMPERATURE)
        real(kind(1d0)), parameter :: dephasing_factor = 1.d-3
        real(kind(1d0)), parameter :: RELAXATION_FACTOR = 1.d-3

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

    end subroutine ComputeLiouvillian0

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
    !        integer :: iPol, iOmega, iLiou, iState
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
    !                do iState = 1, nStates
    !                    DipoleFTplus(iPol, iOmega) = DipoleFTplus(iPol, iOmega) + zmat2(iState, iState)
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
!!$  call LoadDipoleMO( InpDir, nOrb, ivOrb, MuOrb )
!!$
!!$  !.. Check Dipole Matrix Elements between States
!!$  do iStatei = 1, nStates
!!$     do iStatej = 1, nStates
!!$        TME = 0.d0
!!$        do i = 1, nOrb
!!$           do j = 1, nOrb
!!$              TME = TME + TDM(i,j,iStatei,iStatej) * MuOrb(j,i)
!!$           enddo
!!$        enddo
!!$        write(*,*) iStatei, iStatej, TME
!!$     enddo
!!$  enddo


!        if (R_i < 0.36) then
!            r_m = R_i
!        else
!            r_m = R_i * 0.5d0
!        end if