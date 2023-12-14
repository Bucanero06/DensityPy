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

    !.. Run-time parameters
    !..
    character(len = :), allocatable :: InpDir, OutDir, FileGeometry
    real(kind(1d0)) :: StepTime, StepWidth
    character(len = :), allocatable :: Ext_field_file
    logical :: Verbous
    integer :: nOmegas, nTauOmegas
    real   (kind(1d0)) :: OmegaMin, OmegaMax, TauOmegaMin, TauOmegaMax

    !.. Data for FT
    !..
    complex(kind(1d0)), allocatable :: DipoleFTtotal(:, :), DipoleFTminus(:, :), DipoleFTplus (:, :)
    complex(kind(1d0)), allocatable :: AtomicChargeFT(:, :)
    !    complex(kind(1d0)), allocatable :: AtomicChargeFT(:, :, :)
    complex(kind(1d0)), allocatable :: DipoleFTwt(:, :, :), DipoleFTww(:, :, :)
    complex(kind(1d0)), allocatable :: ChargeFTwt(:, :, :), ChargeFTww(:, :, :)
    real   (kind(1d0)), allocatable :: OmegaVec(:), TauOmegaVec(:)
    real   (kind(1d0)), allocatable :: tvec(:)
    integer :: iOmega, ntimes, i
    real(kind(1d0)) :: tmin, tmax

    !.. Local parameters
    !..
    character(len = *), parameter :: DIPOLE_FT_PATH_ALL = "DipoleFT_ALL"
    character(len = *), parameter :: CHARGE_FT_PATH_ALL = "/AtomicCharge/AtomicChargeFT_ALL"
    integer, parameter :: GS_IDX = 1
    integer :: nStates
    integer :: uid_AtomicChargeALL, uid
    real(kind(1d0)), allocatable :: Evec(:)
    real(kind(1d0)), external :: NCD_Phi
    complex(kind(1d0)), external :: zdotu

    !.. Expectation Value of the Dipole Moment (Mu)
    !..
    complex(kind(1d0)), allocatable :: zMuEV(:, :)

    !.. Molecular Geometry
    !..
    integer :: nAtoms
    real(kind(1d0)), allocatable :: AtCoord(:, :) ! 3 x npts
    character(len = 16), allocatable :: AtomName(:)

    !.. Statistical Density Matrix
    !..
    real   (kind(1d0)) :: dt
    real   (kind(1d0)), allocatable :: AtomicChargeEvolution(:, :), AtomicChargeEvolutionComponent(:, :, :)

    !.. Pulse parameters
    !..
    integer :: iSim, N_Simulations
    character(len = 64), pointer :: Simulation_Tagv(:)
    type(pulse_train), pointer :: train(:)

    !.. XUV_Dipole and XUV_Charge
    !..
    complex(kind(1d0)), allocatable :: XUVDipole(:, :), XUVCharge(:, :)!, XUVChargeComponent(:, :, :)
    complex(kind(1d0)), allocatable :: XUVDipoleFT(:, :)

    !..Test
    !..
    complex(kind(1d0)), allocatable :: AtomicChargeFT_new(:, :, :), ChargeFTwt_new(:, :, :, :), ChargeFTww_new(:, :, :, :), XUVCharge_new(:, :, :)
    complex(kind(1d0)), allocatable :: XUVDipole_I(:, :, :, :)
    real   (kind(1d0)), allocatable :: AtomicChargeEvolution_new(:, :, :)
    complex(kind(1d0)), allocatable :: Debug(:, :, :)
    external :: system
    integer :: iPol, iAtom, it, iCoord, i_excitation, icomponent_RI
    character(len = 64) :: PART, excitation_string
    !#######################################################################!
    real   (kind(1d0)) :: Mean, Variance, StdDev, Tau_Time, MASK  ! results
    real   (kind(1d0)), allocatable :: Tau_Time_Vec(:)
    character(len = 256) :: output_file
    integer :: ppos
    character(len = 3) :: new_ext = "csv"
    !#######################################################################!

    call GetRunTimeParameters(InpDir, OutDir, FileGeometry, StepTime, StepWidth, &
            Ext_field_file, verbous, nOmegas, OmegaMin, OmegaMax, nTauOmegas, TauOmegaMin, TauOmegaMax, i_excitation, icomponent_RI)

    allocate(OmegaVec, source = [(OmegaMin + (OmegaMax - OmegaMin) / dble(nOmegas - 1) * dble(iOmega - 1), iOmega = 1, nOmegas)])
    allocate(TauOmegaVec, source = [(TauOmegaMin + (TauOmegaMax - TauOmegaMin) / dble(nTauOmegas - 1) * dble(iOmega - 1), iOmega = 1, nTauOmegas)])

    call Set_CD_IO_Verbous(.FALSE.)

    open(newunit = uid, &
            file = Ext_field_file, &
            form = "formatted", &
            status = "old")
    call Parse_Simulation_File(uid, OUTPUT_UNIT, N_Simulations, Simulation_Tagv, train)
    close(uid)

    !.. Loads Root Ernergies and Molecular Geometry
    call LoadEnergies(InpDir // "/ROOT_ENERGIES", nStates, Evec)
    call LoadGeometry(nAtoms, AtCoord, FileGeometry, AtomName)

    !.. Determines the number of times available, as well as the minimum and maximum time
    write(excitation_string, '(a,i0.2)') OutDir // '/Dipole/Dipole_State_01'
    call DetermineNumberofTimes (trim(excitation_string) // "/Dipole" // trim(Simulation_tagv(1)) // "_Real", tmin, tmax, nTimes)
    dt = (tmax - tmin) / dble(nTimes - 1)

    if (i_excitation ==0 .and. icomponent_RI==0)  then !non_parallel
        write(excitation_string, '(a)') trim('rm ' // OutDir // '/Dipole/Dipole_State_01/DipoleFT_ALL_Real')
        call system(excitation_string)
        write(excitation_string, '(a)') trim('rm ' // OutDir // '/Dipole/Dipole_State_01/DipoleFT_ww_Real')
        call system(excitation_string)
        do i_excitation = 2, nStates
            write(excitation_string, '(a,i0.2, a)') 'rm ' // OutDir // '/Dipole/Dipole_State_', i_excitation, '/DipoleFT_ALL_Real'
            call system(excitation_string)
            !            write(excitation_string, '(a,i0.2, a)') 'rm ' // OutDir // '/Dipole/Dipole_State_', i_excitation, '/DipoleFT_ALL_Imaginary'
            !            call system(excitation_string)
            write(excitation_string, '(a,i0.2, a)') 'rm ' // OutDir // '/Dipole/Dipole_State_', i_excitation, '/DipoleFT_ww_Real'
            call system(excitation_string)
            write(excitation_string, '(a,i0.2, a)') 'rm ' // OutDir // '/Dipole/Dipole_State_', i_excitation, '/DipoleFT_ww_Imaginary'
            call system(excitation_string)
        end do
        i_excitation = 0
        icomponent_RI = 0
    else !parallel
        write(excitation_string, '(a,i0.2)') OutDir // '/Dipole/Dipole_State_', i_excitation
        call system("rm " // excitation_string // "/DipoleFT_ALL")
        call system("rm " // excitation_string // "/DipoleFT_ww")
    end if

    write(*, *) nStates
    write(*, *) i_excitation, icomponent_RI

    !#######################################################################!
    allocate(Tau_Time_Vec(N_Simulations - 1))
    !
    Mean = 0.0
    ppos = scan(trim(Simulation_tagv(1)), "PP", BACK = .true.)! compute mean
    do iSim = 1, N_Simulations - 1
        !
        !..Get Tau_Time and append to Tau_Time_Vec
        output_file = trim(ADJUSTL(Simulation_tagv(iSim)(ppos + 1:)))
        read(output_file, *) Tau_Time !convert string to real number
        Tau_Time_Vec(iSim) = Tau_Time ! append to Tau_Time_Vec
        !
        Mean = Mean + Tau_Time
        !
    end do
    Mean = Mean / (N_Simulations - 1)
    !
    Variance = 0.0                       ! compute variance
    do iSim = 1, (N_Simulations - 1)
        !
        Variance = Variance + (Tau_Time_Vec(iSim) - Mean)**2
        !
    end do
    Variance = Variance / (N_Simulations - 1)
    !
    StdDev = SQRT(Variance)         ! compute standard deviation
    !#######################################################################!
    !
    !    do iSim = 1, (N_Simulations - 1)
    !        MASK = ((Tau_Time_Vec(iSim) - Tau_Time_Vec(1)) / StdDev)
    !        write(*, *)MASK, (Tau_Time_Vec(iSim) - Tau_Time_Vec(1))
    !    end do
    !    write(*, *) "Mean:", Mean, "Var:", Variance, "StdDev:", StdDev
    !    open(newunit = uid, &
    !            file = "MASK_FUNCTION", &
    !            form = "formatted", &
    !            status = "unknown", &
    !            action = "write")
    !    do iSim = 1, (N_Simulations - 1)
    !        MASK = (Tau_Time_Vec(iSim) - Tau_Time_Vec(1)) / StdDev
    !        write(*, *)Tau_Time_Vec(iSim), MASK
    !        write(uid, "(*(x,e24.16))") Tau_Time_Vec(iSim), MASK, Variance, StdDev, (Variance**(1 / 2))
    !    end do
    !    close(uid)
    !    stop
    !
    !###########################################################################################
    !#  Loads and Regularizes the Response to the XUV alone for each molecular excitation (Re and Im)
    allocate(XUVDipole_I(3, ntimes, (nStates - 1) * 2, 2))
    !
    if (i_excitation ==0 .and. icomponent_RI==0)  then !Non_Parallel
        XUV_Excitation_Loop : do i_excitation = 1, nStates
            !
            do icomponent_RI = 1, 2
                !
                if ((i_excitation == GS_IDX) .and. (icomponent_RI == 2)) cycle
                !
                if (icomponent_RI == 1) PART = "Real"
                if (icomponent_RI == 2) PART = "Imaginary"
                !
                write(excitation_string, '(a,i0.2)') OutDir // '/Dipole/Dipole_State_', i_excitation
                !#
                !############################################################################################
                !.. Load XUVDipole from file, Regularize it and Compute FT of XUVDipole
                call Load_XUVDipole(trim(excitation_string) // "/Dipole" // trim(Simulation_tagv(N_Simulations)) // "_" // trim(PART), XUVDipole, nTimes)
                call Regularize_XUVDipole(XUVDipole, tmin, dt, nTimes, StepTime, StepWidth)
                !############################################################################################
                !#
                XUVDipole_I(:, :, i_excitation, icomponent_RI) = XUVDipole
                !
            end do
            !
        end do XUV_Excitation_Loop
        i_excitation = 0
        icomponent_RI = 0
    else
        if (icomponent_RI == 1) PART = "Real"
        if (icomponent_RI == 2) PART = "Imaginary"
        !
        if ((i_excitation == GS_IDX) .and. (icomponent_RI == 2)) stop !stop the rest of the code for imaginary ground state
        !
        write(excitation_string, '(a,i0.2)') OutDir // '/Dipole/Dipole_State_', i_excitation
        !#
        !############################################################################################
        !.. Load XUVDipole from file, Regularize it and Compute FT of XUVDipole
        call Load_XUVDipole(trim(excitation_string) // "/Dipole" // trim(Simulation_tagv(N_Simulations)) // "_" // trim(PART), XUVDipole, nTimes)
        call Regularize_XUVDipole(XUVDipole, tmin, dt, nTimes, StepTime, StepWidth)
        !############################################################################################
        !#
        XUVDipole_I(:, :, i_excitation, icomponent_RI) = XUVDipole
    end if
    !#
    !############################################################################################

    !..Start Sim_loop
    !..
    allocate(DipoleFTwt   (3, nOmegas, N_Simulations))
    allocate(DipoleFTplus (3, nOmegas))
    allocate(DipoleFTminus(3, nOmegas))
    allocate(DipoleFTtotal(3, nOmegas))
    allocate(zMuEV        (3, nTimes))
    !
    N_Simulations = N_Simulations - 1 !Do not Compute Sim = XUV alone
    !
    !############################################################################################
    !#
    if (i_excitation ==0 .and. icomponent_RI==0) then !Non_Parallel
        Single_Excitation_Loop : do i_excitation = 1, nStates
            !
            do icomponent_RI = 1, 2
                !
                if ((i_excitation == GS_IDX) .and. (icomponent_RI == 2)) cycle
                !
                if (icomponent_RI == 1) PART = "Real"
                if (icomponent_RI == 2) PART = "Imaginary"
                !
                write(excitation_string, '(a,i0.2)') OutDir // '/Dipole/Dipole_State_', i_excitation
                write(*, '(a,i0.2,a)') "Starting Sim Loop for Molecular Excitation ", i_excitation, " " // trim(PART)
                !
                !############################################################################################
                !#
                Sim_loop : do iSim = 1, N_Simulations
                    write(*, *) "State", i_excitation, PART
                    write(*, *) "iSimFT =", iSim, "   N_Simulations =", N_Simulations
                    !..
                    call LoadDipoles(trim(excitation_string) // "/Dipole" // trim(Simulation_tagv(iSim)) // "_" // trim(PART), nTimes, zMuEV)
                    !.. Compute the regularized dipole $\mu_-(t)$ and Atomic Charges and FT
                    call Regularize_Dipole(zMuEV, tmin, dt, nTimes, StepTime, StepWidth)
                    zMuEV = zMuEV - XUVDipole_I(:, :, i_excitation, icomponent_RI)
                    call ComputeFT_Dipole(DipoleFTminus, zMuEV, OmegaVec, tmin, dt, nTimes, nOmegas)
                    !..Compute DipoleFTtotal
                    DipoleFTplus = Z0
                    DipoleFTtotal = (DipoleFTminus + DipoleFTplus)
                    !
                    !
                    !#######################################################################!
                    !                    MASK = (Tau_Time_Vec(iSim) - Tau_Time_Vec(1)) / StdDev
                    !                    if (MASK > 1) Mask = 1.
                    call Regularize_FTwrtTime(DipoleFTtotal, Tau_Time_Vec, iSim, StepTime, StepWidth)
                    !
                    DipoleFTtotal = DipoleFTtotal
                    !#######################################################################!
                    !
                    !.. Write the FT of Dipoles to a single file
                    call WriteAllDipoleFTtoSingleFile(trim(excitation_string) // "/DipoleFT_ALL_" // trim(PART), DipoleFTtotal, OmegaVec, &
                            nOmegas, iSim, train)
                end do Sim_loop
                !#
                !############################################################################################
                !
            end do
            !
        end do Single_Excitation_Loop
        i_excitation = 0
        icomponent_RI = 0
    else
        !
        if (icomponent_RI == 1) PART = "Real"
        if (icomponent_RI == 2) PART = "Imaginary"
        !
        write(excitation_string, '(a,i0.2)') OutDir // '/Dipole/Dipole_State_', i_excitation
        write(*, '(a,i0.2,a)') "1Starting Sim Loop for Molecular Excitation ", i_excitation, " " // trim(PART)
        !
        Sim_loop_parallel : do iSim = 1, N_Simulations
            write(*, *) "State", i_excitation, PART
            write(*, *) "iSimFT =", iSim, "   N_Simulations =", N_Simulations
            !..
            call LoadDipoles(trim(excitation_string) // "/Dipole" // trim(Simulation_tagv(iSim)) // "_" // trim(PART), nTimes, zMuEV)
            !.. Compute the regularized dipole $\mu_-(t)$ and Atomic Charges and FT
            call Regularize_Dipole(zMuEV, tmin, dt, nTimes, StepTime, StepWidth)
            zMuEV = zMuEV - XUVDipole_I(:, :, i_excitation, icomponent_RI)
            call ComputeFT_Dipole(DipoleFTminus, zMuEV, OmegaVec, tmin, dt, nTimes, nOmegas)
            !..Compute DipoleFTtotal
            DipoleFTplus = Z0
            DipoleFTtotal = (DipoleFTminus + DipoleFTplus)
            !
            !
            !#######################################################################!
            MASK = (Tau_Time_Vec(iSim) - Tau_Time_Vec(1)) / StdDev
            if (MASK > 1) Mask = 1. ! only place momentarily as Mask kept grow y
            write(*, *) MASK
            DipoleFTtotal = DipoleFTtotal * MASK
            !#######################################################################!
            !
            !.. Write the FT of Dipoles to a single file
            call WriteAllDipoleFTtoSingleFile(trim(excitation_string) // "/DipoleFT_ALL_" // trim(PART), DipoleFTtotal, OmegaVec, &
                    nOmegas, iSim, train)
        end do Sim_loop_parallel
    end if
    !#
    !############################################################################################


    !############################################################################################
    !#
    if (i_excitation ==0 .and. icomponent_RI==0) then !Non_Parallel
        Single_Excitation_Loop_FT : do i_excitation = 1, nStates
            do icomponent_RI = 1, 2
                !
                if ((i_excitation == GS_IDX) .and. (icomponent_RI == 2)) cycle
                !
                if (icomponent_RI == 1) PART = "Real"
                if (icomponent_RI == 2) PART = "Imaginary"
                !
                write(excitation_string, '(a,i0.2)') OutDir // '/Dipole/Dipole_State_', i_excitation
                !############################################################################################
                !#
                !.. COMPUTE 2D SPECTRUM DIPOLE
                !
                !.. Load the FT of the dipole, as a function of the time delay
                call LoadFTDipole_asfuncitonof_TimeDelay(trim(excitation_string) // "/DipoleFT_ALL_" // trim(PART), N_simulations, nOmegas, DipoleFTwt, tvec)
                !
                !.. Regularizes the dipole with respect to the time-delay edges and Compute FT with respect to time
                call Regularize_Dipole_withRespectto_TimeDelay(DipoleFTwt, tvec, N_Simulations)
                call Compute_FTofDipole_withRespectto_Time(DipoleFTwt, DipoleFTww, tvec, TauOmegaVec, N_Simulations, nTauOmegas)
                !
                !.. Write the Bidimensional spectrum to file
                call SaveBidimentioal_Dipole_Spectrum(trim(excitation_string) // "/DipoleFT_ww_" // trim(PART), DipoleFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
                deallocate(tvec)
                !#
                !############################################################################################
            end do
        end do Single_Excitation_Loop_FT
        i_excitation = 0
        icomponent_RI = 0
    else
        !
        if (icomponent_RI == 1) PART = "Real"
        if (icomponent_RI == 2) PART = "Imaginary"
        !
        write(excitation_string, '(a,i0.2)') OutDir // '/Dipole/Dipole_State_', i_excitation
        !############################################################################################
        !#
        !.. COMPUTE 2D SPECTRUM DIPOLE
        !
        !.. Load the FT of the dipole, as a function of the time delay
        call LoadFTDipole_asfuncitonof_TimeDelay(trim(excitation_string) // "/DipoleFT_ALL_" // trim(PART), N_simulations, nOmegas, DipoleFTwt, tvec)
        !
        !.. Regularizes the dipole with respect to the time-delay edges and Compute FT with respect to time
        call Regularize_Dipole_withRespectto_TimeDelay(DipoleFTwt, tvec, N_Simulations)
        call Compute_FTofDipole_withRespectto_Time(DipoleFTwt, DipoleFTww, tvec, TauOmegaVec, N_Simulations, nTauOmegas)
        !
        !.. Write the Bidimensional spectrum to file
        call SaveBidimentioal_Dipole_Spectrum(trim(excitation_string) // "/DipoleFT_ww_" // trim(PART), DipoleFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
        deallocate(tvec)
        !#
        !############################################################################################
    end if
    !#
    !############################################################################################
    stop


contains


    !    subroutine GetAtomicChargeComponent (AtomicChargeEvolution, AtomicChargeEvolutionComponent, AtCoord, nAtoms, nTimes)
    !        real   (kind(1d0)), intent(in) :: AtomicChargeEvolution(:, :)
    !        complex(kind(1d0)) :: XUVCharge(:, :), XUVChargeComponent(:, :, :)
    !        real   (kind(1d0)), allocatable, intent(out) :: AtomicChargeEvolutionComponent(:, :, :)
    !        real(kind(1d0)), intent(in) :: AtCoord(:, :)
    !        integer, intent(in) :: nAtoms, ntimes
    !
    !        integer :: iAtom, iCoord, it, uid
    !
    !        allocate(AtomicChargeEvolutionComponent(3, nAtoms, nTimes))
    !        do iCoord = 1, 3
    !            do iAtom = 1, nAtoms
    !                do it = 1, nTimes
    !                    AtomicChargeEvolutionComponent(iCoord, iAtom, it) = AtomicChargeEvolution(iAtom, it) * AtCoord(iCoord, iAtom)
    !                end do
    !            end do
    !        end do
    !
    !    end subroutine GetAtomicChargeComponent

    subroutine DetermineNumberofTimes (FileName, tmin, tmax, nTimes)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(out) :: tmin, tmax
        integer, intent(out) :: nTimes

        real   (kind(1d0)) :: dBuf
        integer :: uid_dipole, iBuf, iostat

        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")
        !
        read(uid_dipole, *, iostat = iostat) iBuf, dBuf
        if(iostat/=0)then
            write(*, *) FileName // "   IS EMPTY"
            stop
        end if
        tmin = dBuf
        tmax = tmin
        nTimes = 1
        do
            read(uid_dipole, *, iostat = iostat) iBuf, dBuf
            if(iostat/=0)exit
            tmax = dBuf
            nTimes = nTimes + 1
        enddo
        close(uid_dipole)
        write(*, *) "nTimes = ", nTimes
    end subroutine DetermineNumberofTimes


    subroutine Regularize_XUVDipole(XUVDipole, tmin, dt, nTimes, StepTime, StepWidth)
        complex(kind(1d0)), intent(inout) :: XUVDipole(:, :)
        real(kind(1d0)), intent(in) :: tmin, dt, StepTime, StepWidth
        integer, intent(in) :: nTimes

        real(kind(1d0)), external :: NCD_Phi
        integer :: it, uid_NCD_phi
        real(kind(1d0)) :: dstep, t, StepT

        !        open(newunit = uid_NCD_phi, &
        !                file = "NCD_phi", &
        !                form = "formatted", &
        !                status = "unknown", &
        !                action = "write")
        !
        !
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            StepT = StepTime
            if(abs(t - StepTime)<=dt / 2.d0)then
                StepT = t
            endif
            dstep = 1.d0 - NCD_Phi(t, StepT, StepWidth)
            !
            !            write(uid_NCD_phi, "(*(x,e24.16))")t, dstep
            !
            XUVDipole(:, it) = XUVDipole(:, it) * dstep
        enddo
        !
        !
        !        close(uid_NCD_phi)
    end subroutine Regularize_XUVDipole

    subroutine Regularize_FTwrtTime(Dipole, Tau_Time_Vec, iSim, StepTime, StepWidth)
        complex(kind(1d0)), intent(inout) :: Dipole(:, :)
        real(kind(1d0)), intent(in) :: Tau_Time_Vec(:), StepTime, StepWidth
        integer, intent(in) :: iSim

        real(kind(1d0)), external :: NCD_Phi
        integer :: uid_NCD_phi, nTau
        real(kind(1d0)) :: dstep, t, StepT, TAU_CURRENT, TAU_MIN, TAU_MAX, dTau

        nTau = size(Tau_Time_Vec)
        !
        TAU_MIN = Tau_Time_Vec(1)
        TAU_MAX = Tau_Time_Vec(nTau)
        TAU_CURRENT = Tau_Time_Vec(iSim)
        !
        dTau = (tmax - tmin) / dble(nTimes - 1)
        !
        StepT = StepTime
        !        if(abs(TAU_CURRENT - StepTime)<=dTau / 2.d0)then
        !            StepT = TAU_CURRENT
        !        endif
        dstep = 1.d0 - NCD_Phi(TAU_CURRENT, 100.d0, 100.d0)
        !
        Dipole = Dipole * dstep

    end subroutine Regularize_FTwrtTime


    subroutine Regularize_XUVAtomicCharge(XUVCharge, tmin, dt, nTimes, StepTime, StepWidth)
        complex(kind(1d0)), intent(inout) :: XUVCharge(:, :)
        real(kind(1d0)), intent(in) :: tmin, dt, StepTime, StepWidth
        integer, intent(in) :: nTimes

        real(kind(1d0)), external :: NCD_Phi
        integer :: it
        real(kind(1d0)) :: dstep, t, StepT

        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            StepT = StepTime
            if(abs(t - StepTime)<=dt / 2.d0)then
                StepT = t
            endif
            dstep = 1.d0 - NCD_Phi(t, StepT, StepWidth)
            XUVCharge(:, it) = XUVCharge(:, it) * dstep
        enddo
    end subroutine Regularize_XUVAtomicCharge
    !
    !
    subroutine Regularize_XUVAtomicCharge1(Charge, tmin, dt, nTimes, StepTime, StepWidth)
        complex(kind(1d0)), intent(inout) :: Charge(:, :, :)
        real(kind(1d0)), intent(in) :: tmin, dt, StepTime, StepWidth
        integer, intent(in) :: nTimes

        real(kind(1d0)), external :: NCD_Phi
        integer :: it
        real(kind(1d0)) :: dstep, t, StepT

        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            StepT = StepTime
            if(abs(t - StepTime)<=dt / 2.d0)then
                StepT = t
            endif
            dstep = 1.d0 - NCD_Phi(t, StepT, StepWidth)
            do iPol = 1, 3
                Charge(iPol, :, it) = Charge(iPol, :, it) * dstep
            end do
        enddo
    end subroutine Regularize_XUVAtomicCharge1
    !
    !
    subroutine ComputeFT_Dipole(DipoleFTminus, Dipole, OmegaVec, tmin, dt, nTimes, nOmegas)
        complex(kind(1d0)), intent(out) :: DipoleFTminus(:, :)
        complex(kind(1d0)), intent(in) :: Dipole(:, :)
        real   (kind(1d0)), intent(in) :: OmegaVec(:)
        real(kind(1d0)), intent(in) :: tmin, dt
        integer, intent(in) :: nTimes, nOmegas

        complex(kind(1d0)) :: zExpFact
        real(kind(1d0)) :: t, w
        integer :: iOmega, it, iPol
        !
        if (.not.allocated(AtomicChargeFT_new))allocate(AtomicChargeFT_new(3, nOmegas, nAtoms))
        !
        DipoleFTminus = Z0
        !
        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            do it = 1, nTimes
                t = tmin + dt * dble(it - 1)
                zexpFact = exp(Zi * w * t)
                DipoleFTminus(:, iOmega) = DipoleFTminus(:, iOmega) + zExpFact * Dipole(:, it)
            enddo
        end do
        !
        DipoleFTminus = DipoleFTminus * dt / (2.d0 * PI)
    end subroutine ComputeFT_Dipole


    subroutine Regularize_Dipole(Dipole, tmin, dt, nTimes, StepTime, StepWidth)
        complex(kind(1d0)), intent(inout) :: Dipole(:, :)
        real(kind(1d0)), intent(in) :: tmin, dt, StepTime, StepWidth
        integer, intent(in) :: nTimes

        real(kind(1d0)), external :: NCD_Phi
        integer :: it
        real(kind(1d0)) :: dstep, t, StepT

        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            StepT = StepTime
            if(abs(t - StepTime)<=dt / 2.d0)then
                StepT = t
            endif
            dstep = 1.d0 - NCD_Phi(t, StepT, StepWidth)
            Dipole(:, it) = Dipole(:, it) * dstep
        end do
    end subroutine Regularize_Dipole

    !############################################################################################
    !..Bidimentional Spectrum Subroutines
    subroutine Regularize_Dipole_withRespectto_TimeDelay(DipoleFTwt, tvec, N_Simulations)
        complex(kind(1d0)), intent(inout) :: DipoleFTwt(:, :, :)
        real   (kind(1d0)), intent(in) :: tvec(:)
        integer, intent(in) :: N_Simulations

        real(kind(1d0)), external :: NCD_Phi
        integer :: iSim
        real(kind(1d0)) :: dstep, t, tStep1, tStep2, StepWidth

        !        tStep1 = tvec(N_Simulations / 20)
        !        tStep2 = tvec(9 * N_Simulations / 20)
        !        StepWidth = (tvec(N_Simulations) - tvec(1)) / 30.d0

        tStep1 = tvec(max(nint(0.2d0 * dble(N_Simulations)), 1))
        tStep2 = tvec(nint(8.d0 * N_Simulations / 10.d0))
        StepWidth = (tvec(N_Simulations) - tvec(1)) / 15.d0

        do iSim = 1, N_Simulations
            t = tvec(iSim)
            dstep = NCD_Phi(t, tStep1, StepWidth) * (1.d0 - NCD_Phi(t, tStep2, StepWidth))
            DipoleFTwt(:, :, iSim) = dstep * DipoleFTwt(:, :, iSim)
        enddo
    end subroutine Regularize_Dipole_withRespectto_TimeDelay


    subroutine Compute_FTofDipole_withRespectto_Time(DipoleFTwt, DipoleFTww, tvec, TauOmegaVec, N_Simulations, nTauOmegas)
        complex(kind(1d0)), intent(in) :: DipoleFTwt(:, :, :)
        complex(kind(1d0)), allocatable, intent(out) :: DipoleFTww(:, :, :)
        real   (kind(1d0)), intent(in) :: tvec(:)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:)
        integer, intent(in) :: N_Simulations, nTauOmegas

        complex(kind(1d0)) :: zExpFact
        integer :: iOmega, iSim
        real(kind(1d0)) :: t, w, dt

        allocate(DipoleFTww(3, nOmegas, nTauOmegas))

        DipoleFTww = Z0
        dt = tvec(2) - tvec(1)
        do iOmega = 1, nTauOmegas
            w = TauOmegaVec(iOmega)
            do iSim = 1, N_Simulations
                t = tvec(iSim)
                zexpFact = exp(Zi * w * t)
                DipoleFTww(:, :, iOmega) = DipoleFTww(:, :, iOmega) + zExpFact * DipoleFTwt(:, :, iSim)
            end do
        enddo
        DipoleFTww = DipoleFTww * dt / (2.d0 * PI)!((2.d0 * PI)**(1 / 2))
    end subroutine Compute_FTofDipole_withRespectto_Time


    !.. Regularizes the charge with respect to the time-delay edges
    !..

    subroutine Regularize_Charge_withtimedelay (ChargeFTwt, tvec, nAtoms, N_Simulations)
        complex(kind(1d0)), intent(inout) :: ChargeFTwt(:, :, :)
        real   (kind(1d0)), intent(in) :: tvec(:)
        integer, intent(in) :: nAtoms, N_Simulations

        real(kind(1d0)), external :: NCD_Phi
        real(kind(1d0)) :: dstep, t, tStep1, tStep2, StepWidth
        integer :: iAtom, iSim

        tStep1 = tvec(max(nint(0.2d0 * dble(N_Simulations)), 1))
        tStep2 = tvec(nint(8.d0 * N_Simulations / 10.d0))
        StepWidth = (tvec(N_Simulations) - tvec(1)) / 15.d0

        do iAtom = 1, nAtoms
            do iSim = 1, N_Simulations
                t = tvec(iSim)
                dstep = NCD_Phi(t, tStep1, StepWidth) * (1.d0 - NCD_Phi(t, tStep2, StepWidth))
                ChargeFTwt(:, iSim, iAtom) = dstep * ChargeFTwt(:, iSim, iAtom)
            enddo
        enddo

    end subroutine Regularize_Charge_withtimedelay


    subroutine Compute_FTCharge_withTimeDelay(ChargeFTwt, tvec, ChargeFTww, N_simulations, nOmegas, nAtoms, nTauOmegas, TauOmegaVec)
        complex(kind(1d0)), intent(in) :: ChargeFTwt(:, :, :)
        real   (kind(1d0)), intent(in) :: tvec(:), TauOmegaVec(:)
        complex(kind(1d0)), allocatable, intent(out) :: ChargeFTww(:, :, :)
        integer, intent(in) :: N_simulations, nOmegas, nAtoms, nTauOmegas

        complex(kind(1d0)) :: zExpFact
        real(kind(1d0)) :: dt, w, t
        integer :: iAtom, iOmega, iSim

        allocate(ChargeFTww(nOmegas, nTauOmegas, nAtoms))
        dt = (tvec(N_simulations) - tvec(1)) / dble(N_simulations - 1)
        do iAtom = 1, nAtoms
            do iOmega = 1, nTauOmegas
                w = TauOmegaVec(iOmega)
                do iSim = 1, N_Simulations
                    t = tvec(iSim)
                    zexpFact = exp(Zi * w * t)
                    !!$           if(iAtom==1)write(*,*) iSim, iOmega, t, w
                    ChargeFTww(:, iOmega, iAtom) = ChargeFTww(:, iOmega, iAtom) + &
                            zExpFact * ChargeFTwt(:, iSim, iAtom)
                end do
            enddo
        enddo
        ChargeFTww = ChargeFTww * dt / (2.d0 * PI)! ((2.d0 * PI)**(1 / 2))
    end subroutine Compute_FTCharge_withTimeDelay

    !#########################
    !    trick
    subroutine LoadMadeUpDipoles(FileName, N_simulations, nOmegas, DipoleFTwt, tvec, OmegaVec)!, TauOmegaVec)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: N_simulations, nOmegas
        complex(kind(1d0)), allocatable, intent(out) :: DipoleFTwt(:, :, :)
        real   (kind(1d0)), allocatable, intent(out) :: tvec(:), OmegaVec(:)!, TauOmegaVec(:)

        real   (kind(1d0)) :: drx
        integer :: iSim, iOmega, uid_dipoleFT
        integer :: iostat
        character(len = 1000) :: iomsg

        allocate(DipoleFTwt(3, nOmegas, N_Simulations))
        allocate(tvec(N_Simulations))
        allocate(OmegaVec(nOmegas))
        !allocate(TauOmegaVec(nOmegas))

        open(newunit = uid_dipoleFT, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read", &
                iostat = iostat, &
                iomsg = iomsg)
        read(uid_dipoleFT, *)
        do iSim = 1, N_simulations
            do iOmega = 1, nOmegas
                read(uid_dipoleFT, *, iostat = iostat) tvec(iSim), OmegaVec(iOmega), drx
                if(iostat/=0)exit
                DipoleFTwt(1, iOmega, iSim) = Z1 * drx
                DipoleFTwt(2, iOmega, iSim) = Z1 * drx
                DipoleFTwt(3, iOmega, iSim) = Z1 * drx
            end do

            read(uid_dipoleFT, *, iostat = iostat)
            if(iostat/=0)exit
        enddo
        close(uid_dipoleFT)

        !
        !do iOmega = 1, size(OmegaVec)
        !  TauOmegaVec(iOmega) = OmegaVec(iOmega)
        !end do
        close(uid_dipoleFT)
    end subroutine LoadMadeUpDipoles
    !#########################


    subroutine Load_Q_Charge_and_Write2FTALL1(FileName, Charge, nTimes, tmin, dt, nAtoms, iSim, uid_AtomicChargeALL)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), allocatable, intent(out) :: Charge(:, :, :)
        real   (kind(1d0)), intent(in) :: tmin, dt
        integer, intent(in) :: nTimes, nAtoms, iSim, uid_AtomicChargeALL

        real   (kind(1d0)) :: t, dBuf
        integer :: uid_AtomicCharge, iPol, it, iAtom

        if (.not.allocated(Charge))allocate(Charge(3, nAtoms, nTimes))

        !.. Save Q_Charge
        open(newunit = uid_AtomicCharge, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")

        !..New
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            read(uid_AtomicCharge, "(*(x,e24.16))") dBuf, dBuf, ((Charge(iPol, iAtom, it), iPol = 1, 3), iAtom = 1, nAtoms)
            write(uid_AtomicChargeALL, "(*(x,e24.16))") t, dble(iSim), ((Charge(iPol, iAtom, it), iPol = 1, 3), iAtom = 1, nAtoms)
        enddo
        !
        write(uid_AtomicChargeALL, *)
        close(uid_AtomicCharge)
    end subroutine Load_Q_Charge_and_Write2FTALL1

    subroutine WriteAllAtomicChargeFTtoSingleFile1(FileName, AtomicChargeFT_new, OmegaVec, nOmegas, nAtoms, iSim, train)
        complex(kind(1d0)) :: AtomicChargeFT_new(:, :, :)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nOmegas, nAtoms, iSim
        type(pulse_train), pointer, intent(in) :: train(:)

        real   (kind(1d0)) :: OmegaVec(:)
        real   (kind(1d0)) :: w
        integer :: uid_AtomicChargeFT, iOmega, iAtom, iPol

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write", &
                position = "append")
        do iOmega = 1, nOmegas
            call train(iSim)%PrintPulses(uid_AtomicChargeFT)
            w = OmegaVec(iOmega)
            write(uid_AtomicChargeFT, "(i4,*(x,E24.16))") iOmega, w, ((&
                    (dble(AtomicChargeFT_new(iPol, iOmega, iAtom)), &
                    aimag(AtomicChargeFT_new(iPol, iOmega, iAtom))), iPol = 1, 3), &
                    iAtom = 1, nAtoms)
        end do
        write(uid_AtomicChargeFT, *)
        close(uid_AtomicChargeFT)
    end subroutine WriteAllAtomicChargeFTtoSingleFile1

    subroutine LoadFTofChargeasFuncofTimeDelay1(FileName, N_simulations, nOmegas, nAtoms, tvec, ChargeFTwt_new)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: N_simulations, nOmegas, nAtoms
        real   (kind(1d0)), allocatable, intent(out) :: tvec(:)
        complex(kind(1d0)), allocatable, intent(out) :: ChargeFTwt_new(:, :, :, :)

        real   (kind(1d0)), allocatable :: dAtomFTRe(:, :), dAtomFTIm(:, :)
        integer :: uid_AtomicChargeFT, iSim, iOmega, iBuf, iAtom, i, iPol
        real(kind(1d0)) :: dBuf

        allocate(tvec(N_simulations))
        tvec = 0.d0
        allocate(dAtomFTRe(3, nAtoms))
        allocate(dAtomFTIm(3, nAtoms))
        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")
        allocate(ChargeFTwt_new(3, nOmegas, N_Simulations, nAtoms))
        do iSim = 1, N_simulations
            do iOmega = 1, nOmegas
                read(uid_AtomicChargeFT, *)iBuf, (dBuf, i = 1, 7), tvec(iSim), (dBuf, i = 1, 6), iBuf, dBuf, &
                (((dAtomFTRe(iPol, iAtom), dAtomFTIm(iPol, iAtom)), iPol = 1, 3), iAtom = 1, nAtoms)
                do iPol = 1, 3
                    do iAtom = 1, nAtoms
                        ChargeFTwt_new(iPol, iOmega, iSim, iAtom) = Z1 * dAtomFTRe(iPol, iAtom) + Zi * dAtomFTIm(iPol, iAtom)
                    enddo
                end do
            enddo
        end do
        read(uid_AtomicChargeFT, *)
        deallocate(dAtomFTRe, dAtomFTIm)
        close(uid_AtomicChargeFT)
    end subroutine LoadFTofChargeasFuncofTimeDElay1

    subroutine Regularize_Charge_withtimedelay1 (ChargeFTwt_new, tvec, nAtoms, N_Simulations)
        complex(kind(1d0)), intent(inout) :: ChargeFTwt_new(:, :, :, :)
        real   (kind(1d0)), intent(in) :: tvec(:)
        integer, intent(in) :: nAtoms, N_Simulations

        real(kind(1d0)), external :: NCD_Phi
        real(kind(1d0)) :: dstep, t, tStep1, tStep2, StepWidth
        integer :: iAtom, iSim

        tStep1 = tvec(max(nint(0.2d0 * dble(N_Simulations)), 1))
        tStep2 = tvec(nint(8.d0 * N_Simulations / 10.d0))
        StepWidth = (tvec(N_Simulations) - tvec(1)) / 15.d0

        do iAtom = 1, nAtoms
            do iSim = 1, N_Simulations
                t = tvec(iSim)
                dstep = NCD_Phi(t, tStep1, StepWidth) * (1.d0 - NCD_Phi(t, tStep2, StepWidth))
                do iPol = 1, 3
                    ChargeFTwt_new(iPol, :, iSim, iAtom) = dstep * ChargeFTwt_new(iPol, :, iSim, iAtom)
                end do
            enddo
        enddo
    end subroutine Regularize_Charge_withtimedelay1

    subroutine Compute_FTCharge_withTimeDelay1(ChargeFTwt_new, tvec, ChargeFTww_new, N_simulations, nOmegas, nAtoms, nTauOmegas, TauOmegaVec)
        complex(kind(1d0)), intent(in) :: ChargeFTwt_new(:, :, :, :)
        real   (kind(1d0)), intent(in) :: tvec(:), TauOmegaVec(:)
        complex(kind(1d0)), allocatable, intent(out) :: ChargeFTww_new(:, :, :, :)
        integer, intent(in) :: N_simulations, nOmegas, nAtoms, nTauOmegas

        complex(kind(1d0)) :: zExpFact
        real(kind(1d0)) :: dt, w, t
        integer :: iAtom, iOmega, iSim

        allocate(ChargeFTww_new(3, nOmegas, nTauOmegas, nAtoms))
        dt = (tvec(N_simulations) - tvec(1)) / dble(N_simulations - 1)
        do iAtom = 1, nAtoms
            do iOmega = 1, nTauOmegas
                w = TauOmegaVec(iOmega)
                do iSim = 1, N_Simulations
                    t = tvec(iSim)
                    zexpFact = exp(Zi * w * t)
                    !!$           if(iAtom==1)write(*,*) iSim, iOmega, t, w
                    do iPol = 1, 3
                        ChargeFTww_new(iPol, :, iOmega, iAtom) = ChargeFTww_new(iPol, :, iOmega, iAtom) + &
                                zExpFact * ChargeFTwt_new(iPol, :, iSim, iAtom)
                    end do
                end do
            enddo
        enddo
        ChargeFTww_new = ChargeFTww_new * dt / (2.d0 * PI)                  ! ((2.d0 * PI)**(1 / 2))
    end subroutine Compute_FTCharge_withTimeDelay1

    subroutine Write_BidimentionalChargeFTww1(FileName, ChargeFTww_new, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), intent(in) :: ChargeFTww_new(:, :, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas, nAtoms

        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iAtom, iPol

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                write(uid_AtomicChargeFT, "(*(x,e24.16))") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
                        (((dble(ChargeFTww_new(iPol, iOmega, iOmegaTau, iAtom)), &
                        (aimag(ChargeFTww_new(iPol, iOmega, iOmegaTau, iAtom)))), iPol = 1, 3), iAtom = 1, nAtoms)
            end do
            write(uid_AtomicChargeFT, *)
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Write_BidimentionalChargeFTww1
end program ChargeMigration

!!!!!this are notes and aspirations .. :]

!    DipoleFTminus = DipoleFTminus * dt / (2.d0 * PI)
!     AtomicChargeFT=AtomicChargeFT * dt / (2.d0 * PI)
!
!*** MUST GO TO A SEPARATE SUBROUTINE
!!$     !.. Save Dipole FT File
!!$     open(newunit = uid_dipoleFT, &
!!$          file    ="DipoleFT-"//trim(Simulation_tagv(iSim)), &
!!$          form    ="formatted", &
!!$          status  ="unknown"  , &
!!$          action  ="write"    )
!!$     do iOmega = 1, nOmegas
!!$        w = OmegaVec( iOmega )
!!$        write(uid_dipoleFT,"(i4,*(x,E24.16))") iOmega, w, &
!!$             ((dble(DipoleFTminus(iPol,iOmega)),aimag(DipoleFTminus(iPol,iOmega))),iPol=1,3)
!!$     end do
!!$     close( uid_dipoleFT )


!*** THE FT OF THE DIPOLE FOR LARGE TIME STILL DOES NOT WORK
!     DipoleFTplus = Z0
!!$     call ComputeDipoleFTplus( L0_Eval, L0_LEvec, L0_REvec, &
!!$          OmegaVec, zStatRho0, DipoleFTplus, StepTime, StepWidth )
!!$     !
!!$     !*** MUST GO TO A SEPARATE SUBROUTINE
!!$     !.. Save Dipole FT File
!!$     open(newunit = uid_dipoleFT, &
!!$          file    ="DipoleFT+"//trim(Simulation_tagv(iSim)), &
!!$          form    ="formatted", &
!!$          status  ="unknown"  , &
!!$          action  ="write"    )
!!$     do iOmega = 1, nOmegas
!!$        w = OmegaVec( iOmega )
!!$        write(uid_dipoleFT,"(i4,*(x,E24.16))") iOmega, w, &
!!$             ((dble(DipoleFTplus(iPol,iOmega)),aimag(DipoleFTplus(iPol,iOmega))),iPol=1,3)
!!$     end do
!!$     close( uid_dipoleFT )


!     DipoleFTtotal = DipoleFTminus + DipoleFTplus
