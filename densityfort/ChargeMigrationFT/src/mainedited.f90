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
    !.. Local modules
    use ModuleRTP
    use Module_CD_IO
    implicit none
    !.. Run-time parameters
    character(len = :), allocatable :: InpDir, OutDir, FileGeometry, Ext_field_file
    real   (kind(1d0)) :: OmegaMin, OmegaMax, TauOmegaMin, TauOmegaMax, StepTime, StepWidth
    integer :: nOmegas, nTauOmegas
    logical :: Verbous
    !.. Data for FT
    complex(kind(1d0)), allocatable :: DipoleFTtotal(:, :), DipoleFTminus(:, :), DipoleFTplus (:, :)
    complex(kind(1d0)), allocatable :: DipoleFTwt(:, :, :), DipoleFTww(:, :, :), ReplaceDipoleFTwt(:, :, :)
    real   (kind(1d0)), allocatable :: OmegaVec(:), TauOmegaVec(:), tvec(:)
    real(kind(1d0)) :: tmin, tmax
    integer :: iOmega, ntimes, i, iPol
    !.. Local parameters
    character(len = *), parameter :: DIPOLE_FT_PATH_ALL = "/Dipole/DipoleFT_ALL"
    integer, parameter :: GS_IDX = 1
    integer :: nStates
    integer :: uid
    real(kind(1d0)), allocatable :: Evec(:)
    real(kind(1d0)), external :: NCD_Phi
    complex(kind(1d0)), external :: zdotu
    !.. Expectation Value of the Dipole Moment (Mu)
    complex(kind(1d0)), allocatable :: zMuEV(:, :)
    !.. Statistical Density Matrix
    real   (kind(1d0)) :: dt
    !.. Pulse parameters
    integer :: iSim, N_Simulations
    character(len = 64), pointer :: Simulation_Tagv(:)
    type(pulse_train), pointer :: train(:)
    !.. XUV_Dipole
    complex(kind(1d0)), allocatable :: XUVDipole(:, :)
    complex(kind(1d0)), allocatable :: XUVDipoleFT(:, :)

    external :: system
    call GetRunTimeParameters(InpDir, OutDir, FileGeometry, StepTime, StepWidth, &
            Ext_field_file, verbous, nOmegas, OmegaMin, OmegaMax, nTauOmegas, TauOmegaMin, TauOmegaMax)
    allocate(OmegaVec, source = [(OmegaMin + (OmegaMax - OmegaMin) / dble(nOmegas - 1) * dble(iOmega - 1), iOmega = 1, nOmegas)])
    allocate(TauOmegaVec, source = [(TauOmegaMin + (TauOmegaMax - TauOmegaMin) / dble(nTauOmegas - 1) * dble(iOmega - 1), iOmega = 1, nTauOmegas)])
    call Set_CD_IO_Verbous(.FALSE.)

    !.. Parse Simulation File
    open(newunit = uid, &
            file = Ext_field_file, &
            form = "formatted", &
            status = "old")
    call Parse_Simulation_File(uid, OUTPUT_UNIT, N_Simulations, Simulation_Tagv, train)
    close(uid)


    !.. Loads Root Ernergies and Molecular Geometry
    !    call LoadEnergies(InpDir // "/ROOT_ENERGIES", nStates, Evec)

    !.. Determines the number of times available, as well as the minimum and maximum time
    call DetermineNumberofTimes (OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(1)), tmin, tmax, nTimes)

    !    call DetermineNumberofTimes ("GaussTestDipole/Dipole" // trim(Simulation_tagv(1)), tmin, tmax, nTimes)!debug

    dt = (tmax - tmin) / dble(nTimes - 1)

    !.. Load XUVDipole from file, Regularize it and Compute FT of XUVDipole
    call Load_XUVDipole(OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(N_Simulations)), XUVDipole, nTimes)!debug uncomment after
    call Regularize_XUVDipole(XUVDipole, tmin, dt, nTimes, StepTime, StepWidth)
!    call Compute_XUVDipoleFT(XUVDipole, XUVDipoleFT, OmegaVec, tmin, dt, nTimes, nOmegas)


    call system("rm " // OutDir // "/Dipole/DipoleFT_ALL")
    call system("rm " // OutDir // "/Dipole/DipoleFT_ww")


    !..Redefine N_Simulatins (Last simulation is the XUV with no IR)
    N_Simulations = N_Simulations - 1
    !..Start Sim_loop
    !..
    allocate(DipoleFTwt   (3, nOmegas, N_Simulations))
    allocate(DipoleFTplus (3, nOmegas))
    allocate(DipoleFTminus(3, nOmegas))
    allocate(DipoleFTtotal(3, nOmegas))
    allocate(zMuEV        (3, nTimes ))

    Sim_loop : do iSim = 1, N_Simulations
        write(*, *) "iSim =", iSim, "   N_Simulations =", N_Simulations

        !..Load Dipole from File
        call LoadDipoles(OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(iSim)), nTimes, zMuEV)

        !        call LoadDipoles("GaussTestDipole/Dipole" // trim(Simulation_tagv(iSim)), nTimes, zMuEV) !debug

        !.. Compute the regularized dipole $\mu_-(t)$ and Atomic Charges and FT
        call Regularize_Dipole_and_AtomicCharge(zMuEV, tmin, dt, nTimes, StepTime, StepWidth)

        zMuEV=zMuEV - XUVDipole

        call ComputeFT_Dipole_and_AtomicCharges(DipoleFTminus, zMuEV, OmegaVec, tmin, dt, nTimes, nOmegas)

        !..Compute DipoleFTtotal
        DipoleFTplus = Z0
        DipoleFTtotal = (DipoleFTminus + DipoleFTplus) !- XUVDipoleFT
        !.. Save FT of Dipoles
!        call SaveDipoleFTFile (OutDir // "/Dipole/DipoleFT" // trim(Simulation_tagv(iSim)), &
!                DipoleFTtotal, OmegaVec, nOmegas)

        !        !.. Write the FT of Dipoles to a single file
        call WriteAllDipoleFTtoSingleFile(OutDir // DIPOLE_FT_PATH_ALL, DipoleFTtotal, OmegaVec, &
                nOmegas, iSim, train)

    end do Sim_loop



    !    !.. COMPUTE 2D SPECTRUM DIPOLE
    !    !..
    !.. Load the FT of the dipole, as a function of the time delay
    !    allocate(ReplaceDipoleFTwt(3, nOmegas, N_Simulations)) !debug
    call LoadFTDipole_asfuncitonof_TimeDelay(OutDir // DIPOLE_FT_PATH_ALL, N_simulations, nOmegas, DipoleFTwt, tvec)
    !    !        !*** HOTFIX
    !    !    do iSim=1,N_Simulations
    !    !        do iOmega=1,nOmegas
    !    !            do iPol=1,3
    !    !                if (OmegaVec(iOmega) < 0.30d0 .or. OmegaVec(iOmega) > 0.55d0) DipoleFTwt(iPol,iOmega,iSim)=0.d0
    !    !            end do
    !    !        end do
    !    !    end do
    !    !
    !    !
    !    !.. Regularizes the dipole with respect to the time-delay edges and Compute FT with respect to time
    call Regularize_Dipole_withRespectto_TimeDelay(DipoleFTwt, tvec, N_Simulations)

    call Compute_FTofDipole_withRespectto_Time(DipoleFTwt, DipoleFTww, tvec, TauOmegaVec, N_Simulations, nTauOmegas)
    !.. Write the Bidimensional spectrum to file
    call SaveBidimentioal_Dipole_Spectrum(OutDir // "/Dipole/DipoleFT_ww", DipoleFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
    deallocate(tvec)

    stop
contains


    subroutine DetermineNumberofTimes (FileName, tmin, tmax, nTimes)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(out) :: tmin, tmax
        integer, intent(out) :: nTimes

        real   (kind(1d0)) :: dBuf
        integer :: uid_dipole, iBuf, iostat, i

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
        integer :: it
        real(kind(1d0)) :: dstep, t, StepT

        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            StepT = StepTime
            if(abs(t - StepTime)<=dt / 2.d0)then
                StepT = t
            endif
            dstep = 1.d0 - NCD_Phi(t, StepT, StepWidth)
            XUVDipole(:, it) = XUVDipole(:, it) * dstep
        enddo
    end subroutine Regularize_XUVDipole


    subroutine Compute_XUVDipoleFT(XUVDipole, XUVDipoleFT, OmegaVec, tmin, dt, nTimes, nOmegas)
        complex(kind(1d0)), intent(in) :: XUVDipole(:, :)
        complex(kind(1d0)), allocatable, intent(out) :: XUVDipoleFT(:, :)
        real   (kind(1d0)), intent(in) :: OmegaVec(:)
        real(kind(1d0)), intent(in) :: tmin, dt
        integer, intent(in) :: nTimes, nOmegas

        integer :: iOmega, it
        real   (kind(1d0)) :: w, t
        complex(kind(1d0)) :: zExpFact

        allocate(XUVDipoleFT(3, nOmegas))

        XUVDipoleFT = Z0
        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            w = OmegaVec(iOmega)
            do it = 1, nTimes
                t = tmin + dt * dble(it - 1)
                zexpFact = exp(Zi * w * t)
                XUVDipoleFT(:, iOmega) = XUVDipoleFT(:, iOmega) + zExpFact * XUVDipole(:, it)
            enddo
        end do
        XUVDipoleFT = XUVDipoleFT * dt / ((2.d0 * PI)**(1 / 2))
    end subroutine Compute_XUVDipoleFT


    subroutine Regularize_Dipole_and_AtomicCharge(Dipole, tmin, dt, nTimes, StepTime, StepWidth)
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
        enddo
    end subroutine Regularize_Dipole_and_AtomicCharge


    subroutine ComputeFT_Dipole_and_AtomicCharges(DipoleFTminus, Dipole, OmegaVec, &
            tmin, dt, nTimes, nOmegas)
        complex(kind(1d0)), intent(out) :: DipoleFTminus(:, :)
        complex(kind(1d0)), intent(in) :: Dipole(:, :)
        real   (kind(1d0)), intent(in) :: OmegaVec(:)
        real(kind(1d0)), intent(in) :: tmin, dt
        integer, intent(in) :: nTimes, nOmegas

        complex(kind(1d0)) :: zExpFact
        real(kind(1d0)) :: t, w
        integer :: iOmega, it

        DipoleFTminus = Z0
        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            do it = 1, nTimes
                t = tmin + dt * dble(it - 1)
                zexpFact = exp(Zi * w * t)
                DipoleFTminus(:, iOmega) = DipoleFTminus(:, iOmega) + zExpFact * Dipole(:, it)
            enddo
        end do
        DipoleFTminus = DipoleFTminus * dt / ((2.d0 * PI)**(1 / 2))
    end subroutine ComputeFT_Dipole_and_AtomicCharges

    !############################################################################################
    !..Bidimentional Spectrum Subroutines
    subroutine Regularize_Dipole_withRespectto_TimeDelay(DipoleFTwt, tvec, N_Simulations)
        complex(kind(1d0)), intent(inout) :: DipoleFTwt(:, :, :)
        real   (kind(1d0)), intent(in) :: tvec(:)
        integer, intent(in) :: N_Simulations

        real(kind(1d0)), external :: NCD_Phi
        integer :: iSim
        real(kind(1d0)) :: dstep, t, tStep1, tStep2, StepWidth

        tStep1 = tvec(N_Simulations / 20)
        tStep2 = tvec(9 * N_Simulations / 20)
        StepWidth = (tvec(N_Simulations) - tvec(1)) / 30.d0
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
        DipoleFTww = DipoleFTww * dt / ((2.d0 * PI)**(1 / 2))
    end subroutine Compute_FTofDipole_withRespectto_Time

end program ChargeMigration

