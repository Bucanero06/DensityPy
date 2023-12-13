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
program ChargeMigrationFT

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
    complex(kind(1d0)), allocatable :: DipoleFTwt(:, :, :), DipoleFTww(:, :, :)
    complex(kind(1d0)), allocatable :: AtomicChargeFT(:, :, :), ChargeFTwt(:, :, :, :), ChargeFTww(:, :, :, :), XUVCharge(:, :, :)

    real   (kind(1d0)), allocatable :: OmegaVec(:), TauOmegaVec(:)
    real   (kind(1d0)), allocatable :: tvec(:)
    integer :: iOmega, ntimes, i
    real(kind(1d0)) :: tmin, tmax

    !.. Local parameters
    !..
    character(len = *), parameter :: DIPOLE_FT_PATH_ALL = "/Dipole/DipoleFT_ALL.csv"
    character(len = *), parameter :: CHARGE_FT_PATH_ALL = "/AtomicCharge/AtomicChargeFT_ALL.csv"
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
    !.. Pulse parameters
    !..
    integer :: iSim, N_Simulations
    character(len = 64), pointer :: Simulation_Tagv(:)
    type(pulse_train), pointer :: train(:)

    !.. XUV_Dipole and XUV_Charge
    !..
    complex(kind(1d0)), allocatable :: XUVDipole(:, :)
    complex(kind(1d0)), allocatable :: XUVDipoleFT(:, :)



    !..Test
    !..
    real   (kind(1d0)), allocatable :: AtomicChargeEvolution(:, :, :)
    complex(kind(1d0)), allocatable :: Debug(:, :, :)
    external :: system
    integer :: iPol, iAtom, it, iCoord

    call GetRunTimeParameters(InpDir, OutDir, FileGeometry, StepTime, StepWidth, &
            Ext_field_file, verbous, nOmegas, OmegaMin, OmegaMax, nTauOmegas, TauOmegaMin, TauOmegaMax)

    ! These need to be deleted before running the program since they are appended
    call system("rm " // OutDir // DIPOLE_FT_PATH_ALL)
    call system("rm " // OutDir // CHARGE_FT_PATH_ALL)
    ! These are just deleted for completeness
    call system("rm " // OutDir // "/Dipole/DipoleFT_ww.csv")
    call system("rm " // OutDir // "/AtomicCharge/AtomicChargeFT_ww.csv")

    allocate(OmegaVec, source = [(OmegaMin + (OmegaMax - OmegaMin) / dble(nOmegas - 1) * dble(iOmega - 1), iOmega = 1, nOmegas)])
    allocate(TauOmegaVec, source = [(TauOmegaMin + (TauOmegaMax - TauOmegaMin) / dble(nTauOmegas - 1) * dble(iOmega - 1), iOmega = 1, nTauOmegas)])

    call Set_CD_IO_Verbous(.FALSE.)

    open(newunit = uid, &
            file = Ext_field_file, &
            form = "formatted", &
            status = "old")
    call Parse_Simulation_File(uid, OUTPUT_UNIT, N_Simulations, Simulation_Tagv, train)
    close(uid)
    !    write(*, "(a,i3)")"N_Simulations=", N_Simulations

    !.. Loads Root Ernergies and Molecular Geometry
    call LoadEnergies(InpDir // "/ROOT_ENERGIES", nStates, Evec)
    call LoadGeometry(nAtoms, AtCoord, FileGeometry, AtomName)

    !.. Determines the number of times available, as well as the minimum and maximum time
    call DetermineNumberofTimes (OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(1)) // '.csv', tmin, tmax, nTimes)
    dt = (tmax - tmin) / dble(nTimes - 1)

    Write(*, *) "tmin=", tmin, "tmax=", tmax, "nTimes=", nTimes, "dt=", dt

    !#
    !############################################################################################
    !.. Load XUVDipole from file, Regularize it and Compute FT of XUVDipole
    call Load_Dipole(OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(N_Simulations)) // '.csv', XUVDipole, nTimes)
    call Regularize_XUVDipole(XUVDipole, tmin, dt, nTimes, StepTime, StepWidth)
    !    call Compute_XUVDipoleFT(XUVDipole, XUVDipoleFT, OmegaVec, tmin, dt, nTimes, nOmegas)
    !############################################################################################
    !#
    !############################################################################################
    !.. Load XUVCharge from file, Regularize it and Compute FT of XUVCharge
    call Load_XUVAtomicCharge(OutDir // "/AtomicCharge/AtomicCharge" // trim(Simulation_tagv(N_Simulations)) // '.csv', XUVCharge, nTimes, nAtoms)
    call Regularize_XUVAtomicCharge(XUVCharge, tmin, dt, nTimes, StepTime, StepWidth)
    !##########################################################################################
    !#
    !############################################################################################
    !.. Open Q_Charge File
    open(newunit = uid_AtomicChargeALL, &
            file = OutDir // "/AtomicCharge/AtomicCharge_ALL.csv", &
            form = "formatted", &
            status = "unknown", &
            action = "write")

    !..Redefine N_Simulatins (Last simulation is the XUV with no IR)
    N_Simulations = N_Simulations - 1
    !
    !..Start Sim_loop
    !..
    allocate(DipoleFTwt   (3, nOmegas, N_Simulations))
    allocate(DipoleFTplus (3, nOmegas))
    allocate(DipoleFTminus(3, nOmegas))
    allocate(DipoleFTtotal(3, nOmegas))
    allocate(zMuEV        (3, nTimes))
    allocate(AtomicChargeEvolution(3, nAtoms, nTimes)) !delete
    allocate(AtomicChargeFT(3, nOmegas, nAtoms))
    Sim_loop : do iSim = 1, N_Simulations
        write(*, *) "iSim =", iSim, "   N_Simulations =", N_Simulations

        !..
        call Load_Dipole(OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(iSim)) // '.csv', zMuEV, nTimes)

        call Load_Q_Charge_and_Write2ALL(OutDir // "/AtomicCharge/AtomicCharge" // trim(Simulation_tagv(iSim)) // '.csv', &
                AtomicChargeEvolution, nTimes, tmin, dt, nAtoms, iSim, AtomName, uid_AtomicChargeALL)


        !.. Compute the regularized dipole $\mu_-(t)$ and Atomic Charges and FT
        call Regularize_Dipole_and_AtomicCharge1(zMuEV, AtomicChargeEvolution, tmin, dt, nTimes, StepTime, StepWidth)

        zMuEV = zMuEV - XUVDipole
        AtomicChargeEvolution = AtomicChargeEvolution - XUVCharge!!!

        call ComputeFT_Dipole_and_AtomicCharges(DipoleFTminus, AtomicChargeFT, zMuEV, OmegaVec, &
                AtomicChargeEvolution, tmin, dt, nTimes, nOmegas)

        !..Compute DipoleFTtotal
        DipoleFTplus = Z0
        DipoleFTtotal = (DipoleFTminus + DipoleFTplus)

        !.. Save FT of Dipoles
        call WriteDipoleFTFile (OutDir // "/Dipole/DipoleFT" // trim(Simulation_tagv(iSim)) // ".csv", &
                DipoleFTtotal, OmegaVec, nOmegas)

        !.. Write the FT of Dipoles to a single file
        call AppendDipole2FTAllFile(OutDir // DIPOLE_FT_PATH_ALL, DipoleFTtotal, OmegaVec, &
                nOmegas, iSim, train)

        !.. Save Atomic Charge FT
        call WriteAtomicChargeFT(OutDir // "/AtomicCharge/AtomicChargeFT" // trim(Simulation_tagv(iSim)) // ".csv", &
                AtomicChargeFT, OmegaVec, nOmegas, nAtoms, AtomName)
        !
        !.. Save Atomic Charge FT in a single file
        call WriteAllAtomicChargeFTtoSingleFile(OutDir // CHARGE_FT_PATH_ALL, AtomicChargeFT, OmegaVec, &
                nOmegas, nAtoms, iSim, train, AtomName)

    end do Sim_loop
    close(uid_AtomicChargeALL)
    !############################################################################################
    !#

    !#
    !############################################################################################
    !.. COMPUTE 2D SPECTRUM DIPOLE
    !
    !.. Load the FT of the dipole, as a function of the time delay
    call LoadFTDipole_asfuncitonof_TimeDelay(OutDir // DIPOLE_FT_PATH_ALL, N_simulations, nOmegas, DipoleFTwt, tvec)
    !
    !.. Regularizes the dipole with respect to the time-delay edges and Compute FT with respect to time
    call Regularize_Dipole_withRespectto_TimeDelay(DipoleFTwt, tvec, N_Simulations)
    call Compute_FTofDipole_withRespectto_Time(DipoleFTwt, DipoleFTww, tvec, TauOmegaVec, N_Simulations, nTauOmegas)
    !
    !.. Write the Bidimensional spectrum to file
    call SaveBidimentioal_Dipole_Spectrum(OutDir // "/Dipole/DipoleFT_ww.csv", DipoleFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
    deallocate(tvec)
    !############################################################################################
    !#

    !#
    !############################################################################################
    !.. COMPUTE 2D SPECTRUM CHARGE
    !..
    !.. Load the FT of the charge, as a functAtomicChargeion of the time delay
    call LoadFTofChargeasFuncofTimeDelay(OutDir // CHARGE_FT_PATH_ALL, N_simulations, nOmegas, nAtoms, tvec, ChargeFTwt)
    !

    !.. Regularizes the charge with respect to the time-delay edges
    call Regularize_Charge_withtimedelay (ChargeFTwt, tvec, nAtoms, N_Simulations)

    !.. Compute the FT wrt the time delay
    call Compute_FTCharge_withTimeDelay(ChargeFTwt, tvec, ChargeFTww, N_simulations, nOmegas, nAtoms, nTauOmegas, TauOmegaVec)

    !.. Write the Bidimensional spectrum to file
    !    call Write_BidimentionalChargeFTww(OutDir // "/AtomicCharge/AtomicChargeFT_ww", ChargeFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
    call Write_BidimentionalChargeFTww(OutDir // "/AtomicCharge/AtomicChargeFT_ww.csv", ChargeFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms, AtomName)
    !############################################################################################
    !#

    stop


contains

    subroutine DetermineNumberOfTimes(FileName, tmin, tmax, nTimes)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(out) :: tmin, tmax
        integer, intent(out) :: nTimes

        real(kind(1d0)) :: dBuf
        integer :: uid_dipole, iostat, iBuf
        character(200) :: line    ! buffer to read the header

        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")

        ! Read header line and ignore it
        read(uid_dipole, '(a)') line

        read(uid_dipole, *, iostat = iostat) iBuf, dBuf
        if (iostat /= 0) then
            write(*, *) FileName // " IS EMPTY"
            stop
        end if

        tmin = dBuf
        tmax = tmin
        nTimes = 1

        do
            read(uid_dipole, *, iostat = iostat) iBuf, dBuf
            if (iostat /= 0) exit
            tmax = dBuf
            nTimes = nTimes + 1
        end do

        close(uid_dipole)
        write(*, *) "nTimes = ", nTimes
    end subroutine DetermineNumberOfTimes


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

    subroutine Regularize_XUVAtomicCharge(Charge, tmin, dt, nTimes, StepTime, StepWidth)
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
    end subroutine Regularize_XUVAtomicCharge


    subroutine Regularize_Dipole_and_AtomicCharge1(Dipole, Charge, tmin, dt, nTimes, StepTime, StepWidth)
        complex(kind(1d0)), intent(inout) :: Dipole(:, :)
        real(kind(1d0)), intent(inout) :: Charge(:, :, :)
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
            do iPol = 1, 3
                Charge(iPol, :, it) = Charge(iPol, :, it) * dstep
            enddo
        end do
    end subroutine Regularize_Dipole_and_AtomicCharge1


    subroutine ComputeFT_Dipole_and_AtomicCharges(DipoleFTminus, AtomicChargeFT, Dipole, OmegaVec, &
            AtomicChargeEvolution, tmin, dt, nTimes, nOmegas)
        complex(kind(1d0)), intent(out) :: DipoleFTminus(:, :)
        complex(kind(1d0)), allocatable, intent(out) :: AtomicChargeFT(:, :, :)
        complex(kind(1d0)), intent(in) :: Dipole(:, :)
        real   (kind(1d0)), intent(in) :: OmegaVec(:), AtomicChargeEvolution(:, :, :)
        real(kind(1d0)), intent(in) :: tmin, dt
        integer, intent(in) :: nTimes, nOmegas

        complex(kind(1d0)) :: zExpFact
        real(kind(1d0)) :: t, w
        integer :: iOmega, it, iPol

        if (.not.allocated(AtomicChargeFT))allocate(AtomicChargeFT(3, nOmegas, nAtoms))

        DipoleFTminus = Z0
        AtomicChargeFT = Z0

        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            do it = 1, nTimes
                t = tmin + dt * dble(it - 1)
                zexpFact = exp(Zi * w * t)
                DipoleFTminus(:, iOmega) = DipoleFTminus(:, iOmega) + zExpFact * Dipole(:, it)
                do iPol = 1, 3
                    AtomicChargeFT(iPol, iOmega, :) = AtomicChargeFT(iPol, iOmega, :) + &
                            zExpFact * AtomicChargeEvolution(iPol, :, it)
                end do
!            WRITE(*, *) "AtomicChargeEvolution(iPol, :, it)=", AtomicChargeEvolution(1, 1, it)
            enddo
        end do
        DipoleFTminus = DipoleFTminus * dt / (2.d0 * PI)                 !((2.d0 * PI)**(1 / 2))
        AtomicChargeFT = AtomicChargeFT * dt / (2.d0 * PI)       !((2.d0 * PI)**(1 / 2))
    end subroutine ComputeFT_Dipole_and_AtomicCharges
    !
    !

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
        complex(kind(1d0)), intent(inout) :: ChargeFTwt(:, :, :, :)
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
                    ChargeFTwt(iPol, :, iSim, iAtom) = dstep * ChargeFTwt(iPol, :, iSim, iAtom)
                end do
            enddo
        enddo
    end subroutine Regularize_Charge_withtimedelay


    subroutine Compute_FTCharge_withTimeDelay(ChargeFTwt, tvec, ChargeFTww, N_simulations, nOmegas, nAtoms, nTauOmegas, TauOmegaVec)
        complex(kind(1d0)), intent(in) :: ChargeFTwt(:, :, :, :)
        real   (kind(1d0)), intent(in) :: tvec(:), TauOmegaVec(:)
        complex(kind(1d0)), allocatable, intent(out) :: ChargeFTww(:, :, :, :)
        integer, intent(in) :: N_simulations, nOmegas, nAtoms, nTauOmegas

        complex(kind(1d0)) :: zExpFact
        real(kind(1d0)) :: dt, w, t
        integer :: iAtom, iOmega, iSim

        allocate(ChargeFTww(3, nOmegas, nTauOmegas, nAtoms))
        dt = (tvec(N_simulations) - tvec(1)) / dble(N_simulations - 1)
        do iAtom = 1, nAtoms
            do iOmega = 1, nTauOmegas
                w = TauOmegaVec(iOmega)
                do iSim = 1, N_Simulations
                    t = tvec(iSim)
                    zexpFact = exp(Zi * w * t)
                    !!$           if(iAtom==1)write(*,*) iSim, iOmega, t, w
                    do iPol = 1, 3
                        ChargeFTww(iPol, :, iOmega, iAtom) = ChargeFTww(iPol, :, iOmega, iAtom) + &
                                zExpFact * ChargeFTwt(iPol, :, iSim, iAtom)
                    end do
                end do
            enddo
        enddo
        ChargeFTww = ChargeFTww * dt / (2.d0 * PI)                  ! ((2.d0 * PI)**(1 / 2))
    end subroutine Compute_FTCharge_withTimeDelay

end program ChargeMigrationFT

!!!!!this are notes and aspirations .. :]
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



!!*** APPARENTLY, IT IS NOT WORKING YET!
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

