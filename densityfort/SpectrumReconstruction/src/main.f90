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

    !.. Run-time parameters
    !..
    character(len = :), allocatable :: InpDir, OutDir, FileGeometry
    real(kind(1d0)) :: StepTime, StepWidth
    character(len = :), allocatable :: Ext_Field_File
    logical :: Verbous
    integer :: nOmegas, nTauOmegas
    real   (kind(1d0)) :: OmegaMin, OmegaMax, TauOmegaMin, TauOmegaMax

    !.. Data for FT
    !..
    complex(kind(1d0)), allocatable :: DipoleFTtotal(:, :), DipoleFTminus(:, :), DipoleFTplus (:, :)
    complex(kind(1d0)), allocatable :: AtomicChargeFT(:, :)
    complex(kind(1d0)), allocatable :: DipoleFTwt(:, :, :), DipoleFTww(:, :, :), Reconstructed_Dipole(:, :, :)
    complex(kind(1d0)), allocatable :: ChargeFTwt(:, :, :), ChargeFTww(:, :, :), ChargeFTwwComponent(:, :, :, :), ChargeFTww_new(:, :, :, :)
    real   (kind(1d0)), allocatable :: OmegaVec(:), TauOmegaVec(:)
    real   (kind(1d0)), allocatable :: tvec(:)
    integer :: iOmega, iOmegaTau, ntimes, i, iAtom
    real(kind(1d0)) :: tmin, tmax

    !.. Local parameters
    !..
    character(len = *), parameter :: DIPOLE_FT_PATH_ALL = "/Dipole/DipoleFT_ALL"
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
    real   (kind(1d0)), allocatable :: AtomicChargeEvolution(:, :)

    !.. Pulse parameters
    !..
    integer :: iSim, N_Simulations
    character(len = 64), pointer :: Simulation_Tagv(:)
    type(pulse_train), pointer :: train(:)

    !.. XUV_Dipole
    !..
    complex(kind(1d0)), allocatable :: XUVDipole(:, :)
    complex(kind(1d0)), allocatable :: XUVDipoleFT(:, :), DipoleDifference(:, :, :)

    !    trick
    complex(kind(1d0)), allocatable :: Debug(:, :, :)
    external :: system

    call GetRunTimeParameters(InpDir, OutDir, FileGeometry, StepTime, StepWidth, &
            Ext_field_file, verbous, nOmegas, OmegaMin, OmegaMax, nTauOmegas, TauOmegaMin, TauOmegaMax)

    call system("rm " // OutDir // "/Dipole/DipoleFT_ww_reconstructed")

    allocate(OmegaVec, source = [(OmegaMin + (OmegaMax - OmegaMin) / dble(nOmegas - 1) * dble(iOmega - 1), iOmega = 1, nOmegas)])
    allocate(TauOmegaVec, source = [(TauOmegaMin + (TauOmegaMax - TauOmegaMin) / dble(nTauOmegas - 1) * dble(iOmega - 1), iOmega = 1, nTauOmegas)])

    call LoadGeometry(nAtoms, AtCoord, FileGeometry, AtomName)


    !>Load Dipole and Charge ww
    call LoadBidimentioal_Dipole_Spectrum(OutDir // "/Dipole/DipoleFT_ww", DipoleFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)

    !    call Load_BidimentionalChargeFTww(OutDir // "/AtomicCharge/AtomicChargeFT_ww", ChargeFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
    call Load_BidimentionalChargeFTww1(OutDir // "/AtomicCharge/AtomicChargeFT_ww", ChargeFTww_new, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
    write(*, *) sum(ChargeFTww_new(1, 5, 5, :))
    !    call GetBiDimentionalAtomicChargeComponent (ChargeFTww, ChargeFTwwComponent, AtCoord, nAtoms, nOmegas, nTauOmegas)
    !    call Write_BidimentionalChargeFTwwComponent(OutDir // "/AtomicCharge/AtomicChargeFT_ww_", ChargeFTwwComponent, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)

    !>Reconstruct and Save Dipole from Charge
    !    call ReconstructDipole_from_AtomicCharge_times_XYZ (ChargeFTwwComponent, Reconstructed_Dipole, AtCoord, nAtoms, nOmegas, nTauOmegas)
    call ReconstructDipole_from_AtomicCharge_times_XYZ1 (ChargeFTww_new, Reconstructed_Dipole, AtCoord, nAtoms, nOmegas, nTauOmegas)
    write(*, *) sum(ChargeFTww_new(1, 5, 5, :))
    call Write_2DReconstructDipole(OutDir // "/Dipole/DipoleFT_ww_reconstructed", Reconstructed_Dipole, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
    !
    !    allocate(DipoleDifference(3, nOmegas, nTauOmegas))
    !    !    DipoleDifference = DipoleFTww - Reconstructed_Dipole
    !    call Write_DipoleDifference(OutDir // "/Dipole/DipoleFT_ww_difference", DipoleFTww, Reconstructed_Dipole, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
    !    call Compute_Charge_Coefficients
contains
    subroutine Write_BidimentionalChargeFTww(FileName, ChargeFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), intent(in) :: ChargeFTww(:, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas, nAtoms

        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iAtom

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                write(uid_AtomicChargeFT, "(*(x,e24.16))") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
                        ((dble(ChargeFTww(iOmega, iOmegaTau, iAtom)), &
                        (aimag(ChargeFTww(iOmega, iOmegaTau, iAtom)))), iAtom = 1, nAtoms)
            end do
            write(uid_AtomicChargeFT, *)
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Write_BidimentionalChargeFTww

    !    subroutine Write_DipoleDifference(FileName, DipoleFTww, Reconstructed_Dipole, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
    !        character(len = *), intent(in) :: FileName
    !        complex(kind(1d0)), intent(in) :: DipoleFTww(:, :, :), Reconstructed_Dipole(:, :, :)
    !        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
    !        integer, intent(in) :: nTauOmegas, nOmegas
    !
    !        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iAtom, iCoord
    !
    !        open(newunit = uid_AtomicChargeFT, &
    !                file = FileName, &
    !                form = "formatted", &
    !                status = "unknown", &
    !                action = "write")
    !        do iOmegaTau = 1, nTauOmegas
    !            do iOmega = 1, nOmegas
    !                write(uid_AtomicChargeFT, "(*(x,e24.16))") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
    !                        ((dble(DipoleFTww(iCoord, iOmega, iOmegaTau)), &
    !                        aimag(DipoleFTww(iCoord, iOmega, iOmegaTau))), iCoord = 1, 3), ((dble(Reconstructed_Dipole(iCoord, iOmega, iOmegaTau)), &
    !                        aimag(Reconstructed_Dipole(iCoord, iOmega, iOmegaTau))), iCoord = 1, 3)
    !            end do
    !            write(uid_AtomicChargeFT, *)
    !        enddo
    !        close(uid_AtomicChargeFT)
    !    end subroutine Write_DipoleDifference

    subroutine ReconstructDipole_from_AtomicCharge_times_XYZ (ChargeFTwwComponent, Reconstructed_Dipole, AtCoord, nAtoms, nOmegas, nTauOmegas)
        complex(kind(1d0)), intent(in) :: ChargeFTwwComponent(:, :, :, :)
        real(kind(1d0)), intent(in) :: AtCoord(:, :)
        complex(kind(1d0)), allocatable, intent(out) :: Reconstructed_Dipole(:, :, :)
        integer, intent(in) :: nAtoms, nOmegas, nTauOmegas

        integer :: iAtom, iCoord, iOmega, iOmegaTau, uid

        allocate(Reconstructed_Dipole(3, nOmegas, nTauOmegas))
        do iCoord = 1, 3
            do iAtom = 1, nAtoms
                do iOmega = 1, nOmegas
                    do iOmegaTau = 1, nTauOmegas
                        Reconstructed_Dipole(iCoord, iOmega, iOmegaTau) = Reconstructed_Dipole(iCoord, iOmega, iOmegaTau) &
                                + ChargeFTwwComponent(iCoord, iOmega, iOmegaTau, iAtom)
                    end do
                end do
            end do
        end do
    end subroutine ReconstructDipole_from_AtomicCharge_times_XYZ

    subroutine ReconstructDipole_from_AtomicCharge_times_XYZ1 (ChargeFTww_new, Reconstructed_Dipole, AtCoord, nAtoms, nOmegas, nTauOmegas)
        complex(kind(1d0)), intent(in) :: ChargeFTww_new(:, :, :, :)
        real(kind(1d0)), intent(in) :: AtCoord(:, :)
        complex(kind(1d0)), allocatable, intent(out) :: Reconstructed_Dipole(:, :, :)
        integer, intent(in) :: nAtoms, nOmegas, nTauOmegas

        integer :: iAtom, iPol, iOmega, iOmegaTau, uid

        allocate(Reconstructed_Dipole(3, nOmegas, nTauOmegas))
        !        do iAtom = 1, nAtoms
        do iOmega = 1, nOmegas
            do iOmegaTau = 1, nTauOmegas
                do iPol = 1, 3
                    Reconstructed_Dipole(iPol, iOmega, iOmegaTau) = sum(ChargeFTww_new(iPol, iOmega, iOmegaTau, :))
                    !                        Reconstructed_Dipole(iPol, iOmega, iOmegaTau) = Reconstructed_Dipole(iPol, iOmega, iOmegaTau) &
                    !                                + ChargeFTww_new(iPol, iOmega, iOmegaTau, iAtom)
                end do
            end do
        end do
        !        end do

    end subroutine ReconstructDipole_from_AtomicCharge_times_XYZ1

    !    subroutine Compute_Charge_Coefficients (ChargeFTww,ACoeff, nAtoms, nOmegas, nTauOmegas)
    !        complex(kind(1d0)), intent(in) :: ChargeFTww(:, :, :, :)
    !        real(kind(1d0)), intent(in) :: AtCoord(:, :)
    !
    !        integer, intent(in) :: nAtoms, nOmegas, nTauOmegas
    !
    !
    !
    !    end subroutine Compute_Charge_Coefficients

    subroutine GetBiDimentionalAtomicChargeComponent (ChargeFTww, ChargeFTwwComponent, AtCoord, nAtoms, nOmegas, nTauOmegas)
        complex(kind(1d0)), intent(in) :: ChargeFTww(:, :, :)
        real(kind(1d0)), intent(in) :: AtCoord(:, :)
        complex(kind(1d0)), allocatable, intent(out) :: ChargeFTwwComponent(:, :, :, :)
        integer, intent(in) :: nAtoms, nOmegas, nTauOmegas

        integer :: iAtom, iCoord, iOmega, iOmegaTau, uid

        allocate(ChargeFTwwComponent(3, nOmegas, nTauOmegas, nAtoms))
        do iAtom = 1, nAtoms
            do iCoord = 1, 3
                do iOmegaTau = 1, nTauOmegas
                    do iOmega = 1, nOmegas
                        ChargeFTwwComponent(iCoord, iOmega, iOmegaTau, iAtom) = ChargeFTww(iOmega, iOmegaTau, iAtom) * AtCoord(iCoord, iAtom)
                    end do
                end do
            end do
        end do
    end subroutine GetBiDimentionalAtomicChargeComponent


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

    !    subroutine DetermineFrequencies (FileName, OmegaMin, OmegaMax, nOmegas, TauOmegaMin, TauOmegaMax, nTauOmegas)
    !        character(len = *), intent(in) :: FileName
    !        integer, intent(out) :: nTauOmegas, nOmegas
    !        real   (kind(1d0)), intent(out) :: TauOmegaMin, TauOmegaMax, OmegaMin, OmegaMax
    !
    !        real   (kind(1d0)) :: dBuf1, dBuf2, dBuf
    !        integer :: uid, iBuf, iostat
    !
    !        open(newunit = uid, &
    !                file = FileName, &
    !                form = "formatted", &
    !                status = "old", &
    !                action = "read")
    !        !
    !        read(uid, *, iostat = iostat) dBuf1, dBuf2, dBuf, dBuf, dBuf, dBuf, dBuf, dBuf
    !        if(iostat/=0)then
    !            write(*, *) FileName // "   IS EMPTY"
    !            stop
    !        end if
    !        TauOmegaMin = dBuf1
    !        TauOmegaMax = TauOmegaMin
    !        nTauOmegas = 1
    !
    !        OmegaMin = dBuf2
    !        OmegaMax = OmegaMin
    !        nOmegas = 1
    !        do
    !            read(uid, *, iostat = iostat) dBuf1, dBuf2
    !            if(iostat/=0)exit
    !            TauOmegaMax = dBuf1
    !            nTauOmegas = nTauOmegas + 1
    !
    !            OmegaMax = dBuf2
    !            nOmegas = nOmegas + 1
    !        enddo
    !        close(uid)
    !        write(*, *) "nOmegas = ", nOmegas, "nTauOmegas = ", nTauOmegas
    !end subroutine DetermineFrequencies


end program ChargeMigration

