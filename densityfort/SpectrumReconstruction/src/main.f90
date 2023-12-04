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
    character(len = :), allocatable :: Ext_field_file
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

    call system("rm " // OutDir // "/Dipole/DipoleFT_ww_reconstructed.csv")

    allocate(OmegaVec, source = [(OmegaMin + (OmegaMax - OmegaMin) / dble(nOmegas - 1) * dble(iOmega - 1), iOmega = 1, nOmegas)])
    allocate(TauOmegaVec, source = [(TauOmegaMin + (TauOmegaMax - TauOmegaMin) / dble(nTauOmegas - 1) * dble(iOmega - 1), iOmega = 1, nTauOmegas)])

    call LoadGeometry(nAtoms, AtCoord, FileGeometry, AtomName)

    !>Load Dipole and Charge ww
    call LoadBidimentioal_Dipole_Spectrum(OutDir // "/Dipole/DipoleFT_ww.csv", DipoleFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)

    call Load_BidimentionalChargeFTww(OutDir // "/AtomicCharge/AtomicChargeFT_ww.csv", ChargeFTww_new, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)

    !>Reconstruct and Save Dipole from Charge
    call ReconstructDipole_from_AtomicCharge_times_XYZ (ChargeFTww_new, Reconstructed_Dipole, AtCoord, nAtoms, nOmegas, nTauOmegas)

    write(*, *) sum(ChargeFTww_new(1, 5, 5, :))
    call Write_2DReconstructDipole(OutDir // "/Dipole/DipoleFT_ww_reconstructed.csv", Reconstructed_Dipole, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
    !
contains



    subroutine ReconstructDipole_from_AtomicCharge_times_XYZ (ChargeFTww_new, Reconstructed_Dipole, AtCoord, nAtoms, nOmegas, nTauOmegas)
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
                end do
            end do
        end do
        !        end do

    end subroutine ReconstructDipole_from_AtomicCharge_times_XYZ





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

