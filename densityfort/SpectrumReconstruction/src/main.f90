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
program SpectrumReconstruction

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
    use Module_SR_RTP
    use Module_SR_CD_IO

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
    complex(kind(1d0)), allocatable :: AtomicChargeFT(:, :)
    complex(kind(1d0)), allocatable :: DipoleFTww(:, :, :), Reconstructed_Dipole(:, :, :)
    complex(kind(1d0)), allocatable :: ChargeFTww(:, :, :), ChargeFTww_new(:, :, :, :)
    real   (kind(1d0)), allocatable :: OmegaVec(:), TauOmegaVec(:)
    integer :: iOmega, iOmegaTau, ntimes, i, iAtom
    real(kind(1d0)) :: tmin, tmax

    !.. Molecular Geometry
    !..
    integer :: nAtoms
    real(kind(1d0)), allocatable :: AtCoord(:, :) ! 3 x npts
    real(kind(1d0)), allocatable :: R_el_bc(:, :) ! 3 x npts
    character(len = 16), allocatable :: AtomName(:)

    real   (kind(1d0)) :: dt
    integer :: uid
    external :: system

    call GetRunTimeParameters(InpDir, OutDir, FileGeometry, StepTime, StepWidth, &
            Ext_field_file, verbous, nOmegas, OmegaMin, OmegaMax, nTauOmegas, TauOmegaMin, TauOmegaMax)

    call system("rm " // OutDir // "/Dipole/DipoleFT_ww_reconstructed.csv")

    allocate(OmegaVec, source = [(OmegaMin + (OmegaMax - OmegaMin) / dble(nOmegas - 1) * dble(iOmega - 1), iOmega = 1, nOmegas)])
    allocate(TauOmegaVec, source = [(TauOmegaMin + (TauOmegaMax - TauOmegaMin) / dble(nTauOmegas - 1) * dble(iOmega - 1), iOmega = 1, nTauOmegas)])

    call LoadGeometry(nAtoms, AtCoord, FileGeometry, AtomName)

    ! Load R el barycenter
    ! example file
    ! "Atom_Index","Atom_Name","X_Position","Y_Position","Z_Position"
    !    1,                   O,   -0.07126411307121,    2.30457942716860,    0.00000141047634
    !    2,                   N,    1.04185194584435,   -0.61416187727967,   -0.00000245456477
    !    3,                   C,    0.03487782751090,    0.04679493632641,   -0.00000115408771
    !    4,                   C,    2.16700544144725,   -0.62352078448538,   -0.00001000839462
    !    5,                   C,   -0.98099316311318,   -0.47892496257572,   -0.00000319061686
    !    6,                   H,    1.69427031327723,   -1.77124405319568,   -0.00000241500998
    !    7,                   H,    2.77985142247970,   -1.48736577051015,   -0.00000591538630
    !    8,                   H,   -1.65782286861781,    1.26578032027380,   -0.00000324475991
    !    9,                   H,   -1.11377254826595,   -1.12889045879090,   -1.50089541713585
    !   10,                   H,   -1.11377490986092,   -1.12888690420504,    1.50089663009806
    !   11,                   H,    1.93015526019949,    1.13437013574216,   -1.36058305350774
    !   12,                   H,    1.93018964279058,    1.13430325044931,    1.36062049358343

    open(newunit = uid, file = OutDir // "/R_el_bc.csv", status = "old", action = "read")
    allocate(R_el_bc(3, nAtoms))
    read(uid, *)
    do iAtom = 1, nAtoms
        read(uid, *) i, AtomName(iAtom), R_el_bc(1, iAtom), R_el_bc(2, iAtom), R_el_bc(3, iAtom)
    end do
    close(uid)

!    >Load Dipole and Charge ww
!    call LoadBidimentioal_Dipole_Spectrum(OutDir // "/Dipole/DipoleFT_ww.csv", DipoleFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)

    call Load_BidimentionalChargeFTww(OutDir // "/AtomicCharge/AtomicChargeFT_ww.csv", ChargeFTww_new, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)

    !>Reconstruct and Save Dipole from Charge
    call ReconstructDipole_from_AtomicCharge_times_XYZ(ChargeFTww_new, Reconstructed_Dipole, AtCoord, nAtoms, nOmegas, nTauOmegas)

    call Write_2DReconstructDipole(OutDir // "/Dipole/DipoleFT_ww_reconstructed.csv", Reconstructed_Dipole, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
    !
contains


    subroutine ReconstructDipole_from_AtomicCharge_times_XYZ (ChargeFTww_new, Reconstructed_Dipole, R_el_bc, nAtoms, nOmegas, nTauOmegas)
        complex(kind(1d0)), intent(in) :: ChargeFTww_new(:, :, :, :)
        real(kind(1d0)), intent(in) :: R_el_bc(:, :)
        complex(kind(1d0)), allocatable, intent(out) :: Reconstructed_Dipole(:, :, :)
        integer, intent(in) :: nAtoms, nOmegas, nTauOmegas

        integer :: iAtom, iPol, iOmega, iOmegaTau, uid

        allocate(Reconstructed_Dipole(3, nOmegas, nTauOmegas))
        Reconstructed_Dipole = 0d0
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                do iPol = 1, 3
                    do iAtom = 1, nAtoms
                        !                        Reconstructed_Dipole(iPol, iOmega, iOmegaTau) = sum(ChargeFTww_new(iPol, iOmega, iOmegaTau, :))
                        Reconstructed_Dipole(iPol, iOmega, iOmegaTau) = &
                                Reconstructed_Dipole(iPol, iOmega, iOmegaTau) + &
                                        ChargeFTww_new(iPol, iOmega, iOmegaTau, iAtom) !* &
                                        !R_el_bc(iPol, iAtom)
                    end do
                end do
            end do
        end do

        ! Write to screen the sum of theReal part of the X component of the reconstructed dipole
        WRITE(*, *) "Sum of the Real part of the X component of the reconstructed dipole: ", sum(real(Reconstructed_Dipole(1, :, :)))
        WRITE(*, *) "Sum of the Real part of the Y component of the reconstructed dipole: ", sum(real(Reconstructed_Dipole(2, :, :)))
        WRITE(*, *) "Sum of the Real part of the Z component of the reconstructed dipole: ", sum(real(Reconstructed_Dipole(3, :, :)))
    end subroutine ReconstructDipole_from_AtomicCharge_times_XYZ


end program SpectrumReconstruction

