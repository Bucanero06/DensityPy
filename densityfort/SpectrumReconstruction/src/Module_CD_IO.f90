Module Module_CD_IO

    use, intrinsic :: ISO_FORTRAN_ENV
    use ModulePulses_3D
    use ModuleConstants

    implicit none

    private

    logical :: Verbous

    public :: &
            Set_CD_IO_Verbous, &
            LoadGeometry, &
            LoadBidimentioal_Dipole_Spectrum, &
            Write_2DReconstructDipole, &
            Load_BidimentionalChargeFTww


contains

    subroutine Set_CD_IO_Verbous(logi)
        logical, intent(in) :: logi
        Verbous = logi
    end subroutine Set_CD_IO_Verbous

    subroutine LoadGeometry(nAtoms, AtCoord, FileName, atom_names)
        !
        integer, intent(out) :: nAtoms
        real(kind(1d0)), allocatable, intent(out) :: AtCoord(:, :)
        character(len = *), intent(in) :: FileName
        !
        integer :: iAtom, iCoord, uid
        character(len = 16), allocatable, intent(out) :: atom_names(:)
        !
        write(*, "(a)") "Loading Geometry"
        !.. Open file with geometry
        open(newunit = uid, file = trim(FileName), form = "formatted", status = "old")
        !*** Skip to the line specifying the number of atoms
        ! determine the number of atoms
        read(uid, *) nAtoms

        !*** Skip to the line where the coordinates start to be listed
        !allocate the matrix of coordinates
        read(uid, *)
        allocate(AtCoord(3, nAtoms))
        allocate(atom_names(nAtoms))
        do iAtom = 1, nAtoms
            read(uid, *) atom_names(iAtom), (AtCoord(iCoord, iAtom), iCoord = 1, 3)
        enddo
        close(uid)
        write(*, *) "nAtoms", nAtoms
    end subroutine LoadGeometry

    subroutine replace_char(strn, ch1, ch2)
        character(len = *), intent(inout) :: strn
        character, intent(in) :: ch1, ch2
        integer :: i
        do
            i = index(strn, ch1)
            if(i<=0)exit
            strn(i:i) = ch2
        enddo
    end subroutine replace_char

    subroutine LoadBidimentioal_Dipole_Spectrum(FileName, DipoleFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), allocatable, intent(out) :: DipoleFTww(:, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas

        real   (kind(1d0)) :: dBuf
        integer :: iOmegaTau, iOmega, iPol, uid_dipoleFT
        character(200) :: line   ! buffer to read the header

        allocate(DipoleFTww(3, nOmegas, nTauOmegas))
        open(newunit = uid_dipoleFT, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")

        read(uid_dipoleFT, '(a)') line

        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                read(uid_dipoleFT, "(*(x,e24.16))") dBuf, dBuf, &
                        (((DipoleFTww(iPol, iOmega, iOmegaTau))), iPol = 1, 3)
                !                        dble(DipoleFTww(1, iOmega, iOmegaTau)), aimag(DipoleFTww(1, iOmega, iOmegaTau)), &
                !                        dble(DipoleFTww(2, iOmega, iOmegaTau)), aimag(DipoleFTww(2, iOmega, iOmegaTau)), &
                !                        dble(DipoleFTww(3, iOmega, iOmegaTau)), aimag(DipoleFTww(3, iOmega, iOmegaTau))
            end do
        enddo
        close(uid_dipoleFT)
    end subroutine LoadBidimentioal_Dipole_Spectrum


    subroutine Write_2DReconstructDipole(FileName, Dipole, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), intent(in) :: Dipole(:, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas

        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iAtom, iCoord

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")

        ! Write Header
        !        write(uid_AtomicChargeFT, *) "OmegaVec,TauOmegaVec", "2DDipoleX_Re", "2DDipoleX_Im", "2DDipoleY_Re", "2DDipoleY_Im", "2DDipoleZ_Re", "2DDipoleZ_Im"
        write(uid_AtomicChargeFT, "(a)") "OmegaVec,TauOmegaVec,2DDipoleX_Re,2DDipoleX_Im,2DDipoleY_Re,2DDipoleY_Im,2DDipoleZ_Re,2DDipoleZ_Im"
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                write(uid_AtomicChargeFT, "(*(x,e24.16,','))", advance = "no") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
                        dble(Dipole(1, iOmega, iOmegaTau)), aimag(Dipole(1, iOmega, iOmegaTau)), &
                        dble(Dipole(2, iOmega, iOmegaTau)), aimag(Dipole(2, iOmega, iOmegaTau))
                write(uid_AtomicChargeFT, "(E24.16,',',E24.16)") &
                        dble(Dipole(3, iOmega, iOmegaTau)), aimag(Dipole(3, iOmega, iOmegaTau))
            end do
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Write_2DReconstructDipole


    subroutine Write_ReconstructedBidimentionalDipole(FileName, Dipole, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), intent(in) :: Dipole(:, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas

        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iCoord

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                write(uid_AtomicChargeFT, "(*(x,e24.16))") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
                        ((dble(Dipole(iCoord, iOmega, iOmegaTau)), &
                        (aimag(Dipole(iCoord, iOmega, iOmegaTau)))), iCoord = 1, 3)
            end do
            !            write(uid_AtomicChargeFT, *)
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Write_ReconstructedBidimentionalDipole

    subroutine Load_BidimentionalChargeFTww(FileName, ChargeFTww_new, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), allocatable, intent(out) :: ChargeFTww_new(:, :, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas, nAtoms

        real   (kind(1d0)) :: dBuf
        real(kind(1d0)), allocatable :: dvec1(:), dvec2(:)
        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iAtom, iCoord, iPol
        character(200) :: line   ! buffer to read the header
        integer :: ioerr

        allocate(ChargeFTww_new(3, nOmegas, nTauOmegas, nAtoms))

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")

        read(uid_AtomicChargeFT, *) line

        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                !                read(uid_AtomicChargeFT, '(a)') line
                read(uid_AtomicChargeFT, "(*(x,e24.16))", IOSTAT = ioerr) dBuf, dBuf, &
                        ((ChargeFTww_new(iPol, iOmega, iOmegaTau, iAtom), iPol = 1, 3), iAtom = 1, nAtoms)
            end do
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Load_BidimentionalChargeFTww


    subroutine Load_BidimentionalChargeFTww2(FileName, ChargeFTww_new, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), allocatable, intent(out) :: ChargeFTww_new(:, :, :, :)
        real(kind(1d0)), allocatable, intent(inout) :: OmegaVec(:), TauOmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas, nAtoms
        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iAtom, iCoord, iError
        real(kind(1d0)) :: rePart, imPart
        character(len = 200) :: headerLine

        ! Open the file for reading
        open(newunit = uid_AtomicChargeFT, file = FileName, form = "formatted", status = "old", action = "read", iostat = iError)
        if (iError /= 0) then
            print *, "Error opening file:", FileName
            return
        endif

        ! Read the header line
        read(uid_AtomicChargeFT, "(a)", iostat = iError) headerLine
        if (iError /= 0) then
            print *, "Error reading header from file:", FileName
            close(uid_AtomicChargeFT)
            return
        endif

        ! Allocate the ChargeFTww_new array
        allocate(ChargeFTww_new(3, nOmegas, nTauOmegas, nAtoms))

        ! Read the data
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                read(uid_AtomicChargeFT, "(*(x,e24.16,','))", iostat = iError) OmegaVec(iOmega), TauOmegaVec(iOmegaTau)
                if (iError /= 0) exit

                do iAtom = 1, nAtoms
                    do iCoord = 1, 3
                        read(uid_AtomicChargeFT, "(*(x,e24.16,','))", iostat = iError) rePart, imPart
                        if (iError /= 0) exit
                        ChargeFTww_new(iCoord, iOmega, iOmegaTau, iAtom) = cmplx(rePart, imPart)
                    end do
                end do
            end do
            if (iError /= 0) then
                print *, "Error reading data from file:", FileName
                exit
            endif
        end do

        close(uid_AtomicChargeFT)
    end subroutine Load_BidimentionalChargeFTww2


end Module Module_CD_IO
