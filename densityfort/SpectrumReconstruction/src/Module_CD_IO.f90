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
            LoadEnergies, &
            SaveDipoleFTFile, &
            SaveAtomicChargeFT, &
            WriteAllDipoleFTtoSingleFile, &
            WriteAllAtomicChargeFTtoSingleFile, &
            LoadDipoles, &
            ReadAtomicCharges, &
            SaveRegularizedDipole, &
            LoadFTDipole_asfuncitonof_TimeDelay, &
            Load_XUVDipole, &
            SaveBidimentioal_Dipole_Spectrum, &
            LoadFTofChargeasFuncofTimeDelay, &
            LoadBidimentioal_Dipole_Spectrum, &
            Load_BidimentionalChargeFTww, &
            Write_BidimentionalChargeFTwwComponent, &
            Write_2DReconstructDipole, &
            Load_BidimentionalChargeFTww1


contains

    subroutine Set_CD_IO_Verbous(logi)
        logical, intent(in) :: logi
        Verbous = logi
    end subroutine Set_CD_IO_Verbous


    !> Load the position of the atomic nuclei
    subroutine LoadGeometry(nAtoms, AtCoord, FileName, AtomName)
        !
        integer, intent(out) :: nAtoms
        real(kind(1d0)), allocatable, intent(out) :: AtCoord(:, :)
        character(len = *), intent(in) :: FileName
        !
        integer :: iAtom, iCoord, uid
        character(len = 16), allocatable, intent(out) :: AtomName(:)
        !
        !.. Open file with geometry
        open(newunit = uid, file = trim(FileName), form = "formatted", status = "old")
        !*** Skip to the line specifying the number of atoms
        ! determine the number of atoms
        read(uid, *) nAtoms
        !*** Skip to the line where the coordinates start to be listed
        !allocate the matrix of coordinates
        read(uid, *)
        allocate(AtCoord(3, nAtoms))
        allocate(AtomName(nAtoms))
        do iAtom = 1, nAtoms
            read(uid, *) AtomName(iAtom), (AtCoord(iCoord, iAtom), iCoord = 1, 3)
        enddo
        close(uid)
    end subroutine LoadGeometry


    !> Loads the Energies found inside the InpDir
    subroutine LoadEnergies(FileName, nStates, Evec)
        !
        use ModuleErrorHandling
        use ModuleString
        !
        implicit none
        !
        character(len = *), intent(in) :: FileName
        integer, intent(out) :: nStates
        real(kind(1d0)), allocatable, intent(out) :: Evec(:)

        integer :: uid, iostat
        character(len = 1000) :: iomsg
        real(kind(1d0)) :: E
        integer :: i

        open(&
                newunit = uid, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read", &
                iostat = iostat, &
                iomsg = iomsg)
        if(iostat /= 0) then
            call errormessage(trim(iomsg))
            stop
        endif

        nStates = 0
        do
            read(uid, *, iostat = iostat) E
            if(iostat/= 0) exit
            nStates = nStates + 1
        enddo
        rewind(uid)
        allocate(Evec(nStates))
        write(*, "(a)") " State Energies "
        do i = 1, nStates
            read(uid, *)Evec(i)
            !write(*,*) i, Evec(i)
        enddo
        write(*, *)
        close(uid)
        !
    end subroutine LoadEnergies


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


    subroutine Load_XUVDipole(FileName, Dipole, nTimes)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nTimes
        complex(kind(1d0)), allocatable, intent(out) :: Dipole(:, :)

        real   (kind(1d0)) :: dBuf
        real   (kind(1d0)), allocatable :: dvec1(:), dvec2(:)
        real(kind(1d0)), parameter :: IMAG_THRESHOLD = 1.d-12
        integer :: uid_dipole, it, iBuf, iPol

        allocate(Dipole(3, ntimes))
        Dipole = Z0
        allocate(dvec1(3), dvec2(3))
        !..Load dipole from file
        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")

        do it = 1, nTimes
            read(uid_dipole, "(i4,*(x,E24.16))") iBuf, dBuf, ((dvec1(iPol), dvec2(iPol)), iPol = 1, 3)
            ! write(*,*) it, iBuf, dBuf
            if(sum(abs(dvec2))>IMAG_THRESHOLD)then
                write(ERROR_UNIT, "(a,d14.6)") "non-zero imaginary dipole", sum(abs(dvec2))
            endif
            do iPol = 1, 3
                Dipole(iPol, it) = Z1 * dvec1(iPol)
            enddo
        enddo
        close(uid_dipole)
    end subroutine Load_XUVDipole


    subroutine LoadDipoles(FileName, nTimes, Dipole)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nTimes
        complex(kind(1d0)), allocatable, intent(out) :: Dipole(:, :)

        real   (kind(1d0)) :: dBuf
        real   (kind(1d0)), allocatable :: dvec1(:), dvec2(:)
        real(kind(1d0)), parameter :: IMAG_THRESHOLD = 1.d-12
        integer :: uid_dipole, it, iBuf, iPol

        allocate(Dipole(3, nTimes))
        Dipole = Z0
        allocate(dvec1(3), dvec2(3))

        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")
        do it = 1, nTimes
            read(uid_dipole, *) iBuf, dBuf, ((dvec1(iPol), dvec2(iPol)), iPol = 1, 3) !"(i4,*(x,E24.16))"
            ! write(*,*) it, iBuf, dBuf
            if(sum(abs(dvec2))>IMAG_THRESHOLD)then
                write(ERROR_UNIT, "(a,d14.6)") "non-zero imaginary dipole", sum(abs(dvec2))
            endif
            do iPol = 1, 3
                Dipole(iPol, it) = Z1 * dvec1(iPol)
            enddo
        enddo
        close(uid_dipole)
    end subroutine LoadDipoles


    subroutine ReadAtomicCharges(FileName, AtomicChargeEvolution, nTimes, nAtoms, iSim, tmin, dt, uid_AtomicChargeALL)
        real(kind(1d0)), intent(in) :: tmin, dt
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nTimes, iSim, nAtoms, uid_AtomicChargeALL
        real(kind(1d0)), allocatable, intent(out) :: AtomicChargeEvolution(:, :)

        real   (kind(1d0)) :: dBuf, t
        integer :: uid_AtomicCharge, it, iAtom

        allocate(AtomicChargeEvolution(nAtoms, nTimes))

        open(newunit = uid_AtomicCharge, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            read(uid_AtomicCharge, "(*(x,e24.16))") dBuf, dBuf, (AtomicChargeEvolution(iAtom, it), iAtom = 1, nAtoms)
            write(uid_AtomicChargeALL, "(*(x,e24.16))") t, dble(iSim), (AtomicChargeEvolution(iAtom, it), iAtom = 1, nAtoms)
        enddo
        write(uid_AtomicChargeALL, *)
        close(uid_AtomicCharge)
    end subroutine ReadAtomicCharges


    subroutine SaveRegularizedDipole(FileName, zMuEV, nTimes, tmin, dt)
        complex(kind(1d0)), intent(in) :: zMuEV(:, :)
        real(kind(1d0)), intent(in) :: tmin, dt
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nTimes

        integer :: uid_dipole, it, iPol
        real(kind(1d0)) :: t

        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            write(uid_dipole, "(i4,*(x,E24.16))") it, t, ((dble(zMuEV(iPol, it)), aimag(zMuEV(iPol, it))), iPol = 1, 3)
        enddo
        close(uid_dipole)
    end subroutine SaveRegularizedDipole


    subroutine SaveDipoleFTFile (FileName, DipoleFTtotal, OmegaVec, nOmegas)
        complex(kind(1d0)), intent(in) :: DipoleFTtotal(:, :)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nOmegas

        real   (kind(1d0)), allocatable :: OmegaVec(:)
        real   (kind(1d0)) :: w
        integer :: uid_dipoleFT, iOmega, iPol

        open(newunit = uid_dipoleFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            write(uid_dipoleFT, "(i4,*(x,E24.16))") iOmega, w, &
                    ((dble(DipoleFTtotal(iPol, iOmega)), aimag(DipoleFTtotal(iPol, iOmega))), iPol = 1, 3)
        end do
        close(uid_dipoleFT)
    end subroutine SaveDipoleFTFile


    subroutine WriteAllDipoleFTtoSingleFile (FileName, DipoleFTtotal, OmegaVec, nOmegas, iSim, train)
        complex(kind(1d0)), intent(in) :: DipoleFTtotal(:, :)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nOmegas, iSim
        type(pulse_train), pointer, intent(in) :: train(:)

        real   (kind(1d0)) :: OmegaVec(:)
        real   (kind(1d0)) :: w
        integer :: uid_dipoleALLFT, iOmega, iPol

        open(newunit = uid_dipoleALLFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write", &
                position = "append")
        do iOmega = 1, nOmegas
            call train(iSim)%PrintPulses(uid_dipoleALLFT)
            w = OmegaVec(iOmega)
            write(uid_dipoleALLFT, "(i4,*(x,E24.16))") iOmega, w, &
                    ((dble(DipoleFTtotal(iPol, iOmega)), aimag(DipoleFTtotal(iPol, iOmega))), iPol = 1, 3)
        end do
        write(uid_dipoleALLFT, *)
        close(uid_dipoleALLFT)
    end subroutine WriteAllDipoleFTtoSingleFile


    subroutine SaveAtomicChargeFT (FileName, AtomicChargeFT, OmegaVec, nOmegas, nAtoms)
        complex(kind(1d0)), intent(in) :: AtomicChargeFT(:, :)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nOmegas, nAtoms

        real   (kind(1d0)) :: OmegaVec(:)
        real   (kind(1d0)) :: w
        integer :: uid_AtomicChargeFT, iOmega, iAtom

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            write(uid_AtomicChargeFT, "(i4,*(x,E24.16))") iOmega, w, ((&
                    dble (AtomicChargeFT(iOmega, iAtom)), &
                    aimag(AtomicChargeFT(iOmega, iAtom))), &
            iAtom = 1, nAtoms)
        end do
        close(uid_AtomicChargeFT)
    end subroutine SaveAtomicChargeFT


    subroutine WriteAllAtomicChargeFTtoSingleFile(FileName, AtomicChargeFT, OmegaVec, nOmegas, nAtoms, iSim, train)
        complex(kind(1d0)) :: AtomicChargeFT(:, :)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nOmegas, nAtoms, iSim
        type(pulse_train), pointer, intent(in) :: train(:)

        real   (kind(1d0)) :: OmegaVec(:)
        real   (kind(1d0)) :: w
        integer :: uid_AtomicChargeFT, iOmega, iAtom

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
                    dble (AtomicChargeFT(iOmega, iAtom)), &
                    aimag(AtomicChargeFT(iOmega, iAtom))   ), &
            iAtom = 1, nAtoms)
        end do
        write(uid_AtomicChargeFT, *)
        close(uid_AtomicChargeFT)
    end subroutine WriteAllAtomicChargeFTtoSingleFile


    !############################################################################################
    !..Bidimentional Spectrum Subroutines
    subroutine LoadFTDipole_asfuncitonof_TimeDelay(FileName, N_simulations, nOmegas, DipoleFTwt, tvec)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: N_simulations, nOmegas
        complex(kind(1d0)), allocatable, intent(out) :: DipoleFTwt(:, :, :)
        real   (kind(1d0)), allocatable, intent(out) :: tvec(:)

        real   (kind(1d0)) :: dBuf, drx, dix, dry, diy, drz, diz
        integer :: iSim, iOmega, uid_dipoleFT, iBuf, i

        allocate(DipoleFTwt(3, nOmegas, N_Simulations))
        allocate(tvec(N_Simulations))

        open(newunit = uid_dipoleFT, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")
        do iSim = 1, N_simulations
            do iOmega = 1, nOmegas
                read(uid_dipoleFT, *)iBuf, (dBuf, i = 1, 7), tvec(iSim), (dBuf, i = 1, 6), iBuf, dBuf, &
                        drx, dix, dry, diy, drz, diz
                DipoleFTwt(1, iOmega, iSim) = Z1 * drx + Zi * dix
                DipoleFTwt(2, iOmega, iSim) = Z1 * dry + Zi * diy
                DipoleFTwt(3, iOmega, iSim) = Z1 * drz + Zi * diz
            end do
            read(uid_dipoleFT, *)
        enddo
        close(uid_dipoleFT)
    end subroutine LoadFTDipole_asfuncitonof_TimeDelay


    subroutine SaveBidimentioal_Dipole_Spectrum(FileName, DipoleFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), intent(in) :: DipoleFTww(:, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas

        integer :: iOmegaTau, iOmega, iPol, uid_dipoleFT

        open(newunit = uid_dipoleFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                write(uid_dipoleFT, "(*(x,e24.16))") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
                        ((dble(DipoleFTww(iPol, iOmega, iOmegaTau)), (aimag(DipoleFTww(iPol, iOmega, iOmegaTau)))), iPol = 1, 3)
            end do
            write(uid_dipoleFT, *)
        enddo
        close(uid_dipoleFT)
    end subroutine SaveBidimentioal_Dipole_Spectrum

    subroutine LoadBidimentioal_Dipole_Spectrum(FileName, DipoleFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), allocatable, intent(out) :: DipoleFTww(:, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas

        real   (kind(1d0)) :: dBuf
        integer :: iOmegaTau, iOmega, iPol, uid_dipoleFT

        allocate(DipoleFTww(3, nOmegas, nTauOmegas))
        open(newunit = uid_dipoleFT, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                !                read(uid_dipoleFT, "(*(x,e24.16))") dBuf, dBuf, &
                !                ((dble(DipoleFTww(iPol, iOmega, iOmegaTau)), (aimag(DipoleFTww(iPol, iOmega, iOmegaTau)))), iPol = 1, 3)
                read(uid_dipoleFT, "(*(x,e24.16))") dBuf, dBuf, &
                        (((DipoleFTww(iPol, iOmega, iOmegaTau))), iPol = 1, 3)

            end do
            read(uid_dipoleFT, *)
        enddo
        close(uid_dipoleFT)
    end subroutine LoadBidimentioal_Dipole_Spectrum


    !.. Load the FT of the charge, as a function of the time delay
    subroutine LoadFTofChargeasFuncofTimeDelay(FileName, N_simulations, nOmegas, nAtoms, tvec, ChargeFTwt)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: N_simulations, nOmegas, nAtoms
        real   (kind(1d0)), allocatable, intent(out) :: tvec(:)
        complex(kind(1d0)), allocatable, intent(out) :: ChargeFTwt(:, :, :)

        real   (kind(1d0)), allocatable :: dAtomFTRe(:), dAtomFTIm(:)
        integer :: uid_AtomicChargeFT, iSim, iOmega, iBuf, iAtom, i
        real(kind(1d0)) :: dBuf

        allocate(tvec(N_simulations))
        tvec = 0.d0
        allocate(dAtomFTRe(nAtoms))
        allocate(dAtomFTIm(nAtoms))
        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")
        allocate(ChargeFTwt(nOmegas, N_Simulations, nAtoms))
        do iSim = 1, N_simulations
            do iOmega = 1, nOmegas
                read(uid_AtomicChargeFT, *)iBuf, (dBuf, i = 1, 7), tvec(iSim), (dBuf, i = 1, 6), iBuf, dBuf, &
                ((dAtomFTRe(iAtom), dAtomFTIm(iAtom)), iAtom = 1, nAtoms)
                do iAtom = 1, nAtoms
                    ChargeFTwt(iOmega, iSim, iAtom) = Z1 * dAtomFTRe(iAtom) + Zi * dAtomFTIm(iAtom)
                enddo
            enddo
        end do
        read(uid_AtomicChargeFT, *)
        deallocate(dAtomFTRe, dAtomFTIm)
        close(uid_AtomicChargeFT)
    end subroutine LoadFTofChargeasFuncofTimeDElay


    !.. Write the Bidimensional spectrum to file
    !    subroutine Write_BidimentionalChargeFTww(FileName, ChargeFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
    !        character(len = *), intent(in) :: FileName
    !        complex(kind(1d0)), intent(in) :: ChargeFTww(:, :, :)
    !        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
    !        integer, intent(in) :: nTauOmegas, nOmegas, nAtoms
    !
    !        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iAtom
    !
    !        open(newunit = uid_AtomicChargeFT, &
    !                file = FileName, &
    !                form = "formatted", &
    !                status = "unknown", &
    !                action = "write")
    !        do iOmegaTau = 1, nTauOmegas
    !            do iOmega = 1, nOmegas
    !                write(uid_AtomicChargeFT, "(*(x,e24.16))") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
    !                        ((dble(ChargeFTww(iOmega, iOmegaTau, iAtom)), &
    !                        (aimag(ChargeFTww(iOmega, iOmegaTau, iAtom)))), iAtom = 1, nAtoms)
    !            end do
    !            write(uid_AtomicChargeFT, *)
    !        enddo
    !        close(uid_AtomicChargeFT)
    !    end subroutine Write_BidimentionalChargeFTww

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
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                write(uid_AtomicChargeFT, "(*(x,e24.16))") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
                        ((dble(Dipole(iCoord, iOmega, iOmegaTau)), &
                        aimag(Dipole(iCoord, iOmega, iOmegaTau))), iCoord = 1, 3)
            end do
            write(uid_AtomicChargeFT, *)
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Write_2DReconstructDipole

    subroutine Write_BidimentionalChargeFTwwComponent(FileName, Dipole, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), intent(in) :: Dipole(:, :, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas, nAtoms

        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iCoord, iAtom
        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                write(uid_AtomicChargeFT, "(*(x,e24.16))") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
                        (((dble(Dipole(iCoord, iOmega, iOmegaTau, iAtom)), &
                        (aimag(Dipole(iCoord, iOmega, iOmegaTau, iAtom)))), iCoord = 1, 3), iAtom =1, nAtoms)
            end do
            write(uid_AtomicChargeFT, *)
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Write_BidimentionalChargeFTwwComponent

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
            write(uid_AtomicChargeFT, *)
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Write_ReconstructedBidimentionalDipole

    subroutine Load_BidimentionalChargeFTww(FileName, ChargeFTww, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), allocatable, intent(out) :: ChargeFTww(:, :, :)
        real   (kind(1d0)), allocatable, intent(out) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas, nAtoms

        real   (kind(1d0)) :: dBuf
        real(kind(1d0)), allocatable :: dvec1(:), dvec2(:)
        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iAtom, iCoord

        allocate(ChargeFTww(nOmegas, nTauOmegas, nAtoms))
        allocate(TauOmegaVec(nTauOmegas))
        allocate(OmegaVec(nOmegas))

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                read(uid_AtomicChargeFT, "(*(x,e24.16))") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
                        (ChargeFTww(iOmega, iOmegaTau, iAtom), iAtom = 1, nAtoms)
            end do
            read(uid_AtomicChargeFT, *)
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Load_BidimentionalChargeFTww
    !
    !
    subroutine Load_BidimentionalChargeFTww1(FileName, ChargeFTww_new, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), allocatable, intent(out) :: ChargeFTww_new(:, :, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas, nAtoms

        real   (kind(1d0)) :: dBuf
        real(kind(1d0)), allocatable :: dvec1(:), dvec2(:)
        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iAtom, iCoord, iPol

        allocate(ChargeFTww_new(3, nOmegas, nTauOmegas, nAtoms))

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                read(uid_AtomicChargeFT, "(*(x,e24.16))") dBuf, dBuf, &
                        ((ChargeFTww_new(iPol, iOmega, iOmegaTau, iAtom), iPol = 1, 3), iAtom = 1, nAtoms)
            end do
            read(uid_AtomicChargeFT, *)
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Load_BidimentionalChargeFTww1


end Module Module_CD_IO
