Module Module_CD_IO

    use, intrinsic :: ISO_FORTRAN_ENV
    use ModulePulses_3D
    use ModuleConstants
    use ModuleErrorHandling

    implicit none

    private

    logical :: Verbous

    public :: &
            Set_CD_IO_Verbous, &
            LoadGeometry, &
            LoadEnergies, &
            WriteDipoleFTFile, &
            Load_Q_Charge_and_Write2ALL, &
            WriteAtomicChargeFT, &
            AppendDipole2FTAllFile, &
            WriteAllAtomicChargeFTtoSingleFile, &
            Load_Dipole, &
            LoadAtomicCharges, &
            LoadFTDipole_asfuncitonof_TimeDelay, &
            SaveBidimentioal_Dipole_Spectrum, &
            LoadFTofChargeasFuncofTimeDelay, &
            Write_BidimentionalChargeFTww, &
            Load_XUVAtomicCharge


contains

    subroutine Set_CD_IO_Verbous(logi)
        logical, intent(in) :: logi
        Verbous = logi
    end subroutine Set_CD_IO_Verbous
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

    function int2strn(int_num) result(strn)
        integer, intent(in) :: int_num
        character(len = 16) :: strn
        write(strn, "(i16)") int_num
        strn = adjustl(strn)
        strn = trim(strn)
    end function int2strn


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


    !
    subroutine Load_XUVAtomicCharge(FileName, Charge, nTimes, nAtoms)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nTimes, nAtoms
        complex(kind(1d0)), allocatable, intent(out) :: Charge(:, :, :)

        real   (kind(1d0)) :: dBuf
        real(kind(1d0)), parameter :: IMAG_THRESHOLD = 1.d-12
        integer :: uid, it, iBuf, iAtom, iPol
        character(200) :: line   ! buffer to read the header

        allocate(Charge(3, nAtoms, ntimes))
        !..Load dipole from file
        open(newunit = uid, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")


        ! Read and ignore the header
        read(uid, '(a)') line
        do it = 1, nTimes
            read(uid, *) iBuf, dBuf, dBuf, ((Charge(iPol, iAtom, it), iPol = 1, 3), iAtom = 1, nAtoms)
        enddo
        close(uid)
    end subroutine Load_XUVAtomicCharge

    subroutine Load_Dipole(FileName, Dipole, nTimes)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nTimes
        complex(kind(1d0)), allocatable, intent(out) :: Dipole(:, :)

        real(kind(1d0)) :: dBuf
        real(kind(1d0)), allocatable :: dvecReal(:), dvecImag(:)
        real(kind(1d0)), parameter :: IMAG_THRESHOLD = 1.d-12
        integer :: uid_dipole, it, iBuf, iPol
        character(200) :: line   ! buffer to read the header

        allocate(Dipole(3, nTimes))
        Dipole = Z0
        allocate(dvecReal(3), dvecImag(3))

        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")

        ! Read and ignore the header
        read(uid_dipole, '(a)') line

        do it = 1, nTimes
            read(uid_dipole, *) iBuf, dBuf, dvecReal(1), dvecImag(1), dvecReal(2), dvecImag(2), dvecReal(3), dvecImag(3)

            ! Check for significant imaginary components
            if(sum(abs(dvecImag)) > IMAG_THRESHOLD) then
                write(ERROR_UNIT, "(a,d14.6)") "non-zero imaginary dipole", sum(abs(dvecImag))
            endif

            do iPol = 1, 3
                !   Dipole(iPol, it) = cmplx(dvecReal(iPol), 0.0, kind = 1d0)  ! As the imaginary part is assumed to be negligible
                Dipole(iPol, it) = Z1 * dvecReal(iPol)
            enddo
        enddo

        close(uid_dipole)
    end subroutine Load_Dipole
    subroutine Load_Q_Charge_and_Write2ALL(FileName, Charge, nTimes, tmin, dt, nAtoms, iSim, atom_names, uid_AtomicChargeALL)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), allocatable, intent(out) :: Charge(:, :, :)
        real   (kind(1d0)), intent(in) :: tmin, dt
        character(len = 16), intent(in) :: atom_names(:)
        integer, intent(in) :: nTimes, nAtoms, iSim, uid_AtomicChargeALL

        character(200) :: line   ! buffer to read the header

        real   (kind(1d0)) :: t, dBuf
        integer :: uid_AtomicCharge, iPol, it, iAtom, iBuf

        if (.not.allocated(Charge))allocate(Charge(3, nAtoms, nTimes))

        !.. Save Q_Charge
        open(newunit = uid_AtomicCharge, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")

        if (iSim == 1) then
            write(uid_AtomicChargeALL, '(a)', advance = "no") '"itime","iSim",'

            do iAtom = 1, nAtoms - 1
                write(uid_AtomicChargeALL, "(a)", advance = "no") "" &
                        // '"Atom_' // trim(atom_names(iAtom)) // '_ChargeX",' &
                        // '"Atom_' // trim(atom_names(iAtom)) // '_ChargeY",' &
                        // '"Atom_' // trim(atom_names(iAtom)) // '_ChargeZ",'
            end do
            write(uid_AtomicChargeALL, '(a)') '' &
                    // '"Atom_' // trim(atom_names(nAtoms)) // '_ChargeX",' &
                    // '"Atom_' // trim(atom_names(nAtoms)) // '_ChargeY",' &
                    // '"Atom_' // trim(atom_names(nAtoms)) // '_ChargeZ"'
        end if

        ! skip the header
        read(uid_AtomicCharge, "(a)") line
        !..New
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            read(uid_AtomicCharge, *) iBuf, dBuf, dBuf, ((Charge(iPol, iAtom, it), iPol = 1, 3), iAtom = 1, nAtoms)
!            write(uid_AtomicChargeALL, "(*(x,e24.16,','))") t, dble(iSim), ((Charge(iPol, iAtom, it), iPol = 1, 3), iAtom = 1, nAtoms)
            write(uid_AtomicChargeALL, "(*(x,e24.16,','))", advance = "no") t, dble(iSim)
            do iAtom = 1, nAtoms - 1
                write(uid_AtomicChargeALL, "(*(x,e24.16,','))", advance = "no") &
                        Charge(1, iAtom, it), &
                        Charge(2, iAtom, it), &
                        Charge(3, iAtom, it)
            end do
            write(uid_AtomicChargeALL, "(*(x,e24.16,','))", advance = "no") &
                    Charge(1, nAtoms, it), &
                    Charge(2, nAtoms, it)
            write(uid_AtomicChargeALL, "(E24.16)") &
                    Charge(3, nAtoms, it)
        enddo
        !
        write(uid_AtomicChargeALL, *)
        close(uid_AtomicCharge)
    end subroutine Load_Q_Charge_and_Write2ALL
    subroutine LoadAtomicCharges(FileName, AtomicChargeEvolution, nTimes, nAtoms, iSim, tmin, dt, uid_AtomicChargeALL)
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
    end subroutine LoadAtomicCharges

    subroutine WriteDipoleFTFile (FileName, DipoleFTtotal, OmegaVec, nOmegas)
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

        ! Header
        write(uid_dipoleFT, "(a)", advance = "no")  '"iOmega","OmegaVec",'
        write(uid_dipoleFT, "(a)") '' &
                // '"FTDipoleX_Re","FTDipoleX_Im",' &
                // '"FTDipoleY_Re","FTDipoleY_Im",' &
                // '"FTDipoleZ_Re","FTDipoleZ_Im"'

        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            write(uid_dipoleFT, "(i4,',',*(x,E24.16,','))", advance = "no")  iOmega, w, &
                    dble(DipoleFTtotal(1, iOmega)), aimag(DipoleFTtotal(1, iOmega)), &
                    dble(DipoleFTtotal(2, iOmega)), aimag(DipoleFTtotal(2, iOmega))
            write(uid_dipoleFT, "(E24.16,',',E24.16)") &
                    dble(DipoleFTtotal(3, iOmega)), aimag(DipoleFTtotal(3, iOmega))

        end do
        close(uid_dipoleFT)
    end subroutine WriteDipoleFTFile

    subroutine Write_Pulse_Columns(train, uid)
        class(pulse_train), intent(in) :: train
        integer :: uid
        integer :: iPulse

        !# todo bug review, in CM-FT CD_IO Module this returns an array not a scalar e.g. train%n = 2 2 2 2 2 1
        write(uid, "(x,i5,',')", advance = "no") train%n
        do iPulse = 1, train%n
            write(uid, "(*(x,e14.6,','))", advance = "no") &
                    train%p(iPulse)%t, &
                    train%p(iPulse)%o, &
                    train%p(iPulse)%f, &
                    train%p(iPulse)%d, &
                    train%p(iPulse)%i, &
                    train%p(iPulse)%a, &
                    train%p(iPulse)%p
        enddo
    end subroutine Write_Pulse_Columns

    subroutine AppendDipole2FTAllFile (FileName, DipoleFTtotal, OmegaVec, nOmegas, iSim, train)
        complex(kind(1d0)), intent(in) :: DipoleFTtotal(:, :)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nOmegas, iSim
        type(pulse_train), pointer, intent(in) :: train(:)

        real   (kind(1d0)) :: OmegaVec(:)
        real   (kind(1d0)) :: w
        integer :: uid_dipoleALLFT, iOmega, iPol, iPulse
        character(len = 16) :: iPulseStr
        character(len = 32) :: header
        character(len = 32) :: trailing_text

        open(newunit = uid_dipoleALLFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write", &
                position = "append")

        ! Header ... only write the first time, meaning if it is the first line
        if (iSim == 1) then
            !    Attribute Symbols  |  Description
            !    ------------------------------------------
            !    n                 |  Not directly found in 'pulse', but in 'pulse_train'. Represents the number of pulses.
            !    t                 |  Central Time (in atomic units) - The time at which the pulse is centered or reaches its peak.
            !    o                 |  Carrier Frequency (in atomic units) - The frequency of the carrier wave of the pulse.
            !    f                 |  Full Width at Half Maximum (FWHM) - The width of the pulse at half its maximum amplitude.
            !    d                 |  Carrier Envelope Phase (in degrees) - Phase difference between the pulse's carrier frequency and its envelope.
            !    i                 |  Intensity (in PW/cm^2) - Power of the pulse per unit area.
            !    a                 |  Amplitude (in atomic units) - Maximum amplitude or height of the pulse.
            !    p                 |  Period (in atomic units) - Time for one complete cycle of the wave.
            write(uid_dipoleALLFT, "(a)", advance = "no")  '"number_of_pulses",'

            do iPulse = 1, train(iSim)%n
                iPulseStr = int2strn(iPulse)
                write(uid_dipoleALLFT, "(a)", advance = "no") '"central_time_' // trim(iPulseStr) // '",'
                write(uid_dipoleALLFT, "(a)", advance = "no") '"carrier_frequency_' // trim(iPulseStr) // '",'
                write(uid_dipoleALLFT, "(a)", advance = "no") '"fwhm_' // trim(iPulseStr) // '",'
                write(uid_dipoleALLFT, "(a)", advance = "no") '"carrier_envelope_phase_' // trim(iPulseStr) // '",'
                write(uid_dipoleALLFT, "(a)", advance = "no") '"intensity_' // trim(iPulseStr) // '",'
                write(uid_dipoleALLFT, "(a)", advance = "no") '"amplitude_' // trim(iPulseStr) // '",'
                write(uid_dipoleALLFT, "(a)", advance = "no") '"period_' // trim(iPulseStr) // '",'

            end do
            write(uid_dipoleALLFT, "(a)", advance = "no") '"iOmega","OmegaVec",'
            write(uid_dipoleALLFT, '(a)') '' &
                    // '"FTDipoleX_Re","FTDipoleX_Im",' &
                    // '"FTDipoleY_Re","FTDipoleY_Im",' &
                    // '"FTDipoleZ_Re","FTDipoleZ_Im"'
        end if

        do iOmega = 1, nOmegas
            call Write_Pulse_Columns(train(iSim), uid_dipoleALLFT)
            w = OmegaVec(iOmega)
            !
            write(uid_dipoleALLFT, "(i4,',',*(x,E24.16,','))", advance = "no")  iOmega, w, &
                    dble(DipoleFTtotal(1, iOmega)), aimag(DipoleFTtotal(1, iOmega)), &
                    dble(DipoleFTtotal(2, iOmega)), aimag(DipoleFTtotal(2, iOmega))
            write(uid_dipoleALLFT, "(E24.16,',',E24.16)") &
                    dble(DipoleFTtotal(3, iOmega)), aimag(DipoleFTtotal(3, iOmega))
            !
        end do
        close(uid_dipoleALLFT)
    end subroutine AppendDipole2FTAllFile

    subroutine WriteAtomicChargeFT (FileName, AtomicChargeFT, OmegaVec, nOmegas, nAtoms, AtomName)
        complex(kind(1d0)), intent(in) :: AtomicChargeFT(:, :)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nOmegas, nAtoms
        character(len = 16), intent(in) :: AtomName(:)

        real   (kind(1d0)) :: OmegaVec(:)
        real   (kind(1d0)) :: w
        integer :: uid_AtomicChargeFT, iOmega, iAtom

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")

        write(uid_AtomicChargeFT, "(a)", advance = "no") '"iOmega","OmegaVec",'
        do iAtom = 1, nAtoms - 1
            write(uid_AtomicChargeFT, "(a)", advance = "no") '' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeX_Re",' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeX_Im",' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeY_Re",' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeY_Im",' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeZ_Re",' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeZ_Im",'

        end do
        write(uid_AtomicChargeFT, '(a)') '' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeX_Re",' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeX_Im",' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeY_Re",' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeY_Im",' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeZ_Re",' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeZ_Im"'

        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
!            write(uid_AtomicChargeFT, "(i4,*(x,E24.16,','))") iOmega, w, ((&
!                    dble (AtomicChargeFT(iOmega, iAtom)), &
!                    aimag(AtomicChargeFT(iOmega, iAtom))), &
!            iAtom = 1, nAtoms)

            write(uid_AtomicChargeFT, "(i4,',',*(x,E24.16,','))", advance = "no") iOmega, w
            do iAtom = 1, nAtoms - 1
                write(uid_AtomicChargeFT, "(*(x,e24.16,','))", advance = "no") &
                        dble(AtomicChargeFT(iOmega, iAtom)), &
                        aimag(AtomicChargeFT(iOmega, iAtom)), &
                        dble(AtomicChargeFT(iOmega, iAtom)), &
                        aimag(AtomicChargeFT(iOmega, iAtom)), &
                        dble(AtomicChargeFT(iOmega, iAtom)), &
                        aimag(AtomicChargeFT(iOmega, iAtom))
            end do
            write(uid_AtomicChargeFT, "(*(x,e24.16,','))", advance = "no") &
                    dble(AtomicChargeFT(iOmega, nAtoms)), &
                    aimag(AtomicChargeFT(iOmega, nAtoms)), &
                    dble(AtomicChargeFT(iOmega, nAtoms)), &
                    aimag(AtomicChargeFT(iOmega, nAtoms))
            write(uid_AtomicChargeFT, "(E24.16,',',E24.16)") &
                    dble(AtomicChargeFT(iOmega, nAtoms)), &
                    aimag(AtomicChargeFT(iOmega, nAtoms))
            !
        end do
        close(uid_AtomicChargeFT)
    end subroutine WriteAtomicChargeFT

    subroutine WriteAllAtomicChargeFTtoSingleFile(FileName, AtomicChargeFT_new, OmegaVec, nOmegas, nAtoms, iSim, train, AtomName)
        use ModuleErrorHandling
        complex(kind(1d0)) :: AtomicChargeFT_new(:, :, :)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nOmegas, nAtoms, iSim
        type(pulse_train), pointer, intent(in) :: train(:)
        character(len = 16), intent(in) :: AtomName(:)
        real   (kind(1d0)) :: OmegaVec(:)
        real   (kind(1d0)) :: w
        integer :: uid_AtomicChargeFT, iOmega, iAtom, iPol, iPulse
        character(len = 16) :: iPulseStr

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write", &
                position = "append")
        if (iSim == 1) then
            write(uid_AtomicChargeFT, "(a)", advance = "no")  '"number_of_pulses",'
            do iPulse = 1, train(iSim)%n
                iPulseStr = int2strn(iPulse)
                write(uid_AtomicChargeFT, "(a)", advance = "no") '"central_time_' // trim(iPulseStr) // '",'
                write(uid_AtomicChargeFT, "(a)", advance = "no") '"carrier_frequency_' // trim(iPulseStr) // '",'
                write(uid_AtomicChargeFT, "(a)", advance = "no") '"fwhm_' // trim(iPulseStr) // '",'
                write(uid_AtomicChargeFT, "(a)", advance = "no") '"carrier_envelope_phase_' // trim(iPulseStr) // '",'
                write(uid_AtomicChargeFT, "(a)", advance = "no") '"intensity_' // trim(iPulseStr) // '",'
                write(uid_AtomicChargeFT, "(a)", advance = "no") '"amplitude_' // trim(iPulseStr) // '",'
                write(uid_AtomicChargeFT, "(a)", advance = "no") '"period_' // trim(iPulseStr) // '",'
            end do
            write(uid_AtomicChargeFT, "(a)", advance = "no") '"iOmega","OmegaVec",'
            do iAtom = 1, nAtoms - 1
                write(uid_AtomicChargeFT, "(a)", advance = "no") '' &
                        // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeX_Re",' &
                        // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeX_Im",' &
                        // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeY_Re",' &
                        // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeY_Im",' &
                        // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeZ_Re",' &
                        // '"Atom_' // trim(AtomName(iAtom)) // '_FTChargeZ_Im",'

            end do
            write(uid_AtomicChargeFT, '(a)') '' &
                    // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeX_Re",' &
                    // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeX_Im",' &
                    // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeY_Re",' &
                    // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeY_Im",' &
                    // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeZ_Re",' &
                    // '"Atom_' // trim(AtomName(nAtoms)) // '_FTChargeZ_Im"'
        end if

        do iOmega = 1, nOmegas
            call Write_Pulse_Columns(train(iSim), uid_AtomicChargeFT)
            w = OmegaVec(iOmega)
            write(uid_AtomicChargeFT, "(i4,',',*(x,E24.16,','))", advance = "no") iOmega, w
            do iAtom = 1, nAtoms - 1
                write(uid_AtomicChargeFT, "(*(x,e24.16,','))", advance = "no") &
                        dble(AtomicChargeFT_new(1, iOmega, iAtom)), &
                        aimag(AtomicChargeFT_new(1, iOmega, iAtom)), &
                        dble(AtomicChargeFT_new(2, iOmega, iAtom)), &
                        aimag(AtomicChargeFT_new(2, iOmega, iAtom)), &
                        dble(AtomicChargeFT_new(3, iOmega, iAtom)), &
                        aimag(AtomicChargeFT_new(3, iOmega, iAtom))
            end do
            write(uid_AtomicChargeFT, "(*(x,e24.16,','))", advance = "no") &
                    dble(AtomicChargeFT_new(1, iOmega, nAtoms)), &
                    aimag(AtomicChargeFT_new(1, iOmega, nAtoms)), &
                    dble(AtomicChargeFT_new(2, iOmega, nAtoms)), &
                    aimag(AtomicChargeFT_new(2, iOmega, nAtoms))
            write(uid_AtomicChargeFT, "(E24.16,',',E24.16)") &
                    dble(AtomicChargeFT_new(3, iOmega, nAtoms)), &
                    aimag(AtomicChargeFT_new(3, iOmega, nAtoms))
        end do
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

        ! Skip the header
        read(uid_dipoleFT, *)

        do iSim = 1, N_simulations
            do iOmega = 1, nOmegas
                read(uid_dipoleFT, *) iBuf, (dBuf, i = 1, 7), tvec(iSim), (dBuf, i = 1, 6), iBuf, dBuf, &
                        drx, dix, dry, diy, drz, diz
                DipoleFTwt(1, iOmega, iSim) = Z1 * drx + Zi * dix
                DipoleFTwt(2, iOmega, iSim) = Z1 * dry + Zi * diy
                DipoleFTwt(3, iOmega, iSim) = Z1 * drz + Zi * diz
            end do
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

        ! Headers
        write(uid_dipoleFT, "(a)") "OmegaVec,TauOmegaVec,2DDipoleX_Re,2DDipoleX_Im,2DDipoleY_Re,2DDipoleY_Im,2DDipoleZ_Re,2DDipoleZ_Im"
        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                write(uid_dipoleFT, "(*(x,E24.16,','))", advance = "no") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
                        dble(DipoleFTww(1, iOmega, iOmegaTau)), aimag(DipoleFTww(1, iOmega, iOmegaTau)), &
                        dble(DipoleFTww(2, iOmega, iOmegaTau)), aimag(DipoleFTww(2, iOmega, iOmegaTau))
                write(uid_dipoleFT, "(E24.16,',',E24.16)") &
                        dble(DipoleFTww(3, iOmega, iOmegaTau)), aimag(DipoleFTww(3, iOmega, iOmegaTau))
            end do
        enddo
        close(uid_dipoleFT)
    end subroutine SaveBidimentioal_Dipole_Spectrum

    !.. Load the FT of the charge, as a function of the time delay
    subroutine LoadFTofChargeasFuncofTimeDelay(FileName, N_simulations, nOmegas, nAtoms, tvec, ChargeFTwt_new)
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


        ! Skip the header
        read(uid_AtomicChargeFT, *)

        do iSim = 1, N_simulations
            do iOmega = 1, nOmegas
                read(uid_AtomicChargeFT, *) iBuf, (dBuf, i = 1, 7), tvec(iSim), (dBuf, i = 1, 6), iBuf, dBuf, &
                (((dAtomFTRe(iPol, iAtom), dAtomFTIm(iPol, iAtom)), iPol = 1, 3), iAtom = 1, nAtoms)
                do iPol = 1, 3
                    do iAtom = 1, nAtoms
                        ChargeFTwt_new(iPol, iOmega, iSim, iAtom) = Z1 * dAtomFTRe(iPol, iAtom) + Zi * dAtomFTIm(iPol, iAtom)
                    enddo
                end do
            enddo
        end do
        deallocate(dAtomFTRe, dAtomFTIm)
        close(uid_AtomicChargeFT)
    end subroutine LoadFTofChargeasFuncofTimeDElay

    !.. Write the Bidimensional spectrum to file
    subroutine Write_BidimentionalChargeFTww(FileName, ChargeFTww_new, TauOmegaVec, OmegaVec, nTauOmegas, nOmegas, nAtoms, AtomName)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), intent(in) :: ChargeFTww_new(:, :, :, :)
        real   (kind(1d0)), intent(in) :: TauOmegaVec(:), OmegaVec(:)
        integer, intent(in) :: nTauOmegas, nOmegas, nAtoms
        character(len = 16), intent(in) :: AtomName(:)

        integer :: uid_AtomicChargeFT, iOmegaTau, iOmega, iAtom, iPol

        open(newunit = uid_AtomicChargeFT, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")

        write(uid_AtomicChargeFT, "(a)", advance = "no") "OmegaVec,TauOmegaVec,"
        do iAtom = 1, nAtoms - 1
            write(uid_AtomicChargeFT, "(a)", advance = "no") '' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_2DChargeX_Re",' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_2DChargeX_Im",' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_2DChargeY_Re",' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_2DChargeY_Im",' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_2DChargeZ_Re",' &
                    // '"Atom_' // trim(AtomName(iAtom)) // '_2DChargeZ_Im",'

        end do
        write(uid_AtomicChargeFT, "(a)") "" &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_2DChargeX_Re",' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_2DChargeX_Im",' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_2DChargeY_Re",' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_2DChargeY_Im",' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_2DChargeZ_Re",' &
                // '"Atom_' // trim(AtomName(nAtoms)) // '_2DChargeZ_Im"'

        do iOmegaTau = 1, nTauOmegas
            do iOmega = 1, nOmegas
                write(uid_AtomicChargeFT, "(*(x,e24.16,','))", advance = "no") OmegaVec(iOmega), TauOmegaVec(iOmegaTau)
                do iAtom = 1, nAtoms - 1
                    write(uid_AtomicChargeFT, "(*(x,e24.16,','))", advance = "no") &
                            dble(ChargeFTww_new(1, iOmega, iOmegaTau, iAtom)), &
                            aimag(ChargeFTww_new(1, iOmega, iOmegaTau, iAtom)), &
                            dble(ChargeFTww_new(2, iOmega, iOmegaTau, iAtom)), &
                            aimag(ChargeFTww_new(2, iOmega, iOmegaTau, iAtom)), &
                            dble(ChargeFTww_new(3, iOmega, iOmegaTau, iAtom)), &
                            aimag(ChargeFTww_new(3, iOmega, iOmegaTau, iAtom))

                end do
                write(uid_AtomicChargeFT, "(*(x,e24.16,','))", advance = "no") &
                        dble(ChargeFTww_new(1, iOmega, iOmegaTau, nAtoms)), &
                        aimag(ChargeFTww_new(1, iOmega, iOmegaTau, nAtoms)), &
                        dble(ChargeFTww_new(2, iOmega, iOmegaTau, nAtoms)), &
                        aimag(ChargeFTww_new(2, iOmega, iOmegaTau, nAtoms))
                write(uid_AtomicChargeFT, "(e24.16,',',e24.16)") &
                        dble(ChargeFTww_new(3, iOmega, iOmegaTau, nAtoms)), &
                        aimag(ChargeFTww_new(3, iOmega, iOmegaTau, nAtoms))

            end do
        enddo
        close(uid_AtomicChargeFT)
    end subroutine Write_BidimentionalChargeFTww

end Module Module_CD_IO
