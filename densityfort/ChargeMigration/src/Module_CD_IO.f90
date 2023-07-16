Module Module_CD_IO

    use, intrinsic :: ISO_FORTRAN_ENV
    use ModulePulses_3D
    implicit none

    private

    logical :: Verbous

    public :: &
            Set_CD_IO_Verbous, &
            LoadGeometry, &
            LoadEnergies, &
            LoadGrid, &
            LoadOrbitals, &
            LoadDipoleME, &
            LoadTDMs, &
            LoadDipoleMO, &
            Save_Dipole, &
            Save_Q_Charge, &
            Save_Dipole1, &
            Save_Q_Charge1, &
            SaveChDen, &
            Write_Weights, &
            Read_Weights, &
            Write_Summary

contains

    !>> Commons <<!
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
    subroutine Write_Summary(FileName, nPts, nAtoms, Volume, Computed_Volume, nTimes, tmin, tmax, AtomName, Radius_BS, nOrb, OrbTab)
        integer, intent(in) :: npts, nAtoms, nTimes, nOrb
        real(kind(1d0)), intent(in) :: Volume, Computed_Volume, tmin, tmax, Radius_BS(:), OrbTab(:, :)
        character(len = 16), intent(in) :: AtomName(:)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)) :: Orbital_overlap_self, Orbital_overlap_other
        integer :: uid, iAtom, iOrb, jOrb, iPts

        write(*, "(a)") "Finished Sim Loop"
        write(*, *)
        write(*, *)
        write(*, "(a)")                     "SUMMARY"
        write(*, "(a,i0)")                  "    # of Points=", nPts
        write(*, "(a,i0)")                  "    # of Atoms=", nAtoms
        write(*, "(a,e14.6,a,e14.6)")       "    Volume =", Volume, " Computed Volume =", Computed_Volume
        write(*, "(a,i0,e14.6,e14.6)")      "    nTimes, tmin, tmax=", nTimes, tmin, tmax
        write(*, "(a)")                     "    Bragg-Slater radii fo Becke's Weights"
        do iAtom = 1, nAtoms
            if (iAtom .ne. 1)then
                if (AtomName(iAtom)==AtomName(iAtom - 1)) cycle
            end if
            write(*, "(a,a,a,e14.4)")                 "        ", AtomName(iAtom), "=", Radius_BS(iAtom)!write Bragg-Slater radii
        end do

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
            Orbital_overlap_other = Orbital_overlap_other * Volume

            !Write Results to Screen
            write(*, "(a,i0,a,e14.6,a,e14.6)") "    Norm(", iOrb, ") =", Orbital_overlap_self, "    Overlap other =", Orbital_overlap_other
        end do
        !__________________________________________________________________________________
        !>>Write Summary to File
        !..
        open(newunit = uid, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        write(uid, "(a)")                     "SUMMARY"
        write(uid, "(a,i0)")                  "    # of Points=", nPts
        write(uid, "(a,i0)")                  "    # of Atoms=", nAtoms
        write(uid, "(a,e14.6,a,e14.6)")       "    Volume =", Volume, " Computed Volume =", Computed_Volume
        write(uid, "(a,i0,e14.6,e14.6)")      "    nTimes, tmin, tmax=", nTimes, tmin, tmax
        write(uid, "(a)")                     "    Bragg-Slater radii fo Becke's Weights"
        do iAtom = 1, nAtoms
            write(uid, "(a,a,a,e14.4)")                 "        ", AtomName(iAtom), "=", Radius_BS(iAtom)!write Bragg-Slater radii
        end do

        write(uid, *)
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
            Orbital_overlap_other = Orbital_overlap_other * Volume

            !Write Results to Screen
            write(uid, "(a,i0,a,e14.6,a,e14.6)") "    Norm(", iOrb, ") =", Orbital_overlap_self, "    Overlap other =", Orbital_overlap_other
        end do
        close(uid)
        !_________________________________________
    end subroutine Write_Summary

    !>> Load Subroutines <<!
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
        allocate(AtomName(nAtoms))
        do iAtom = 1, nAtoms
            read(uid, *) AtomName(iAtom), (AtCoord(iCoord, iAtom), iCoord = 1, 3)
        enddo
        close(uid)
        write(*, *) "nAtoms", nAtoms
    end subroutine LoadGeometry
    !
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

        write(*, "(a)") "Loading Energies"
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
            write(*, *) i, Evec(i)
        enddo
        write(*, *)
        close(uid)
        !
    end subroutine LoadEnergies
    !
    !> Loads the Grid found inside the InpDir
    subroutine LoadGrid(FileName, npts, gridv)
        !
        use ModuleErrorHandling
        use ModuleString
        !
        implicit none
        !
        character(len = *), intent(in) :: FileName
        integer, intent(out) :: npts
        real(kind(1d0)), allocatable, intent(out) :: gridv(:, :)

        integer :: uid, iostat
        character(len = 1000) :: iomsg
        real(kind(1d0)) :: x, y, z
        integer :: i, j

        write(*, "(a)") "Loading Grid"
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

        npts = 0
        do
            read(uid, *, iostat = iostat) x, y, z
            if(iostat/= 0) exit
            npts = npts + 1
        enddo
        rewind(uid)
        allocate(gridv(3, npts))
        do i = 1, npts
            read(uid, *) (gridv(j, i), j = 1, 3)
            !!$      write(*,"(i5,*(x,e14.6))") i,(gridv(j,i),j=1,3)
        enddo
        close(uid)
        !
    end subroutine LoadGrid
    !
    !> Loads the Orbitals found inside the InpDir
    subroutine LoadOrbitals(Dir, nOrb, ivOrb, npts, OrbTab)
        !
        use ModuleErrorHandling
        use ModuleString
        !
        implicit none
        !
        character(len = *), intent(in) :: Dir
        integer, intent(in) :: nOrb, npts, ivOrb(:)
        real(kind(1d0)), allocatable, intent(out) :: OrbTab(:, :)

        integer :: uid, iostat, i, iOrb
        character(len = 1000) :: iomsg
        character(len = 6) :: istrn
        !
        write(*, "(a)") "Loading Orbitals"
        !
        allocate(OrbTab(npts, nOrb))
        do iOrb = 1, nOrb
            write(istrn, "(i0)") ivOrb(iOrb)
            istrn = adjustl(istrn)
            open(&
                    newunit = uid, &
                    file = Dir // "/grid" // trim(istrn), &
                    form = "formatted", &
                    status = "old", &
                    action = "read", &
                    iostat = iostat, &
                    iomsg = iomsg)
            if(iostat /= 0) then
                call errormessage(trim(iomsg))
                stop
            endif
            do i = 1, npts
                read(uid, *) OrbTab(i, iOrb)
            enddo
            close(uid)
            !
        enddo
    end subroutine LoadOrbitals
    !
    !> Loads all the Dipoles found inside the InpDir
    subroutine LoadDipoleME(Dmat, InpDir, nStates)
        use ModuleErrorHandling
        implicit none
        real(kind(1d0)), allocatable, intent(out) :: Dmat(:, :, :)
        character(len = *), intent(in) :: InpDir
        integer, intent(in) :: nStates
        !
        real(kind(1d0)), allocatable :: dBufm(:, :)
        !
        write(*, "(a)") "Loading Density Matrix"
        !
        if(nStates <= 0)then
            call ErrorMessage("Invalid nStates in LoadDipoleME")
            stop
        endif
        if(allocated(Dmat)) deallocate(Dmat)
        allocate(Dmat(nStates, nStates, 3), dBufM(nStates, nStates))
        call LoadDipoles (InpDir // "/X_DIPOLE", nStates, dBufM)
        Dmat(:, :, 1) = dBufM
        call LoadDipoles (InpDir // "/Y_DIPOLE", nStates, dBufM)
        Dmat(:, :, 2) = dBufM
        call LoadDipoles (InpDir // "/Z_DIPOLE", nStates, dBufM)
        Dmat(:, :, 3) = dBufM
    end subroutine LoadDipoleME
    !
    !> Loads the Dipoles found inside the InpDir
    subroutine LoadDipoles(FileName, nStates, Dmat)
        !
        use ModuleErrorHandling
        use ModuleString
        !
        implicit none
        !
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nStates
        real(kind(1d0)), intent(inout) :: Dmat(:, :)

        integer :: uid, iostat
        character(len = 1000) :: iomsg
        integer :: i, j

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

        write(*, *)
        write(*, "(a)") trim(FileName)
        read(uid, *)
        do i = 1, nStates
            read(uid, *) (Dmat(i, j), j = 1, nStates)
            write(*, "(*(x,e14.6))") (Dmat(i, j), j = 1, nStates)
        enddo
        close(uid)
        !
    end subroutine LoadDipoles
    !
    !> Loads the TDM found inside the InpDir
    subroutine LoadTDMs(FileNameDM, FileNameTDM, nStates, nOrb, TDM)
        !
        use ModuleErrorHandling
        use ModuleString
        !
        implicit none
        !
        character(len = *), intent(in) :: FileNameDM
        character(len = *), intent(in) :: FileNameTDM
        integer, intent(in) :: nStates
        integer, intent(out) :: nOrb
        real(kind(1d0)), allocatable, intent(out) :: TDM(:, :, :, :)

        integer :: uid, iostat
        character(len = 1000) :: iomsg
        character(len = 10000) :: line

        integer :: iState, iStatei, iStatej, iOrbi, iOrbj, nlines

        write(*, "(a)") "Loading Transition Density Matrix"
        !
        !.. Open the file with the Diagonal part of TDM, i.e., the "DM"
        open(&
                newunit = uid, &
                file = FileNameDM, &
                form = "formatted", &
                status = "old", &
                action = "read", &
                iostat = iostat, &
                iomsg = iomsg)
        if(iostat /= 0) then
            call errormessage(trim(iomsg))
            stop
        endif

        !.. Determines the number of orbitals
        nlines = 0
        do
            read(uid, "(a)", iostat = iostat) line
            if(iostat/= 0) exit
            if(len_trim(line)==0)cycle
            nlines = nlines + 1
        enddo
        if(mod(nlines, nStates)/=0)then
            call ErrorMessage("Inconsistent number of lines in " // FileNameDM)
        endif
        nOrb = nlines / nStates

        !.. Allocate the whole TDM array
        allocate(TDM(nOrb, nOrb, nStates, nStates))
        TDM = 0.d0

        !.. Reads the diagonal blocks of TDM
        rewind(uid)
        iState = 0
        iOrbi = 0
        do
            read(uid, "(a)", iostat = iostat) line
            if(iostat/= 0) exit
            if(len_trim(line)==0)cycle
            if(iOrbi == 0)then
                iState = iState + 1
                iOrbi = 1
            else
                iOrbi = iOrbi + 1
            endif
            call replace_char(line, ",", " ")
            read(line, *) (TDM(iOrbi, iOrbj, iState, iState), iOrbj = 1, nOrb)
            if(iOrbi==nOrb) iOrbi = 0
        enddo
        close(uid)
        !


        !.. Open the file with the Off-Diagonal part of TDM, i.e., the "TDM"
        open(&
                newunit = uid, &
                file = FileNameTDM, &
                form = "formatted", &
                status = "old", &
                action = "read", &
                iostat = iostat, &
                iomsg = iomsg)
        if(iostat /= 0) then
            call errormessage(trim(iomsg))
            stop
        endif

        !.. Reads the diagonal blocks of TDM
        read(uid, *)
        do iStatej = 2, nStates
            do iStatei = 1, iStatej - 1
                do iOrbi = 1, nOrb
                    read(uid, "(a)") line
                    call replace_char(line, ",", " ")
                    read(line, *) (TDM(iOrbi, iOrbj, iStatei, iStatej), iOrbj = 1, nOrb)
                enddo
                TDM(:, :, iStatej, iStatei) = transpose(TDM(:, :, iStatei, iStatej))
            enddo
        enddo
        close(uid)
        !
        if(Verbous)then
            write(*, *) "nOrb = ", nOrb
            do iStatei = 1, nStates
                do iStatej = 1, nStates
                    write(*, *) "STATES ", iStatei, " ", iStatej
                    do iOrbi = 1, nOrb
                        write(*, "(*(x,e14.6))") (TDM(iOrbi, iOrbj, iStatei, iStatej), iOrbj = 1, nOrb)
                    enddo
                    if(iState < nStates) write(*, *)
                enddo
            enddo
        endif

    end subroutine LoadTDMs
    !
    !> Loads Dipole Matrix Elements between molecular orbitals
    subroutine LoadDipoleMO(InpDir, nOrb, ivOrb, MuOrb) !*** IDEALLY, SHOULD COMPUTE THE DIPOLE FROM THE AO - AO DIPOLES.
        !
        use ModuleErrorHandling
        use ModuleString
        !
        implicit none
        !
        character(len = *), intent(in) :: InpDir
        integer, intent(in) :: nOrb
        integer, intent(in) :: ivOrb(:)
        real(kind(1d0)), allocatable, intent(out) :: MuOrb(:, :)

        character(len = *), parameter :: FILE_AO_DIPOLE_X = "AO_MLTPL_X"
        character(len = *), parameter :: FILE_MO_VECTORS = "MO_VECTORS"

        integer :: uid, iostat
        character(len = 1000) :: iomsg
        character(len = 10000) :: line

        integer :: iAOi, iAOj, iAO, nAO
        integer :: iMO
        integer :: iOrbi, iOrbj
        real(kind(1d0)), allocatable :: MuAO(:, :), dBufv(:)
        real(kind(1d0)), allocatable :: MoVec(:, :)

        !.. Read file dipole between AOs
        open(&
                newunit = uid, &
                file = InpDir // "/" // FILE_AO_DIPOLE_X, &
                form = "formatted", &
                status = "old", &
                action = "read", &
                iostat = iostat, &
                iomsg = iomsg)
        if(iostat /= 0) then
            call errormessage(trim(iomsg))
            stop
        endif

        !.. Determines the number of AOs
        nAO = 0
        do
            read(uid, "(a)", iostat = iostat)line
            if(iostat/= 0) exit
            nAO = nAO + 1
        enddo
        write(*, *) "nAO : ", nAO

        !.. Reads the Dipole between AOs
        allocate(MuAO(nAO, nAO))
        MuAO = 0.d0
        rewind(uid)
        do iAOi = 1, nAO
            read(uid, "(a)", iostat = iostat)line
            call replace_char(line, ",", " ")
            read(line, *) (MuAO(iAOi, iAOj), iAOj = 1, nAO)
        enddo
        close(uid)

        !.. Read file MO coefficients
        open(&
                newunit = uid, &
                file = InpDir // "/" // FILE_MO_VECTORS, &
                form = "formatted", &
                status = "old", &
                action = "read", &
                iostat = iostat, &
                iomsg = iomsg)
        if(iostat /= 0) then
            call errormessage(trim(iomsg))
            stop
        endif

        read(uid, *)
        allocate(MOVec(nAO, nAO))
        MOVec = 0.d0
        do iMO = 1, nAO
            do iAO = 1, nAO
                read(uid, *) MOVec(iAO, iMO)
            enddo
        enddo
        close(uid)
        write(*, *) "MOVec check:", sum(abs(MOVec))
        write(*, "(a,*(x,i4))") "ivOrb:", ivOrb

        allocate(dBufv(nAO))
        allocate(MuOrb(nOrb, nOrb))
        dBufv = 0.d0
        MuOrb = 0.d0
        do iOrbj = 1, nOrb
            dBufv = matmul(MuAO, MOVec(:, ivOrb(iOrbj)))
            do iOrbi = 1, nOrb
                MuOrb(iOrbi, iOrbj) = dot_product(MOVec(:, ivOrb(iOrbi)), dBufv)
            enddo
        enddo

        deallocate(MOVec, dBufv, MuAO)

        write(*, *)
        write(*, *) "MO DIPOLE X"
        do iOrbi = 1, nOrb
            write(*, "(*(x,e14.6))") (MuOrb(iOrbi, iOrbj), iOrbj = 1, nOrb)
        enddo

    end subroutine LoadDipoleMO

    !>> Save Subroutines
    subroutine Save_Dipole(FileName, Dipole, nTimes, tmin, dt)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), intent(in) :: Dipole(:, :)
        real   (kind(1d0)), intent(in) :: tmin, dt
        integer, intent(in) :: nTimes

        real   (kind(1d0)) :: t
        integer :: uid_dipole, iPol, it

        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            write(uid_dipole, "(i4,*(x,E24.16))") it, t, ((dble(Dipole(iPol, it)), aimag(Dipole(iPol, it))), iPol = 1, 3)
        enddo
        close(uid_dipole)
    end subroutine Save_Dipole
    subroutine Save_Q_Charge(FileName, Charge, nTimes, tmin, dt, nAtoms)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(in) :: Charge(:, :)
        real   (kind(1d0)), intent(in) :: tmin, dt
        integer, intent(in) :: nTimes, nAtoms

        real   (kind(1d0)) :: t
        integer :: uid_AtomicCharge, iPol, it, iAtom

        !.. Save Q_Charge
        open(newunit = uid_AtomicCharge, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        !..Regular
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            write(uid_AtomicCharge, "(*(x,e24.16))") t, sum(Charge(:, it)), (Charge(iAtom, it), iAtom = 1, nAtoms)
        enddo
        !
        !
        !..New
        !        do it = 1, nTimes
        !            t = tmin + dt * dble(it - 1)
        !            write(uid_AtomicCharge, "(*(x,e24.16))") t, sum(Charge(:, :, it)), ((Charge(iPol, iAtom, it), iPol = 1, 3), iAtom = 1, nAtoms)
        !        enddo
        !        !
        close(uid_AtomicCharge)
    end subroutine Save_Q_Charge
    subroutine Save_Dipole1(FileName, Dipole, nTimes, tmin, dt)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(in) :: Dipole(:, :)
        !        complex(kind(1d0)), intent(in) :: Dipole(:, :)
        real   (kind(1d0)), intent(in) :: tmin, dt
        integer, intent(in) :: nTimes

        real   (kind(1d0)) :: t
        integer :: uid_dipole, iPol, it

        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            !            write(uid_dipole, "(i4,*(x,E24.16))") it, t, ((dble(Dipole(iPol, it)), aimag(Dipole(iPol, it))), iPol = 1, 3)
            write(uid_dipole, "(i4,*(x,E24.16))") it, t, ((Dipole(iPol, it), 0.d0), iPol = 1, 3)
        enddo
        close(uid_dipole)
    end subroutine Save_Dipole1
    subroutine Save_Q_Charge1(FileName, Charge, nTimes, tmin, dt, nAtoms)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(in) :: Charge(:, :, :)
        real   (kind(1d0)), intent(in) :: tmin, dt
        integer, intent(in) :: nTimes, nAtoms

        real   (kind(1d0)) :: t
        integer :: uid_AtomicCharge, iPol, it, iAtom

        !.. Save Q_Charge
        open(newunit = uid_AtomicCharge, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")

        !..New
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            write(uid_AtomicCharge, "(*(x,e24.16))") t, sum(Charge(:, :, it)), ((Charge(iPol, iAtom, it), iPol = 1, 3), iAtom = 1, nAtoms)
        enddo
        !
        close(uid_AtomicCharge)
    end subroutine Save_Q_Charge1
    subroutine SaveChDen(FileName, npts, gridv, ChDen, Weightv, nAtoms)
        !
        use ModuleErrorHandling
        use ModuleString
        !
        implicit none
        !
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nPts, nAtoms
        real(kind(1d0)), allocatable, intent(in) :: gridv(:, :)
        real(kind(1d0)), allocatable, intent(in) :: ChDen(:)
        real(kind(1d0)), allocatable, intent(in) :: Weightv(:, :)
        integer :: uid, iostat
        character(len = 1000) :: iomsg
        integer :: iPts, iAtom, j

        open(&
                newunit = uid, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write", &
                iostat = iostat, &
                iomsg = iomsg)
        if(iostat /= 0) then
            call errormessage(trim(iomsg))
            stop
        endif

        do iPts = 1, nPts
            !            if (gridv(1, iPts) .ne. gridv(1, iPts-1)) then
            !                    write(uid, *)
            !            end if

            write(uid, "(*(x,e24.14e3))") (gridv(j, iPts), j = 1, 3), ChDen(iPts)!, (ChDen(iPts) * Weightv(iPts, iAtom), iAtom=1, nAtoms),(Weightv(iPts, iAtom), iAtom=1, nAtoms)

        enddo

        close(uid)
        !
    end subroutine SaveChDen

    !..Write and Read Weights Subroutines
    subroutine Write_Weights(FileName, WEIGHTV, gridv, nAtoms, nPts)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(in) :: gridv(:, :)
        integer, intent(in) :: nPts, nAtoms
        real   (kind(1d0)), intent(in) :: WEIGHTV(:, :)
        integer :: uid, iPts, iAtom, i
        !
        write(*, "(a)") "Writing Weights to File"
        !
        !..Write Weights to File
        open(newunit = uid, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do iPts = 1, nPts
            write(uid, "(*(x,e24.14e3))") (gridv(i, iPts), i = 1, 3), (WEIGHTV(iPts, iAtom), iAtom = 1, nAtoms)
        end do
        close(uid)
    end subroutine Write_Weights
    subroutine Read_Weights(FileName, WEIGHTV, nAtoms, nPts)
        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nPts, nAtoms
        real   (kind(1d0)), allocatable, intent(out) :: WEIGHTV(:, :)
        real   (kind(1d0)) :: dBuf(3)
        integer :: uid, iPts, iAtom, i
        allocate(WEIGHTV(nPts, nAtoms))
        !
        write(*, "(a)") "Reading Weights from File"
        !
        !..Read Weights from File
        open(newunit = uid, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read")
        do iPts = 1, nPts
            read(uid, "(*(x,e24.14e3))") (dBuf(i), i = 1, 3), (WEIGHTV(iPts, iAtom), iAtom = 1, nAtoms)
        end do
        close(uid)
    end subroutine Read_Weights


end Module Module_CD_IO
