Module Module_CD_IO

    use, intrinsic :: ISO_FORTRAN_ENV
    use ModuleErrorHandling
    use ModuleString

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
            Write_R_el_bc, &
            Write_Dipole, &
            Write_Q_Charge, &
            Write_Charge_Density, &
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
    subroutine Write_Summary(FileName, nPts, nAtoms, volume, Computed_volume, n_times, t_min, t_max, atom_names, Radius_BS, number_of_orbitals, OrbTab)
        integer, intent(in) :: nPts, nAtoms, n_times, number_of_orbitals
        real(kind(1d0)), intent(in) :: volume, Computed_volume, t_min, t_max, Radius_BS(:), OrbTab(:, :)
        character(len = 16), intent(in) :: atom_names(:)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)) :: self_overlap, other_overlap
        integer :: uid, atom_idx, orbital_idx1, orbital_idx2, point_idx, ierr

        ! Call the new subroutine to write the summary to the screen and the file

        open(newunit = uid, file = FileName, form = "formatted", status = "unknown", action = "write", iostat = ierr)
        if(ierr /= 0) then
            call assert("Write_Summary: Failed to open file for writing: " // FileName, Assertion%LEVEL_SEVERE)
            return
        end if

        write(uid, "(a)")                     "SUMMARY"
        write(uid, "(a,i0)")                  "    Number of Points=", nPts
        write(uid, "(a,i0)")                  "    Number of Atoms=", nAtoms
        write(uid, "(a,e14.6,a,e14.6)")       "    volume =", volume, " Computed volume =", Computed_volume
        write(uid, "(a,i0,e14.6,e14.6)")      "    n_times, t_min, t_max=", n_times, t_min, t_max
        write(uid, "(a)")                     "    Bragg-Slater radii for Becke's Weights"
        do atom_idx = 1, nAtoms
            if (atom_idx .ne. 1) then
                if (atom_names(atom_idx) == atom_names(atom_idx - 1)) cycle
            end if
            write(uid, "(a,a,a,e14.4)")        "        ", atom_names(atom_idx), "=", Radius_BS(atom_idx) !write Bragg-Slater radii
        end do
        write(uid, *)
        do orbital_idx1 = 1, number_of_orbitals
            self_overlap = calculate_overlap(OrbTab, nPts, orbital_idx1, orbital_idx1) * Computed_volume
            other_overlap = 0.d0
            do orbital_idx2 = 1, number_of_orbitals
                if (orbital_idx2 == orbital_idx1) cycle
                other_overlap = other_overlap + calculate_overlap(OrbTab, nPts, orbital_idx1, orbital_idx2)
            end do
            other_overlap = other_overlap * volume
            write(uid, "(a,i0,a,e14.6,a,e14.6)") "    Norm(", orbital_idx1, ") =", self_overlap, "    Overlap other =", other_overlap
        end do
        close(uid)

        ! Write to screen as well - calculte_overlap is only computed once, returning a value after 1st run
        write(*, "(a)")                     "SUMMARY"
        write(*, "(a,i0)")                  "    Number of Points=", nPts
        write(*, "(a,i0)")                  "    Number of Atoms=", nAtoms
        write(*, "(a,e14.6,a,e14.6)")       "    volume =", volume, " Computed volume =", Computed_volume
        write(*, "(a,i0,e14.6,e14.6)")      "    n_times, t_min, t_max=", n_times, t_min, t_max
        write(*, "(a)")                     "    Bragg-Slater radii for Becke's Weights"
        do atom_idx = 1, nAtoms
            if (atom_idx .ne. 1) then
                if (atom_names(atom_idx) == atom_names(atom_idx - 1)) cycle
            end if
            write(*, "(a,a,a,e14.4)")        "        ", atom_names(atom_idx), "=", Radius_BS(atom_idx) !write Bragg-Slater radii
        end do
        write(*, *)
        do orbital_idx1 = 1, number_of_orbitals
            self_overlap = calculate_overlap(OrbTab, nPts, orbital_idx1, orbital_idx1) * Computed_volume
            other_overlap = 0.d0
            do orbital_idx2 = 1, number_of_orbitals
                if (orbital_idx2 == orbital_idx1) cycle
                other_overlap = other_overlap + calculate_overlap(OrbTab, nPts, orbital_idx1, orbital_idx2)
            end do
            other_overlap = other_overlap * volume
            write(*, "(a,i0,a,e14.6,a,e14.6)") "    Norm(", orbital_idx1, ") =", self_overlap, "    Overlap other =", other_overlap
        end do

    contains


        real function calculate_overlap(OrbTab, nPts, orbital_idx1, orbital_idx2)
            real(kind(1d0)), intent(in) :: OrbTab(:, :) ! 3 x npts
            integer, intent(in) :: nPts, orbital_idx1, orbital_idx2
            integer :: point_idx
            calculate_overlap = 0.d0
            do point_idx = 1, nPts
                calculate_overlap = calculate_overlap + OrbTab(point_idx, orbital_idx1) * OrbTab(point_idx, orbital_idx2)
            end do
        end function calculate_overlap

    end subroutine Write_Summary


    !>> Load Subroutines <<!
    !> Load the position of the atomic nuclei
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
    !
    !> Loads the Energies found inside the input_directory
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
    !> Loads the Grid found inside the input_directory
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
        character(len = 100) :: line
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

        ! Skip the header line
        read(uid, '(A)') line

        npts = 0
        do
            read(uid, '(A)', iostat = iostat) line
            if(iostat /= 0) exit
            read(line, *, iostat = iostat) x, y, z
            if(iostat /= 0) exit
            npts = npts + 1
        enddo
        rewind(uid)

        ! Skip the header line again
        read(uid, '(A)') line

        allocate(gridv(3, npts))
        do i = 1, npts
            read(uid, '(A)', iostat = iostat) line
            read(line, *) x, y, z
            gridv(:, i) = [x, y, z]
            ! Display the gridv values
        enddo
        close(uid)
    end subroutine LoadGrid

    !
    !> Loads the Orbitals found inside the input_directory
    !> Loads the Orbitals from the grid_density.csv file inside the input_directory
    subroutine LoadOrbitals(Dir, number_of_orbitals, npts, OrbTab)
        !
        use ModuleErrorHandling
        use ModuleString
        !
        implicit none
        !
        character(len = *), intent(in) :: Dir
        integer, intent(in) :: number_of_orbitals, npts
        real(kind(1d0)), allocatable, intent(out) :: OrbTab(:, :)

        integer :: uid, iostat, i, iOrb
        character(len = 1000) :: iomsg
        character(len = 1000) :: headerLine, dataLine
        !
        write(*, "(a)") "Loading Orbitals"
        !
        allocate(OrbTab(npts, number_of_orbitals))
        ! Make sure that OrbTab is filled
        open(&
                newunit = uid, &
                file = Dir // "/grid_density.csv", &
                form = "formatted", &
                status = "old", &
                action = "read", &
                iostat = iostat, &
                iomsg = iomsg)

        if(iostat /= 0) then
            call Assert(trim(iomsg), Assertion%LEVEL_SEVERE)
            call StopExecution()
        endif

        ! Skip the header line
        read(uid, '(a)') headerLine

        do i = 1, npts
            read(uid, '(a)') dataLine
            read(dataLine, *, iostat = iostat) (OrbTab(i, iOrb), iOrb = 1, number_of_orbitals)
            if (iostat /= 0) then
                call errormessage("Error reading row " // ConvertToStrn(i) // " from grid_density.csv.")
                stop
            endif
        enddo

        close(uid)

        ! Make sure that OrbTab is filled
        do iOrb = 1, number_of_orbitals
            do i = 1, npts
                if (OrbTab(i, iOrb) == 0) then
                    call errormessage("Error in row " // ConvertToStrn(i) // " of OrbTab. Not filled properly.")
                    stop
                endif
            enddo

        end do
        !
    end subroutine LoadOrbitals

    !    subroutine LoadOrbitals(Dir, number_of_orbitals, ivOrb, npts, OrbTab)
    !        !
    !        use ModuleErrorHandling
    !        use ModuleString
    !        !
    !        implicit none
    !        !
    !        character(len = *), intent(in) :: Dir
    !        integer, intent(in) :: number_of_orbitals, npts, ivOrb(:)
    !        real(kind(1d0)), allocatable, intent(out) :: OrbTab(:, :)
    !
    !        integer :: uid, iostat, i, iOrb
    !        character(len = 1000) :: iomsg
    !        character(len = 6) :: istrn
    !        !
    !        write(*, "(a)") "Loading Orbitals"
    !        !
    !        allocate(OrbTab(npts, number_of_orbitals))
    !        do iOrb = 1, number_of_orbitals
    !            write(istrn, "(i0)") ivOrb(iOrb)
    !            istrn = adjustl(istrn)
    !            open(&
    !                    newunit = uid, &
    !                    file = Dir // "/grid" // trim(istrn), &
    !                    form = "formatted", &
    !                    status = "old", &
    !                    action = "read", &
    !                    iostat = iostat, &
    !                    iomsg = iomsg)
    !            if(iostat /= 0) then
    !                call Assert(trim(iomsg), Assertion%LEVEL_SEVERE)
    !                call StopExecution()
    !            endif
    !            do i = 1, npts
    !                read(uid, *) OrbTab(i, iOrb)
    !            enddo
    !            close(uid)
    !            !
    !        enddo
    !    end subroutine LoadOrbitals
    !
    !> Loads all the Dipoles found inside the input_directory
    subroutine LoadDipoleME(Dmat, input_directory, nStates)
        use ModuleErrorHandling
        implicit none
        real(kind(1d0)), allocatable, intent(out) :: Dmat(:, :, :)
        character(len = *), intent(in) :: input_directory
        integer, intent(in) :: nStates

        real(kind(1d0)), allocatable :: dBufM(:, :)
        integer :: i

        write(*, "(a)") "Loading Dipole Matrices for XYZ components"

        if(nStates <= 0)then
            call ErrorMessage("Invalid nStates in LoadDipoleME")
            stop
        endif

        if(allocated(Dmat)) deallocate(Dmat)
        allocate(Dmat(nStates, nStates, 3), dBufM(nStates, nStates))
        call LoadDipoles (input_directory // "/X_DIPOLE.csv", nStates, dBufM)
        Dmat(:, :, 1) = dBufM
        call LoadDipoles (input_directory // "/Y_DIPOLE.csv", nStates, dBufM)
        Dmat(:, :, 2) = dBufM
        call LoadDipoles (input_directory // "/Z_DIPOLE.csv", nStates, dBufM)
        Dmat(:, :, 3) = dBufM
    end subroutine LoadDipoleME


    !> Loads the Dipoles found inside the input_directory
    subroutine LoadDipoles(FileName, nStates, Dmat)
        use ModuleErrorHandling
        use ModuleString
        implicit none

        character(len = *), intent(in) :: FileName
        integer, intent(in) :: nStates
        real(kind(1d0)), intent(inout) :: Dmat(:, :)

        integer :: uid, iostat
        character(len = 1000) :: iomsg
        integer :: i, j
        character(len = 10) :: StateStr
        character(len = 20) :: iStr
        open(&
                newunit = uid, &
                file = FileName, &
                form = "formatted", &
                status = "old", &
                action = "read", &
                iostat = iostat, &
                iomsg = iomsg)
        if(iostat /= 0) then
            call ErrorMessage(trim(iomsg))
            stop
        endif

        write(*, *)
        write(*, "(a)") trim(FileName)

        ! Read and ignore header
        read(uid, *, iostat = iostat, iomsg = iomsg) StateStr
        if (iostat /= 0) then
            call ErrorMessage("Error reading header: " // trim(iomsg))
            close(uid)
            return
        endif

        do i = 1, nStates
            read(uid, *, iostat = iostat, iomsg = iomsg) StateStr, (Dmat(i, j), j = 1, nStates)

            if (iostat /= 0) then

                write(iStr, '(I20)') i
                call ErrorMessage("Error reading row " // trim(iStr) // ": " // trim(iomsg))

                close(uid)
                return
            endif

            write(*, "(*(x,e14.6))") (Dmat(i, j), j = 1, nStates)
        enddo

        close(uid)
    end subroutine LoadDipoles

    !
    !> Loads the TDM found inside the input_directory
    !---------------------------------------------------------------------
    ! LoadTDMs Subroutine
    !---------------------------------------------------------------------
    !
    ! Description:
    ! This subroutine reads in two files, one containing the density matrix (DM)
    ! and the other containing off-diagonal elements of the transition density matrix (TDM).
    ! The DM values are used for the diagonal blocks of the TDM, representing transitions
    ! from a state to itself.
    !
    ! Linear Algebra Context:
    ! The density matrix, rho, for a quantum state |psi> is defined as:
    ! rho = |psi><psi|
    ! The diagonal elements of rho provide the probability of the system being in
    ! each of the basis states.
    !
    ! The transition density matrix, D, for an initial state |psi_i> and a final state |psi_f> is:
    ! D = |psi_f><psi_i|
    ! When the initial and final states are the same (i.e., |psi_i> = |psi_f>),
    ! the TDM reduces to the DM.
    !
    ! In this subroutine:
    ! - Diagonal blocks of TDM (from the file FileNameDM) correspond to the DM
    !   (i.e., transitions from a state to itself).
    ! - Off-diagonal blocks of TDM (from the file FileNameTDM) represent transitions
    !   between different states.
    !
    ! Inputs:
    ! FileNameDM     - File containing the density matrix (DM) values.
    ! FileNameTDM    - File containing the off-diagonal elements of the TDM.
    ! nStates        - Number of quantum states.
    !
    ! Outputs:
    ! number_of_orbitals - Total number of active orbitals read from the DM file.
    ! TDM               - 4D array containing the transition density matrix.
    !
    ! Usage:
    ! call LoadTDMs(FileNameDM, FileNameTDM, nStates, number_of_orbitals, TDM)
    !
    !---------------------------------------------------------------------
    subroutine LoadTDMs(FileNameDM, FileNameTDM, nStates, number_of_orbitals, TDM)
        ! The density matrix (DM) is used as the diagonal blocks of the transition density matrix (TDM). This actually
        !   makes sense when you think about it. The diagonal blocks of the TDM refer to transitions from a state to
        !   itself. For electronic systems, the transition density matrix's diagonal blocks can be equivalent to the
        !   ground state density matrix, especially when considering transitions between ground states or reference
        !   states to excited states.
        ! During External Perturbations: If the system is under the influence of an external field or perturbation,
        !   the diagonal elements of the TDM might contain additional terms or corrections that differ from the DM???
        !
        use ModuleErrorHandling
        use ModuleString
        !
        implicit none
        !
        character(len = *), intent(in) :: FileNameDM
        character(len = *), intent(in) :: FileNameTDM
        integer, intent(in) :: nStates
        integer, intent(out) :: number_of_orbitals
        real(kind(1d0)), allocatable, intent(out) :: TDM(:, :, :, :)

        integer :: uid, iostat
        character(len = 1000) :: iomsg
        character(len = 10000) :: line

        integer :: i_state, ii_state, ij_state, ii_orb, ij_orb, nlines

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
            line = adjustl(trim(line))  ! This removes leading and trailing spaces
            if(len_trim(line)==0)cycle
            nlines = nlines + 1
        enddo
        if(mod(nlines, nStates)/=0)then
            call ErrorMessage("Inconsistent number of lines in " // FileNameDM)
        endif
        number_of_orbitals = nlines / nStates
        write(*, "(a, i0)") "Number of Active Orbitals: ", number_of_orbitals

        !.. Allocate the whole TDM array
        allocate(TDM(number_of_orbitals, number_of_orbitals, nStates, nStates))
        TDM = 0.d0

        !.. Reads the diagonal blocks of TDM
        rewind(uid)
        i_state = 0
        ii_orb = 0
        do
            read(uid, "(a)", iostat = iostat) line
            if(iostat/= 0) exit
            if(len_trim(line)==0)cycle
            if(ii_orb == 0)then
                i_state = i_state + 1
                ii_orb = 1
            else
                ii_orb = ii_orb + 1
            endif
            call replace_char(line, ",", " ")
            read(line, *) (TDM(ii_orb, ij_orb, i_state, i_state), ij_orb = 1, number_of_orbitals)
            if(ii_orb==number_of_orbitals) ii_orb = 0
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
        do ij_state = 2, nStates
            do ii_state = 1, ij_state - 1
                do ii_orb = 1, number_of_orbitals
                    read(uid, "(a)") line
                    call replace_char(line, ",", " ")
                    read(line, *) (TDM(ii_orb, ij_orb, ii_state, ij_state), ij_orb = 1, number_of_orbitals)
                enddo
                TDM(:, :, ij_state, ii_state) = transpose(TDM(:, :, ii_state, ij_state))
            enddo
        enddo
        close(uid)
        !
        if(Verbous)then
            write(*, *) "number_of_orbitals = ", number_of_orbitals
            do ii_state = 1, nStates
                do ij_state = 1, nStates
                    write(*, *) "STATES ", ii_state, " ", ij_state
                    do ii_orb = 1, number_of_orbitals
                        write(*, "(*(x,e14.6))") (TDM(ii_orb, ij_orb, ii_state, ij_state), ij_orb = 1, number_of_orbitals)
                    enddo
                    if(i_state < nStates) write(*, *)
                enddo
            enddo
        endif
    end subroutine LoadTDMs
    !
    !> Loads Dipole Matrix Elements between molecular orbitals
    subroutine LoadDipoleMO(input_directory, number_of_orbitals, ivOrb, MuOrb) !*** IDEALLY, SHOULD COMPUTE THE DIPOLE FROM THE AO - AO DIPOLES.
        !
        use ModuleErrorHandling
        use ModuleString
        !
        implicit none
        !
        character(len = *), intent(in) :: input_directory
        integer, intent(in) :: number_of_orbitals
        integer, intent(in) :: ivOrb(:)
        real(kind(1d0)), allocatable, intent(out) :: MuOrb(:, :)

        character(len = *), parameter :: FILE_AO_DIPOLE_X = "AO_MLTPL_X"
        character(len = *), parameter :: FILE_MO_VECTORS = "MO_VECTORS"

        integer :: uid, iostat
        character(len = 1000) :: iomsg
        character(len = 10000) :: line

        integer :: iAOi, iAOj, iAO, nAO
        integer :: iMO
        integer :: ii_orb, ij_orb
        real(kind(1d0)), allocatable :: MuAO(:, :), dBufv(:)
        real(kind(1d0)), allocatable :: MoVec(:, :)

        !.. Read file dipole between AOs
        open(&
                newunit = uid, &
                file = input_directory // "/" // FILE_AO_DIPOLE_X, &
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
                file = input_directory // "/" // FILE_MO_VECTORS, &
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
        allocate(MuOrb(number_of_orbitals, number_of_orbitals))
        dBufv = 0.d0
        MuOrb = 0.d0
        do ij_orb = 1, number_of_orbitals
            dBufv = matmul(MuAO, MOVec(:, ivOrb(ij_orb)))
            do ii_orb = 1, number_of_orbitals
                MuOrb(ii_orb, ij_orb) = dot_product(MOVec(:, ivOrb(ii_orb)), dBufv)
            enddo
        enddo

        deallocate(MOVec, dBufv, MuAO)

        write(*, *)
        write(*, *) "MO DIPOLE X"
        do ii_orb = 1, number_of_orbitals
            write(*, "(*(x,e14.6))") (MuOrb(ii_orb, ij_orb), ij_orb = 1, number_of_orbitals)
        enddo

    end subroutine LoadDipoleMO


    !>> Save Subroutines
    subroutine Write_R_el_bc(output_directory, atom_names, nAtoms, R_el)
        ! Inputs:
        character(len = *), intent(in) :: output_directory ! Directory path to save the file
        character(len = 16), allocatable, intent(in) :: atom_names(:)
        integer, intent(in) :: nAtoms ! Number of atoms
        real(kind(1d0)), allocatable, intent(in) :: R_el(:, :)

        ! Local variables:
        integer :: uid ! File unit identifier
        integer :: iAtom, iPol ! Loop counters

        ! Attempt to open the file for writing
        open(newunit = uid, file = output_directory // "/" // "R_el_bc.csv", form = "formatted", status = "unknown", action = "write")

        ! Write headers with fixed-width format
        write(uid, "(A14,',',A20,',',A20,',',A20,',',A20)") "Atom Index", "Atom Name", "X Position", "Y Position", "Z Position"

        ! Loop over atoms and write their positions to the file
        do iAtom = 1, nAtoms
            write(uid, "(i5,',',A20,',',3(f20.14,','))") iAtom, trim(atom_names(iAtom)), (R_el(iPol, iAtom), iPol = 1, 3)
        end do


        ! Close the file after writing
        close(uid)

    end subroutine Write_R_el_bc


    subroutine Write_Dipole(FileName, Dipole, n_times, t_min, dt)
        character(len = *), intent(in) :: FileName
        complex(kind(1d0)), intent(in) :: Dipole(:, :)
        real   (kind(1d0)), intent(in) :: t_min, dt
        integer, intent(in) :: n_times

        real   (kind(1d0)) :: t
        integer :: uid_dipole, iPol, it

        open(newunit = uid_dipole, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")

        write(uid_dipole, '(a)') '"itime","Time","DipoleX_Re","DipoleX_Im","DipoleY_Re","DipoleY_Im","DipoleZ_Re","DipoleZ_Im"'

        do it = 1, n_times
            t = t_min + dt * dble(it - 1)
            write(uid_dipole, "(i4,',',E24.16,',')", advance = 'no') it, t
            do iPol = 1, 2
                write(uid_dipole, "(E24.16,',',E24.16,',')", advance = 'no')  dble(Dipole(iPol, it)), aimag(Dipole(iPol, it))
            enddo
            write(uid_dipole, "(E24.16,',',E24.16)") dble(Dipole(3, it)), aimag(Dipole(3, it))
        enddo
        close(uid_dipole)
    end subroutine Write_Dipole
    subroutine Write_Q_Charge(FileName, Charge, n_times, t_min, dt, nAtoms, atom_names)
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), intent(in) :: Charge(:, :, :)
        real   (kind(1d0)), intent(in) :: t_min, dt
        integer, intent(in) :: n_times, nAtoms
        character(len = 16), intent(in) :: atom_names(:)

        real   (kind(1d0)) :: t
        integer :: uid_AtomicCharge, iPol, it, iAtom

        open(newunit = uid_AtomicCharge, &
                file = FileName, &
                form = "formatted", &
                status = "unknown", &
                action = "write")

        write(uid_AtomicCharge, '(a)', advance = "no") '"itime","Time","TotalCharge",'
        do iAtom = 1, nAtoms - 1
            write(uid_AtomicCharge, "(a)", advance = "no") "" &
                    // '"Atom_' // trim(atom_names(iAtom)) // '_ChargeX",' &
                    // '"Atom_' // trim(atom_names(iAtom)) // '_ChargeY",' &
                    // '"Atom_' // trim(atom_names(iAtom)) // '_ChargeZ",'

        end do
        write(uid_AtomicCharge, '(a)') '' &
                // '"Atom_' // trim(atom_names(nAtoms)) // '_ChargeX",' &
                // '"Atom_' // trim(atom_names(nAtoms)) // '_ChargeY",' &
                // '"Atom_' // trim(atom_names(nAtoms)) // '_ChargeZ",'

        do it = 1, n_times
            t = t_min + dt * dble(it - 1)
            write(uid_AtomicCharge, "(i4,',',E24.16,',')", advance = 'no') it, t
            write(uid_AtomicCharge, "(*(x,e24.16,','))", advance = 'no')  sum(Charge(:, :, it))
            do iAtom = 1, nAtoms - 1
                write(uid_AtomicCharge, "(*(x,e24.16,','))", advance = "no") &
                        Charge(1, iAtom, it), &
                        Charge(2, iAtom, it), &
                        Charge(3, iAtom, it)
            end do
            write(uid_AtomicCharge, "(*(x,e24.16,','))") &
                    Charge(1, nAtoms, it), &
                    Charge(2, nAtoms, it), &
                    Charge(3, nAtoms, it)
        enddo

        close(uid_AtomicCharge)
    end subroutine Write_Q_Charge

    subroutine Write_Charge_Density(FileName, npts, gridv, ChDen, Weightv, nAtoms)
        use ModuleErrorHandling
        use ModuleString

        implicit none

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
            write(uid, "(*(x,e24.14e3,','))") (gridv(j, iPts), j = 1, 3), ChDen(iPts)!, (ChDen(iPts) * Weightv(iPts, iAtom), iAtom=1, nAtoms),(Weightv(iPts, iAtom), iAtom=1, nAtoms)
        enddo

        close(uid)
    end subroutine Write_Charge_Density

    !> Write Weights Subroutine
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
        open(newunit = uid, file = FileName, form = "formatted", status = "unknown", action = "write")
        do iPts = 1, nPts
            write(uid, "(*(e24.14e3,:,','))") (gridv(i, iPts), i = 1, 3), (WEIGHTV(iPts, iAtom), iAtom = 1, nAtoms)
        end do
        close(uid)
    end subroutine Write_Weights

    !> Read Weights Subroutine
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
        open(newunit = uid, file = FileName, form = "formatted", status = "old", action = "read")
        do iPts = 1, nPts
            read(uid, "(*(e24.14e3,:,','))") (dBuf(i), i = 1, 3), (WEIGHTV(iPts, iAtom), iAtom = 1, nAtoms)
        end do
        close(uid)
    end subroutine Read_Weights


    !    subroutine Read_Weights(FileName, WEIGHTV, nAtoms, nPts)
    !        character(len = *), intent(in) :: FileName
    !        integer, intent(in) :: nPts, nAtoms
    !        real   (kind(1d0)), allocatable, intent(out) :: WEIGHTV(:, :)
    !        real   (kind(1d0)) :: dBuf(3)
    !        integer :: uid, iPts, iAtom, i
    !        allocate(WEIGHTV(nPts, nAtoms))
    !        !
    !        write(*, "(a)") "Reading Weights from File"
    !        !
    !        !..Read Weights from File
    !        open(newunit = uid, &
    !                file = FileName, &
    !                form = "formatted", &
    !                status = "old", &
    !                action = "read")
    !        do iPts = 1, nPts
    !            read(uid, "(*(x,e24.14e3))") (dBuf(i), i = 1, 3), (WEIGHTV(iPts, iAtom), iAtom = 1, nAtoms)
    !        end do
    !        close(uid)
    !    end subroutine Read_Weights


end Module Module_CD_IO
