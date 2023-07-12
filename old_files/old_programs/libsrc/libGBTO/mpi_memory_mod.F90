! Copyright 2019
!
! Zdenek Masin with contributions from others (see the UK-AMOR website)
!
! This file is part of GBTOlib.
!
!     GBTOlib is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     GBTOlib is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with  GBTOlib (in trunk/COPYING). Alternatively, you can also visit
!     <https://www.gnu.org/licenses/>.
!

!> \brief   Local and shared memory allocation and deallocation routines
!> \authors Z Masin, J Benda
!> \date    2017 - 2019
!>
!> This module provides convenience wrappers around memory allocation calls. All interfaces support three scenarios:
!>  1. The parameter is an allocatable array. The (de)allocation will then use the standard Fortran statement (DE)ALLOCATE.
!>     The (de)allocated array is local to every process that calls the interface.
!>  2. The parameter is a pointer to array and `shared_enabled` is `.false.`. This will have the same effect as above.
!>  3. The parameter is a pointer to array and `shared_enabled` is `.true.`. This will (de)allocate a single chunk of memory
!>     shared by all processes that call the interface.
!>
!> \warning This module is also used by SCATCI but always in DOUBLE PRECISION. However, SCATCI must always be
!>          linked to the integral library of precision used to calculate the integrals file.
!>          Therefore, in order to ensure correct functionality of SCATCI we need to ensure this module supports
!>          calls to both DOUBLE and QUAD precision MPI routines within the same program run.
!>          This is achieved by interfacing all routines manipulating floating point variables for both wp and ep precisions.
!>
!> \note 05/02/2019 - Jakub Benda: Split subroutines into generic pairs (ALLOCATABLE/POINTER). This is to allow allocation
!>       of both types with this module and hence support both local (ALLOCATABLE) and shared (POINTER) arrays.
!> \note 24/02/2019 - Jakub Benda: Allow operation on user-provided MPI communicators.
!>
module mpi_memory_mod

    use precisn,       only: wp, ep, wp_bytes, ep1_bytes, longint
    use mpi_mod,       only: mpiint, mpiaddr, local_rank, shared_enabled, shared_communicator, mpi_xermsg, mpi_mod_rank
    use const,         only: stdout
    use iso_c_binding, only: c_ptr, c_f_pointer, c_null_ptr

#ifdef usempi
    ! This list should also include MPI_Win_allocate_shared and MPI_Win_shared_query, but these are omitted due to a bug
    ! in Intel MPI (2017.2, 2019.3), which does not include them in its MPI module.
    use mpi,           only: MPI_INFO_NULL, MPI_SUCCESS, MPI_Win_fence, MPI_Win_free, MPI_Barrier
#endif

    implicit none

    integer(kind=mpiint), parameter :: local_master = 0

    interface mpi_memory_allocate_integer
        module procedure mpi_memory_allocate_integer_alc          ! accepts allocatable array
        module procedure mpi_memory_allocate_integer_ptr          ! accepts pointer to array
    end interface

    interface mpi_memory_deallocate_integer
        module procedure mpi_memory_deallocate_integer_alc        ! accepts allocatable array
        module procedure mpi_memory_deallocate_integer_ptr        ! accepts pointer to array
    end interface

    interface mpi_memory_allocate_integer_2dim
        module procedure mpi_memory_allocate_integer_2dim_alc     ! accepts allocatable array
        module procedure mpi_memory_allocate_integer_2dim_ptr     ! accepts pointer to array
    end interface

    interface mpi_memory_deallocate_integer_2dim
        module procedure mpi_memory_deallocate_integer_2dim_alc   ! accepts allocatable array
        module procedure mpi_memory_deallocate_integer_2dim_ptr   ! accepts pointer to array
    end interface

    interface mpi_memory_allocate_real
        module procedure mpi_memory_allocate_real_wp_alc          ! accepts allocatable array
        module procedure mpi_memory_allocate_real_ep_alc          ! accepts allocatable array
        module procedure mpi_memory_allocate_real_wp_ptr          ! accepts pointer to array
        module procedure mpi_memory_allocate_real_ep_ptr          ! accepts pointer to array
    end interface

    interface mpi_memory_deallocate_real
        module procedure mpi_memory_deallocate_real_wp_alc        ! accepts allocatable array
        module procedure mpi_memory_deallocate_real_ep_alc        ! accepts allocatable array
        module procedure mpi_memory_deallocate_real_wp_ptr        ! accepts pointer to array
        module procedure mpi_memory_deallocate_real_ep_ptr        ! accepts pointer to array
    end interface

    interface mpi_memory_allocate_real_2dim
        module procedure mpi_memory_allocate_real_2dim_wp_alc     ! accepts allocatable array
        module procedure mpi_memory_allocate_real_2dim_ep_alc     ! accepts allocatable array
        module procedure mpi_memory_allocate_real_2dim_wp_ptr     ! accepts pointer to array
        module procedure mpi_memory_allocate_real_2dim_ep_ptr     ! accepts pointer to array
    end interface

    interface mpi_memory_deallocate_real_2dim
        module procedure mpi_memory_deallocate_real_2dim_wp_alc   ! accepts allocatable array
        module procedure mpi_memory_deallocate_real_2dim_ep_alc   ! accepts allocatable array
        module procedure mpi_memory_deallocate_real_2dim_wp_ptr   ! accepts pointer to array
        module procedure mpi_memory_deallocate_real_2dim_ep_ptr   ! accepts pointer to array
    end interface

contains

    !> \brief  Dummy initialization routine
    !> \author Z Masin
    !> \date   2017
    !>
    !> Has no effect at the moment.
    !>
    subroutine mpi_memory_setup
    end subroutine mpi_memory_setup


    !> \brief  Allocate 1D allocatable integer array
    !> \author Z Masin
    !> \date   2017
    !>
    integer function mpi_memory_allocate_integer_alc (array, nelem, comm)

        integer,      allocatable, intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem
        integer(mpiint), optional, intent(in)    :: comm

        integer :: error, n_bytes

        allocate(array(nelem), stat = error)

        if (error /= 0) then
            call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_integer', 'Memory allocation failed.', error, 1)
        end if

        mpi_memory_allocate_integer_alc = -1

    end function mpi_memory_allocate_integer_alc


    !> \brief  Allocate 1D allocatable default integer array
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    integer function mpi_memory_allocate_integer_2dim_alc (array, nelem1, nelem2, comm)

        integer,      allocatable, intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem1, nelem2
        integer(mpiint), optional, intent(in)    :: comm

        integer :: error, n_bytes

        allocate(array(nelem1, nelem2), stat = error)

        if (error /= 0) then
            call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_integer_2dim', 'Memory allocation failed.', error, 1)
        end if

        mpi_memory_allocate_integer_2dim_alc = -1

    end function mpi_memory_allocate_integer_2dim_alc


    !> \brief  Allocate 1D allocatable real(wp) array
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    integer function mpi_memory_allocate_real_wp_alc (array, nelem, comm)

        real(wp),     allocatable, intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem
        integer(mpiint), optional, intent(in)    :: comm

        integer :: error

        allocate(array(nelem), stat = error)

        if (error /= 0) then
            call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_real', 'Memory allocation failed.', error, 1)
        end if

        mpi_memory_allocate_real_wp_alc = -1

    end function mpi_memory_allocate_real_wp_alc


    !> \brief  Allocate 1D allocatable real(ep) array
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    integer function mpi_memory_allocate_real_ep_alc (array, nelem, comm)

        real(ep),     allocatable, intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem
        integer(mpiint), optional, intent(in)    :: comm

        integer :: error

        allocate(array(nelem), stat = error)

        if (error /= 0) then
            call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_real', 'Memory allocation failed.', error, 1)
        end if

        mpi_memory_allocate_real_ep_alc = -1

    end function mpi_memory_allocate_real_ep_alc


    !> \brief  Allocate 2D allocatable real(wp) array
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    integer function mpi_memory_allocate_real_2dim_wp_alc (array, nelem1, nelem2, comm)

        real(wp),     allocatable, intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem1, nelem2
        integer(mpiint), optional, intent(in)    :: comm

        integer :: error

        allocate(array(nelem1, nelem2), stat = error)

        if (error /= 0) then
            call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_real_2dim_wp', 'Memory allocation failed.', error, 1)
        end if

        mpi_memory_allocate_real_2dim_wp_alc = -1

    end function mpi_memory_allocate_real_2dim_wp_alc


    !> \brief  Allocate 2D allocatable real(ep) array
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    integer function mpi_memory_allocate_real_2dim_ep_alc (array, nelem1, nelem2, comm)

        real(ep),      allocatable, intent(inout) :: array(:,:)
        integer,                    intent(in)    :: nelem1, nelem2
        integer(mpiint),  optional, intent(in)    :: comm

        integer :: error

        allocate(array(nelem1, nelem2), stat = error)

        if (error /= 0) then
            call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_real_2dim_ep', 'Memory allocation failed.', error, 1)
        end if

        mpi_memory_allocate_real_2dim_ep_alc = -1

    end function mpi_memory_allocate_real_2dim_ep_alc


    !> \brief  Deallocate 1D allocatable real(wp) array
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_real_wp_alc (array, nelem, window, comm)

        real(wp),     allocatable, intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        deallocate(array)

    end subroutine mpi_memory_deallocate_real_wp_alc


    !> \brief  Deallocate 1D allocatable real(wp) array
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_real_ep_alc (array, nelem, window, comm)

        real(ep),     allocatable, intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        deallocate(array)

    end subroutine mpi_memory_deallocate_real_ep_alc


    !> \brief  Deallocate 1D allocatable real(wp) array
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_integer_2dim_alc (array, nelem, window, comm)

        integer,      allocatable, intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        deallocate(array)

    end subroutine mpi_memory_deallocate_integer_2dim_alc


    !> \brief  Deallocate 1D allocatable default integer array
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_integer_alc (array, nelem, window, comm)

        integer,      allocatable, intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        deallocate(array)

    end subroutine mpi_memory_deallocate_integer_alc


    !> \brief  Deallocate 2D allocatable real(wp) array
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_real_2dim_wp_alc (array, nelem, window, comm)

        real(wp),     allocatable, intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        deallocate(array)

    end subroutine mpi_memory_deallocate_real_2dim_wp_alc


    !> \brief  Deallocate 2D allocatable real(ep) array
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_real_2dim_ep_alc (array, nelem, window, comm)

        real(ep),     allocatable, intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        deallocate(array)

    end subroutine mpi_memory_deallocate_real_2dim_ep_alc

    !> \brief   Shared part of the shared memory allocation routines
    !> \authors Z Masin, J Benda
    !> \date    2017 - 2019
    !>
    !> Allocates the given amount of bytes in the MPI shared memory and retrieves the (local master's) pointer to the
    !> beginning of the allocated chunk for all processes (rather than pointers to sections associated with them).
    !>
    subroutine mpi_memory_allocate_shared_bytes (alloc_size, groupcomm, baseptr, win)

        integer(mpiaddr), intent(in)  :: alloc_size
        integer(mpiint),  intent(in)  :: groupcomm
        integer(mpiint),  intent(out) :: win
        type(c_ptr),      intent(out) :: baseptr

#if defined(usempi) && defined(mpithree)
        integer(mpiaddr) :: asize
        integer(mpiint)  :: disp_unit, ierror

        win = -1
        disp_unit = 1
        ierror = MPI_SUCCESS

        call MPI_Win_allocate_shared(alloc_size, disp_unit, MPI_INFO_NULL, groupcomm, baseptr, win, ierror)

        if (ierror /= MPI_SUCCESS) then
            call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_shared_bytes', &
                            'Error when calling MPI_WIN_ALLOCATE_SHARED', int(ierror), 1)
        end if

        call MPI_Win_shared_query(win, local_master, asize, disp_unit, baseptr, ierror)

        if (ierror /= MPI_SUCCESS) then
            call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_shared_bytes', &
                            'Error when calling MPI_WIN_SHARED_QUERY', int(ierror), 1)
        end if
#else
        baseptr = c_null_ptr
        win = -1
#endif

    end subroutine mpi_memory_allocate_shared_bytes


    !> \brief   Deallocate 1D shared MPI default integer array pointer
    !> \authors Z Masin, J Benda
    !> \date    2017 - 2019
    !>
    !> If the sub-group communicator 'comm' is not given, then the memory will be allocated on 'shared_communicator',
    !> which is set up by mpi_mod.
    !>
    integer function mpi_memory_allocate_integer_ptr (array, nelem, comm)

        integer,         pointer,  intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiaddr) :: alloc_size
        integer(mpiint)  :: disp_unit, info, win, ierror, groupcomm, grouprank
        integer          :: error, n_bytes

        type(c_ptr) :: baseptr

        if (shared_enabled) then
            if (present(comm)) then; groupcomm = comm; else; groupcomm = shared_communicator; end if
            call mpi_mod_rank(grouprank, groupcomm)

            n_bytes    = bit_size(array) / 8
            alloc_size = int(merge(n_bytes * nelem, 0, grouprank == local_master), mpiaddr)

            call mpi_memory_allocate_shared_bytes(alloc_size, groupcomm, baseptr, win)
            call c_f_pointer(baseptr, array, [nelem])

            mpi_memory_allocate_integer_ptr = win
        else
            allocate (array(nelem), stat = error)

            if (error /= 0) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_integer', 'Memory allocation failed.', error, 1)
            end if

            mpi_memory_allocate_integer_ptr = -1
        end if

    end function mpi_memory_allocate_integer_ptr


    !> \brief  Allocate 2D shared MPI default integer array pointer
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    !> If the sub-group communicator 'comm' is not given, then the memory will be allocated on 'shared_communicator',
    !> which is set up by mpi_mod.
    !>
    integer function mpi_memory_allocate_integer_2dim_ptr (array, nelem1, nelem2, comm)

        integer,         pointer,  intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem1, nelem2
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiaddr) :: alloc_size
        integer(mpiint)  :: length, disp_unit, info, win, ierror, groupcomm, grouprank
        integer          :: error, n_bytes

        type(c_ptr) :: baseptr

        if (shared_enabled) then
            if (present(comm)) then; groupcomm = comm; else; groupcomm = shared_communicator; end if
            call mpi_mod_rank(grouprank, groupcomm)

            n_bytes    = bit_size(array) / 8
            alloc_size = int(merge(n_bytes * nelem1 * nelem2, 0, grouprank == local_master), mpiaddr)

            call mpi_memory_allocate_shared_bytes(alloc_size, groupcomm, baseptr, win)
            call c_f_pointer(baseptr, array, [nelem1, nelem2])

            mpi_memory_allocate_integer_2dim_ptr = win
        else
            allocate (array(nelem1, nelem2), stat = error)

            if (error /= 0) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_integer_2dim', 'Memory allocation failed.', error, 1)
            end if

            mpi_memory_allocate_integer_2dim_ptr = -1
        end if

    end function mpi_memory_allocate_integer_2dim_ptr


    !> \brief  Allocate 1D shared MPI real(wp) array pointer
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    !> If the sub-group communicator 'comm' is not given, then the memory will be allocated on 'shared_communicator',
    !> which is set up by mpi_mod.
    !>
    integer function mpi_memory_allocate_real_wp_ptr (array, nelem, comm)

        real(wp),        pointer,  intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiaddr) :: alloc_size
        integer(mpiint)  :: disp_unit, info, win, ierror, groupcomm, grouprank
        integer          :: error

        type(c_ptr) :: baseptr

        if (shared_enabled) then
            if (present(comm)) then; groupcomm = comm; else; groupcomm = shared_communicator; end if
            call mpi_mod_rank(grouprank, groupcomm)

            alloc_size = int(merge(wp_bytes * nelem, 0, grouprank == local_master), mpiaddr)

            call mpi_memory_allocate_shared_bytes(alloc_size, groupcomm, baseptr, win)
            call c_f_pointer(baseptr, array, [nelem])

            mpi_memory_allocate_real_wp_ptr = win
        else
            allocate (array(nelem), stat = error)

            if (error /= 0) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_real', 'Memory allocation failed.', error, 1)
            end if

            mpi_memory_allocate_real_wp_ptr = -1
        end if

    end function mpi_memory_allocate_real_wp_ptr


    !> \brief  Allocate 1D shared MPI real(ep) array pointer
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    !> If the sub-group communicator 'comm' is not given, then the memory will be allocated on 'shared_communicator',
    !> which is set up by mpi_mod.
    !>
    integer function mpi_memory_allocate_real_ep_ptr (array, nelem, comm)

        real(ep),        pointer,  intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiaddr) :: alloc_size
        integer(mpiint)  :: info, win, ierror, groupcomm, grouprank
        integer          :: error

        type(c_ptr) :: baseptr

        if (shared_enabled) then
            if (present(comm)) then; groupcomm = comm; else; groupcomm = shared_communicator; end if
            call mpi_mod_rank(grouprank, groupcomm)

            alloc_size = int(merge(ep1_bytes * nelem, 0, grouprank == local_master), mpiaddr)

            call mpi_memory_allocate_shared_bytes(alloc_size, groupcomm, baseptr, win)
            call c_f_pointer(baseptr, array, [nelem])

            mpi_memory_allocate_real_ep_ptr = win
        else
            allocate (array(nelem), stat = error)

            if (error /= 0) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_real', 'Memory allocation failed.', error, 1)
            end if

            mpi_memory_allocate_real_ep_ptr = -1
        end if

    end function mpi_memory_allocate_real_ep_ptr


    !> \brief  Allocate 2D shared MPI real(wp) array pointer
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    !> If the sub-group communicator 'comm' is not given, then the memory will be allocated on 'shared_communicator',
    !> which is set up by mpi_mod.
    !>
    integer function mpi_memory_allocate_real_2dim_wp_ptr (array, nelem1, nelem2, comm)

        real(wp),        pointer,  intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem1, nelem2
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiaddr) :: alloc_size
        integer(mpiint)  :: win, ierror, grouprank, groupcomm
        integer          :: error

        type(c_ptr) :: baseptr

        write (stdout, '("Allocating memory of size ",2I18)') nelem1, nelem2

        if (shared_enabled) then
            if (present(comm)) then; groupcomm = comm; else; groupcomm = shared_communicator; end if
            call mpi_mod_rank(grouprank, groupcomm)

            alloc_size = int(merge(wp_bytes * nelem1 * nelem2, 0, grouprank == local_master), mpiaddr)

            call mpi_memory_allocate_shared_bytes(alloc_size, groupcomm, baseptr, win)
            call c_f_pointer(baseptr, array, [nelem1, nelem2])

            mpi_memory_allocate_real_2dim_wp_ptr = win
        else
            allocate (array(nelem1, nelem2), stat = error)

            if (error /= 0) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_real_2dim_wp', 'Memory allocation failed.', error, 1)
            end if

            mpi_memory_allocate_real_2dim_wp_ptr = -1
        end if

    end function mpi_memory_allocate_real_2dim_wp_ptr


    !> \brief  Allocate 2D shared MPI real(ep) array pointer
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    !> If the sub-group communicator 'comm' is not given, then the memory will be allocated on 'shared_communicator',
    !> which is set up by mpi_mod.
    !>
    integer function mpi_memory_allocate_real_2dim_ep_ptr (array, nelem1, nelem2, comm)

        real(ep),        pointer,  intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem1, nelem2
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiaddr) :: alloc_size
        integer(mpiint)  :: win, ierror, groupcomm, grouprank
        integer          :: error

        type(c_ptr) :: baseptr

        write (stdout, '("Allocating memory of size ",2I18)') nelem1, nelem2

        if (shared_enabled) then
            if (present(comm)) then; groupcomm = comm; else; groupcomm = shared_communicator; end if
            call mpi_mod_rank(grouprank, groupcomm)

            alloc_size = int(merge(ep1_bytes * nelem1 * nelem2, 0, local_rank == local_master), mpiaddr)

            call mpi_memory_allocate_shared_bytes(alloc_size, groupcomm, baseptr, win)
            call c_f_pointer(baseptr, array, [nelem1, nelem2])

            mpi_memory_allocate_real_2dim_ep_ptr = win
        else
            allocate (array(nelem1, nelem2), stat = error)

            if (error /= 0) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_allocate_real_2dim_ep', 'Memory allocation failed.', error, 1)
            end if

            mpi_memory_allocate_real_2dim_ep_ptr = -1
        end if

    end function mpi_memory_allocate_real_2dim_ep_ptr


    !> \brief  Deallocate 1D shared MPI real(wp) array pointer
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_real_wp_ptr (array, nelem, window, comm)

        real(wp),        pointer,  intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiint) :: ierror

        if (shared_enabled) then
            if (window == -1) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_deallocate_real_wp_ptr', &
                                'Illegal shared memory window declared', 1, 1)
            else
                call mpi_memory_synchronize(window, comm)
                call mpi_memory_win_free(window, ierror)
            end if
        else
            deallocate (array)
        end if

    end subroutine mpi_memory_deallocate_real_wp_ptr


    !> \brief  Deallocate 1D shared MPI real(ep) array pointer
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_real_ep_ptr (array, nelem, window, comm)

        real(ep),        pointer,  intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiint) :: ierror

        if (shared_enabled) then
            if (window == -1) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_deallocate_real_ep_ptr', &
                                'Illegal shared memory window declared', 1, 1)
            else
                call mpi_memory_synchronize(window)
                call mpi_memory_win_free(window, ierror)
            end if
        else
            deallocate (array)
        end if

    end subroutine mpi_memory_deallocate_real_ep_ptr


    !> \brief  Deallocate 2D shared MPI integer array pointer
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_integer_2dim_ptr (array, nelem, window, comm)

        integer,         pointer,  intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiint) :: ierror

        if (shared_enabled) then
            if (window == -1) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_deallocate_integer_2dim_ptr', &
                                'Illegal shared memory window declared', 1, 1)
            else
                call mpi_memory_synchronize(window, comm)
                call mpi_memory_win_free(window, ierror)
            end if
        else
            deallocate (array)
        end if

    end subroutine mpi_memory_deallocate_integer_2dim_ptr


    !> \brief  Deallocate 1D shared MPI integer array pointer
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_integer_ptr (array, nelem, window, comm)

        integer,         pointer,  intent(inout) :: array(:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiint) :: ierror

        if (shared_enabled) then
            if (window == -1) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_deallocate_integer_ptr', &
                                'Illegal shared memory window declared', 1, 1)
            else
                call mpi_memory_synchronize(window, comm)
                call mpi_memory_win_free(window, ierror)
            end if
        else
            deallocate (array)
        end if

    end subroutine mpi_memory_deallocate_integer_ptr


    !> \brief  Deallocate 2D shared MPI real(wp) array pointer
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_real_2dim_wp_ptr (array, nelem, window, comm)

        real(wp),        pointer,  intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiint) :: ierror

        if (shared_enabled) then
            if (window == -1) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_deallocate_real_2dim_wp_ptr', &
                                'Illegal shared memory window declared', 1, 1)
            else
                call mpi_memory_synchronize(window, comm)
                call mpi_memory_win_free(window, ierror)
            end if
        else
            deallocate (array)
        end if

    end subroutine mpi_memory_deallocate_real_2dim_wp_ptr


    !> \brief  Deallocate 2D shared MPI real(ep) array pointer
    !> \author Z Masin
    !> \date   2017
    !>
    subroutine mpi_memory_deallocate_real_2dim_ep_ptr (array, nelem, window, comm)

        real(ep),        pointer,  intent(inout) :: array(:,:)
        integer,                   intent(in)    :: nelem, window
        integer(mpiint), optional, intent(in)    :: comm

        integer(mpiint) :: ierror

        if (shared_enabled) then
            if (window == -1) then
                call mpi_xermsg('mpi_memory_mod', 'mpi_memory_deallocate_real_2dim_ep_ptr', &
                                'Illegal shared memory window declared', 1, 1)
            else
                call mpi_memory_synchronize(window, comm)
                call mpi_memory_win_free(window, ierror)
            end if
        else
            deallocate (array)
        end if

    end subroutine mpi_memory_deallocate_real_2dim_ep_ptr


    !> \brief  MPI shared memory barrier
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    !> If the sub-group communicator 'comm' is not given, then the memory will be synchronized on 'shared_communicator',
    !> which is set up by mpi_mod.
    !>
    subroutine mpi_memory_synchronize (window, comm)

        integer,                   intent(in) :: window
        integer(mpiint), optional, intent(in) :: comm

#if (defined(usempi) && defined(mpithree))
        integer(mpiint) :: win, ierror

        win = window

        if (win /= -1) then
            call MPI_Win_fence(0_mpiint, win, ierror)
        end if

        if (present(comm)) then
            call MPI_Barrier(comm, ierror)
        else
            call MPI_Barrier(shared_communicator, ierror)
        end if
#endif

    end subroutine mpi_memory_synchronize


    !> \brief  MPI shared memory deallocation
    !> \author Z Masin, J Benda
    !> \date   2017 - 2019
    !>
    subroutine mpi_memory_win_free (window, ierror)

        integer,         intent(in)    :: window
        integer(mpiint), intent(inout) :: ierror

#if (defined(usempi) && defined(mpithree))
        integer(mpiint) :: win

        win = window
        call MPI_Win_free(win, ierror)
#endif

    end subroutine mpi_memory_win_free

end module mpi_memory_mod
