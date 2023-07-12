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
module orbital_routines
   use precisn
   use mpi_mod

   private

   public GS_ortho_routine, SYM_ortho_routine, check_orbital, write_orbital, read_orbital

contains

   !> Checks that data for orbital_data_obj are internally consistent.
   function check_orbital(number_of_functions,number_of_coefficients,coefficients,spin,energy,occup,irr,point_group)
      implicit none
      integer, intent(in) :: number_of_functions,number_of_coefficients,irr,point_group
      real(kind=cfp), allocatable :: coefficients(:,:),energy(:),occup(:)
      integer, allocatable :: spin(:)
      integer :: check_orbital
 
         check_orbital = 0

         if (number_of_functions < 0) then
            print *,number_of_functions,irr
            check_orbital = 1
            return
         endif

         if (number_of_coefficients < 0) then
            print *,number_of_coefficients
            check_orbital = 2
            return
         endif

         if (irr .le. 0 .or. irr > 8) then
            check_orbital = 7
            return
         endif

         if (point_group .le. 0 .or. point_group > 8) then
            check_orbital = 8
            return
         endif

         if (number_of_coefficients > 0 .and. number_of_functions > 0) then
            if (size(coefficients,1) .ne. number_of_coefficients .or. size(coefficients,2) .ne. number_of_functions) then
               print *,size(coefficients,1),size(coefficients,2),number_of_coefficients,number_of_functions
               check_orbital = 3
               return
            endif
   
            if (size(spin,1) .ne. number_of_functions) then
               check_orbital = 4
               return
            endif
   
            if (size(energy,1) .ne. number_of_functions) then
               check_orbital = 5
               return
            endif
   
            if (size(occup,1) .ne. number_of_functions) then
               check_orbital = 6
               return
            endif
         endif

   end function check_orbital

   subroutine write_orbital(number_of_functions,number_of_coefficients,coefficients,spin,energy,&
                            occup,irr,point_group,norm,lunit,posit,pos_after_rw)
      use utils
      implicit none
      integer, intent(in) :: number_of_functions,number_of_coefficients,irr,point_group
      real(kind=cfp), allocatable :: coefficients(:,:),energy(:),occup(:)
      integer, allocatable :: spin(:)
      real(kind=cfp), intent(in) :: norm
      integer, intent(in) :: lunit, posit
      integer, intent(out) :: pos_after_rw

         pos_after_rw = 0
         if (myrank .eq. master) then
            write(lunit,pos=posit,err=10) number_of_functions, number_of_coefficients
            write(lunit,err=10) point_group, irr

            if (number_of_functions > 0 .and. number_of_coefficients > 0) then
               write(lunit,err=10) coefficients(1:number_of_coefficients,1:number_of_functions)
               write(lunit,err=10) energy(1:number_of_functions)
               write(lunit,err=10) spin(1:number_of_functions)
               write(lunit,err=10) occup(1:number_of_functions)
               write(lunit,err=10) norm
            endif

            inquire(lunit,pos=pos_after_rw)
         endif

         !master ensures all processes know where the record ends
         call mpi_mod_bcast(pos_after_rw,master)

         return

 10      call xermsg ('orbital_routines', 'write_orbital', 'Error writing the orbital data into the file and position given.', 1, 1)

   end subroutine write_orbital

   subroutine read_orbital(number_of_functions,number_of_coefficients,coefficients,spin,&
                           energy,occup,irr,point_group,norm,lunit,posit,pos_after_rw)
      use utils
      implicit none
      integer, intent(out) :: number_of_functions,number_of_coefficients,irr,point_group
      real(kind=cfp), allocatable :: coefficients(:,:),energy(:),occup(:)
      integer, allocatable :: spin(:)
      real(kind=cfp), intent(out) :: norm
      integer, intent(in) :: lunit, posit
      integer, intent(out) :: pos_after_rw

      integer :: err, i

         if (allocated(coefficients)) deallocate(coefficients)
         if (allocated(energy)) deallocate(energy)
         if (allocated(spin)) deallocate(spin)
         if (allocated(occup)) deallocate(occup)

         pos_after_rw = 0
         if (myrank .eq. master) then
            read(lunit,pos=posit,err=10) number_of_functions, number_of_coefficients
            read(lunit,err=10) point_group, irr

            if (number_of_functions > 0 .and. number_of_coefficients > 0) then
               allocate(coefficients(number_of_coefficients,number_of_functions),&
                        energy(number_of_functions),spin(number_of_functions),occup(number_of_functions),stat=err)
               if (err /= 0) then
                  call xermsg ('orbital_routines', 'read_orbital', 'Memory allocation error on the master task.', err, 1)
               end if
   
               read(lunit,err=10) coefficients(1:number_of_coefficients,1:number_of_functions)
               read(lunit,err=10) energy(1:number_of_functions)
               read(lunit,err=10) spin(1:number_of_functions)
               read(lunit,err=10) occup(1:number_of_functions)
               read(lunit,err=10) norm
            endif
            inquire(lunit,pos=pos_after_rw)
         endif

         !master ensures all processes know where the record ends
         call mpi_mod_bcast(pos_after_rw,master)

         !master broadcasts all its data to the other processes: start with scalars
         call mpi_mod_bcast(number_of_coefficients,master)
         call mpi_mod_bcast(number_of_functions,master)
         call mpi_mod_bcast(point_group,master)
         call mpi_mod_bcast(irr,master)

         if (number_of_functions > 0 .and. number_of_coefficients > 0) then
            call mpi_mod_bcast(pos_after_rw,master)
            call mpi_mod_bcast(norm,master)
   
            !all other processes allocate space for the coefficients which will be broadcast last
            if (myrank .ne. master) then
               allocate(coefficients(number_of_coefficients,number_of_functions), &
                        energy(number_of_functions),spin(number_of_functions),occup(number_of_functions),stat=err)
               if (err /= 0) then
                  call xermsg ('orbital_data_obj', 'read_orbital', 'Memory allocation error on the child process(es).', err, 1)
               end if
            endif
   
            !finally arrays:
            !todo this loop is lame the bcast should include bcast of 2d arrays!
            do i=1,number_of_functions
               call mpi_mod_bcast(coefficients(1:number_of_coefficients,i),master)
            enddo
            call mpi_mod_bcast(energy,master)
            call mpi_mod_bcast(occup,master)
            call mpi_mod_bcast(spin(1:number_of_functions),master)
         else
            if (allocated(coefficients)) deallocate(coefficients)
            if (allocated(energy)) deallocate(energy)
            if (allocated(spin)) deallocate(spin)
            if (allocated(occup)) deallocate(occup)
         endif

         return

 10      call xermsg ('orbital_routines', 'read_orbital', &
                      'Error reading the orbital data from the file at the position given.', 1, 1)

   end subroutine read_orbital

end module orbital_routines
