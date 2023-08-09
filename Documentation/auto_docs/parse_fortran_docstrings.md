### Gets This:

```
============================================================================
subroutine Write_R_el_bc
============================================================================

:Synopsis:
   This subroutine writes atomic positions (R_el) to a file in the specified output directory.

:Parameters:
   :output_directory: (character(len = *), intent(in))
       Directory path to save the file.
   :atom_names: (character(len = 16), allocatable, intent(in))
       Names of atoms.
   :nAtoms: (integer, intent(in))
       Number of atoms.
   :R_el: (real(kind(1d0)), allocatable, intent(in))
       Atomic positions.

:Local Variables:
   :uid: (integer)
       File unit identifier.
   :iAtom, iPol: (integer)
       Loop counters.

============================================================================
```

### From this:

```fortran
subroutine Write_R_el_bc(output_directory, atom_names, nAtoms, R_el)
    !s* Synopsis: This subroutine writes atomic positions (R_el) to a file in the specified
    !*  output directory.

    ! Inputs:
    character(len = *), intent(in) :: output_directory !* Directory path to save the file
    character(len = 16), allocatable, intent(in) :: atom_names(:) !* Names of atoms
    integer, intent(in) :: nAtoms !* Number of atoms
    real(kind(1d0)), allocatable, intent(in) :: R_el(:, :) !* Electronic Barycenter positions for each atom

    ! Local variables:
    integer :: uid !* File unit identifier 
    integer :: iAtom, iPol !* Loop counters

    !* Attempt to open the file for writing
    open(newunit = uid, file = output_directory // "/" // "R_el_bc", form = "formatted", status = "unknown")

    ! Write headers with fixed-width format
    write(uid, "(A14, 4(A20))") "Atom Index", "Atom Name", "X Position", "Y Position", "Z Position"

    ! Loop over atoms and write their positions to the file
    do iAtom = 1, nAtoms
        write(uid, "(i5, A20, 3(f20.14))") iAtom, trim(atom_names(iAtom)), (R_el(iPol, iAtom), iPol = 1, 3)
    end do

    ! Close the file after writing
    close(uid)

end subroutine Write_R_el_bc
```

1.What needs to happen is all subroutines/functions/etc... need to be parsed and categorized on their type of
procedure/type
2.Then, each procedure needs to be parsed for its parameters and their types, intents, allocatability, etc...,
subroutine/function/module/program will be refered to as "procedure" from now on

* Rules for parsing:
    * Docstring Comments contain "!" followed by its [optional] argument "{x}" followed by "*"
    * Docstring Comments are the lines we care about
    * Ignore comments that start with "!" but do not contain "*"
    * The procedure name is the first word after the procedure keyword
    * Procedure arguements are parsed and labeled based on the declaration of the procedure, e.g. inputs, outputs, local
      variables, etc...
    * Actions - Links docstring with a given line of code. Depends on the contents of the line. If variable declaration
      then it acts as the description of the variable, else acts as a step in the procedure:
        * ">" — Links to the right meaning it can only be a comment, so this indicates a step, to be used when writing
          the documentation to indicate a step in the procedure. These will be numbered automatically
        * "<" — Links to the left
        * "^" - Links up
        * "v" - Links down
        * "Synopsis" or "S" - The synopsis/description of the procedure. This acts as "S" + "v" which will link to the
          next lines of !* as a multiline synopsis/description until an empty line is reached, a native comment, line of
          code, a stop link like "-", or another docstring action is reached. If not set, the first line of docstring
          comments found in a subroutine/function/module/program is the synopsis/description
    * If !Section {name}* is found outside of a procedure, it is still used as documentation, but it is not parsed for
      arguments and
      acts as a native comment which will appear in the location it is found in the file. !Section* can be used but the
      first Word will be used as the section name which might lead to errors
