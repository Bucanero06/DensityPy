Module Module_Becke

    use, intrinsic :: ISO_FORTRAN_ENV

    implicit none

    private

    logical :: Verbous

    public :: &
        EuclDist, &
        EllipticalCoord, &
        Radius_Table, &
        get_param_a, &
        new_mu_transformation, &
        fkfun, &
        skfun, &
        skfunab, &
        Pkfuna, &
        PkFunTot, &
        wkfun

contains
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! FUNCTIONS ... **notes can be found right below, while  documentation is present right above each function.**
    !Radius_Table: Given the name of an atom, returns its atomic radius. Only supports H, C, N, O, and Mg.
    !get_param_a: Given two radii R_i and R_j, returns a transformed parameter a_ij. If the absolute value of a_ij is greater than 0.5, it's set to 0.5.
    !new_mu_transformation: Returns a new mu parameter, given the original mu and a_ij.
    !wkfun: Calculates a weighted function for atom i.
    !PkFunTot: Sums up the function Pkfuna over all atoms.
    !Pkfuna: Calculates a product of skfunab over all atoms excluding atom i.
    !skfunab: Transforms an elliptical coordinate to a new mu, then returns the result of skfun with this new mu.
    !skfun: Transforms a mu parameter into a step function.
    !fkfun: Recursively calculates a function of mu.
    !EllipticalCoord: Given coordinates rvec, avec, and bvec, calculates an elliptical coordinate.
    !EuclDist: Given two points in space, calculates their Euclidean distance.
    !The main considerations:
    !This code does not handle the case where the input atomic name does not match one of the provided options in Radius_Table. This could lead to an error or unexpected behavior.
    !The recursion in fkfun does not seem to have a base case when k is less than 1, which could potentially lead to a stack overflow error.
    !There are no checks for zero or negative distances in EuclDist which might cause issues if not handled properly upstream.
    !The input values are not validated for possible errors. Adding error handling and input validation could make the code more robust.
    !The functions are tightly coupled, which might make modifications and debugging more difficult. Consider making them more modular.
    ! ----------------------------------------------------
    !Result =P_a(\vec{r})/\sum_bP_b(\vec{r}) $
    ! ----------------------------------------------------
    !.. Euclidean distance
    ! This function calculates the Euclidean distance between two points in 3D space. The inputs are two 3-component vectors representing the coordinates of the points, and the output is a scalar representing the distance between them.
    real(kind(1d0)) function EuclDist(avec, bvec) result(res)
        real(kind(1d0)), intent(in) :: avec(3), bvec(3)
        real(kind(1d0)) :: cvec(3)
        cvec = (avec - bvec)
        res = sqrt(sum(cvec * cvec))
    end function EuclDist
    ! ----------------------------------------------------
    !.. Eq. 1 $\mu_{ab}(\vec{r})$
    ! This function calculates the elliptical coordinate ($\mu_{ab}$) of a point relative to two other points (defined as atoms 'a' and 'b'). This coordinate is a measure of how much closer the point is to atom 'a' compared to atom 'b'.
    real(kind(1d0)) function EllipticalCoord(rvec, avec, bvec) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        real(kind(1d0)) :: r_a, r_b, R_ab ! $r_a$, $r_b$, $R_ab$
        r_a = EuclDist(rvec, avec)
        r_b = EuclDist(rvec, bvec)
        R_ab = EuclDist(avec, bvec)
        res = (r_a - r_b) / R_ab
    end function EllipticalCoord
    ! ----------------------------------------------------
    !..Radius_Table for atomic radius based on atomic name
    ! This function provides the atomic radius based on the atomic name. The input is a string representing the atomic name and the output is the corresponding atomic radius. This is useful for defining the size of atoms in molecular modelling.
    real(kind(1d0)) function Radius_Table(Atomic_Name) result(Atomic_Radius)
        character(len = 16), intent(in) :: Atomic_Name

        if (Atomic_Name == "H") then
            Atomic_Radius = 0.35

        elseif (Atomic_Name == "C") then
            Atomic_Radius = 0.70

        elseif (Atomic_Name == "N") then
            Atomic_Radius = 0.65

        elseif (Atomic_Name == "O") then
            Atomic_Radius = 0.60

        elseif (Atomic_Name == "Mg") then
            Atomic_Radius = 1.50
        end if

    end function Radius_Table
    ! ----------------------------------------------------
    !..get_param_aij
    ! This function calculates the parameter 'a' (a_ij) which depends on the radii of two atoms. This parameter is used in the transformation of the elliptical coordinate.
    real(kind(1d0)) function get_param_a(R_i, R_j) result(a_ij)
        real(kind(1d0)), intent(in) :: R_i, R_j
        real(kind(1d0)) :: Chi
        Chi = R_i / R_j
        a_ij = (Chi - 1) / (Chi + 1)
        if (abs(a_ij)>(0.5d0)) then
            a_ij = 0.5d0
        endif
    end function get_param_a
    ! ----------------------------------------------------
    !..new_mu_transformation
    ! This function transforms the elliptical coordinate using the parameter 'a' (a_ij). This transformation is used to handle cases where the point is much closer to one atom than the other.
    real(kind(1d0)) function new_mu_transformation(mu, a_ij) result(res)
        real(kind(1d0)), intent(in) :: mu, a_ij
        res = mu + a_ij * (1 - mu**2)
    end function new_mu_transformation
    ! ----------------------------------------------------
    !.. Eq. 4 Recursive $f_k(\mu)$
    ! This function calculates a recursive function of the transformed elliptical coordinate. The output of this function will be used in the calculation of the step function skfun.
    real(kind(1d0)) recursive function fkfun(mu, k) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu
        integer, intent(in) :: k
        res = 0.5d0 * mu * (3.d0 - mu**2)
        if(k==1)return
        res = fkfun(res, k - 1)
    end function fkfun
    ! ----------------------------------------------------
    !.. $s(\mu): [-1,1]  \rightarrow [0,1]$  (Step Function)
    ! This function calculates a step function of the transformed elliptical coordinate. This function maps the transformed elliptical coordinate from the range [-1,1] to [0,1].
    real(kind(1d0)) function skfun(mu, k) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu
        integer, intent(in) :: k
        real(kind(1d0)) :: r_m

        res = 0.5d0 * (1.d0 - fkfun(mu, k))
    end function skfun
    ! ----------------------------------------------------
    !.. This function combines the calculations of the transformed elliptical coordinate and the step function. The output of this function will be used in the calculation of the partial partition function Pkfuna.
    real(kind(1d0)) function skfunab(rvec, avec, bvec, k, a_ij) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        integer, intent(in) :: k
        real(kind(1d0)) :: mu_old, mu_new, a_ij
        mu_old = EllipticalCoord(rvec, avec, bvec)
        mu_new = new_mu_transformation(mu_old, a_ij)
        res = skfun(mu_new, k)
    end function skfunab
    ! ----------------------------------------------------
    !.. Eq. 2 P_a(\vec{r}) in and nominator of Eq. 3
    ! This function calculates the partial partition function for a specific atom. The output of this function represents the contribution of a specific atom to the total partition function.
    real(kind(1d0)) function Pkfuna(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS) result(res)
        integer, intent(in) :: nAtoms, iAtom
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms), Radius_BS(:)
        integer, intent(in) :: k
        integer :: jAtom
        real(kind(1d0)) :: a_ij, R_i, R_j
        res = 1.d0
        do jAtom = 1, nAtoms
            if(jAtom == iAtom) cycle
            R_i = Radius_BS(iAtom)
            R_j = Radius_BS(jAtom)
            a_ij = get_param_a(R_i, R_j)
            res = res * skfunab(rvec, AtCoord(:, iAtom), AtCoord(:, jAtom), k, a_ij)
        enddo
    end function Pkfuna
    ! ----------------------------------------------------
    !..Eq. 3 \sum_bP_b(\vec{r}) in the denominator
    ! This function calculates the total partition function. This is done by summing the partial partition functions for all atoms.
    real(kind(1d0)) function PkFunTot(rvec, AtCoord, nAtoms, k, Radius_BS) result(res)
        integer, intent(in) :: nAtoms
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms), Radius_BS(:)
        integer, intent(in) :: k
        integer :: iAtom
        res = 0.d0
        do iAtom = 1, nAtoms
            res = res + Pkfuna(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS)
        enddo
    end function PkFunTot
    ! ----------------------------------------------------
    !.. Eq. 3 $w_a \vec(r)$  (WEIGHTS)
    ! This function calculates the weight of an atom. The weight is defined as the ratio of the partial partition function of the atom to the total partition function. This weight represents the relative importance of the atom.
    real(kind(1d0)) function wkfun(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS) result(res)
        !
        integer, intent(in) :: iAtom, k, nAtoms
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms), Radius_BS(:)
        !
        res = Pkfuna(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS) / PkfunTot(rvec, AtCoord, nAtoms, k, Radius_BS)
    end function wkfun
    ! ----------------------------------------------------
    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


End Module Module_Becke


















