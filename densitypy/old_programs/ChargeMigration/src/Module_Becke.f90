Module Module_Becke

    use, intrinsic :: ISO_FORTRAN_ENV

    implicit none

    private

    logical :: Verbous

    public :: &
            ComputeAtomicWeights


contains
    !>Determines 1 Volume Unit of the Grid
    !
    subroutine ComputeVolume(nPts, Volume, gridv)
        integer, intent(in) :: nPts
        real(kind(1d0)), intent(out) :: Volume
        real(kind(1d0)), allocatable, intent(in) :: gridv(:, :)
        real(kind(1d0)) :: coord1, coord2
        integer :: iPts, iCoord
        Volume = 1.d0
        do iCoord = 1, 3
            coord1 = gridv(iCoord, 1)
            do iPts = 2, nPts
                coord2 = gridv(iCoord, iPts)
                if (coord1 .ne. coord2) then
                    Volume = Volume * abs(coord2 - coord1)
                    exit
                endif
            enddo
        enddo
    end subroutine ComputeVolume

    subroutine ComputeAtomicCharges(OrbitalDensity, BeckeMatrix, QchargeVec)
        real(kind(1d0)), intent(in) :: OrbitalDensity(:, :)
        real(kind(1d0)), intent(in) :: BeckeMatrix   (:, :, :)
        real(kind(1d0)), allocatable, intent(out) :: QchargeVec    (:)
        integer :: iAtom, nAtoms, nOrbs
        real(kind(1d0)), external :: DDOT
        nOrbs = size(BeckeMatrix, 1)
        nAtoms = size(BeckeMatrix, 3)
        if(.not.allocated(QchargeVec))allocate(QchargeVec(nAtoms))
        QchargeVec = 0.d0
        do iAtom = 1, nAtoms
            QchargeVec(iAtom) = ddot(nOrbs * nOrbs, OrbitalDensity, 1, BeckeMatrix(1, 1, iAtom), 1)
        enddo
    end subroutine ComputeAtomicCharges

    subroutine ComputeBeckeMatrix(WeightV, OrbTab, Becke)

        real(kind(1d0)), allocatable, intent(in) :: WeightV(:, :)
        real(kind(1d0)), allocatable, intent(in) :: OrbTab (:, :)
        real(kind(1d0)), allocatable, intent(out) :: Becke(:, :, :)

        integer :: nPts, nOrbs, nAtoms
        integer :: iPts, iOrb1, iOrb2, iAtom
        real(kind(1d0)) :: dsum

        nAtoms = size(WeightV, 2)
        nPts = size(WeightV, 1)
        nOrbs = size(OrbTab, 2)

        if(allocated(Becke))deallocate(Becke)
        allocate(Becke(nOrbs, nOrbs, nAtoms))

        Becke = 0.d0
        do iAtom = 1, nAtoms
            do iOrb2 = 1, nOrbs
                do iOrb1 = 1, nOrbs
                    dsum = 0.d0
                    do iPts = 1, nPts
                        dsum = dsum + WeightV(iPts, iAtom) * OrbTab(iPts, iOrb1) * OrbTab(iPts, iOrb2)
                    enddo
                    Becke(iOrb1, iOrb2, iAtom) = dsum
                enddo
            enddo
        enddo

    end subroutine ComputeBeckeMatrix


    subroutine ComputeAtomicWeights(nPts, gridv, nAtoms, AtCoord, WeightV, Radius_BS)
        integer, intent(in) :: nPts, nAtoms
        real(kind(1d0)), intent(in) :: gridv(:, :), AtCoord(:, :), Radius_BS(:)
        real(kind(1d0)), allocatable, intent(out) :: WeightV(:, :)

        integer, parameter :: k = 4
        integer :: iPts, iAtom
        real(kind(1d0)) :: rvec(3)

        allocate(WEIGHTV(nPts, nAtoms))
        do iPts = 1, nPts
            rvec = gridv(:, iPts)
            do iAtom = 1, nAtoms
                WeightV(iPts, iAtom) = wkfun(rvec, iAtom, AtCoord, nAtoms, k, Radius_BS)
            enddo
        enddo
    end subroutine ComputeAtomicWeights
    !..FUNCTIONS
    !
    !.. Eq. 3 $w_a \vec(r)$  (WEIGHTS)
    !Result =P_a(\vec{r})/\sum_bP_b(\vec{r}) $
    real(kind(1d0)) function wkfun(rvec, iAtom, AtCoord, &
            nAtoms, k, Radius_BS) result(res)
        !
        integer, intent(in) :: iAtom, k, nAtoms
        real(kind(1d0)), intent(in) :: rvec(3), AtCoord(3, nAtoms), Radius_BS(:)
        !
        res = Pkfuna(rvec, iAtom, AtCoord, nAtoms, k, R_i, R_j) / PkfunTot(rvec, AtCoord, nAtoms, k, Radius_BS)
    end function wkfun

    !..Eq. 3 \sum_bP_b(\vec{r}) in the denominator
    !
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

    !.. Eq. 2 P_a(\vec{r}) in and nominator of Eq. 3
    !
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
            res = res * skfunab(rvec, AtCoord(:, iAtom), AtCoord(:, jAtom), k, a_ij, R_i)
        enddo
    end function Pkfuna


    !..
    !
    real(kind(1d0)) function skfunab(rvec, avec, bvec, k, a_ij, R_i) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        integer, intent(in) :: k
        real(kind(1d0)) :: mu, a_ij, R_i
        mu = EllipticalCoord(rvec, avec, bvec)
        mu = new_mu_transformation(mu, a_ij)

        res = skfun(mu, k, R_i)
    end function skfunab

    !.. $s(\mu): [-1,1]  \rightarrow [0,1]$  (Step Function)
    !
    real(kind(1d0)) function skfun(mu, k, R_i) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu, R_i
        integer, intent(in) :: k

        res = 0.5d0 * (1.d0 - fkfun(mu, k)) * (0.5d0 * R_i)
    end function skfun

    !.. Eq. 4 Recursive $f_k(\mu)$
    !
    real(kind(1d0)) recursive function fkfun(mu, k) result(res)
        implicit none
        real(kind(1d0)), intent(in) :: mu
        integer, intent(in) :: k
        res = 0.5d0 * mu * (3.d0 - mu**2)
        if(k==1)return
        res = fkfun(res, k - 1)
    end function fkfun

    !.. Eq. 1 $\mu_{ab}(\vec{r})$
    !
    real(kind(1d0)) function EllipticalCoord(rvec, avec, bvec) result(res)
        real(kind(1d0)), intent(in) :: rvec(3), avec(3), bvec(3)
        real(kind(1d0)) :: r_a, r_b, R_ab ! $r_a$, $r_b$, $R_ab$
        r_a = EuclDist(rvec, avec)
        r_b = EuclDist(rvec, bvec)
        R_ab = EuclDist(avec, bvec)
        res = (r_a - r_b) / R_ab
    end function EllipticalCoord

    !.. Distance
    !
    real(kind(1d0)) function EuclDist(avec, bvec) result(res)
        real(kind(1d0)), intent(in) :: avec(3), bvec(3)
        real(kind(1d0)) :: cvec(3)
        cvec = (avec - bvec)
        res = sqrt(sum(cvec * cvec))
    end function EuclDist
    !
    !END OF FUNCTIONS


    !..get_param_aij
    !..
    real(kind(1d0)) function get_param_a(R_i, R_j) result(a_ij)
        real(kind(1d0)), intent(in) :: R_i, R_j
        real(kind(1d0)) :: Chi
        Chi = R_i / R_j
        a_ij = (Chi - 1) / (Chi + 1)
        if (abs(a_ij)<=1 / 2) then
            return
        elseif (abs(a_ij)>1 / 2) then
            a_ij = 1 / 2
        end if
    end function get_param_a

    real(kind(1d0)) function new_mu_transformation(a_ij) result(res)
        real(kind(1d0)), intent(in) :: mu, a_ij
        res = mu + a_ij * (1 - mu**2)
    end function new_mu_transformation


End Module Module_Becke

subroutine AtomicRadius_Bragg_Slater_Becke(AtomName, nAtom, Radius_BS)
    real(kind((1d0))), allocatable, intent(out) :: Radius_BS
    integer, intent(in) :: nAtom
    character(len = 16), intent(in) :: AtomName(:)
    integer :: iAtom
    allocate(Radius_BS(nAtom))

    do iAtom = 1, nAtom
        Radius_BS(iAtom) = Radius_Table(AtomName(iAtom))
    end do
end subroutine AtomicRadius_Bragg_Slater_Becke


function Radius_Table(Atomic_Name) result(res)
    character(len = 16), intent(in) :: Atomic_Name
    real(kind(1d0)), intent(out) :: Atomic_Radius

    if (Atomic_Name == "H") then
        Atomic_Radius = 0.35

    elseif (Atomic_Name == "C") then
        Atomic_Radius = 0.70

    elseif (Atomic_Name == "N") then
        Atomic_Radius = 0.65

    elseif (Atomic_Name == "O") then
        Atomic_Radius = 0.60
    end if

    res = Atomic_Radius
end function Radius_Table
















