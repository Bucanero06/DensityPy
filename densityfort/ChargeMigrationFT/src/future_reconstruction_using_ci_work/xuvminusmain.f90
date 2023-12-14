! {{{ Detailed description

!> \mainpage Program <ProgramName> <Insert here what the program does>
!!
!! Synopsis:
!! ---------
!!
!!     <Program Name> <mandatory run-time parameters (RTP)> [<optional RTP>]
!!
!! ___
!! Description:
!! ------------
!!
!! Input parameters:      {#Input_Parameters}
!! =================
!! [...Input](@ref ...Input) as specified in the command line.
!!
!> \file
!!
!!
! }}}
program ChargeMigration

    use, intrinsic :: ISO_FORTRAN_ENV
    use ModuleErrorHandling
    use ModuleSystemUtils
    use ModuleString
    use ModuleIO
    use ModuleConstants
    use ModuleDiagonalize
    use ModulePulses_3D
    !
    !.. Local modules
    use ModuleRTP
    use Module_CD_IO

    implicit none

    !.. Run-time parameters
    !..
    character(len = :), allocatable :: InpDir
    character(len = :), allocatable :: OutDir
    character(len = :), allocatable :: FileGeometry
    real(kind(1d0)) :: StepTime
    real(kind(1d0)) :: StepWidth
    character(len = :), allocatable :: Ext_field_file
    logical :: Verbous
    integer :: nOmegas
    real   (kind(1d0)) :: OmegaMin
    real   (kind(1d0)) :: OmegaMax
    integer :: nTauOmegas
    real   (kind(1d0)) :: TauOmegaMin
    real   (kind(1d0)) :: TauOmegaMax

    integer, parameter :: GS_IDX = 1

    character(len = *), parameter :: DIPOLE_FT_PATH_ALL = "/Dipole/DipoleFT_ALL"
    character(len = *), parameter :: CHARGE_FT_PATH_ALL = "/AtomicCharge/AtomicChargeFT_ALL"

    real(kind(1d0)), external :: NCD_Phi

    !.. Data for FT
    !..
    integer :: iOmega
    real   (kind(1d0)), allocatable :: OmegaVec(:), TauOmegaVec(:)
    complex(kind(1d0)) :: zExpFact
    complex(kind(1d0)), allocatable :: DipoleFTtotal(:, :)
    complex(kind(1d0)), allocatable :: DipoleFTminus(:, :)
    complex(kind(1d0)), allocatable :: DipoleFTplus (:, :)
    complex(kind(1d0)), allocatable :: AtomicChargeFT(:, :)
    complex(kind(1d0)), allocatable :: DipoleFTwt(:, :, :)
    complex(kind(1d0)), allocatable :: DipoleFTww(:, :, :)
    complex(kind(1d0)), allocatable :: ChargeFTwt(:, :, :)
    complex(kind(1d0)), allocatable :: ChargeFTww(:, :, :)
    real   (kind(1d0)) :: dBuf, drx, dix, dry, diy, drz, diz, tStep1, tStep2
    real   (kind(1d0)), allocatable :: tvec(:)
    real   (kind(1d0)), allocatable :: dAtomFTRe(:), dAtomFTIm(:)
    integer :: iOmegaTau, iBuf, ntimes
    real(kind(1d0)) :: tmin, tmax


    !.. Local parameters
    integer :: nStates, i
    real(kind(1d0)), allocatable :: Evec(:)
    !.. Expectation Value of the Dipole Moment (Mu)
    complex(kind(1d0)), allocatable :: zMuEV(:, :)

    integer :: nAtoms
    real(kind(1d0)), allocatable :: AtCoord(:, :) ! 3 x npts
    real(kind(1d0)) :: dstep

    !.. Statistical Density Matrix
    real   (kind(1d0)) :: t, dt, w
    real   (kind(1d0)), allocatable :: dvec1(:), dvec2(:)
    real   (kind(1d0)), allocatable :: AtomicChargeEvolution(:, :)

    character(len = 16), allocatable :: AtomName(:)

    integer :: it, iPol, iAtom
    integer :: uid_AtomicCharge, uid_AtomicChargeALL, uid_AtomicChargeFT, uid_dipole, uid_dipoleFT

    !.. Pulse parameters
    integer :: iSim, N_Simulations, uid, iostat
    character(len = 64), pointer :: Simulation_Tagv(:)
    type(pulse_train), pointer :: train(:)

    !..
    !    B_{ij}^{\alpha} = \int d^3r \phi_i(\vec{r})\phi_j(\vec{r}) w_\alpha(\vec{r})
    !                    = BeckeMatrix(i,j,alpha)
    !..

    complex(kind(1d0)), external :: zdotu

    !!!!!!!!!!!!!!!!!!!!!!Rubenjan21!!!!!!!!!
    complex(kind(1d0)), allocatable :: XUVDipole(:, :)
    complex(kind(1d0)), allocatable :: XUVDipoleFT(:, :)
    complex(kind(1d0)), allocatable :: DipoleFTtotalminXUV(:, :)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call GetRunTimeParameters(InpDir, OutDir, FileGeometry, StepTime, StepWidth, &
            Ext_field_file, verbous, nOmegas, OmegaMin, OmegaMax, nTauOmegas, TauOmegaMin, TauOmegaMax)

    allocate(OmegaVec, source = [(OmegaMin + (OmegaMax - OmegaMin) / dble(nOmegas - 1) * iOmega, iOmega = 1, nOmegas)])
    allocate(TauOmegaVec, source = [(TauOmegaMin + (TauOmegaMax - TauOmegaMin) / dble(nTauOmegas - 1) * iOmega, iOmega = 1, nTauOmegas)])

    call Set_CD_IO_Verbous(.FALSE.)

    open(newunit = uid, &
            file = Ext_field_file, &
            form = "formatted", &
            status = "old")
    call Parse_Simulation_File(uid, OUTPUT_UNIT, N_Simulations, Simulation_Tagv, train)
    close(uid)
    write(*, "(a,i3)")"N_Simulations=", N_Simulations

    call LoadEnergies(InpDir // "/ROOT_ENERGIES", nStates, Evec)
    call LoadGeometry(nAtoms, AtCoord, FileGeometry, AtomName)

    !.. Determines the number of times available, as well as the minimum
    !   and maximum time
    !..
    open(newunit = uid_dipole, &
            file = OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(1)), &
            form = "formatted", &
            status = "old", &
            action = "read")
    !
    read(uid_dipole, *, iostat = iostat) iBuf, dBuf
    if(iostat/=0)then
        write(*, *) OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(1)) // "   IS EMPTY"
        stop
    end if
    tmin = dBuf
    tmax = tmin
    nTimes = 1
    do
        read(uid_dipole, *, iostat = iostat) iBuf, dBuf
        if(iostat/=0)exit
        tmax = dBuf
        nTimes = nTimes + 1
    enddo
    close(uid_dipole)
    dt = (tmax - tmin) / dble(nTimes - 1)
    write(*, *) "nTimes = ", nTimes

    allocate(zMuEV(3, nTimes))
    zMuEV = Z0
    allocate(AtomicChargeEvolution(nAtoms, nTimes))
    allocate(AtomicChargeFT(nOmegas, nAtoms))
    AtomicChargeFT = Z0
    allocate(dvec1(3), dvec2(3))

    allocate(DipoleFTtotal(3, nOmegas))
    allocate(DipoleFTminus(3, nOmegas))
    allocate(DipoleFTplus (3, nOmegas))

    !!!!!!!!!!!!!!Rubenjan21!!!!!!!!!!!!
    allocate(XUVDipole(3, ntimes))
    XUVDipole = Z0
    allocate(XUVDipoleFT(3, nOmegas))
    allocate(DipoleFTtotalminXUV(3, nOmegas))

    !.. Load XUVDipole from file
    call Load_XUVDipole(OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(N_Simulations)), XUVDipole, nTimes)

    !..Regularize XUVDipole
    do it = 1, nTimes
        t = tmin + dt * dble(it - 1)
        dstep = 1.d0 - NCD_Phi(t, StepTime, StepWidth)
        XUVDipole(:, it) = XUVDipole(:, it) * dstep
    enddo


    !..Compute FT of XUV
    XUVDipoleFT = Z0
    do iOmega = 1, nOmegas
        w = OmegaVec(iOmega)
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            zexpFact = exp(Zi * w * t)
            XUVDipoleFT(:, iOmega) = XUVDipoleFT(:, iOmega) + zExpFact * XUVDipole(:, it)
        enddo
    end do
    XUVDipoleFT = XUVDipoleFT * dt / (2.d0 * PI)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !.. Open Q_Charge File
    open(newunit = uid_AtomicChargeALL, &
            file = OutDir // "/AtomicCharge/AtomicCharge_ALL", &
            form = "formatted", &
            status = "unknown", &
            action = "write")

    N_Simulations = N_Simulations - 1 !Last simulation is the XUV no probe
    Sim_loop : do iSim = 1, N_Simulations

        !.. Load Dipole from file
        !..
        open(newunit = uid_dipole, &
                file = OutDir // "/Dipole/Dipole" // trim(Simulation_tagv(iSim)), &
                form = "formatted", &
                status = "old", &
                action = "read")
        write(*, *) "nTimes=", nTimes
        do it = 1, nTimes
            read(uid_dipole, "(i4,*(x,E24.16))") iBuf, dBuf, ((dvec1(iPol), dvec2(iPol)), iPol = 1, 3)
            write(*, *) it, iBuf, dBuf
            do iPol = 1, 3
                zMuEV(iPol, it) = Z1 * dvec1(iPol) + Zi * dvec2(iPol)
            enddo
        enddo
        close(uid_dipole)

        !.. Open Q_Charge File
        open(newunit = uid_AtomicCharge, &
                file = OutDir // "/AtomicCharge/AtomicCharge" // trim(Simulation_tagv(iSim)), &
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

        !.. Compute the regularized dipole $\mu_-(t)$ and Atomic Charges
        !..
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            dstep = 1.d0 - NCD_Phi(t, StepTime, StepWidth)
            zMuEV(:, it) = zMuEV(:, it) * dstep
            AtomicChargeEvolution(:, it) = AtomicChargeEvolution(:, it) * dstep
        enddo
        !
        !.. Save Regularized Dipole to file
        open(newunit = uid_dipole, &
                file = OutDir // "/Dipole/Dipole-" // trim(Simulation_tagv(iSim)), &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do it = 1, nTimes
            t = tmin + dt * dble(it - 1)
            write(uid_dipole, "(i4,*(x,E24.16))") it, t, ((dble(zMuEV(iPol, it)), aimag(zMuEV(iPol, it))), iPol = 1, 3)
        enddo
        close(uid_dipole)


        !.. Compute the FT of the dipole $\mu_-(t)$ and of the Atomic Charges
        DipoleFTminus = Z0
        AtomicChargeFT = Z0
        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            do it = 1, nTimes
                t = tmin + dt * dble(it - 1)
                zexpFact = exp(Zi * w * t)
                DipoleFTminus(:, iOmega) = DipoleFTminus(:, iOmega) + zExpFact * zMuEV(:, it)
                AtomicChargeFT(iOmega, :) = AtomicChargeFT(iOmega, :) + &
                        zExpFact * AtomicChargeEvolution(:, it)
            enddo
        end do
        DipoleFTminus = DipoleFTminus * dt / (2.d0 * PI)
        AtomicChargeFT = AtomicChargeFT * dt / (2.d0 * PI)
        !
        !*** MUST GO TO A SEPARATE SUBROUTINE
        !!$     !.. Save Dipole FT File
        !!$     open(newunit = uid_dipoleFT, &
        !!$          file    ="DipoleFT-"//trim(Simulation_tagv(iSim)), &
        !!$          form    ="formatted", &
        !!$          status  ="unknown"  , &
        !!$          action  ="write"    )
        !!$     do iOmega = 1, nOmegas
        !!$        w = OmegaVec( iOmega )
        !!$        write(uid_dipoleFT,"(i4,*(x,E24.16))") iOmega, w, &
        !!$             ((dble(DipoleFTminus(iPol,iOmega)),aimag(DipoleFTminus(iPol,iOmega))),iPol=1,3)
        !!$     end do
        !!$     close( uid_dipoleFT )


        !*** THE FT OF THE DIPOLE FOR LARGE TIME STILL DOES NOT WORK
        DipoleFTplus = Z0
        !!$     call ComputeDipoleFTplus( L0_Eval, L0_LEvec, L0_REvec, &
        !!$          OmegaVec, zStatRho0, DipoleFTplus, StepTime, StepWidth )
        !!$     !
        !!$     !*** MUST GO TO A SEPARATE SUBROUTINE
        !!$     !.. Save Dipole FT File
        !!$     open(newunit = uid_dipoleFT, &
        !!$          file    ="DipoleFT+"//trim(Simulation_tagv(iSim)), &
        !!$          form    ="formatted", &
        !!$          status  ="unknown"  , &
        !!$          action  ="write"    )
        !!$     do iOmega = 1, nOmegas
        !!$        w = OmegaVec( iOmega )
        !!$        write(uid_dipoleFT,"(i4,*(x,E24.16))") iOmega, w, &
        !!$             ((dble(DipoleFTplus(iPol,iOmega)),aimag(DipoleFTplus(iPol,iOmega))),iPol=1,3)
        !!$     end do
        !!$     close( uid_dipoleFT )

        DipoleFTtotal = DipoleFTminus + DipoleFTplus


        !!!!jan21Ruben!!!!!!!!!!!!!!!!!!!!!!1!
        XUVDipoleFT = XUVDipoleFT + DipoleFTplus
        do iPol = 1, 3
            do iOmega = 1, nOmegas !DipoleFTtotalminXUV
                DipoleFTtotal(iPol, iOmega) = DipoleFTtotal(iPol, iOmega) - XUVDipoleFT(iPol, iOmega)
            enddo
        end do

        !.. Write the FT in a single file
        open(newunit = uid_dipoleFT, &
                file = OutDir // DIPOLE_FT_PATH_ALL // "minXUV", &
                form = "formatted", &
                status = "unknown", &
                action = "write", &
                position = "append")
        do iOmega = 1, nOmegas
            call train(iSim)%PrintPulses(uid_dipoleFT)
            w = OmegaVec(iOmega)
            write(uid_dipoleFT, "(i4,*(x,E24.16))") iOmega, w, &
                    ((dble(DipoleFTtotalminXUV(iPol, iOmega)), aimag(DipoleFTtotalminXUV(iPol, iOmega))), iPol = 1, 3)
        end do
        write(uid_dipoleFT, *)
        close(uid_dipoleFT)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        !
        !*** MUST GO TO A SEPARATE SUBROUTINE
        !.. Save Dipole FT File
        open(newunit = uid_dipoleFT, &
                file = OutDir // "/Dipole/DipoleFT" // trim(Simulation_tagv(iSim)), &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            write(uid_dipoleFT, "(i4,*(x,E24.16))") iOmega, w, &
                    ((dble(DipoleFTtotal(iPol, iOmega)), aimag(DipoleFTtotal(iPol, iOmega))), iPol = 1, 3)
        end do
        close(uid_dipoleFT)

        !.. Write the FT in a single file
        open(newunit = uid_dipoleFT, &
                file = OutDir // DIPOLE_FT_PATH_ALL, &
                form = "formatted", &
                status = "unknown", &
                action = "write", &
                position = "append")
        do iOmega = 1, nOmegas
            call train(iSim)%PrintPulses(uid_dipoleFT)
            w = OmegaVec(iOmega)
            write(uid_dipoleFT, "(i4,*(x,E24.16))") iOmega, w, &
                    ((dble(DipoleFTtotal(iPol, iOmega)), aimag(DipoleFTtotal(iPol, iOmega))), iPol = 1, 3)
        end do
        write(uid_dipoleFT, *)
        close(uid_dipoleFT)

        !.. Save Atomic Charge FT
        open(newunit = uid_AtomicChargeFT, &
                file = OutDir // "/AtomicCharge/AtomicChargeFT" // trim(Simulation_tagv(iSim)), &
                form = "formatted", &
                status = "unknown", &
                action = "write")
        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            write(uid_AtomicChargeFT, "(i4,*(x,E24.16))") iOmega, w, ((&
                    dble (AtomicChargeFT(iOmega, iAtom)), &
                    aimag(AtomicChargeFT(iOmega, iAtom))   ), &
            iAtom = 1, nAtoms)
        end do
        close(uid_AtomicChargeFT)

        !.. Save Atomic Charge FT in a single file
        open(newunit = uid_AtomicChargeFT, &
                file = OutDir // CHARGE_FT_PATH_ALL, &
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

    end do Sim_loop
    close(uid_AtomicChargeALL)


    !.. COMPUTE 2D SPECTRUM DIPOLE
    !
    !.. Load the FT of the dipole, as a function of the time delay
    open(newunit = uid_dipoleFT, &
            file = OutDir // DIPOLE_FT_PATH_ALL, &
            form = "formatted", &
            status = "old", &
            action = "read")
    allocate(DipoleFTwt(3, nOmegas, N_Simulations))
    DipoleFTwt = Z0
    allocate(tvec(N_Simulations))
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
    !
    !.. Regularizes the dipole with respect to the time-delay edges
    !..
    tStep1 = tvec(N_Simulations / 10)
    tStep2 = tvec(9 * N_Simulations / 10)
    StepWidth = (tvec(N_Simulations) - tvec(1)) / 20.d0
    do iSim = 1, N_Simulations
        t = tvec(iSim)
        dstep = NCD_Phi(t, tStep1, StepWidth) * (1.d0 - NCD_Phi(t, tStep2, StepWidth))
        DipoleFTwt(:, :, iSim) = dstep * DipoleFTwt(:, :, iSim)
    enddo
    !
    !.. Compute the FT wrt the time delay
    !..
    allocate(DipoleFTww(3, nOmegas, nTauOmegas))
    DipoleFTww = Z0
    dt = tvec(2) - tvec(1)
    do iOmega = 1, nTauOmegas
        w = TauOmegaVec(iOmega)
        do iSim = 1, N_Simulations
            t = tvec(iSim)
            zexpFact = exp(Zi * w * t)
            DipoleFTww(:, :, iOmega) = DipoleFTww(:, :, iOmega) + zExpFact * DipoleFTwt(:, :, iSim)
        end do
    enddo
    DipoleFTww = DipoleFTww * dt / (2.d0 * PI)
    !
    !.. Write the Bidimensional spectrum to file
    open(newunit = uid_dipoleFT, &
            file = OutDir // "/Dipole/DipoleFT_ww", &
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
    deallocate(tvec)


    !.. COMPUTE 2D SPECTRUM CHARGE
    !
    !.. Load the FT of the charge, as a function of the time delay
    allocate(tvec(N_simulations))
    tvec = 0.d0
    allocate(dAtomFTRe(nAtoms))
    allocate(dAtomFTIm(nAtoms))
    open(newunit = uid_AtomicChargeFT, &
            file = OutDir // CHARGE_FT_PATH_ALL, &
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
    !
    !.. Regularizes the charge with respect to the time-delay edges
    !   *** THIS SEEMS TO WORK ***
    !..
    open(newunit = uid_AtomicChargeFT, &
            file = OutDir // CHARGE_FT_PATH_ALL // "reg", &
            form = "formatted", &
            status = "unknown", &
            action = "write")
    tStep1 = tvec(max(nint(0.2d0 * dble(N_Simulations)), 1))
    tStep2 = tvec(nint(8.d0 * N_Simulations / 10.d0))
    StepWidth = (tvec(N_Simulations) - tvec(1)) / 15.d0
    do iAtom = 1, nAtoms
        do iSim = 1, N_Simulations
            t = tvec(iSim)
            dstep = NCD_Phi(t, tStep1, StepWidth) * (1.d0 - NCD_Phi(t, tStep2, StepWidth))
            ChargeFTwt(:, iSim, iAtom) = dstep * ChargeFTwt(:, iSim, iAtom)
        enddo
    enddo
    do iSim = 1, N_simulations
        t = tvec(iSim)
        do iOmega = 1, nOmegas
            w = OmegaVec(iOmega)
            write(uid_AtomicChargeFT, "(*(x,e24.16))") t, w, &
                    ((dble(ChargeFTwt(iOmega, iSim, iAtom)), &
                    aimag(ChargeFTwt(iOmega, iSim, iAtom))), iAtom = 1, nAtoms)

        enddo
        write(uid_AtomicChargeFT, *)
    enddo
    !
    !.. Compute the FT wrt the time delay
    !   ****** POSSIBLY A MISTAKE IN WHAT FOLLOWS *******
    !..
    allocate(ChargeFTww(nOmegas, nTauOmegas, nAtoms))
    ChargeFTww = Z0
    dt = (tvec(N_simulations) - tvec(1)) / dble(N_simulations - 1)
    do iAtom = 1, nAtoms
        do iOmega = 1, nTauOmegas
            w = TauOmegaVec(iOmega)
            do iSim = 1, N_Simulations
                t = tvec(iSim)
                zexpFact = exp(Zi * w * t)
                !!$           if(iAtom==1)write(*,*) iSim, iOmega, t, w
                ChargeFTww(:, iOmega, iAtom) = ChargeFTww(:, iOmega, iAtom) + &
                        zExpFact * ChargeFTwt(:, iSim, iAtom)
            end do
        enddo
    enddo
    ChargeFTww = ChargeFTww * dt / (2.d0 * PI)
    !
    !.. Write the Bidimensional spectrum to file
    open(newunit = uid_AtomicChargeFT, &
            file = OutDir // "/AtomicCharge/AtomicChargeFT_ww", &
            form = "formatted", &
            status = "unknown", &
            action = "write")
    do iOmegaTau = 1, nTauOmegas
        do iOmega = 1, nOmegas
            write(uid_AtomicChargeFT, "(*(x,e24.16))") OmegaVec(iOmega), TauOmegaVec(iOmegaTau), &
                    ((dble(ChargeFTww(iOmega, iOmegaTau, iAtom)), &
                    (aimag(ChargeFTww(iOmega, iOmegaTau, iAtom)))), iAtom = 1, nAtoms)
        end do
        write(uid_AtomicChargeFT, *)
    enddo
    close(uid_AtomicChargeFT)

    stop

end program ChargeMigration

