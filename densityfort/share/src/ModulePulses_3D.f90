!! CONFIDENTIAL
module ModulePulses_3D

    use ModuleConstants

    implicit none

    private

    !!$  real(kind(1d0)), private, parameter :: CAU   = 137.03599911d0
    !!$  real(kind(1d0)), private, parameter :: ALPHA = 1.d0/CAU
    integer, parameter :: MAXNPULSES = 64

    real(kind(1d0)), external :: NCD_Phi
    real(kind(1d0)), external :: dNCD_Phi

    integer, parameter :: GAUSSIAN = 0
    integer, parameter :: COS_SQUARE = 1
    integer, parameter :: TRANSVERSE_RAMP_UP = 2
    integer, parameter :: LONGITUDINAL_RAMP_UP = 3

    type pulse
        !
        !.. Shape
        !..
        integer :: envelope! gaussian, cos2, transverse ramp-up, longitudinal ramp-up
        !
        !.. Central Time (a.u.)
        !..
        real(kind(1d0)) :: t
        !
        !.. Carrier Frequency (a.u.)
        !..
        real(kind(1d0)) :: o
        !
        !.. Full Width at Half Maximum,
        !   in units of carrier periods
        !..
        real(kind(1d0)) :: f
        !
        !.. Carrier Envelope Phase (degs)
        !..
        real(kind(1d0)) :: d
        !
        !.. Intensity (PW/cm^2)
        !..
        real(kind(1d0)) :: i
        !
        !.. Amplitude (a.u.)
        !..
        real(kind(1d0)) :: a
        !
        !.. Period (a.u.)
        !..
        real(kind(1d0)) :: p
        !
        !.. Direction of the linear-polarization axis
        !   in terms of spherical coordinates in degrees
        real(kind(1d0)) :: theta
        real(kind(1d0)) :: phi
        !
    end type pulse

    type pulse_train
        !
        integer :: n
        type(pulse) :: p(MAXNPULSES)
        !
    contains
        !
        generic, public :: PrintPulses => pulse_train_PrintPulses
        generic, public :: Write => pulse_trainWrite
        generic, public :: WriteFTA => pulse_trainWriteFTA
        generic, public :: GetMaxFieldStrength => pulse_trainGetMaxFieldStrength
        generic, public :: free => FreePulseTrain
        generic, public :: add => AddPulseStrn, AddPulseParams
        generic, public :: A => AField
        generic, public :: E => EField
        generic, public :: FTA => FTAField
        generic, public :: FTE => FTEField
        generic, public :: GetFrequencySpace => pulse_trainGetFrequencySpace
        !
        procedure, private :: pulse_trainWrite
        procedure, private :: pulse_trainWriteFTA
        procedure, private :: pulse_trainGetMaxFieldStrength
        procedure, private :: freePulseTrain
        procedure, private :: AddPulseStrn
        procedure, private :: AddPulseParams
        procedure, private :: AField
        procedure, private :: EField
        procedure, private :: FTAField
        procedure, private :: FTEField
        procedure, private :: pulse_trainGetFrequencySpace
        procedure, private :: pulse_train_PrintPulses
        !
    end type pulse_train

    public :: pulse_train
    public :: Parse_Simulation_File, ExternalVectorPotential, ExternalElectricField
    public :: FTVectorPotential, FTElectricField
    public :: ExternalElectricFieldCart

contains

    subroutine pulse_train_PrintPulses(train, uid)
        class(pulse_train), intent(in) :: train
        integer, intent(in) :: uid
        integer :: iPulse
        write(uid, "(x,i5)", advance = "no") train%n
        do iPulse = 1, train%n
            write(uid, "(*(x,e14.6))", advance = "no") &
                    train%p(iPulse)%t, &
                    train%p(iPulse)%o, &
                    train%p(iPulse)%f, &
                    train%p(iPulse)%d, &
                    train%p(iPulse)%i, &
                    train%p(iPulse)%a, &
                    train%p(iPulse)%p
        enddo
    end subroutine pulse_train_PrintPulses


    real(kind(1d0)) function pulse_trainGetMaxFieldStrength(train, tmin, tmax, dt) result(Strength)
        class(pulse_train), intent(in) :: train
        real(kind(1d0)), intent(in) :: tmin, tmax, dt
        real(kind(1d0)) :: time, dw
        complex(kind(1d0)) :: zw
        integer :: iMu

        Strength = 0.d0
        time = tmin
        do while(time<tmax)
            dw = 0.d0
            do iMu = -1, 1
                zw = ExternalVectorPotential(time, train, iMu)
                dw = dw + abs(zw**2)
            enddo
            Strength = max(Strength, sqrt(dw))
            time = time + dt
        enddo

    end function pulse_trainGetMaxFieldStrength


    subroutine pulse_trainWrite(train, fileName, tmin, tmax, dt)
        class(pulse_train), intent(in) :: train
        character(len = *), intent(in) :: fileName
        real(kind(1d0)), intent(in) :: tmin, tmax, dt
        integer :: uid
        real(kind(1d0)) :: time, Ax, Ay, Az
        complex(kind(1d0)) :: zw_m1, zw_p1, zw_p0

        open(newunit = uid, &
                file = trim(fileName), &
                form = "formatted", &
                status = "unknown")
        time = tmin
        do while(time<tmax)
            zw_m1 = ExternalVectorPotential(time, train, -1)
            zw_p0 = ExternalVectorPotential(time, train, 0)
            zw_p1 = ExternalVectorPotential(time, train, 1)
            Ax = dble (zw_m1 - zw_p1) / sqrt(2.d0)
            Ay = aimag(zw_m1 + zw_p1) / sqrt(2.d0)
            Az = dble (zw_p0)
            write(uid, "(*(x,e22.14))")time, Ax, Ay, Az, &
                    dble(zw_m1), aimag(zw_m1), &
                    dble(zw_p0), aimag(zw_p0), &
                    dble(zw_p1), aimag(zw_p1)
            time = time + dt
        enddo
        close(uid)

    end subroutine pulse_trainWrite


    !.. Returns a list of frequencies where the fourierstransform of the
    !   pulse is of relevant size.
    !   Currently returns a linearly spaced grid.
    !..
    subroutine pulse_trainGetFrequencySpace(Self, w)
        class(pulse_train), intent(in) :: Self
        real(kind(1d0)), allocatable, intent(out) :: w(:)
        integer :: iPulse, Nsigma, nw, iw
        real(kind(1d0)) :: dwi, dw, maxw

        maxw = 0.d0
        dw = huge(1.d0)
        Nsigma = 4

        do iPulse = 1, Self%n
            dwi = Self%p(iPulse)%o / (2 * PI * Self%p(iPulse)%f)
            dw = min(dw, dwi)
            maxw = max(maxw, Self%p(iPulse)%o + Nsigma * dwi)
        enddo

        nw = nint(maxw / dw)
        if(allocated(w)) deallocate(w)
        allocate(w(nw + 1))
        do iw = 0, nw
            w(iw + 1) = iw * dw
        enddo

    end subroutine pulse_trainGetFrequencySpace

    !.. Writes the Fourier transform of the vector potential
    !   to the file "FileName".
    !..
    subroutine pulse_trainWriteFTA(Self, FileName)
        implicit none
        class(pulse_train), intent(inout) :: Self
        character(len = *), intent(in) :: FileName
        real(kind(1d0)), allocatable :: Freq(:)
        integer :: uid, iFreq, mu
        real   (kind(1d0)) :: FTAx, FTAy, FTAz
        complex(kind(1d0)) :: FT(3)

        call Self%GetFrequencySpace(Freq)

        open(newunit = uid, &
                file = FileName, &
                form = "formatted", &
                status = "unknown")

!        write(uid, "(a)") "# w, FT(Ax), FT(Ay), FT(Az), Re(FT(A^-1(w))), Im(FT(A^-1(w))), Re(FT(A^0(w))), ..."
        !
        do iFreq = 1, size(Freq)
            do mu = -1, 1
                FT(mu + 2) = Self%FTA(Freq(iFreq), mu)
            enddo
            FTAx = dble (FT(1) - FT(3)) / sqrt(2.d0)
            FTAy = aimag(FT(1) + FT(3)) / sqrt(2.d0)
            FTAz = dble (FT(2))
            write(uid, "(*(x,e22.14))") Freq(iFreq), FTAx, FTAy, FTAz, &
                    dble(FT(1)), aimag(FT(1)), &
                    dble(FT(2)), aimag(FT(2)), &
                    dble(FT(3)), aimag(FT(3))
        enddo

        close(uid)
        deallocate(Freq)
    end subroutine pulse_trainWriteFTA


    subroutine Parse_Simulation_File(fp, fp_out, NSIM, SIM_TAGS, SIM_PTR)
        !
        implicit none
        !
        integer, intent(in) :: fp, fp_out
        integer, intent(out) :: NSIM
        character(len = *), pointer :: SIM_TAGS(:)
        type(pulse_train), pointer :: SIM_PTR(:)
        !
        integer :: i, i1, i2, j, stat
        type(pulse_train), pointer :: local_trains(:)
        character(len = 1000000) :: line
        character(len = 1000000) :: text
        integer :: itrain, ntrains
        character(len = 256), allocatable :: train_label(:)
        character(len = 10000), allocatable :: train_strn(:)
        character(len = 10000) :: tmp_train_strn
        !
        !.. Copies the file on a single line,
        !   eliminating comments and indentations.
        !..
        text = " "
        do
            line = " "
            read(fp, "(a)", IOSTAT = stat)line
            if(stat/=0)exit
            i = index(line, "#")
            if(i>0)line(i:) = " "
            if(len_trim(line)==0)cycle
            text = trim(text) // adjustl(line)
        enddo
        !
        !write(*,"(a)") "text :"//trim(text)
        !
        !.. Counts the number of trains
        !..
        ntrains = 0
        j = 0
        i = 0
        do
            i = index(text(j + 1:), "[")
            if(i<=0)exit
            ntrains = ntrains + 1
            j = j + i
        enddo
        !
        !
        write(fp_out, "(a,i4)") "Number of trains:", ntrains
        !
        !.. Create a library of trains
        !..
        allocate(local_trains(ntrains))
        allocate(train_label(ntrains))
        allocate(train_strn(ntrains))
        !
        !.. Parse the simulations
        !..
        i2 = 0
        do itrain = 1, ntrains
            !
            i1 = index(text(i2 + 1:), "[") + i2
            i2 = index(text(i1 + 1:), "]") + i1
            !
            train_label(itrain) = adjustl(text(i1 + 1:i2 - 1))
            !
            i1 = index(text(i2 + 1:), "{") + i2
            i2 = index(text(i1 + 1:), "}") + i1
            train_strn(itrain) = adjustl(text(i1 + 1:i2 - 1))
            !
            write(fp_out, "(a,x,i4,x,a,x,a)") "A", itrain, &
                    trim(train_label(itrain)), trim(train_strn(itrain))
            !
        enddo
        !
        !
        !.. Replace the references to other trains with
        !   the actual list of pulses those trains contain
        !..
        do itrain = 1, ntrains
            !
            i1 = 0
            inner : do
                !
                i2 = index(train_strn(itrain)(i1 + 1:), ";") + i1
                !
                !.. There are no more train specifications
                !..
                if(i2 - i1<=0)exit inner
                !
                if(index(train_strn(itrain)(i1 + 1:i2), "(")<=0)then
                    !
                    !.. Train Label Found
                    !..
                    tmp_train_strn = train_strn(itrain)(i2 + 1:)
                    do i = 1, itrain - 1
                        !
                        if(trim(adjustl(train_strn(itrain)(i1 + 1:i2 - 1)))==&
                                trim(adjustl(train_label(i))))then
                            !
                            train_strn(itrain)(i1 + 1:) = adjustl(train_strn(i))
                            train_strn(itrain) = trim(train_strn(itrain)) // &
                                    adjustl(tmp_train_strn)
                            !
                            i1 = i1 + len_trim(adjustl(train_strn(i)))
                            exit
                            !
                        endif
                        !
                    enddo
                    !
                else
                    !
                    i1 = i2
                    !
                endif
                !
            enddo inner
            !
            write(fp_out, "(a,x,i4,x,a)")"B", itrain, trim(train_strn(itrain))
            !
        enddo
        !
        !
        !.. Now find which trains should actually be used
        !   for the set of simulations
        !..
        i2 = index(text, "EXECUTE")
        i1 = index(text(i2:), "{") + i2 - 1
        i2 = index(text(i1:), "}") + i1 - 1
        i = i1
        NSIM = 0
        do
            j = index(text(i + 1:i2 - 1), ";") + i
            if(j - i<=0)exit
            i = j
            NSIM = NSIM + 1
        enddo
        !
        write(fp_out, "(a,i4)") "NSIM=", NSIM
        !
        allocate(SIM_TAGS(NSIM))
        allocate(SIM_PTR(NSIM))
        i = i1
        NSIM = 0
        do
            j = index(text(i + 1:i2 - 1), ";") + i
            if(j - i<=0)exit
            i = j
            NSIM = NSIM + 1
            SIM_TAGS(NSIM) = adjustl(text(i1 + 1:i - 1))
            i1 = i
        enddo
        !
        !
        do i = 1, NSIM
            !
            do j = 1, ntrains
                !
                if(trim(SIM_TAGS(i))==trim(train_label(j)))then
                    !
                    write(fp_out, "(i4,x,a)")i, trim(SIM_TAGS(i)) // " " // trim(train_label(j))
                    !
                    call Parse_Train(train_strn(j), SIM_PTR(i))
                    !
                    exit
                    !
                endif
                !
            enddo
            !
        enddo
        !
        !
        return
        !
        !
    end subroutine Parse_Simulation_File


    subroutine Parse_Train(strn, train)
        !
        !.. Extracts the parameters of the train
        !   from a string with the format
        !
        !   (<pulse_strn_1>); (<pulse_strn_1>); ...; (<pulse_strn_N>);
        !
        !..
        character(len = *), intent(in) :: strn
        type(pulse_train), intent(out) :: train
        !
        integer :: i1, i2
        character(len = 10000) :: line
        !
        train%n = 0
        i2 = 0
        do
            !
            i1 = index(strn(i2 + 1:), "(") + i2
            if(i1 - i2<=0)exit
            train%n = train%n + 1
            !
            i2 = index(strn(i1 + 1:), ")") + i1
            line = adjustl(strn(i1 + 1:i2 - 1))
            call Parse_Pulse_Strn(line, train%p(train%n))
            i2 = index(strn(i1 + 1:), ";") + i1
            !
        enddo
        !
        !
        return
        !
    end subroutine Parse_Train


    subroutine AddPulseParams(train, &
            Envelope, t0, omega0, nt, cep, I0)
        !
        class(pulse_train), intent(inout) :: train
        character, intent(in) :: Envelope
        real(kind(1d0)), intent(in) :: t0, omega0, nt, cep, I0
        !
        real(kind(1d0)), parameter :: PI = 3.1415926535897932d0
        !
        integer :: i, n
        character(len = 10000) :: line
        !
        if(train%n>=MAXNPULSES)return
        train%n = train%n + 1
        n = train%n
        !
        if(Envelope=="G")train%p(n)%envelope = GAUSSIAN
        if(Envelope=="C")train%p(n)%envelope = COS_SQUARE
        if(Envelope=="R")train%p(n)%envelope = LONGITUDINAL_RAMP_UP
        if(Envelope=="l")train%p(n)%envelope = TRANSVERSE_RAMP_UP
        train%p(n)%t = t0
        train%p(n)%o = omega0
        train%p(n)%f = nt
        train%p(n)%d = cep
        train%p(n)%d = train%p(n)%d * PI / 180.d0
        train%p(n)%i = I0
        line = adjustl(line(i + 1:))
        train%p(n)%p = 2.d0 * PI / train%p(n)%o
        train%p(n)%a = CAU / 18.73d0 / train%p(n)%o * &
                sqrt(10.d0 * train%p(n)%i)
        !
    end subroutine AddPulseParams


    subroutine AddPulseStrn(train, strn)
        !
        !.. Extracts from a string, e.g.
        !   C  0.0 0.057 10 0 0.01
        !   the type of pulse (C/G for cos2 or gaussian,
        !   respectively), its central time (in a.u.),
        !   its energy (in a.u.), its full width at half
        !   maximum in terms of laser periods (pure number),
        !   the carrier-envelope phase (degs), and
        !   the intensity (PW/cm^2)
        !..
        class(pulse_train), intent(inout) :: train
        character(len = *), intent(in) :: strn
        !
        real(kind(1d0)), parameter :: PI = 3.1415926535897932d0
        !
        integer :: i, n
        character(len = 10000) :: line
        logical :: binary
        !
        if(train%n>=MAXNPULSES)return
        train%n = train%n + 1
        n = train%n
        !
        line = adjustl(strn)
        !
        if(line(1:1)=="G")train%p(n)%envelope = GAUSSIAN
        if(line(1:1)=="C")train%p(n)%envelope = COS_SQUARE
        if(line(1:1)=="R")train%p(n)%envelope = LONGITUDINAL_RAMP_UP
        if(line(1:1)=="l")train%p(n)%envelope = TRANSVERSE_RAMP_UP
        binary = .FALSE.
        if(line(1:1)=="b")then
            train%p(n)%envelope = COS_SQUARE
            binary = .TRUE.
        endif
        line = adjustl(line(2:))
        i = index(line, " ")
        if(binary)then
            !
            !.. Convert the time expressed
            !   in femtoseconds and in binary
            !   system to atomic units in
            !   decimal system
            !..
            train%p(n)%t = binfsstrn2decau(line(1:i))
            !
        else
            read(line(1:i), *)train%p(n)%t
        endif
        line = adjustl(line(i + 1:))
        i = index(line, " ")
        read(line(1:i), *)train%p(n)%o
        line = adjustl(line(i + 1:))
        i = index(line, " ")
        read(line(1:i), *)train%p(n)%f
        line = adjustl(line(i + 1:))
        i = index(line, " ")
        read(line(1:i), *)train%p(n)%d
        train%p(n)%d = train%p(n)%d * PI / 180.d0
        line = adjustl(line(i + 1:))
        i = index(line, " ")
        read(line(1:i), *)train%p(n)%i
        line = adjustl(line(i + 1:))
        train%p(n)%p = 2.d0 * PI / train%p(n)%o
        train%p(n)%a = CAU / 18.73d0 / train%p(n)%o * &
                sqrt(10.d0 * train%p(n)%i)
        !
        return
        !
    end subroutine AddPulseStrn

    subroutine FreePulseTrain(train)
        class(pulse_train), intent(inout) :: train
        train%n = 0
    end subroutine FreePulseTrain

    real(kind(1d0)) function AField(train, Time) result(A)
        class(pulse_train), intent(inout) :: train
        real(kind(1d0)), intent(in) :: Time
        A = ExternalVectorPotential(Time, train)
    end function AField

    real(kind(1d0)) function EField(train, Time) result(E)
        class(pulse_train), intent(inout) :: train
        real(kind(1d0)), intent(in) :: Time
        E = ExternalElectricField(Time, train)
    end function EField

    complex(kind(1d0)) function FTAField(train, omega, mu) result(zFTA)
        class(pulse_train), intent(inout) :: train
        real(kind(1d0)), intent(in) :: Omega
        integer, optional, intent(in) :: mu
        zFTA = FTVectorPotential(omega, train, mu)
    end function FTAField

    complex(kind(1d0)) function FTEField(train, omega, mu) result(zFTE)
        class(pulse_train), intent(inout) :: train
        real(kind(1d0)), intent(in) :: Omega
        integer, optional, intent(in) :: mu
        zFTE = FTElectricField(omega, train, mu)
    end function FTEField

    subroutine Parse_Pulse_Strn(strn, OutPulse)
        !
        !.. Extracts from a string, e.g.
        !   C  0.0 0.057 10 0 0.01
        !   the type of pulse (C/G for cos2 or gaussian,
        !   respectively), its central time (in a.u.),
        !   its energy (in a.u.), its full width at half
        !   maximum in terms of laser periods (pure number),
        !   the carrier-envelope phase (degs), and
        !   the intensity (PW/cm^2)
        !..
        character(len = *), intent(in) :: strn
        type(pulse), intent(out) :: OutPulse
        !
        real(kind(1d0)), parameter :: PI = 3.1415926535897932d0
        !
        integer :: i
        character(len = 10000) :: line
        logical :: binary
        !
        !.. With ifort, the reading of the string would
        !   be completely straightforward. Unfortunately,
        !   xlf90 seems to pose more problems, therefore
        !   I'll use a more involved algorithm to parse
        !   the string
        !
        !   Longitudinal ramp-up is a zero-frequency ramp-up
        !   of the external electric field.
        !
        !   Transverse ramp-up is a finite-frequency ramp-up
        !   of the external vector potential.
        !..
        line = adjustl(strn)
        !
        if(line(1:1)=="G")OutPulse%envelope = GAUSSIAN
        if(line(1:1)=="C")OutPulse%envelope = COS_SQUARE
        if(line(1:1)=="R")OutPulse%envelope = LONGITUDINAL_RAMP_UP
        if(line(1:1)=="l")OutPulse%envelope = TRANSVERSE_RAMP_UP
        binary = .FALSE.
        if(line(1:1)=="b")then
            OutPulse%envelope = COS_SQUARE
            binary = .TRUE.
        endif
        line = adjustl(line(2:))
        i = index(line, " ")
        if(binary)then
            !
            !.. Convert the time expressed
            !   in femtoseconds and in binary
            !   system to atomic units in
            !   decimal system
            !..
            OutPulse%t = binfsstrn2decau(line(1:i))
            !
        else
            read(line(1:i), *)OutPulse%t
        endif
        line = adjustl(line(i + 1:))
        i = index(line, " ")
        read(line(1:i), *)OutPulse%o
        line = adjustl(line(i + 1:))
        i = index(line, " ")
        read(line(1:i), *)OutPulse%f
        line = adjustl(line(i + 1:))
        i = index(line, " ")
        read(line(1:i), *)OutPulse%d
        OutPulse%d = OutPulse%d * PI / 180.d0
        line = adjustl(line(i + 1:))
        i = index(line, " ")
        read(line(1:i), *)OutPulse%i
        line = adjustl(line(i + 1:))
        OutPulse%p = 2.d0 * PI / OutPulse%o
        OutPulse%a = CAU / 18.73d0 / OutPulse%o * &
                sqrt(10.d0 * OutPulse%i)
        OutPulse%Theta = 0.d0
        OutPulse%Phi = 0.d0
        if(len_trim(line)/=0)then
            i = index(line, " ")
            read(line(1:i), *)OutPulse%Theta
            OutPulse%Theta = OutPulse%Theta / 180.d0 * PI
            line = adjustl(line(i + 1:))
            i = index(line, " ")
            read(line(1:i), *)OutPulse%Phi
            OutPulse%Phi = OutPulse%Phi / 180.d0 * PI
            line = adjustl(line(i + 1:))
        endif
        !
        return
        !
    end subroutine Parse_Pulse_Strn

    real(kind(1d0)) function binfsstrn2decau(nstrn) result(res)
        !
        ! BINary FemtoSecond STRiNg TO DECimal Atomic Units
        !
        implicit none
        character(len = *), intent(in) :: nstrn
        character(len = 100) :: lstrn, lstrn1, lstrn2
        !
        real(kind(1d0)) :: w
        integer :: comma, i
        !
        lstrn = " "
        lstrn = adjustl(nstrn)
        do i = len(lstrn), 1, -1
            if(lstrn(i:i)/="0")exit
            lstrn(i:i) = " "
        enddo
        do i = 1, len_trim(lstrn)
            if(lstrn(i:i)/="0")exit
            lstrn(i:i) = " "
        enddo
        comma = index(nstrn, ".")
        lstrn1 = adjustl(lstrn(1:comma - 1))
        lstrn2 = adjustl(lstrn(comma + 1:))
        !
        w = 1.d0
        res = 0.d0
        do i = len_trim(lstrn1), 1, -1
            if(lstrn(i:i)=="1")res = res + w
            w = 2.d0 * w
        enddo
        w = 1.d0
        do i = 1, len_trim(lstrn2)
            w = w / 2.d0
            if(lstrn(i:i)=="1")res = res + w
        enddo
        !
        return
        !
    end function binfsstrn2decau


    real(kind(1d0)) function ExternalVectorPotentialAbsVal(time, train) result(res)
        real(kind(1d0)), intent(in) :: time
        type(pulse_train), intent(in) :: train
        integer :: imu
        complex(kind(1d0)) :: zw
        res = 0.d0
        do imu = -1, 1
            zw = ExternalVectorPotential(time, train, imu)
            res = res + abs(zw)**2
        enddo
        res = sqrt(res)
    end function ExternalVectorPotentialAbsVal

    !.. Disregards contributions from possible quasi-static external fields
    !   (i.e., not due to a laser pulse, but to a charged/charging condenser)
    !   Returns the spherical contravariant vector component A^\mu of the vector potential
    !..
    complex(kind(1d0)) function ExternalVectorPotential(time, train, mu_in) result(res)
        !
        use, intrinsic :: ISO_FORTRAN_ENV
        implicit none
        real(kind(1d0)), intent(in) :: time
        type(pulse_train), intent(in) :: train
        integer, optional, intent(in) :: mu_in
        !
        integer :: i, mu
        real   (kind(1d0)) :: x, w1, wc, tau, t0, A0, w0, phi, pulseA
        complex(kind(1d0)) :: geometrical_factor
        !
        real(kind(1d0)), parameter :: SIGMA_THRESHOLD = 10
        res = 0.d0
        mu = 0
        if(present(mu_in))mu = mu_in
        do i = 1, train%n

            pulseA = 0.d0

            x = train%p(i)%o * (time - train%p(i)%t)
            if(train%p(i)%envelope==GAUSSIAN)then
                !.. For an individual pulse,
                !
                !   A(t) = A_0 * \exp\left[ - \ln 2 \left(\frac{\omega (t-t_i)}{\pi F}\right)^2 \right] *
                !                \cos( \omega(t-t_i) + \delta )
                !
                w1 = x / (PI * train%p(i)%f)
                !
                if (abs(w1) < SIGMA_THRESHOLD) then !HOTFIX for residual noise in pulse(Ruben-8th-2021)
                    pulseA = train%p(i)%a * exp(-LN2 * w1 * w1) * cos(x + train%p(i)%d)
                end if
                !
            elseif(train%p(i)%envelope==COS_SQUARE)then
                w1 = x / (4 * train%p(i)%f)
                if(abs(w1)<0.5d0 * PI)then
                    wc = cos(w1)
                    pulseA = train%p(i)%a * wc * wc * cos(x + train%p(i)%d)
                endif
            elseif(train%p(i)%envelope==TRANSVERSE_RAMP_UP)then
                !
                w0 = train%p(i)%o
                t0 = train%p(i)%t
                A0 = train%p(i)%a
                phi = train%p(i)%d
                tau = 2 * PI / w0 * train%p(i)%f
                !
                pulseA = A0 * NCD_Phi(time, t0, tau) * cos(w0 * (time - t0) + phi)
                !
            endif

            select case(mu)
            case(0)
                geometrical_factor = cos(train%p(i)%theta)
            case(+1)
                geometrical_factor = - sin(train%p(i)%theta) / sqrt(2.d0) * exp(-Zi * train%p(i)%phi)
            case(-1)
                geometrical_factor = sin(train%p(i)%theta) / sqrt(2.d0) * exp(Zi * train%p(i)%phi)
            case default
                write(ERROR_UNIT, "(a)") "Invalid magnetic quantum number for vector potential"
                stop
            end select

            res = res + pulseA * geometrical_factor

        enddo
        return
        !
    end function ExternalVectorPotential


    !.. Disregards contributions from possible quasi-static external fields
    !   (i.e., not due to a laser pulse, but to a charged/charging condenser)
    !   WARNING: Does not compute the fourier transform of mu!=0 components correctly!
    !..
    complex(kind(1d0)) function FTVectorPotential(omega, train, mu) result(zres)
        !
        use, intrinsic :: ISO_FORTRAN_ENV
        implicit none
        real   (kind(1d0)), intent(in) :: omega
        type(pulse_train), intent(in) :: train
        integer, optional, intent(in) :: mu
        !
        complex(kind(1d0)), parameter :: Z0 = (0.d0, 0.d0)
        complex(kind(1d0)), parameter :: Z1 = (1.d0, 0.d0)
        complex(kind(1d0)), parameter :: Zi = (0.d0, 1.d0)
        real   (kind(1d0)), parameter :: LN2 = 0.6931471805599453d0
        real   (kind(1d0)), parameter :: PI = 3.1415926535897932d0
        real   (kind(1d0)), parameter :: MAX_POWER = 40.d0
        real   (kind(1d0)) :: t0, A0, w0, f0, d0
        complex(kind(1d0)) :: zw0, zwp, zwm, geometrical_factor, pulseFT
        integer :: i, mu_

        mu_ = 0
        if(present(mu)) mu_ = mu
        zres = Z0
        do i = 1, train%n
            !
            w0 = train%p(i)%o
            t0 = train%p(i)%t
            A0 = train%p(i)%a
            f0 = train%p(i)%f
            d0 = train%p(i)%d
            pulseFT = Z0
            !
            if(train%p(i)%envelope==GAUSSIAN)then
                !
                !.. For an individual pulse,
                !
                !   A_i(t) = A_i * \exp\left[ - \ln 2 \left(\frac{\omega_i (t-t_i)}{\pi F_i}\right)^2 \right]
                !                \cos( \omega(t-t_i) + \delta_i )
                !
                !
                !   A_i(w) = \frac{A_i}{2} \frac{\pi F_i}{\sqrt{2 \ln 2} \omega_i} \exp( i \omega t_i )
                !            {
                !               \exp[ - \frac{(\omega-\omega_i)^2}{4\ln2 \omega_i^2}\pi^2 F_i^2 + i \delta_i ] -
                !               \exp[ - \frac{(\omega+\omega_i)^2}{4\ln2 \omega_i^2}\pi^2 F_i^2 - i \delta_i ]
                !            }
                !
                zw0 = 0.5d0 * A0 * PI * f0 / w0 / sqrt(2.d0 * LN2) * exp(Zi * omega * t0)
                zwp = exp(- Z1 * (1.d0 - omega / w0)**2 / (4.d0 * LN2) * (PI * f0)**2 + Zi * d0)
                zwm = exp(- Z1 * (1.d0 + omega / w0)**2 / (4.d0 * LN2) * (PI * f0)**2 - Zi * d0)
                pulseFT = zw0 * (zwp - zwm)
                !
                !
            elseif(train%p(i)%envelope==COS_SQUARE)then
                !
                zw0 = A0 * 0.5d0 / sqrt(2.d0 * PI) * exp(Zi * (t0 * omega - d0)) * w0**2
                zwp = sin(2.d0 * PI * f0 * (1.d0 - omega / w0)) / &
                        ((omega - w0) * (2.d0 * f0 * (omega - w0) - w0) * (2.d0 * f0 * (omega - w0) + w0))
                zwm = sin(2.d0 * PI * f0 * (1.d0 + omega / w0)) / &
                        ((omega + w0) * (2.d0 * f0 * (omega + w0) - w0) * (2.d0 * f0 * (omega + w0) + w0))
                zwm = zwm * exp(2.d0 * Zi * d0)
                !
                pulseFT = zw0 * (zwp - zwm)
                !
            elseif(train%p(i)%envelope==TRANSVERSE_RAMP_UP)then
                !
                !pulseFT=Z0
                !
            endif

            if(abs(pulseFT)**2 < epsilon(0.d0)) cycle

            select case(mu_)
            case(0)
                geometrical_factor = cos(train%p(i)%theta)
            case(+1)
                geometrical_factor = - sin(train%p(i)%theta) / sqrt(2.d0) * exp(-Zi * train%p(i)%phi)
            case(-1)
                geometrical_factor = sin(train%p(i)%theta) / sqrt(2.d0) * exp(Zi * train%p(i)%phi)
            case default
                write(ERROR_UNIT, "(a)") "Invalid magnetic quantum number for the FT of the vector potential"
                stop
            end select
            !!$      write(*,"(e22.14)",advance="no") atan2(aimag(zres),dble(zres))
            zres = zres + pulseFT * geometrical_factor
            !!$       write(*,"(3(x,e22.14))") atan2(aimag(geometrical_factor),dble(geometrical_factor)),&
            !!$            atan2(aimag(pulseFT),dble(pulseFT)), atan2(aimag(zres),dble(zres))
            !!$       if(atan2(aimag(zres),dble(zres))>1.d0.and.atan2(aimag(zres),dble(zres))<3.d0.and.&
            !!$            abs(zres)>1d-14)  write(*,*) "Entire zres:",zres
            !!$       if(mu_ == +1) write(*,*) "m=+1",zres
            !!$       if(mu_ == -1) write(*,*) "m=-1",zres

            !!$       write(*,"(i2,6(x,e22.14))") mu_, geometrical_factor, pulseFT, zres

        enddo
        return
        !
    end function FTVectorPotential


    complex(kind(1d0)) function FTElectricField(omega, train, mu) result(zres)
        implicit none
        real   (kind(1d0)), intent(in) :: omega
        type(pulse_train), intent(in) :: train
        integer, optional, intent(in) :: mu
        zres = - (0.d0, 1.d0) * omega * ALPHA * FTVectorPotential(omega, train, mu)
    end function FTElectricField


    !.. Returns the amplitude of the electric field computed as
    !   $$
    !      E(t) = - ALPHA * dA(t)/dt + Quasi-static non transversal
    !             external fields
    !   $$
    complex(kind(1d0)) function ExternalElectricField(time, train, mu_in) result(zres)
        !..
        use, intrinsic :: ISO_FORTRAN_ENV
        implicit none
        real(kind(1d0)), intent(in) :: time
        type(pulse_train), intent(in) :: train
        integer, optional, intent(in) :: mu_in
        real(kind(1d0)) :: x, w1, wc, ws, A0, w0, t0, phi, tau, pulseE
        integer :: iPulse, mu
        complex(kind(1d0)) :: geometrical_factor
        !
        real(kind(1d0)), parameter :: SIGMA_THRESHOLD = 10

        mu = 0
        if(present(mu_in)) mu = mu_in

        zres = Z0
        do iPulse = 1, train%n

            pulseE = 0.d0
            !
            x = train%p(iPulse)%o * (time - train%p(iPulse)%t)
            if(train%p(iPulse)%envelope==GAUSSIAN)then
                !
                w1 = x / (PI * train%p(iPulse)%f)
                if (abs(w1) < SIGMA_THRESHOLD) then !$$Luca's Fix to noise at the end (Ruben aug-8th-2021)
                    pulseE = ALPHA * train%p(iPulse)%a * train%p(iPulse)%o * exp(-LN2 * w1 * w1) * &
                            (2.d0 * LN2 / (PI * train%p(iPulse)%f)**2 * x * cos(x + train%p(iPulse)%d) + &
                                    sin(x + train%p(iPulse)%d))
                end if
                !
            elseif(train%p(iPulse)%envelope==COS_SQUARE)then
                !
                w1 = x / (4 * train%p(iPulse)%f)
                if(abs(w1)<0.5d0 * PI)then
                    wc = cos(w1)
                    ws = sin(w1)
                    pulseE = ALPHA * train%p(iPulse)%a * train%p(iPulse)%o * (&
                            wc * ws * cos(x + train%p(iPulse)%d) * 0.5d0 / train%p(iPulse)%f + &
                                    wc * wc * sin(x + train%p(iPulse)%d))
                endif
                !
            elseif(train%p(iPulse)%envelope==TRANSVERSE_RAMP_UP)then
                !
                w0 = train%p(iPulse)%o
                t0 = train%p(iPulse)%t
                A0 = train%p(iPulse)%a
                phi = train%p(iPulse)%d
                tau = 2 * PI / w0 * train%p(iPulse)%f
                !
                pulseE = - ALPHA * A0 * dNCD_Phi(time, t0, tau) * cos(w0 * (time - t0) + phi)
                pulseE = pulseE + ALPHA * A0 * w0 * NCD_Phi(time, t0, tau) * sin(w0 * (time - t0) + phi)
                !
            elseif(train%p(iPulse)%envelope==LONGITUDINAL_RAMP_UP)then
                !
                w0 = train%p(iPulse)%o
                t0 = train%p(iPulse)%t
                A0 = train%p(iPulse)%a
                phi = train%p(iPulse)%d
                tau = 2 * PI / w0 * train%p(iPulse)%f
                !
                pulseE = ALPHA * A0 * NCD_Phi(time, t0, tau)
                !
            endif
            !
            select case(mu)
            case(0)
                geometrical_factor = cos(train%p(iPulse)%theta)
            case(+1)
                geometrical_factor = - sin(train%p(iPulse)%theta) / sqrt(2.d0) * exp(-Zi * train%p(iPulse)%phi)
            case(-1)
                geometrical_factor = sin(train%p(iPulse)%theta) / sqrt(2.d0) * exp(Zi * train%p(iPulse)%phi)
            case default
                write(ERROR_UNIT, "(a)") "Invalid magnetic quantum number for the electric field"
                stop
            end select

            zres = zres + pulseE * geometrical_factor
            !
        enddo
        return
        !
    end function ExternalElectricField


    !.. Returns the amplitude of the electric field computed as
    !   $$
    !      E(t) = - ALPHA * dA(t)/dt + Quasi-static non transversal
    !             external fields
    !   $$
    real(kind(1d0)) function ExternalElectricFieldCart(time, train, iCoord) result(res)
        !..
        use, intrinsic :: ISO_FORTRAN_ENV
        implicit none
        real(kind(1d0)), intent(in) :: time
        type(pulse_train), intent(in) :: train
        integer, intent(in) :: iCoord
        real(kind(1d0)) :: x, w1, wc, ws, A0, w0, t0, phi, tau, pulseE
        integer :: iPulse, mu
        real(kind(1d0)) :: geometrical_factor

        res = 0.d0
        do iPulse = 1, train%n

            pulseE = 0.d0
            !
            x = train%p(iPulse)%o * (time - train%p(iPulse)%t)
            if(train%p(iPulse)%envelope==GAUSSIAN)then
                !
                w1 = x / (PI * train%p(iPulse)%f)
                pulseE = ALPHA * train%p(iPulse)%a * train%p(iPulse)%o * exp(-LN2 * w1 * w1) * &
                        (2.d0 * LN2 / (PI * train%p(iPulse)%f)**2 * x * cos(x + train%p(iPulse)%d) + &
                                sin(x + train%p(iPulse)%d))
                !
            elseif(train%p(iPulse)%envelope==COS_SQUARE)then
                !
                w1 = x / (4 * train%p(iPulse)%f)
                if(abs(w1)<0.5d0 * PI)then
                    wc = cos(w1)
                    ws = sin(w1)
                    pulseE = ALPHA * train%p(iPulse)%a * train%p(iPulse)%o * (&
                            wc * ws * cos(x + train%p(iPulse)%d) * 0.5d0 / train%p(iPulse)%f + &
                                    wc * wc * sin(x + train%p(iPulse)%d))
                endif
                !
            elseif(train%p(iPulse)%envelope==TRANSVERSE_RAMP_UP)then
                !
                w0 = train%p(iPulse)%o
                t0 = train%p(iPulse)%t
                A0 = train%p(iPulse)%a
                phi = train%p(iPulse)%d
                tau = 2 * PI / w0 * train%p(iPulse)%f
                !
                pulseE = - ALPHA * A0 * dNCD_Phi(time, t0, tau) * cos(w0 * (time - t0) + phi)
                pulseE = pulseE + ALPHA * A0 * w0 * NCD_Phi(time, t0, tau) * sin(w0 * (time - t0) + phi)
                !
            elseif(train%p(iPulse)%envelope==LONGITUDINAL_RAMP_UP)then
                !
                w0 = train%p(iPulse)%o
                t0 = train%p(iPulse)%t
                A0 = train%p(iPulse)%a
                phi = train%p(iPulse)%d
                tau = 2 * PI / w0 * train%p(iPulse)%f
                !
                pulseE = ALPHA * A0 * NCD_Phi(time, t0, tau)
                !
            endif
            !
            select case(iCoord)
            case(1)
                geometrical_factor = sin(train%p(iPulse)%theta) * cos(train%p(iPulse)%phi)
            case(2)
                geometrical_factor = sin(train%p(iPulse)%theta) * sin(train%p(iPulse)%phi)
            case(3)
                geometrical_factor = cos(train%p(iPulse)%theta)
            case default
                write(ERROR_UNIT, "(a)") "Invalid cartesian coordinate"
                stop
            end select

            res = res + pulseE * geometrical_factor
            !
        enddo
        return
        !
    end function ExternalElectricFieldCart

    real(kind(1d0)) function ExternalElectricFieldAbsVal(time, train) result(res)
        implicit none
        type(pulse_train), intent(in) :: train
        real(kind(1d0)), intent(in) :: time
        integer :: imu
        complex(kind(1d0)) :: zfield

        res = 0.d0
        do imu = -1, 1
            zfield = ExternalElectricField(time, train, imu)
            res = res + abs(zfield)**2
        enddo
        res = sqrt(res)
    end function ExternalElectricFieldAbsVal


end module ModulePulses_3D
