SUBROUTINE READIN( Freq, MaxMedia, NMedia, &
     TopOpt, CPT, CST, rhoT, BumDen, eta, xi, NG, sigma, Depth, &
     BotOpt, CPB, CSB, rhoB, nc, ssp, sspHS )

  ! Reads in the info in ENVFIL

  IMPLICIT NONE
  INTEGER, PARAMETER :: MaxSSP = 2001
  INTEGER            :: NMedia, NG( * ), MaxMedia, NElts, Medium, Nneeded, nc
  REAL    (KIND=8)   :: alphaR, betaR, alphaI, betaI, rhoR, Freq, rhoT, rhoB, &
                        sigma( * ), Depth( * ), BumDen, eta, xi, rho( 1 ), C, deltaz, ssp(nc,6), sspHS(2,6)
  COMPLEX (KIND=8)   :: CPT, CST, CPB, CSB, CP( MaxSSP ), CS( MaxSSP )
  CHARACTER          :: TopOpt*( * ), BotOpt*( * ), SSPType*1, AttenUnit*2, BCType*1, Task*8
  INTEGER            :: IWRITE
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  IWRITE=0

  alphaR = 1500.0
  betaR  = 0.0
  rhoR   = 1.0
  alphaI = 0.0
  betaI  = 0.0
  NElts  = 0         ! this is a dummy variable, passed to profil during read of SSP

  if (IWRITE == 1) WRITE( 6, "( ' Frequency = ', G11.4, 'Hz   NMedia = ', I3, // )" ) Freq, NMedia

  IF ( NMedia > MaxMedia ) THEN
     if (IWRITE == 1) WRITE( 6, * ) 'MaxMedia = ', MaxMedia
     CALL ERROUT( 6, 'F', 'READIN', 'Too many Media' )
  ENDIF

  ! TOP OPTIONS
  SSPType         = TopOpt(1:1)
  BCType          = TopOpt(2:2)
  AttenUnit       = TopOpt(3:4)

  if (IWRITE == 1) then
    SELECT CASE ( SSPType )        ! SSP approximation options
    CASE ( 'N' )
       WRITE( 6, * ) '    N2-LINEAR approximation to SSP'
    CASE ( 'C' )
       WRITE( 6, * ) '    C-LINEAR approximation to SSP'
    CASE ( 'S' )
       WRITE( 6, * ) '    SPLINE approximation to SSP'
    CASE ( 'A' )
       WRITE( 6, * ) '    ANALYTIC SSP option'
    CASE DEFAULT
       CALL ERROUT( 6, 'F', 'READIN', 'Unknown option for SSP approximation' )
    END SELECT
  end if

  if (IWRITE == 1) then 
    SELECT CASE ( AttenUnit(1:1) ) ! Attenuation options
    CASE ( 'N' )
       WRITE( 6, * ) '    Attenuation units: nepers/m'
    CASE ( 'F' )
       WRITE( 6, * ) '    Attenuation units: dB/mkHz'
    CASE ( 'M' ) 
       WRITE( 6, * ) '    Attenuation units: dB/m'
    CASE ( 'W' )
       WRITE( 6, * ) '    Attenuation units: dB/wavelength'
    CASE ( 'Q' )
       WRITE( 6, * ) '    Attenuation units: Q'
    CASE ( 'L' )
       WRITE( 6, * ) '    Attenuation units: Loss parameter'
    CASE DEFAULT
       CALL ERROUT( 6, 'F', 'READIN', 'Unknown attenuation units' )
    END SELECT
  end if

  if (IWRITE == 1) then
    SELECT CASE ( AttenUnit(2:2) ) !  Added volume attenuation
    CASE ( 'T' )
       WRITE( 6, * ) '    THORP attenuation added'
    END SELECT
  end if

  ! CALL TOPBOT  to read top BC 
  IF ( BCType == 'A' .and. IWRITE==1) &
        WRITE( 6, "( //, '   Z (m)     alphaR (m/s)     betaR   rho (g/cm^3)  alphaI     betaI', / )" )

  !  Internal media 
  IF ( BCType /= 'A' .and. IWRITE==1) &
        WRITE( 6, "( //, '   Z (m)     alphaR (m/s)     betaR   rho (g/cm^3)  alphaI     betaI', / )" )

  CALL TOPBOT( Freq, BCType, AttenUnit, CPT, CST, rhoT, BumDen, eta, xi, sspHS(1,1:6), nc )

  DO Medium = 1, NMedia
     if (IWRITE == 1) WRITE( 6, "( /, '          ( Number of pts = ', I5, '  RMS roughness = ', G10.3, ' )')" ) &
          NG( Medium ),            sigma( Medium )

     !  Call PROFIL to read in SSP 
     Task = 'INIT'
     CALL PROFIL( Depth, CP, CS, rho, Medium, NElts, Freq, SSPType, AttenUnit, Task, nc, ssp)

     ! Estimate number of points needed
     C = alphar
     IF ( betar > 0.0 ) C = betar     ! shear?
     deltaz = C / Freq / 20          ! default sampling: 20 points per wavelength
     Nneeded = ( Depth( Medium + 1 ) - Depth( Medium ) ) / deltaz
     Nneeded = MAX( Nneeded, 10 )     ! require a minimum of 10 points
  
     IF ( NG( Medium ) == 0 ) THEN         ! automatic calculation of f.d. mesh
        NG( Medium ) = Nneeded
       if (IWRITE == 1)  WRITE( 6, * ) 'Number of pts = ', NG( Medium )
     ELSEIF ( NG( Medium ) < Nneeded/2 ) THEN
       print *,'Nneeded/2 = ',Nneeded/2
       CALL ERROUT( 6, 'F', 'READIN', 'Mesh is too coarse' )
     END IF

  END DO   ! next Medium

  ! Bottom properties 
  BCType = BotOpt(1:1)
  if (IWRITE == 1) WRITE( 6, * )
  if (IWRITE == 1) WRITE( 6, "( 33X, '( RMS roughness = ', G10.3, ' )' )" ) sigma( NMedia + 1 )

  ! CALL TOPBOT  to read bottom BC 
  CALL TOPBOT( Freq, BCType, AttenUnit, CPB, CSB, rhoB, BumDen, eta, xi, sspHS(2,1:6), nc )
  RETURN

END SUBROUTINE READIN
!**********************************************************************!
SUBROUTINE TOPBOT( Freq, BCType, AttenUnit, CPHS, CSHS, rhoHS, BumDen, eta, xi, sspTMP, nc )

  ! Handles top and bottom boundary conditions

  ! Input:
  !     Freq:   Frequency
  !     BCType: Boundary condition type
  !
  ! Output:
  !    CPHS:    P-wave speed in halfspace
  !    CSHS:    S-wave speed in halfspace
  !    rhoHS:   density in halfspace

  !    BumDen:  Bump density
  !    eta:     Principal radius 1
  !    xi:      Principal radius 2

  IMPLICIT NONE
  INTEGER          :: IWRITE, nc
  REAL    (KIND=8) :: alphaR, betaR, alphaI, betaI, rhoR, Freq, rhoHS, BumDen, eta, xi, ZTMP, sspTMP(1,6)
  COMPLEX (KIND=8) :: CPHS, CSHS, CRCI
  CHARACTER        :: BCType*1, AttenUnit*2
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  IWRITE=0

  ! Echo user's choice of boundary condition 
if (IWRITE == 1) then
    SELECT CASE ( BCType )
    CASE ( 'S' )
       WRITE( 6, * ) '    Twersky SOFT BOSS scatter model'
    CASE ( 'H' )
       WRITE( 6, * ) '    Twersky HARD BOSS scatter model'
    CASE ( 'T' )
       WRITE( 6, * ) '    Twersky (amplitude only) SOFT BOSS scatter model'
    CASE ( 'I' )
       WRITE( 6, * ) '    Twersky (amplitude only) HARD BOSS scatter model'
    CASE ( 'V' )
       WRITE( 6, * ) '    VACUUM'
    CASE ( 'R' )
       WRITE( 6, * ) '    Perfectly RIGID'
    CASE ( 'A' )
       WRITE( 6, * ) '    ACOUSTO-ELASTIC half-space'
    CASE ( 'F' )
       WRITE( 6, * ) '    FILE used for reflection loss'
    CASE ( 'W' )
       WRITE( 6, * ) '    Writing an IRC file'
    CASE ( 'P' )
       WRITE( 6, * ) '    reading PRECALCULATED IRC'
    CASE DEFAULT
       CALL ERROUT( 6, 'F', 'TOPBOT', 'Unknown boundary condition type' )
    END SELECT
endif

  ! Read in BC parameters depending on particular choice 
  CPHS  = 0.0
  CSHS  = 0.0
  rhoHS = 0.0

  ! Twersky ice model parameters 
  IF ( BCType == 'S' .OR. BCType == 'H' .OR. BCType == 'T' .OR. BCType == 'I' ) THEN
     !READ(  ENVFIL, *    ) BumDen, eta, xi
    if (IWRITE == 1)  WRITE( 6, 1000 ) BumDen, eta, xi
1000 FORMAT( /, ' Twersky ice model parameters:', /, &
          ' Bumden = ', G15.6, '  Eta = ', G11.3, '  Xi = ', G11.3, /)
  ENDIF

  !  Half-space properties 
  IF ( BCType == 'A' ) THEN

    ZTMP=sspTMP(1,1)
    alphaR=sspTMP(1,2);
    betaR=sspTMP(1,3);
    rhoR=sspTMP(1,4);
    alphaI=sspTMP(1,5);
    betaI=sspTMP(1,6);
		
     !READ(  ENVFIL, *    ) ZTMP, alphaR, betaR, rhoR, alphaI, betaI
     if (IWRITE == 1) WRITE( 6, 2000 ) ZTMP, alphaR, betaR, rhoR, alphaI, betaI
2000 FORMAT( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )
     CPHS = CRCI( alphaR, alphaI, Freq, AttenUnit )
     CSHS = CRCI( betaR,  betaI,  Freq, AttenUnit )
     rhoHS = rhoR
     IF ( CPHS == 0.0 .OR. rhoHS == 0.0 ) &
          CALL ERROUT( 6, 'F', 'TOPBOT', 'Sound speed or density vanishes in halfspace' )
  ENDIF

  RETURN
END SUBROUTINE TOPBOT
!**********************************************************************!
SUBROUTINE PROFIL( Depth, CP, CS, rhoT, Medium, N1, Freq, SSPType, AttenUnit, Task, nc, ssp )

  ! Call the particular SSP routine specified by SSPType
  ! PROFIL is expected to perform two Tasks:
  !    Task = 'TAB'  then tabulate CP, CS, rhoT
  !    Task = 'INIT' then initialize
  ! Note that Freq is only need if Task = 'INIT'

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  REAL    (KIND=8) :: rhoT( * ), Depth( * )
  INTEGER          :: nc
  REAL    (KIND=8) :: ssp(nc,6)
!  COMPLEX (KIND=8) :: CPT, CST
  COMPLEX (KIND=8) :: CP( * ), CS( * )
  CHARACTER        :: SSPType*1, AttenUnit*2, Task*8
  INTEGER          :: IWRITE

  IWRITE=0

  SELECT CASE ( SSPType )
  CASE ( 'N' )  !  N2-linear profile option 
     CALL N2LIN(  Depth, CP, CS, rhoT, Medium, N1, Freq, AttenUnit, Task, nc, ssp )
  CASE ( 'C' )  !  C-linear profile option 
     CALL CLIN(   Depth, CP, CS, rhoT, Medium, N1, Freq, AttenUnit, Task, nc, ssp )
  CASE ( 'S' )  !  Cubic spline profile option 
     CALL CCUBIC( Depth, CP, CS, rhoT, Medium, N1, Freq, AttenUnit, Task, nc, ssp )
  CASE DEFAULT  !  Non-existent profile option 
     if (IWRITE == 1) WRITE( 6, * ) 'Profile option: ', SSPType(1:1)
     CALL ERROUT( 6, 'F', 'PROFIL', 'Unknown profile option' )
  END SELECT

  RETURN
END SUBROUTINE PROFIL
!**********************************************************************!
FUNCTION CRCI( C, alpha, Freq, AttenUnit )

  ! Converts real wave speed and attenuation to a single
  !  complex wave speed (with positive imaginary part)

  ! 6 CASES:    N for Nepers/meter
  !             M for dB/meter      (M for Meters)
  !             F for dB/m-kHZ      (F for Frequency dependent)
  !             W for dB/wavelength (W for Wavelength)
  !             Q for Q
  !             L for Loss parameter
  !
  ! second letter adds volume attenuation according to standard laws:
  !             T for Thorp

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  PARAMETER ( PI = 3.1415926535897932D0 )
  COMPLEX (KIND=8) :: CRCI
  CHARACTER           AttenUnit*2

  Omega = 2.0 * PI * Freq

  !  Convert to Nepers/m 
  alphaT = 0.0
  SELECT CASE ( AttenUnit(1:1) )
  CASE ( 'N' )
     alphaT = alpha
  CASE ( 'M' )
     alphaT = alpha / 8.6858896D0
  CASE ( 'F' )
     alphaT = alpha * Freq / 8685.8896D0
  CASE ( 'W' )
     IF ( C /= 0.0 ) alphaT = alpha * Freq / ( 8.6858896D0 * C )
     !        The following lines give f^1.25 Frequency dependence
     !        FAC = SQRT( SQRT( Freq / 50.0 ) )
     !        IF ( C /= 0.0 ) alphaT = FAC * alpha * Freq / ( 8.6858896D0 * C )
  CASE ( 'Q' )
     IF( C * alpha /= 0.0 ) alphaT = Omega / ( 2.0 * C * alpha )
  CASE ( 'L' )   ! loss parameter
     IF ( C /= 0.0 ) alphaT = alpha * Omega / C
  END SELECT

  ! added volume attenuation
  SELECT CASE ( AttenUnit(2:2) )
  CASE ( 'T' )
     F2 = ( Freq / 1000.0 ) **2
     Thorpe = 40.0 * F2 / ( 4100.0 + F2 ) + 0.1 * F2 / ( 1.0 + F2 )
     Thorpe = Thorpe / 914.4D0                 ! dB / m
     Thorpe = Thorpe / 8.6858896D0             ! Nepers / m
     alphaT = alphaT + Thorpe
  END SELECT

  ! Convert Nepers/m to equivalent imaginary sound speed 
  alphaT = alphaT * C * C / Omega
  CRCI = CMPLX( C, alphaT, KIND=8 )

  RETURN
END FUNCTION CRCI
!**********************************************************************!
SUBROUTINE N2LIN( Depth, CP, CS, rhoT, Medium, N1, Freq, AttenUnit, Task, nc, ssp )

  ! Tabulate CP, CS, rho for specified Medium
  ! Uses N2-linear segments for P and S-wave speeds
  ! Uses rho-linear segments for density

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  PARAMETER ( MaxMedia = 501, MaxSSP = 2001 )
  INTEGER          :: Loc( MaxMedia ), NSSPPts( MaxMedia )
  REAL (KIND=8)    :: Depth( * ), rhoT( * ), Z( MaxSSP ), rho( MaxSSP )
  INTEGER          :: nc, ind
  REAL (KIND=8)    :: ssp(nc,6)
  COMPLEX (KIND=8) :: CP( * ), CS( * ), alpha( MaxSSP ), beta( MaxSSP ), N2BOT, N2TOP, CRCI
  CHARACTER        :: AttenUnit*2, Task*8
  INTEGER          :: IWRITE
  SAVE Z, alpha, beta, rho, Loc, NSSPPts
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  IWRITE=0
  ! If Task = 'INIT' then this is the first call and SSP is read.
  ! Any other call is a request for SSP subtabulation.

  IF ( Task(1:4) == 'INIT' ) THEN   ! Task 'INIT' for initialization

     ! The variable Loc( Medium ) points to the starting point for the
     ! data in the arrays Z, alpha, beta and rho
     IF ( Medium == 1 ) THEN
        Loc( Medium ) = 0
     ELSE
        Loc( Medium ) = Loc( Medium - 1 ) + NSSPPts( Medium - 1 )
     ENDIF
     ILoc = Loc( Medium )

     !  Read in data and convert attenuation to Nepers/m 
     N1 = 1
     DO I = 1, MaxSSP
        ind=ILoc+I
        Z(ind)=ssp(ind,1)
        alphaR=ssp(ind,2);
        betaR=ssp(ind,3);
        rhoR=ssp(ind,4);
        alphaI=ssp(ind,5);
        betaI=ssp(ind,6);
        if (IWRITE == 1) WRITE( 6, 2000 ) Z( ILoc + I ), alphaR, betaR, rhoR, alphaI, betaI
2000    FORMAT( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )

        alpha(ILoc + I) = CRCI( alphaR, alphaI, Freq, AttenUnit )
        beta( ILoc + I) = CRCI( betaR,  betaI,  Freq, AttenUnit )
        rho(  ILoc + I) = rhoR

        ! Did we read the last point?
        IF ( ABS( Z( ILoc + I ) - Depth( Medium+1 ) ) < EPSILON( 1.0e0 ) * Depth( Medium+1 ) ) THEN
           NSSPPts( Medium ) = N1
           IF ( Medium == 1 ) Depth( 1 ) = Z( 1 )
           RETURN
        ENDIF

        N1 = N1 + 1
     END DO

     ! Fall through means too many points in the profile
     if (IWRITE == 1) WRITE( 6, * ) 'Max. #SSP points: ', MaxSSP
     CALL ERROUT( 6, 'F', 'N2LIN', 'Number of SSP points exceeds limit' )

  ELSE   ! Task = 'TABULATE'
     ILoc = Loc( Medium )
     N    = N1 - 1
     H    = ( Z( ILoc + NSSPPts( Medium ) ) - Z( ILoc + 1 ) ) / N
     Lay  = 1

     DO I = 1, N1
        ZT = Z( ILoc + 1 ) + ( I - 1 ) * H
        IF ( I == N1 ) ZT = Z( ILoc + NSSPPts( Medium ) )   ! Make sure no overshoot

        DO WHILE ( ZT > Z( ILoc + Lay + 1 ) )
           Lay = Lay + 1
        END DO

        R = ( ZT - Z( ILoc + Lay ) ) / ( Z( ILoc + Lay+1 ) - Z( ILoc + Lay ) )

        ! P-wave
        N2TOP   = 1.0 / alpha( ILoc + Lay     )**2
        N2BOT   = 1.0 / alpha( ILoc + Lay + 1 )**2
        CP( I ) = 1.0 / SQRT( ( 1.0 - R ) * N2TOP + R * N2BOT )

        ! S-wave
        IF ( beta(ILoc + Lay) /= 0.0 ) THEN
           N2TOP   = 1.0 / beta( ILoc + Lay     )**2
           N2BOT   = 1.0 / beta( ILoc + Lay + 1 )**2
           CS( I ) = 1.0 / SQRT( ( 1.0 - R ) * N2TOP + R * N2BOT )
        ELSE
           CS( I ) = 0.0
        ENDIF

        rhoT( I ) = ( 1.0 - R ) * rho( ILoc + Lay ) + R * rho( ILoc + Lay + 1 )
     END DO

  ENDIF

  RETURN
END SUBROUTINE N2LIN
!**********************************************************************!
SUBROUTINE CLIN( Depth, CP, CS, rhoT, Medium, N1, Freq, AttenUnit, Task, nc, ssp )

  ! Tabulate CP, CS, rho for specified Medium

  ! Uses c-linear segments for P and S-wave speeds
  ! Uses rho-linear segments for density

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  PARAMETER ( MaxMedia = 501, MaxSSP = 2001 )
  INTEGER          :: Loc( MaxMedia ), NSSPPts( MaxMedia )
  REAL    (KIND=8) :: Depth( * ), rhoT( * ), Z( MaxSSP ), rho( MaxSSP )
  INTEGER          :: nc,ind
  REAL (KIND=8)    :: ssp(nc,6)
  COMPLEX (KIND=8) :: CP( * ), CS( * ), alpha( MaxSSP ), beta( MaxSSP ), CRCI
  CHARACTER        :: AttenUnit*2, Task*8
  INTEGER          :: IWRITE
  SAVE Z, alpha, beta, rho, Loc, NSSPPts
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  IWRITE=0

  ! If Task = 'INIT' then this is the first call and SSP is read.
  ! Any other call is a request for SSP subtabulation.

  IF ( Task(1:4) == 'INIT' ) THEN   ! Task 'INIT' FOR INITIALIZATION
     NSSPPts( Medium ) = N1

     ! The variable Loc(Medium) points to the starting point for the
     ! data in the arrays Z, alpha, beta and rho

     IF ( Medium == 1 ) THEN
        Loc( Medium ) = 0
     ELSE
        Loc( Medium ) = Loc( Medium - 1 ) + NSSPPts( Medium - 1 )
     ENDIF
     ILoc = Loc( Medium )

     !  Read in data and convert attenuation to Nepers/m 
     N1 = 1
     DO I = 1, MaxSSP
        !READ(  ENVFIL, *    ) Z( ILoc + I ), alphaR, betaR, rhoR, alphaI, betaI
        ind=ILoc+I
        Z(ind)=ssp(ind,1)
        alphaR=ssp(ind,2);
        betaR=ssp(ind,3);
        rhoR=ssp(ind,4);
        alphaI=ssp(ind,5);
        betaI=ssp(ind,6);

        if (IWRITE == 1) WRITE( 6, 2000 ) Z( ILoc + I ), alphaR, betaR, rhoR, alphaI, betaI
2000    FORMAT( F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4 )

        alpha( ILoc + I ) = CRCI( alphaR, alphaI, Freq, AttenUnit )
        beta(  ILoc + I ) = CRCI( betaR,  betaI,  Freq, AttenUnit )
        rho(   ILoc + I ) = rhoR

        ! Did we read the last point?
        IF ( ABS( Z( ILoc + I ) - Depth( Medium+1 ) ) <  EPSILON( 1.0e0 ) * Depth( Medium+1 ) ) THEN
           NSSPPts( Medium ) = N1
           IF ( Medium == 1 ) Depth( 1 ) = Z( 1 )
           RETURN
        ENDIF

        N1 = N1 + 1
     END DO

     ! Fall through means too many points in the profile

     if (IWRITE == 1) WRITE( 6, * ) 'Max. #SSP points: ', MaxSSP
     CALL ERROUT( 6, 'F', 'CLIN', 'Number of SSP points exceeds limit' )

  ELSE   ! Task = 'TABULATE'
     ILoc = Loc( Medium )
     N    = N1 - 1
     H    = ( Z( ILoc + NSSPPts( Medium ) ) - Z( ILoc + 1 ) ) / N
     Lay  = 1

     DO I = 1, N1
        ZT = Z( ILoc + 1 ) + ( I - 1 ) * H
        IF ( I == N1 ) ZT = Z( ILoc + NSSPPts( Medium ) )   ! Make sure no overshoot

        DO WHILE ( ZT > Z( ILoc + Lay + 1 ) )
           Lay = Lay + 1
        END DO

        R = ( ZT - Z( ILoc + Lay ) ) / ( Z( ILoc + Lay + 1 ) - Z( ILoc + Lay ) )
        CP(   I ) = ( 1.0 - R ) * alpha( ILoc + Lay ) + R * alpha( ILoc + Lay+1 )
        CS(   I ) = ( 1.0 - R ) *  beta( ILoc + Lay ) + R *  beta( ILoc + Lay+1 )
        rhoT( I ) = ( 1.0 - R ) *   rho( ILoc + Lay ) + R *   rho( ILoc + Lay+1 )
     END DO
  ENDIF

  RETURN
END SUBROUTINE CLIN
!**********************************************************************!
SUBROUTINE CCUBIC( Depth, CP, CS, rhoT, Medium, N1, Freq, AttenUnit, Task, nc, ssp )

  ! Tabulate CP, CS, rho for specified Medium
  ! using cubic spline interpolation

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  PARAMETER ( MaxMedia = 501, MaxSSP = 2001 )
  INTEGER          :: Loc( MaxMedia ), NSSPPts( MaxMedia )
  REAL    (KIND=8) :: Depth( * ), Z( MaxSSP ), rhoT( * )
  INTEGER          :: nc,ind
  REAL (KIND=8)    :: ssp(nc,6)
  COMPLEX (KIND=8) :: CP( * ), CS( * ), ESPLINE, CRCI, alpha( 4, MaxSSP ), beta( 4, MaxSSP ), rho( 4, MaxSSP )
  CHARACTER        :: AttenUnit*2, Task*8
  INTEGER          :: IWRITE
  SAVE Z, alpha, beta, rho, Loc, NSSPPts
  COMMON /CPREV/ alphaR, betaR, rhoR, alphaI, betaI

  IWRITE=0

  ! If Task = 'INIT' then this is the first call and SSP is read.
  ! Any other call is a request for SSP subtabulation.

  IF ( Task(1:4) == 'INIT' ) THEN   ! --- Task 'INIT' for initialization
     NSSPPts( Medium ) = N1

     ! The variable Loc(Medium) points to the starting point for the
     ! data in the arrays Z, alpha, beta and rho

     IF ( Medium == 1 ) THEN
        Loc( Medium ) = 0
     ELSE
        Loc( Medium ) = Loc( Medium - 1 ) + NSSPPts( Medium - 1 )
     ENDIF
     ILoc = Loc( Medium )

     !  Read in data and convert attenuation to Nepers/m 
     N1 = 1

     DO I = 1, MaxSSP
        ind=ILoc+I
        Z(ind)=ssp(ind,1)
        alphaR=ssp(ind,2);
        betaR=ssp(ind,3);
        rhoR=ssp(ind,4);
        alphaI=ssp(ind,5);
        betaI=ssp(ind,6);

        if (IWRITE == 1) WRITE( 6, 2000 ) Z( ILoc + I ), alphaR, betaR, rhoR, alphaI, betaI
2000    FORMAT(F10.2, 3X, 2F10.2, 3X, F6.2, 3X, 2F10.4)

        alpha(1, ILoc + I) = CRCI( alphaR, alphaI, Freq, AttenUnit )
        beta( 1, ILoc + I) = CRCI( betaR,  betaI,  Freq, AttenUnit )
        rho(  1, ILoc + I) = rhoR

        ! Did we read the last point?
        IF ( ABS( Z( ILoc + I ) - Depth(Medium+1) ) <  EPSILON( 1.0e0 ) * Depth( Medium+1 ) ) THEN
           NSSPPts( Medium ) = N1
           IF ( Medium == 1 ) Depth( 1 ) = Z( 1 )

           !   Compute spline coefs 
           IBCBEG = 0
           IBCEND = 0
           CALL CSPLINE( Z( ILoc + 1 ), alpha( 1, ILoc + 1 ), NSSPPts( Medium ), IBCBEG, IBCEND, NSSPPts( Medium ) )
           CALL CSPLINE( Z( ILoc + 1 ),  beta( 1, ILoc + 1 ), NSSPPts( Medium ), IBCBEG, IBCEND, NSSPPts( Medium ) )
           CALL CSPLINE( Z( ILoc + 1 ),   rho( 1, ILoc + 1 ), NSSPPts( Medium ), IBCBEG, IBCEND, NSSPPts( Medium ) )

           RETURN
        ENDIF

        N1 = N1 + 1
     END DO

     ! Fall through means too many points in the profile

     if (IWRITE == 1) WRITE( 6, * ) 'Max. #SSP points: ', MaxSSP
     CALL ERROUT( 6, 'F', 'CCUBIC', 'Number of SSP points exceeds limit' )

  ELSE   ! Task = 'TABULATE'
     ILoc = Loc( Medium )
     N    = N1 - 1
     H    = ( Z( ILoc + NSSPPts( Medium ) ) - Z( ILoc + 1 ) ) / N
     Lay  = 1
     DO I = 1, N1
        ZT = Z( ILoc + 1 ) + ( I - 1 ) * H
        IF ( I == N1 ) ZT = Z( ILoc + NSSPPts( Medium ) )   ! Make sure no overshoot
        DO WHILE ( ZT > Z( ILoc + Lay + 1 ) )
           Lay = Lay + 1
        END DO

        HSPLNE = ZT - Z( ILoc + Lay )

        CP(   I ) =       ESPLINE( alpha( 1, ILoc + Lay ), HSPLNE )
        CS(   I ) =       ESPLINE(  beta( 1, ILoc + Lay ), HSPLNE )
        rhoT( I ) = DBLE( ESPLINE(   rho( 1, ILoc + Lay ), HSPLNE ) )

     END DO
  ENDIF

  RETURN
END SUBROUTINE CCUBIC
