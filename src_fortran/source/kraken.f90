  subroutine kraken(nm,frq,nl,note1,bb,nc,ssp,note2,bsig,sspHS,clh,rng,nsr,zsr,nrc,zrc,nz,cg,cp,kr,att,zm,modes)

  ! Program for solving for ocean acoustic normal modes
  ! Michael B. Porter
  !
  ! Adapted to be employed as a Matlab mex file by B. Dushaw  3/2010

  USE krakmod
  USE sdrdrmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  INTEGER :: Min_LOC( 1 ), IWRITE
  REAL (KIND=8) :: ZMin, ZMax

  INTEGER       :: nm,nl,nc,nsr,nrc,nz
  DOUBLE PRECISION :: frq,bsig,rng
  REAL (KIND=8) :: bb(nl,3), ssp (nc,6), sspHS(2,6), clh(2), zsr(1), zrc(2)
  CHARACTER     :: note1(3), note2(1)
  REAL (KIND=8) :: cg(nm),cp(nm),kr(nm),att(nm),zm(nz),modes(nz,nm)
  INTEGER       :: il, MMM

interface
  subroutine SOLVE(ERROR, nm,nz,zm,modes)

      USE krakmod

      IMPLICIT REAL (KIND=8) (A-H, O-Z)
      INTEGER :: IWRITE, NOMODES

      INTEGER, intent(in)       :: nm,nz
      REAL (KIND=8),intent(out) :: zm(nz), modes(nz,nm)

  end subroutine SOLVE

end interface

! Zero the return variables
  cg=0.0*cg
  cp=0.0*cp
  kr=0.0*kr
  att=0.0*att
  modes=0.0*modes

  !  Change IWRITE to 1 and recompile to get sundry outputs written to STDOUT
  IWRITE=0

  NV( 1 : 5 ) = (/ 1, 2, 4, 8, 16 /)

  FREQ=frq
  NMedia=nl
  TopOpt=note1(1)//note1(2)//note1(3)
  BotOpt=note2(1)

  Depth(1)=0
  do il=1,nl
    NG(il)=bb(il,1)
    sigma(il)=bb(il,2)
    Depth (il+1)=bb(il,3)
  enddo
  SIGMA(nl+1)=bsig ! Changement Rémi SIGMA(il+1)=bsig

  ! Read in environmental info
  CALL READIN( FREQ, MaxMedium, NMedia, &
       TopOpt, CPT, CST, rhoT, BumDen, eta, xi, NG, sigma, Depth, & 
       BotOpt, CPB, CSB, rhoB, nc, ssp, sspHS )


  CLow=clh(1)
  CHigh=clh(2)
  if (IWRITE == 1) WRITE( 6, '( /, '' CLow = '', G12.5, ''  CHigh = '', G12.5 )'  ) CLow, CHigh

  RMax=rng
  if (IWRITE == 1)  WRITE( 6, * ) 'RMax = ', RMax

  ! Read source/receiver depths
  ZMin = Depth( 1 )
  ZMax = Depth( NMedia + 1 )
  Nsd=nsr
  Nrd=nrc
  CALL SDRD( ZMin, ZMax, Nsd, sd, Nrd, rd, zsr, zrc)

  omega2 = ( 2.0 * PI * FREQ ) ** 2
  if (IWRITE == 1) WRITE( 6, * )
  if (IWRITE == 1) WRITE( 6, * ) 'Mesh multiplier   CPU seconds'

  DO ISet = 1, NSets ! Main loop: solve the problem for a sequence of meshes
     N( 1 : NMedia ) = NG( 1 : NMedia ) * NV( ISet )
     H( 1 : NMedia ) = ( Depth( 2 : NMedia + 1 ) - Depth( 1 : NMedia ) ) / N( 1 : NMedia )
     HV( ISet )    = H( 1 )
     CALL SOLVE( ERROR,nm,nz,zm,modes)
     IF ( ERROR * 1000.0 * RMax < 1.0 ) GOTO 3000
  END DO

  ! Fall through indicates failure to converge
  CALL ERROUT( 6, 'W', 'KRAKEN', 'Too many meshes needed: check convergence' )

3000 OMEGA = SQRT( omega2 ) ! Solution complete: discard modes with phase velocity above CHigh

  Min_LOC = MINLOC( Extrap( 1, 1 : M ), Extrap( 1, 1 : M ) > omega2 / CHigh ** 2 )
  M       = Min_LOC( 1 )

  if (IWRITE == 1) WRITE( 6, * )
  if (IWRITE == 1) WRITE( 6, * ) '   I          K             ALPHA          PHASE SPEED       GROUP SPEED'

  k( 1 : M ) = SQRT( Extrap( 1, 1 : M ) + k( 1 : M ) )
  MMM=min(M,nm)

  DO Mode = 1, MMM
     if (IWRITE == 1) WRITE( 6, "( I5, 4G18.10 )" ) Mode, k( Mode ), OMEGA / DBLE( k( Mode ) ), VG( Mode )
     cp(Mode)=OMEGA/DBLE(k(Mode))
     cg(Mode)=VG(Mode)
	 kr(Mode)=DBLE(k(Mode))
     att(Mode)=AIMAG(k(Mode))
  END DO
   
  return
  END SUBROUTINE KRAKEN

! **********************************************************************!
SUBROUTINE INIT

  ! Initializes arrays defining difference equations

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)

  LOGICAL :: ELFLAG = .FALSE.
  REAL    (KIND=8) :: CP2, CS2
  REAL    (KIND=8) :: dummy(6)
  COMPLEX (KIND=8), ALLOCATABLE :: CP( : ), CS( : )
  CHARACTER :: TASK*8

  CMin     = 1.0E6
  FirstAcoustic    = 0
  LOC( 1 ) = 0

  ! Allocate storage for finite-difference coefficients

  NPTS = SUM( N( 1 : NMedia ) ) + NMedia

  IF ( ALLOCATED( B1 ) ) DEALLOCATE( B1, B1C, B2, B3, B4, rho )
  ALLOCATE ( B1( NPTS ), B1C( NPTS ), B2( NPTS ), B3( NPTS ), B4( NPTS ), rho( NPTS ), &
       CP( NPTS ), CS( NPTS ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) &
       CALL ERROUT( 6, 'F', 'KRAKEN - INIT', 'Insufficient memory: Reduce mesh.' )

  DO Medium = 1, NMedia   ! Loop over media
     IF ( Medium /= 1 ) LOC( Medium ) = LOC( Medium-1 ) + N( Medium-1 ) + 1
     N1  = N(   Medium ) + 1
     II  = LOC( Medium ) + 1

     ! PROFIL reads in the data for a medium
     TASK = 'TAB'
      ! Two variables 1 and dummy at the end to ensure the proper number of variables are sent to PROFIL
      ! These are meant to be nc and ssp in the initial call to PROFIL when the data are read.
      ! They are not needed in this call.
      CALL PROFIL( Depth, CP( II ), CS( II ), rho( II ), Medium, N1, FREQ, TopOpt(1 : 1), TopOpt(3 : 4), TASK, 1, dummy )


     ! Load diagonals of the finite-difference equations
     IF ( REAL( CS( II ) ) == 0.0 ) THEN ! Case of an acoustic medium

        Mater( Medium ) = 'ACOUSTIC'
        IF ( FirstAcoustic == 0 ) FirstAcoustic = Medium
        LastAcoustic = Medium

        CMin = MIN( CMin, MINVAL( DBLE( CP( II : II + N( Medium ) ) ) ) )

        B1(  II : II+N(Medium) ) = -2.0 + H( Medium ) ** 2 * DBLE( omega2 / CP( II : II+N(Medium) ) ** 2 )
        B1C( II : II+N(Medium) ) = AIMAG( omega2 / CP( II : II+N(Medium) ) ** 2 )
     ELSE                               ! Case of an elastic medium

        IF ( SIGMA( Medium ) /= 0.0 )&
             CALL ERROUT( 6, 'F', 'KRAKEN', 'Rough elastic interfaces are not allowed' )

        Mater( Medium ) = 'ELASTIC'
        ELFLAG = .TRUE.
        TWOH   = 2.0 * H( Medium )

        DO J = II, II + N( Medium )
           CMin = MIN( DBLE( CS( J ) ), CMin )

           CP2 = CP( J ) ** 2
           CS2 = CS( J ) ** 2

           B1( J )  = TWOH / ( rho( J ) * CS2 )
           B2( J )  = TWOH / ( rho( J ) * CP2 )
           B3( J )  = 4.0 * TWOH * rho( J ) * CS2 * ( CP2 - CS2 ) / CP2
           B4( J )  = TWOH * ( CP2 - 2.0 * CS2 ) / CP2
           rho( J ) = TWOH * omega2 * rho( J )
        END DO

     ENDIF
  END DO   ! next Medium

  ! (CLow, CHigh) = phase speed interval for the mode search
  ! user specified interval is reduced if it exceeds domain
  ! of possible mode phase speeds

  ! Bottom properties
  IF ( BotOpt(1 : 1) == 'A' ) THEN
     IF ( REAL( CSB ) > 0.0 ) THEN        ! Elastic bottom
        ELFLAG = .TRUE.
        CMin  = MIN( CMin,  DBLE( CSB ) )
        CHigh = MIN( CHigh, DBLE( CSB ) )
     ELSE                                 ! Acoustic bottom
        CMin  = MIN( CMin,  DBLE( CPB ) )
        CHigh = MIN( CHigh, DBLE( CPB ) )
     ENDIF
  ENDIF

  ! Top properties
  IF ( TopOpt(2 : 2) == 'A' ) THEN
     IF ( REAL( CST ) > 0.0 ) THEN        ! Elastic  top half-space
        ELFLAG = .TRUE.
        CMin =  MIN( CMin,  DBLE( CST ) )
        CHigh = MIN( CHigh, DBLE( CST ) )
     ELSE                                 ! Acoustic top half-space
        CMin  = MIN( CMin,  DBLE( CPT ) )
        CHigh = MIN( CHigh, DBLE( CPT ) )
     ENDIF
  ENDIF

  ! If elastic medium then reduce CMin for Scholte wave
  IF ( ELFLAG ) CMin = 0.85 * CMin
  CLow = MAX( CLow, CMin )

  RETURN
END SUBROUTINE INIT
! **********************************************************************!
SUBROUTINE SOLVE( ERROR,nm,nz,zm,modes )

  ! Solves the eigenvalue problem at the current mesh
  ! and produces a new extrapolation

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  INTEGER :: IWRITE, NOMODES

  INTEGER, intent(in)          :: nm,nz
  REAL (KIND=8),intent(out)    :: zm(nz), modes(nz,nm)

interface
  subroutine VECTOR(nm,nz,zm,modes)
    USE krakmod
    USE sdrdrmod

    IMPLICIT REAL (KIND=8) (A-H, O-Z)
    INTEGER, ALLOCATABLE :: IZTAB( : )
    REAL :: ZTAB( NSD + NRD )
    REAL, ALLOCATABLE    :: Z( : ), WTS( : ), modesave(:,:)

    INTEGER,intent(in)           :: nm,nz
    REAL (KIND=8),intent(out)    :: zm(nz), modes(nz,nm)

    REAL (KIND=8), ALLOCATABLE :: PHI( : ), D( : ), E( : ), RV1( : ), RV2( : ), RV3( : ), RV4( : )
    COMPLEX, ALLOCATABLE :: PHITAB( : )  ! this could be made real to save space
    CHARACTER            :: BCTop*1, BCBot*1
    INTEGER :: IWRITE

  end subroutine VECTOR

  subroutine SOLVE1(nm)
    USE krakmod

    IMPLICIT REAL (KIND=8) (A-H, O-Z)

    INTEGER,intent(in)         :: nm
    REAL (KIND=8), ALLOCATABLE :: XL( : ), XR( : )
    CHARACTER :: ERRMSG*80
    INTEGER :: IWRITE
  
  end subroutine SOLVE1


end interface

  IWRITE=0

  ! Change NOMODES to 1 to disable calculating the mode eigenfunctions
  ! This also turns off the group phase speeds, which need the eigenfunctions
  NOMODES=0

  CALL CPU_TIME( Tstart )
  CALL INIT

  ! Choose a solver ...
  IF ( IProf > 1 .AND. ISet <= 2 .AND. TopOpt(4 : 4) == 'C' ) THEN
     ! use continuation from last profile if option selected and we're doing the first or second mesh
     CALL SOLVE3
  ELSE IF ( ( ISet <= 2 ) .AND. ( NMedia <= LastAcoustic-FirstAcoustic+1 ) ) THEN
     ! use bisection for first two sets if possible (applicable if elasticity is limited to halfspaces)
     CALL SOLVE1(nm)
  ELSE
     ! use extrapolation from first two meshes
     CALL SOLVE2
  ENDIF

  Extrap( ISet, 1 : M ) = EVMat( ISet, 1 : M )

  IF ( ISet == 1 .and. NOMODES == 0 ) CALL VECTOR(nm,nz,zm,modes)   ! If this is the first mesh, compute the eigenvectors
!  IF ( ISet == 1 ) CALL VECTOR(nm,nz,zm,modes)   ! If this is the first mesh, compute the eigenvectors

  ! Now do the Richardson extrapolation

  ERROR = 10            ! initialize error to a large number
  KEY = 2 * M / 3 + 1   ! KEY value used to check convergence

  IF ( ISet > 1 ) THEN
     T1 = Extrap( 1, KEY )

     DO J = ISet - 1, 1, -1
        DO Mode = 1, M
           X1 = NV( J    ) ** 2
           X2 = NV( ISet ) ** 2
           F1 = Extrap( J,     Mode )
           F2 = Extrap( J + 1, Mode )
           Extrap( J, Mode ) = F2 - ( F1 - F2 ) / ( X2 / X1 - 1.0 )
        END DO
     END DO

     T2 = Extrap( 1, KEY )
     ERROR = ABS( T2 - T1 )
  ENDIF
	 
  CALL CPU_TIME( Tend )
  ET( ISet ) = Tend - Tstart
  if (IWRITE == 1) WRITE( 6, '( 1X, I8, 6X, G15.3, ''s'' )' ) NV( ISet ), ET( ISet )

  RETURN
END SUBROUTINE SOLVE
! **********************************************************************!
SUBROUTINE SOLVE2

  ! Provides initial guess to root finder for each EVMat(I)

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  REAL (KIND=8):: P( 10 )
  CHARACTER    :: ERRMSG*80
  INTEGER :: IWRITE

  IWRITE=0

  X = omega2 / CLow ** 2
  MaxIT = 500
 
  ! solve1 has already allocated space for the following unless the problem has shear
  IF ( .NOT. ALLOCATED( k ) ) THEN
     M = 3000   ! this sets the upper limit on how many modes can be calculated
     ALLOCATE( EVMat( NSets, M ), Extrap( NSets, M ), k( M ), VG( M ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 6, 'F', 'KRAKEN - SOLVE2', 'Insufficient memory (too many modes).' )
  END IF
	

 
  DO Mode = 1, M

     ! For first or second meshes, use a high guess
     ! Otherwise use extrapolation to produce an initial guess

     X = 1.00001 * X

     IF ( ISet >= 2 ) THEN

        P( 1 : ISet - 1 ) = EVMat( 1 : ISet - 1, Mode )

        IF ( ISet >= 3 ) THEN
           DO II = 1, ISet - 2
              DO J = 1, ISet - II - 1
                 X1 = HV( J      ) ** 2
                 X2 = HV( J + II ) ** 2

                 P( J ) = ( ( HV( ISet ) ** 2 - X2 ) * P( J     ) - &
                            ( HV( ISet ) ** 2 - X1 ) * P( J + 1 ) ) &
                          / ( X1 - X2 )
              END DO
           END DO
           X = P( 1 )
        ENDIF

     ENDIF

     !  Use the secant method to refine the eigenvalue
     TOL = ABS( X ) * 10.0 ** ( 4.0 - PRECISION( X ) )
     CALL ZSECX( X, TOL, IT, MaxIT, ERRMSG )

     IF ( ERRMSG /= ' ' ) THEN
        if (IWRITE == 1) WRITE( 6, * ) 'ISet, Mode = ', ISet, Mode
        CALL ERROUT( 6, 'W', 'KRAKEN-ZSECX', ERRMSG )
        X = TINY( X )   ! make sure value is discarded
     ENDIF

     EVMat( ISet, Mode ) = X

     ! Toss out modes outside user specified spectrum
     IF ( omega2 / X > CHigh ** 2 ) THEN
        M = Mode - 1
        RETURN
     ENDIF

  END DO

  RETURN
END SUBROUTINE SOLVE2
! **********************************************************************!
SUBROUTINE SOLVE3

  ! Provides initial guess to root finder for each EVMat(I)
  ! This solver tries to use eigenvalues from a previous profile
  ! as initial guesses

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)

  CHARACTER :: ERRMSG*80
  INTEGER :: IWRITE
  
  IWRITE=0

  MaxIT = 500

  ! Determine number of modes

  XMin = 1.00001D0 * omega2 / CHigh ** 2

  CALL FUNCT( XMin, DELTA, IPower )
  M = ModeCount

  DO Mode = 1, M

     X = EVMat( ISet, Mode )
     TOL = ABS( X ) * 10.0 ** ( 2.0 - PRECISION( X ) )
     CALL ZSECX( X, TOL, IT, MaxIT, ERRMSG )  ! Use the secant method to refine the eigenvalue
	
     IF ( ERRMSG /= ' ' ) THEN
        if (IWRITE == 1) WRITE( 6, * ) 'ISet, Mode = ', ISet, Mode
        CALL ERROUT( 6, 'W', 'KRAKEN-ZSECX', ERRMSG )
        X = TINY( X )   ! make sure value is discarded
     ENDIF

     EVMat( ISet, Mode ) = X

     IF ( omega2 / X > CHigh ** 2 ) THEN  ! Toss out modes outside user specified spectrum
        M = Mode - 1
        RETURN
     ENDIF

  END DO

  RETURN
END SUBROUTINE SOLVE3
! **********************************************************************!
SUBROUTINE FUNCT( X, DELTA, IPower )

  ! FUNCT( X ) = 0 is the dispersion relation

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)

  PARAMETER  ( Roof = 1.0E5, Floor = 1.0E-5, IPowerR = 5, IPowerF = -5 )
  CHARACTER :: BCType*1

  IF ( X <= omega2 / CHigh ** 2 ) THEN    ! For a k below the cts spectrum limit, force a zero
     DELTA  = 0.0
     IPower = 0
     RETURN
  ENDIF

  ModeCount = 0

  BCType(1 : 1) = BotOpt(1 : 1)
  CALL BCIMP( X, BCType, 'BOT', CPB, CSB, rhoB, F,  G,  IPower  )   ! Bottom impedance
  CALL ACOUST( X, F, G, IPower  )                                   ! Shoot through acoustic layers
  BCType(1 : 1) = TopOpt(2 : 2)
  CALL BCIMP( X, BCType, 'TOP', CPT, CST, rhoT, F1, G1, IPower1 )   ! Top impedance

  DELTA = F * G1 - G * F1
  IF ( G * DELTA > 0.0 ) ModeCount = ModeCount + 1
  IPower = IPower + IPower1

  ! Deflate previous roots
  IF ( ( Mode > 1 ) .AND. ( NMedia > LastAcoustic - FirstAcoustic+1 ) ) THEN
     DELTA = DELTA / ( X - EVMat( ISet, Mode-1 ) )

     IF ( Mode > 2 ) THEN
        DO J = 1, Mode - 2
           DELTA = DELTA / ( X - EVMat( ISet, J ) )

           ! Scale if necessary
           DO WHILE ( ABS( DELTA ) < Floor .AND. ABS( DELTA ) > 0.0 )
              DELTA  = Roof * DELTA
              IPower = IPower - IPowerR
           END DO

           DO WHILE ( ABS( DELTA ) > Roof )
              DELTA  = Floor * DELTA
              IPower = IPower - IPowerF
           END DO

        END DO
     ENDIF
  ENDIF
  RETURN
END SUBROUTINE FUNCT
! **********************************************************************!
SUBROUTINE ACOUST( X, F, G, IPower )

  ! Shoot through acoustic layers

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  PARAMETER  ( Roof = 1.0E20, Floor = 1.0E-20, IPowerF = -20 )

  IF ( FirstAcoustic == 0 ) RETURN

  DO Medium = LastAcoustic, FirstAcoustic, -1    ! Loop over successive acoustic media
     H2K2  = H( Medium ) ** 2 * X
     ii    = LOC( Medium ) + N( Medium ) + 1
     rhoM  = rho(  LOC( Medium ) + 1  )   ! density is made homogeneous using value at top of each medium
     P1    = -2.0 * G
     P2    = ( B1( ii ) - H2K2 ) * G - 2.0 * H( Medium ) * F * rhoM

     ! Shoot through a single medium
     DO ii = LOC( Medium ) + N( Medium ), LOC( Medium ) + 1, -1

        P0 = P1
        P1 = P2
        P2 = ( H2K2 - B1( ii ) ) * P1 - P0

        IF ( P0 * P1 <= 0.0 ) ModeCount = ModeCount + 1
        DO WHILE ( ABS( P2 ) > Roof )   ! Scale if necessary
           P0 = Floor * P0
           P1 = Floor * P1
           P2 = Floor * P2
           IPower = IPower - IPowerF
        END DO

     END DO

     ! F = P' / rho and G = -P since F P + G P' / rho = 0
     rhoM = rho( LOC( Medium ) + 1 )
     F = -( P2 - P0 ) / ( 2.0 * H( Medium ) ) / rhoM
     G = -P1
  END DO

  RETURN
END SUBROUTINE ACOUST

! **********************************************************************!

SUBROUTINE VECTOR(nm,nz,zm,modes)

  ! Do inverse iteration to compute each of the eigenvectors and write these to the disk file

  USE krakmod
  USE sdrdrmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  INTEGER, ALLOCATABLE :: IZTAB( : )
  REAL :: ZTAB( NSD + NRD )
  REAL, ALLOCATABLE    :: Z( : ), WTS( : ), modesave(:,:)

  INTEGER,intent(in)           :: nm,nz
  REAL (KIND=8),intent(out)    :: zm(nz), modes(nz,nm)

  REAL (KIND=8), ALLOCATABLE :: PHI( : ), D( : ), E( : ), RV1( : ), RV2( : ), RV3( : ), RV4( : )
  COMPLEX, ALLOCATABLE :: PHITAB( : )  ! this could be made real to save space
  CHARACTER            :: BCTop*1, BCBot*1
  INTEGER :: IWRITE,MMM

  IWRITE=0
   
  BCTop(1 : 1) = TopOpt(2 : 2)
  BCBot(1 : 1) = BotOpt(1 : 1)

  ! Tabulate z-coordinates and off-diagonals of matrix

  NTot  = SUM( N( FirstAcoustic : LastAcoustic ) )
  NTot1 = NTot + 1

  ALLOCATE( Z( NTot1 ), E( NTot1 + 1 ), D( NTot1 ), PHI( NTot1 ), &
       RV1( NTot1 ), RV2( NTot1 ), RV3( NTot1 ), RV4( NTot1 ), Stat = IAllocStat )
  IF ( IAllocStat /= 0 ) CALL ERROUT( 6, 'F', 'KRAKEN - VECTOR', 'Insufficient memory: Reduce mesh.' )

  J      = 1
  Z( 1 ) = Depth( FirstAcoustic )

  DO Medium = FirstAcoustic, LastAcoustic

     Hrho = H( Medium ) * rho( LOC( Medium ) + 1 )

     E( J+1 : J + N( Medium ) ) = 1.0 / Hrho
     Z( J+1 : J + N( Medium ) ) = Z( J ) + H( Medium ) * (/ (JJ, JJ = 1, N(Medium) ) /)

     J = J + N( Medium )
  END DO

  E( NTot1 + 1 ) = 1.0 / Hrho       ! dummy value; never used

  ! Calculate the indices, weights, ... for mode interpolation
  CALL MERGEV( SD, NSD, RD, NRD, ZTAB, NZTAB )
  ALLOCATE( WTS( NZTAB ), IZTAB( NZTAB ), PHITAB( NZTAB ) )
  CALL WEIGHT( Z, NTot1, ZTAB, NZTAB, WTS, IZTAB )

  ! Main loop: for each eigenvalue call SINVIT to get eigenvector

  allocate(modesave(NZTAB,M))

  DO Mode = 1, M
     X = EVMat( 1, Mode )

     ! Corner elt requires top impedance
     CALL BCIMP( X, BCTop, 'TOP', CPT, CST, rhoT, F, G, IPower )
     IF ( G == 0.0 ) THEN
        D( 1 ) = 1.0
        E( 2 ) = EPSILON( D( 1 ) )
     ELSE
        L      = LOC( FirstAcoustic ) + 1
        XH2    = X * H( FirstAcoustic ) * H( FirstAcoustic )
        Hrho   = H( FirstAcoustic ) * rho( L )
        D( 1 ) = ( B1( L ) - XH2 ) / Hrho / 2.0 + F / G
     ENDIF

     ! Set up the diagonal
     ITP = NTot
     J   = 1
     L   = LOC( FirstAcoustic ) + 1

     DO Medium = FirstAcoustic, LastAcoustic
        XH2  = X * H( Medium ) ** 2
        Hrho = H( Medium ) * rho( LOC( Medium ) + 1 )

        IF ( Medium >= FirstAcoustic + 1 ) THEN
           L      = L + 1
           D( J ) = ( D( J ) + ( B1( L ) - XH2 ) / Hrho ) / 2.0
        ENDIF

        DO ii = 1, N( Medium )
           J      = J + 1
           L      = L + 1
           D( J ) = ( B1( L ) - XH2 ) / Hrho

           IF ( B1( L ) - XH2 + 2.0 > 0.0 ) THEN    ! Find index of turning point nearest top
              ITP = MIN( J, ITP )
           ENDIF
        END DO

     END DO

     ! Corner elt requires bottom impedance
     CALL BCIMP( X, BCBot, 'BOT', CPB, CSB, rhoB, F, G, IPower )
     IF ( G == 0.0 ) THEN
        D( NTot1 ) = 1.0
        E( NTot1 ) = EPSILON( D( NTot1 ) )
     ELSE
        D( NTot1 ) =  D( NTot1 ) / 2.0 - F / G
     ENDIF

     CALL SINVIT( NTot1, D, E, IERR, RV1, RV2, RV3, RV4, PHI )   ! Inverse iteration to compute eigenvector

     IF ( IERR /= 0 ) THEN
        if (IWRITE == 1) WRITE( 6, * ) 'Mode = ', Mode
        CALL ERROUT( 6, 'W', 'KRAKEN-SINVIT', 'Inverse iteration failed to converge' )
     ENDIF

     CALL NORMIZ( PHI, ITP, NTot1, X )  ! Normalize the eigenvector

     ! Tabulate the modes at the source/rcvr depths and write to disK
     PHITAB = PHI( IZTAB ) + WTS * ( PHI( IZTAB + 1 ) - PHI( IZTAB ) )
     modesave(:,Mode)=real(PHITAB,8)
  END DO

     MMM=min(M,nm)

     do MZ=1,NZTAB
         do MODE=1,MMM
           modes(MZ,MODE)=modesave(MZ,MODE)
         end do
        zm(MZ)=ZTAB(MZ)
     enddo

  DEALLOCATE( Z, E, D, PHI, RV1, RV2, RV3, RV4 )

  RETURN
END SUBROUTINE VECTOR
! **********************************************************************!
SUBROUTINE NORMIZ( PHI, ITP, NTot1, X )

  ! Normalize the eigenvector:
  !    SQNRM = Integral(PHI ** 2) by the trapezoidal rule:
  !    Integral(F) = H*( F(1)+...+F(N-1) + 0.5*(F(0)+F(N)) )

  ! Compute perturbation due to material absorption
  ! Compute the group velocity
  ! Call SCAT to figure interfacial scatter loss

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  REAL    (KIND=8) :: PHI( NTot1 ), SLOW
  COMPLEX (KIND=8) :: PERK, DEL
  CHARACTER        :: BCType*1
  INTEGER :: IWRITE
  
  IWRITE=0
 
  SQNRM = 0.0
  PERK  = 0.0
  SLOW  = 0.0

 
  ! Compute perturbation due to loss in top half-space
  IF ( TopOpt(2 : 2) == 'A' ) THEN
     DEL = -0.5*( omega2 / CPT ** 2 - DBLE( omega2 / CPT ** 2 ) ) / &
            SQRT( X                 - DBLE( omega2 / CPT ** 2 ) )
     PERK = PERK - DEL * PHI( 1 ) ** 2 / rhoT
     SLOW = SLOW + PHI( 1 ) ** 2 / ( 2 * SQRT( X - DBLE( omega2 / CPT ** 2 ) ) ) / ( rhoT * DBLE( CPT ) ** 2 )

  ENDIF

  ! Compute norm and pertubation due to material absorption
  L = LOC( FirstAcoustic )
  J = 1

  DO Medium = FirstAcoustic, LastAcoustic
     L       = L + 1
     rhoM    = rho( L )
     rhoOMH2 = rhoM * omega2 * H( Medium ) ** 2

     ! top interface
     SQNRM = SQNRM + 0.5 * H( Medium ) *                  PHI( J ) ** 2 / rhoM
     PERK  = PERK  + 0.5 * H( Medium ) * i * B1C(L)     * PHI( J ) ** 2 / rhoM
     SLOW  = SLOW  + 0.5 * H( Medium ) * ( B1(L) + 2. ) * PHI( J ) ** 2 / rhoOMH2

     ! medium
     L1 = L + 1   ;   L = L + N( Medium ) - 1
     J1 = J + 1   ;   J = J + N( Medium ) - 1
     SQNRM = SQNRM +  H( Medium ) *     SUM(                     PHI( J1 : J ) ** 2 ) / rhoM
     PERK  = PERK  +  H( Medium ) * i * SUM(  B1C(L1 : L) *        PHI( J1 : J ) ** 2 ) / rhoM
     SLOW  = SLOW  +  H( Medium ) *     SUM( ( B1(L1 : L) + 2. ) * PHI( J1 : J ) ** 2 ) / rhoOMH2
  
     ! bottom interface
     L = L + 1   ;   J = J + 1
     SQNRM = SQNRM + 0.5 * H( Medium ) *                  PHI( J ) ** 2 / rhoM
     PERK  = PERK  + 0.5 * H( Medium ) * i * B1C(L)     * PHI( J ) ** 2 / rhoM
     SLOW  = SLOW  + 0.5 * H( Medium ) * ( B1(L) + 2. ) * PHI( J ) ** 2 / rhoOMH2
	
  END DO
  ! Compute perturbation due to loss in bottom half-space
  IF ( BotOpt(1 : 1) == 'A' ) THEN
     DEL  = -0.5*( omega2 / CPB ** 2 - DBLE( omega2 / CPB ** 2 ) ) / &
             SQRT( X                 - DBLE( omega2 / CPB ** 2 ) )
     PERK = PERK - DEL * PHI( J ) ** 2 / rhoB
     SLOW = SLOW + PHI( J ) ** 2 / ( 2 * SQRT( X - DBLE( omega2 / CPB ** 2 ) ) ) / ( rhoB * DBLE( CPB ) ** 2 )
  ENDIF
  
!   WRITE(6,*)  'PERK1=', PERK
  ! Compute derivative of top admitance
  X1 = 0.9999999D0 * X
  X2 = 1.0000001D0 * X

  BCType(1 : 1) = TopOpt(2 : 2)
  CALL BCIMP( X1, BCType, 'TOP', CPT, CST, rhoT, F1, G1, IPower )
  CALL BCIMP( X2, BCType, 'TOP', CPT, CST, rhoT, F2, G2, IPower )
  DrhoDX = 0.0
  IF ( G1 /= 0.0 ) DrhoDX = -( F2 / G2 - F1 / G1 ) / ( X2 - X1 )

  ! Compute derivative of bottom admitance
  BCType(1 : 1) = BotOpt(1 : 1)
  CALL BCIMP( X1, BCType, 'BOT', CPB, CSB, rhoB, F1, G1, IPower )
  CALL BCIMP( X2, BCType, 'BOT', CPB, CSB, rhoB, F2, G2, IPower )
  DetaDX = 0.0
  IF ( G1 /= 0.0 ) DetaDX = -( F2 / G2 - F1 / G1 ) / ( X2 - X1 )

  ! Scale the mode
  RN = SQNRM + DrhoDX * PHI( 1 ) ** 2 - DetaDX * PHI( NTot1 ) ** 2

  IF ( RN <= 0.0 ) THEN
     RN = -RN
     if (IWRITE == 1) WRITE( 6, * ) 'Mode = ', Mode
     CALL ERROUT( 6, 'W', 'KRAKEN', 'Normalization constant non-positive' )
     WRITE( *, * ) PHI
  ENDIF

  SCALEF = 1.0 / SQRT( RN )
  IF ( PHI( ITP ) < 0.0 ) SCALEF = -SCALEF

  PHI( 1 : NTot1 ) = SCALEF * PHI( 1 : NTot1 )
  PERK           = SCALEF ** 2 * PERK
  SLOW           = SCALEF ** 2 * SLOW * SQRT( omega2 / X )
  VG( Mode )     = 1 / SLOW

  CALL SCAT( PERK, PHI, X ) ! Compute interfacial scatter loss
  	
  RETURN
END SUBROUTINE NORMIZ
! **********************************************************************!
SUBROUTINE SCAT( PERK, PHI, X )

  ! Figure scatter loss

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  REAL    (KIND=8) :: PHI( * )
  COMPLEX (KIND=8) :: CImped, KX, TWERSK, PERK, eta1SQ, eta2SQ, KUPING, U, PHIC
  CHARACTER        :: BCType*1

  OMEGA = SQRT( omega2 )
  KX    = SQRT( X )

  ! Top Twersky roughness
  BCType(1 : 1) = TopOpt(2 : 2)
  IF ( BCType(1 : 1) == 'S' .OR. BCType(1 : 1) == 'H' .OR. BCType(1 : 1) == 'T' .OR. BCType(1 : 1) == 'I' ) THEN

     Itop   = LOC( FirstAcoustic ) + N( FirstAcoustic ) + 1
     rhoINS = rho( Itop )
     CINS   = SQRT( omega2 * H( FirstAcoustic ) ** 2 / ( 2.0 + B1( FirstAcoustic ) ) )
     CImped = TWERSK( BCType(1 : 1), OMEGA, BumDen, xi, eta, KX, rhoINS, CINS )
     CImped = CImped / ( -i * OMEGA * rhoINS )
     DPHIDZ = PHI( 2 ) / H( FirstAcoustic )
     PERK   = PERK - CImped * DPHIDZ ** 2
  ENDIF

  ! Bottom Twersky roughness
  BCType(1 : 1) = BotOpt(1 : 1)
  IF ( BCType(1 : 1) == 'S' .OR. BCType(1 : 1) == 'H' .OR. BCType(1 : 1) == 'T' .OR. BCType(1 : 1) == 'I' ) THEN

     Ibot   = LOC( LastAcoustic ) + N( LastAcoustic ) + 1
     rhoINS = rho( Ibot )
     CINS   = SQRT( omega2 * H(LastAcoustic) ** 2 / ( 2.0 + B1(LastAcoustic) ) )
     CImped = TWERSK( BCType(1 : 1), OMEGA, BumDen, xi, eta, KX, rhoINS, CINS )
     CImped = CImped / ( -i * OMEGA * rhoINS )
     DPHIDZ = PHI( 2 ) / H( FirstAcoustic )
     PERK   = PERK - CImped * DPHIDZ ** 2
  ENDIF

  J = 1
  L = LOC( FirstAcoustic )

  DO Medium = FirstAcoustic - 1, LastAcoustic  ! Loop over media

     ! Calculate rho1, eta1SQ, PHI, U

     IF ( Medium == FirstAcoustic - 1 ) THEN
        ! Top properties
        BCType(1 : 1) = TopOpt(2 : 2)
        SELECT CASE ( BCType(1 : 1) )
        CASE ( 'A' )          ! Acousto-elastic
           rho1   = rhoT
           eta1SQ = X - omega2 / CPT ** 2
           U = SQRT( eta1SQ ) * PHI( 1 ) / rhoT
        CASE ( 'V' )          ! Vacuum
           rho1   = 1.0E-9
           eta1SQ = 1.0
           rhoINS = rho( LOC( FirstAcoustic ) + 1 )
           U      = PHI( 2 ) / H( FirstAcoustic ) / rhoINS
        CASE ( 'R' )          ! Rigid
           rho1   = 1.0E+9
           eta1SQ = 1.0
           U      = 0.0
        END SELECT
     ELSE
        H2 = H( Medium ) ** 2
        J  = J + N( Medium )
        L  = LOC( Medium ) + N( Medium ) + 1

        rho1   = rho( L )
        eta1SQ = ( 2.0 + B1( L ) ) / H2 - X
        U = ( -PHI( J-1 ) - 0.5 * ( B1(L) - H2*X ) * PHI( J ) ) / ( H( Medium ) * rho1 )
     ENDIF

     ! Calculate rho2, eta2

     IF ( Medium == LastAcoustic ) THEN
        ! Bottom properties

        BCType(1 : 1) = BotOpt(1 : 1)
        SELECT CASE ( BCType(1 : 1) )
        CASE ( 'A' )          ! Acousto-elastic
           rho2   = rhoB
           eta2SQ = omega2 / CPB ** 2 - X
        CASE ( 'V' )          ! Vacuum
           rho2   = 1.0E-9
           eta2SQ = 1.0
        CASE ( 'R' )          ! Rigid
           rho2   = 1.0E+9
           eta2SQ = 1.0
        END SELECT
     ELSE
        rho2   = rho( L + 1 )
        eta2SQ = ( 2.0 + B1( L + 1 ) ) / H( Medium + 1 ) ** 2 - X
     ENDIF

     PHIC = PHI( J )   ! convert to complex*16
     PERK = PERK + KUPING( SIGMA( Medium + 1 ), eta1SQ, rho1, eta2SQ, rho2, PHIC, U )

  END DO

  k( Mode ) = PERK

  RETURN
END SUBROUTINE SCAT
! **********************************************************************!
SUBROUTINE SOLVE1(nm)

  ! Uses Sturm sequences to isolate the eigenvalues
  ! and Brent's method to refine them

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)

  INTEGER,intent(in)         :: nm
  REAL (KIND=8), ALLOCATABLE :: XL( : ), XR( : )
  CHARACTER :: ERRMSG*80
  INTEGER :: IWRITE

  IWRITE=0

  ! Determine number of modes

  XMin = 1.00001D0 * omega2 / CHigh ** 2

  CALL FUNCT( XMin, DELTA, IPower )
  M = ModeCount
  if (IWRITE == 1) WRITE( 6, * ) ' --- Number of modes = ', M

  IF ( ALLOCATED( XL ) ) DEALLOCATE( XL, XR )
  ALLOCATE( XL( M + 1 ), XR( M + 1 ) )

  IF ( ISet == 1 ) THEN
     IF ( ALLOCATED ( EVMat ) ) DEALLOCATE( EVMat, Extrap, k, VG )
     ALLOCATE( EVMat( NSets, M ), Extrap( NSets, M ), k( M ), VG( M ), Stat = IAllocStat )
     IF ( IAllocStat /= 0 ) &
          CALL ERROUT( 6, 'F', 'KRAKEN - SOLVE1', 'Insufficient memory (too many modes).' )
  END IF

  XMax = omega2 / CLow ** 2
  CALL FUNCT( XMax, DELTA, IPower )
  M = M - ModeCount
  M = min(M,nm+1 )        ! Modification by Porter to save just the first nm modes.
                          ! Use nm+1 to have BISEC properly bracket the eigenvalues

  NTot = SUM( N( FirstAcoustic : LastAcoustic ) )

  IF ( M > NTot / 5 ) THEN
     if (IWRITE == 1) WRITE( 6, * ) 'Approximate number of modes =', M
     CALL ERROUT( 6, 'W', 'KRAKEN', 'Mesh too coarse to sample the modes adequately' )
  ENDIF

  CALL BISECT( XMin, XMax, XL, XR )   ! Initialize upper and lower bounds

  M = min(M,nm)           ! Correct for nm+1 above

  ! Call ZBRENT to refine each eigenvalue in turn
  DO Mode = 1, M
     X1  = XL( Mode )
     X2  = XR( Mode )
     EPS = ABS( X2 ) * 10.0 ** ( 2.0 - PRECISION( X2 ) )
     CALL ZBRENTX( X, X1, X2, EPS, ERRMSG )

     IF ( ERRMSG /= ' ' ) THEN
        if (IWRITE == 1) WRITE( 6, * ) 'ISet, Mode = ', ISet, Mode
        CALL ERROUT( 6, 'W', 'KRAKEN-ZBRENTX', ERRMSG )
     ENDIF

     EVMat( ISet, Mode ) = X
  END DO

  DEALLOCATE( XL, XR )

  RETURN
END SUBROUTINE SOLVE1
! **********************************************************************!
SUBROUTINE BISECT( XMin, XMax, XL, XR )

  ! Returns an isolating interval (XL, XR) for each eigenvalue

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  INTEGER, PARAMETER :: MaxBIS = 50
  REAL (KIND=8)      :: XL( M + 1 ), XR( M + 1 )

  XL = XMin   ! initial left  boundary
  XR = XMax   ! initial right boundary

  CALL FUNCT( XMax, DELTA, IPower )
  NZER1 = ModeCount

  IF ( M == 1 ) RETURN   ! quick exit if only one mode is sought

   MODELOOP: DO Mode = 1, M - 1  ! loop over eigenvalue

     ! Obtain initial guesses for X1 and X2
     IF ( XL( Mode ) == XMin ) THEN
        X2 = XR( Mode )
        X1 = MAX( MAXVAL( XL( Mode + 1 : M ) ), XMin )

        ! Begin bisection (allowing no more than MaxBIS bisections per mode)
        DO J = 1, MaxBIS
           X = X1 + ( X2 - X1 ) / 2
           CALL FUNCT( X, DELTA, IPower )
           NZeros = ModeCount - NZER1

           IF ( NZeros < Mode ) THEN   ! not too many zeros, this is a new right bdry
              X2 = X
              XR( Mode ) = X
           ELSE                        ! this is a new left bdry
              X1 = X
              IF ( XR( NZeros + 1 ) >= X ) XR( NZeros + 1 ) = X
              IF ( XL( NZeros     ) <= X ) XL( NZeros     ) = X
           ENDIF

           ! when we have replaced the default, initial values, we are done
           IF ( XL( Mode ) /= XMin ) CYCLE MODELOOP
        END DO
     ENDIF
  END DO MODELOOP

  RETURN
END SUBROUTINE BISECT
