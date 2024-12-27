MODULE sdrdrmod
  ! Reads in source depths, receiver depths, and receiver ranges

  IMPLICIT NONE
  SAVE
  INTEGER :: Nsd, Nrd, NR, Ntheta
  INTEGER, ALLOCATABLE :: Isd( : ), Ird( : )
  REAL,    ALLOCATABLE ::  sd( : ),  rd( : ), WS( : ), WR( : ), R( : ), theta( : )

CONTAINS

  SUBROUTINE SDRD(ZMIN,ZMAX,Nsd,sd,Nrd,rd,zsr,zrc)

    IMPLICIT NONE
    INTEGER       ::  IS, IR, IAllocStat
    REAL (KIND=8) ::  ZMIN, ZMAX
    INTEGER       :: Nsd,Nrd
    REAL, ALLOCATABLE :: sd ( : ), rd( : )
    REAL (KIND=8)   :: zsr(1),zrc(2)
    INTEGER :: IWRITE
    INTEGER :: IQ

    IWRITE=0

    ! *** Read source depths ***

    if (IWRITE == 1) WRITE( 6, * )
    if (IWRITE == 1) WRITE( 6, * ) 'Number of sources   = ', Nsd
    if (IWRITE == 1) WRITE( 6, * ) 'Source depths (m)'

    IF ( Nsd <= 0 ) THEN
       if (IWRITE == 1) WRITE( 6, * ) ' Nsd = ', Nsd
       CALL ERROUT( 6, 'F', 'SDRD', 'Number of sources must be positive'  )
    ENDIF

    IF ( ALLOCATED( sd ) ) DEALLOCATE( sd, WS, Isd )
    ALLOCATE( sd( MAX( 3, Nsd ) ), WS( Nsd ), Isd( Nsd ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) THEN
       if (IWRITE == 1) WRITE( 6, * ) 'Nsd = ', Nsd
       CALL ERROUT( 6, 'F', 'SDRD', 'Too many sources'  )
    END IF

    sd( 3 ) = -999.9
    IQ=size(zsr)
    sd(1:IQ)=zsr(1:IQ)

    CALL SUBTAB( sd, Nsd )
    !CALL SORT(   sd, Nsd )

    if (IWRITE == 1) then 
      IF ( Nsd >= 1  ) WRITE( 6, "( 5G14.6 )" ) ( sd( IS ), IS = 1, MIN( Nsd, 51 ) )
      IF ( Nsd > 51 ) WRITE( 6, * ) ' ... ', sd( Nsd )
    end if

    ! *** Read receiver depths ***

    if (IWRITE == 1) then
         WRITE( 6, * )
         WRITE( 6, * ) 'Number of receivers = ', Nrd
         WRITE( 6, * ) 'Receiver depths (m)'
    end if 

    IF ( Nrd <= 0 ) THEN
       if (IWRITE == 1) WRITE( 6, * ) ' Nrd = ', Nrd
       CALL ERROUT( 6, 'F', 'SDRD', 'Number of receivers must be positive'  )
    ENDIF

    IF ( ALLOCATED( rd ) ) DEALLOCATE( rd, WR, Ird )
    ALLOCATE( rd( MAX( 3, Nrd ) ), WR( Nrd ), Ird( Nrd ), Stat = IAllocStat  )
    IF ( IAllocStat /= 0 ) THEN
       if (IWRITE == 1) WRITE( 6, * ) 'Nrd = ', Nrd
       CALL ERROUT( 6, 'F', 'SDRD', 'Too many receivers'  )
    END IF

    rd( 3 ) = -999.9
    IQ=size(zrc)
    rd(1:IQ)=zrc(1:IQ)

    CALL SUBTAB( rd, Nrd )
    !CALL SORT(   rd, Nrd )
 
if (IWRITE == 1) then
    IF ( Nrd >= 1 ) WRITE( 6, "( 5G14.6 )" ) ( rd( IR ), IR = 1, MIN( Nrd, 51 ) )
    IF ( Nrd > 51 ) WRITE( 6, * ) ' ... ', rd( Nrd )
end if

    ! *** Check for sd/rd in upper or lower halfspace ***

    DO IS = 1, Nsd
       IF      ( sd( IS ) < ZMIN ) THEN
          sd( IS ) = ZMIN
         ! CALL ERROUT( 6, 'W', 'SdRdRMod', 'Source above the top bdry has been moved down' ) 
       ELSE IF ( sd( IS ) > ZMAX ) THEN
          sd( IS ) = ZMAX
         ! CALL ERROUT( 6, 'W', 'SdRdRMod', 'Source below the bottom bdry has been moved up' ) 
       ENDIF
    END DO

    DO IR = 1, Nrd
       IF      ( rd( IR ) < ZMIN ) THEN
          rd( IR ) = ZMIN
       ELSE IF ( rd( IR ) > ZMAX ) THEN
          rd( IR ) = ZMAX
       ENDIF
    END DO

    RETURN
  END SUBROUTINE SDRD

  !********************************************************************

  SUBROUTINE RANGES(NR,R)

    ! Read receiver ranges

    IMPLICIT NONE
    INTEGER   :: IR, IAllocStat
    INTEGER   :: NR
    REAL, ALLOCATABLE  :: R( : )
    INTEGER   :: IWRITE

    IWRITE=0

    if (IWRITE == 1) WRITE( 6, * )
    if (IWRITE == 1) WRITE( 6, * ) 'Number of ranges   = ', NR
    if (IWRITE == 1) WRITE( 6, * ) 'Receiver ranges (km)'

    IF ( ALLOCATED( R ) ) DEALLOCATE( R )
    ALLOCATE( R( MAX( 3, NR ) ), Stat = IAllocStat )
    IF ( IAllocStat /= 0 ) CALL ERROUT( 6, 'F', 'RANGES', 'Too many range points' )

    R( 3 ) = -999.9

    CALL SUBTAB( R, NR )
    CALL SORT(   R, NR )

    if (IWRITE == 1) then
      WRITE( 6, "( 5G14.6 )" ) ( R( IR ), IR = 1, MIN( NR, 51 ) )
      IF ( NR > 51 ) WRITE( 6, * ) ' ... ', R( NR )
    end if

    R( 1:NR ) = 1000.0 * R( 1:NR )   ! Convert ranges to meters

    ! For a point source can't have receiver at origin
    ! IF ( OPT(1:1) == 'R' .AND. R( 1 ) <= 0.0 ) 
    !IF ( R( 1 ) <= 0.0 ) R( 1 ) = MIN( 1.0, R( 2 ) )

    RETURN
  END SUBROUTINE ranges

END MODULE SdRdRMod
