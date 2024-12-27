SUBROUTINE BCIMP( x, BCType, BOTTOP, CPHS, CSHS, rhoHS, F, G, IPow )

  ! Compute Boundary Condition IMPedance

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  INTEGER, INTENT( OUT ) :: IPow
  REAL (KIND=8) yV( 5 )
  REAL (KIND=8), INTENT( OUT ) :: F, G
  COMPLEX (KIND=8) CPHS, CSHS
  CHARACTER BCType*1, BOTTOP*3

  IPow = 0

  ! Vacuum or Twersky
  IF ( BCType(1:1) == 'V' .OR. &
       BCType(1:1) == 'S' .OR. BCType(1:1) == 'H' .OR. &
       BCType(1:1) == 'T' .OR. BCType(1:1) == 'I' ) THEN
     F = 1.0
     G = 0.0
     yV( 1 : 5 ) = (/ F, G, 0.D0, 0.D0, 0.D0 /)
  ENDIF

  ! Rigid
  IF ( BCType(1:1) == 'R' ) THEN
     F = 0.0
     G = 1.0
     yV( 1 : 5 ) = (/ F, G, 0.D0, 0.D0, 0.D0 /)
  ENDIF

  ! Acousto-elastic half-space
  IF ( BCType(1:1) == 'A' ) THEN
     IF ( REAL( CSHS ) > 0.0 ) THEN
        gammaS2 = x - Omega2 / DBLE( CSHS ) ** 2
        gammaP2 = x - Omega2 / DBLE( CPHS ) ** 2
        gammaS  = SQRT( gammaS2 )   ;   gammaP = SQRT( gammaP2 )
        RMU     = rhoHS * DBLE( CSHS ) ** 2

        yV( 1 ) = ( gammaS*gammaP - x ) / RMU
        yV( 2 ) = ( ( gammaS2 + x ) ** 2 - 4.0*gammaS * gammaP * x ) * RMU
        yV( 3 ) = 2.0 * gammaS * gammaP - gammaS2 - x
        yV( 4 ) = gammaP * ( x - gammaS2 )
        yV( 5 ) = gammaS * ( gammaS2 - x )

        F = Omega2 * yV( 4 )
        G = yV( 2 )
        IF ( G > 0.0 ) ModeCount = ModeCount + 1
     ELSE
        gammaP = SQRT( x - DBLE( Omega2 / CPHS ** 2 ) )
        F = 1.0
        G = rhoHS / gammaP
     ENDIF
  ENDIF

  IF ( BOTTOP(1:3) == 'TOP' ) G = -G

  ! Shoot through elastic layers
  IF ( BOTTOP(1:3) == 'TOP' ) THEN
     ! Shoot down from top
     IF ( FirstAcoustic > 1 ) THEN
        DO Medium = 1, FirstAcoustic - 1
           CALL ELASDN( x, yV, IPow, Medium )
        END DO
        F = Omega2 * yV( 4 )   ;   G = yV( 2 )
     ENDIF
  ELSE
     ! Shoot up from bottom
     IF ( LastAcoustic < NMedia ) THEN
        DO Medium = NMedia, LastAcoustic + 1, -1
           CALL ELASUP( x, yV, IPow, Medium )
        END DO
        F = Omega2 * yV( 4 )
        G = yV( 2 )
     ENDIF
  ENDIF
  RETURN
END SUBROUTINE BCIMP
!**********************************************************************!
SUBROUTINE ELASUP( x, yV, IPow, Medium )

  ! Propagates through an elastic layer using compound matrix formulation

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)

  REAL (KIND=8) :: xV( 5 ), yV( 5 ), zV( 5 )
  PARAMETER ( Roof = 1.0E5, Floor = 1.0E-5, IPowR = 5, IPowF = -5 )

  ! Euler's method for first step
  TWOx   = 2.0 * x
  TWOH   = 2.0 * H( Medium )
  FOURHx = 4.0 * H( Medium ) * x
  J      = LOC( Medium ) + N( Medium ) + 1
  xB3    = x * B3( J ) - rho( J )

  zV(1) = yV(1) - 0.5*(   B1( J ) * yV( 4 ) - B2( J ) * yV( 5 ) )
  zV(2) = yV(2) - 0.5*( -rho( J ) * yV( 4 ) -     xB3 * yV( 5 ) )
  zV(3) = yV(3) - 0.5*(      TWOH * yV( 4 ) + B4( J ) * yV( 5 ) )
  zV(4) = yV(4) - 0.5*(   xB3*yV(1) + B2(J)*yV(2) -TWOx*B4(J)*yV(3))
  zV(5) = yV(5) - 0.5*(rho(J)*yV(1) - B1(J)*yV(2) -    FOURHx*yV(3))

  ! Modified midpoint method
  DO II = N( Medium ), 1, -1
     J = J - 1

     xV = yV
     yV = zV

     xB3 = x * B3( J ) - rho( J )

     zV(1) = xV(1) - (   B1( J ) * yV( 4 ) - B2( J ) * yV( 5 ) )
     zV(2) = xV(2) - ( -rho( J ) * yV( 4 ) -     xB3 * yV( 5 ) )
     zV(3) = xV(3) - (      TWOH * yV( 4 ) + B4( J ) * yV( 5 ) )
     zV(4) = xV(4) - (   xB3*yV(1) + B2(J)*yV(2) - TWOx*B4(J)*yV(3))
     zV(5) = xV(5) - (rho(J)*yV(1) - B1(J)*yV(2) -     FOURHx*yV(3))

     !         P1 = yV(2) * ( yV(4) - yV(5) ) 
     !         P0 = xV(2) * ( xV(4) - xV(5) ) 
     !         IF ( ( P0 > 0.0 ) .AND. (P1 < 0.0) ) ModeCount = ModeCount+1
     !         IF ( ( P0 < 0.0 ) .AND. (P1 > 0.0) ) ModeCount = ModeCount+1

     ! Scale if necessary
     IF ( II /= 1 ) THEN
        IF ( ABS( zV( 2 ) ) < Floor ) THEN
           zV = Roof * zV
           yV = Roof * yV
           IPow = IPow - IPowR
        ENDIF

        IF ( ABS( zV( 2 ) ) > Roof ) THEN
           zV = Floor * zV
           yV = Floor * yV
           IPow = IPow - IPowF
        ENDIF
     ENDIF
  END DO

  yV = ( xV + 2.0 * yV + zV ) / 4.0 ! Apply the standard filter at the terminal point

  RETURN
END SUBROUTINE ELASUP
!**********************************************************************!
SUBROUTINE ELASDN( x, yV, IPow, Medium )

  ! Propagates through an elastic layer using compound matrix formulation

  USE krakmod

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  REAL (KIND=8) :: xV( 5 ), yV( 5 ), zV( 5 )
  PARAMETER ( Roof = 1.0E5, Floor = 1.0E-5, IPowR = 5, IPowF = -5 )

  ! Euler's method for first step
  TWOx   = 2.0 * x
  TWOH   = 2.0 * H( Medium )
  FOURHx = 4.0 * H( Medium ) * x
  J      = LOC( Medium ) + 1
  xB3    = x * B3( J ) - rho( J )

  zV(1) = yV(1) + 0.5*(   B1( J ) * yV( 4 ) - B2( J ) * yV( 5 ) )
  zV(2) = yV(2) + 0.5*( -rho( J ) * yV( 4 ) -     xB3 * yV( 5 ) )
  zV(3) = yV(3) + 0.5*(      TWOH * yV( 4 ) + B4( J ) * yV( 5 ) )
  zV(4) = yV(4) + 0.5*(   xB3*yV(1) + B2(J)*yV(2) -TWOx*B4(J)*yV(3))
  zV(5) = yV(5) + 0.5*(rho(J)*yV(1) - B1(J)*yV(2) -    FOURHx*yV(3))

  ! Modified midpoint method
  DO II = 1, N( Medium )
     J = J + 1

     xV = yV
     yV = zV

     xB3 = x * B3( J ) - rho( J )

     zV(1) = xV(1) + (   B1( J ) * yV( 4 ) - B2( J ) * yV( 5 ) )
     zV(2) = xV(2) + ( -rho( J ) * yV( 4 ) -     xB3 * yV( 5 ) )
     zV(3) = xV(3) + (      TWOH * yV( 4 ) + B4( J ) * yV( 5 ) )
     zV(4) = xV(4) + (   xB3*yV(1) + B2(J)*yV(2) - TWOx*B4(J)*yV(3))
     zV(5) = xV(5) + (rho(J)*yV(1) - B1(J)*yV(2) -     FOURHx*yV(3))

     ! Scale if necessary
     IF ( II /= N( Medium ) ) THEN
        IF ( ABS( zV( 2 ) ) < Floor ) THEN
           zV   = Roof * zV
           yV   = Roof * yV
           IPow = IPow - IPowR
        ENDIF

        IF ( ABS( zV( 2 ) ) > Roof ) THEN
           zV   = Floor * zV
           yV   = Floor * yV
           IPow = IPow - IPowF
        ENDIF
     ENDIF
  END DO

  yV = ( xV + 2.0 * yV + zV ) / 4.0 ! Apply the standard filter at the terminal point

  RETURN
END SUBROUTINE ELASDN
