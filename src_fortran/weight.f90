SUBROUTINE WEIGHT( x, Nx, xTab, NxTab, w, Ix )

  ! Given 
  !    x(*)    abscissas
  !    xTab(*) points for tabulation
  !    Nx      number of x    points
  !    NxTab   number of xtab points

  ! Compute
  !    w(*)    weights for linear interpolation
  !    Ix(*)   indices for    "         "

  IMPLICIT NONE
  INTEGER Nx, NxTab, L, IxTab
  INTEGER Ix( NxTab )
  REAL    x( Nx ), xTab( NxTab ), w( NxTab )

  ! Quick return if just one X value for interpolation ***
  IF ( Nx == 1 ) THEN
     w(  1 ) = 0.0
     Ix( 1 ) = 1
     RETURN
  ENDIF

  L = 1

  DO IxTab = 1, NxTab   ! Loop over each point for which the weights are needed

     ! search for index, L, such that [X(L), X(L+1)] brackets rcvr depth
     DO WHILE ( xTab( IxTab ) > x( L + 1 ) .AND. L < Nx - 1 )
        L = L + 1
     END DO

     ! make note of index, L, and associated weight for interpolation
     Ix( IxTab ) = L
     w(  IxTab ) = ( xTab( IxTab ) - x( L ) ) / ( x( L+1 ) - x( L ) )

  END DO

  RETURN
END SUBROUTINE WEIGHT
