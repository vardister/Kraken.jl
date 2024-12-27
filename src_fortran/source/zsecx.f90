SUBROUTINE ZSECX( x2, TOL, Iteration, MAXIteration, ErrorMessage ) 

  IMPLICIT NONE
  INTEGER       :: Iteration, MAXIteration, IPower0, IPower1
  REAL (KIND=8) :: x0, x1, x2, TOL, shift, F0, F1
  CHARACTER*80     ErrorMessage 

  ! Secant method                                                     

  ErrorMessage = ' ' 
  x1 = x2 + 10.0 * TOL

  CALL FUNCT( x1, F1, IPower1 )
  !WRITE( *, * )
  !WRITE( *, FMT="( 2G24.16, I5 )" ) SQRT( x1 ), F1, IPower1

  DO Iteration = 1, MAXIteration 
     x0      = x1
     F0      = F1 
     IPower0 = IPower1 
     x1      = x2 

     CALL FUNCT( x1, F1, IPower1 ) 

     IF ( F1 == 0.0 ) THEN 
        shift = 0.0 
     ELSE 
        shift = ( x1 - x0 ) / ( 1.0 - F0 / F1 * 10.0 ** ( IPower0 - IPower1 ) )
     ENDIF

     x2 = x1 - shift 
     !WRITE( *, FMT="( 2G24.16, I5 )" ) SQRT( x1 ), F1, IPower1
     IF ( ABS( x2 - x1 ) < TOL .OR. ABS( x2 - x0 ) < TOL ) RETURN
  END DO

  ErrorMessage = ' *** FAILURE TO CONVERGE IN SECANT' 

  RETURN 
END SUBROUTINE ZSECX
