SUBROUTINE ZBRENTX( X, A, B, T, ERRMSG ) 

  ! MICHAEL B. PORTER  8/84                                               

  ! FORTRAN CONVERSION OF ALGOL PROGRAM PUBLISHED IN                      
  !    THE COMPUTER JOURNAL 14(4):422-425 (1971)                          
  !    BY R. P. BRENT                                                     

  ! RETURNS A ZERO X OF THE FUNCTION F IN THE GIVEN INTERVAL              
  !    [A,B], TO WITHIN A TOLERANCE 6*MACHEP*ABS(X)+2*T, WHERE            
  !    MACHEP IS THE RELATIVE MACHINE PRECISION AND T IS A POSITIVE       
  !    TOLERANCE.  THE PROCEDURE ASSUMES THAT FUNCT(A) AND FUNCT(B) HAVE  
  !    DIFFERENT SIGNS.                                                   

  !    THIS IS THE EXTENDED RANGE VERSION WHICH EXPECTS THE FUNCTION       
  !    TO HAVE THE FORM                                                   
  !         SUBROUTINE FUNCT( X, G, IPOW )
  !    WHERE G * 10 ** IPOW GIVES FUNCT( X )

  IMPLICIT REAL (KIND=8) (A-H, O-Z)
  REAL (KIND=8) :: MACHEP, M
  CHARACTER        ERRMSG*80

  ERRMSG = ' ' 
  MACHEP = 1.0E-16 
  TEN = 10.0 

  CALL FUNCT( A, FA, IEXPA ) 
  CALL FUNCT( B, FB, IEXPB ) 

  IF ( ( (FA > 0.0) .AND. (FB > 0.0) ) .OR.                   &
       ( (FA < 0.0) .AND. (FB < 0.0) ) ) THEN                 
     ERRMSG = ' *** ZBRENT ERROR: FUNCT SGN SAME AT INTRVL ENDPTS' 
     RETURN 
  ENDIF

  ! INTERNAL ROOT                                                     

2000 C  = A 
  FC    = FA 
  IEXPC = IEXPA 
  E     = B-A 
  D     = E 

  ! EXTERNAL ROOT                                                     

  IF ( IEXPA < IEXPB ) THEN 
     F1 = FC * TEN ** ( IEXPC - IEXPB ) 
     F2 = FB 
  ELSE 
     F1 = FC 
     F2 = FB * TEN ** ( IEXPB - IEXPC ) 
  ENDIF

3000 IF ( ABS( F1 ) < ABS( F2 ) ) THEN 
     A     = B 
     B     = C 
     C     = A 
     FA    = FB 
     IEXPA = IEXPB 
     FB    = FC 
     IEXPB = IEXPC 
     FC    = FA 
     IEXPC = IEXPA 
  ENDIF

  TOL = 2.0 * MACHEP * ABS( B ) + T 
  M   = 0.5 * ( C - B ) 
  IF ( ( ABS( M ) > TOL) .AND. ( FB /= 0.0 ) ) THEN 

     ! SEE IF A BISECTION IS FORCED                                  
     IF ( IEXPA < IEXPB ) THEN 
        F1 = FA * TEN ** ( IEXPA - IEXPB ) 
        F2 = FB 
     ELSE 
        F1 = FA 
        F2 = FB * TEN ** ( IEXPB - IEXPA ) 
     ENDIF

     IF ( (ABS(E) < TOL) .OR. ( ABS( F1 ) <= ABS( F2 ) ) ) THEN                            
        E = M 
        D = E 
     ELSE 
        S = FB / FA * TEN ** ( IEXPB - IEXPA ) 
        IF ( A == C ) THEN 
           ! LINEAR INTERPOLATION                                 
           P = 2.0 * M * S 
           Q = 1.0 - S 
        ELSE 
           ! INVERSE QUADRATIC INTERPOLATION                      
           Q = FA / FC * TEN ** ( IEXPA - IEXPC ) 
           R = FB / FC * TEN ** ( IEXPB - IEXPC ) 
           P = S * ( 2.0 * M * Q * ( Q - R ) - ( B - A ) * ( R - 1.0 ) ) 
           Q = ( Q - 1.0 ) * ( R - 1.0 ) * ( S - 1.0 ) 
        ENDIF
        IF ( P > 0.0 ) THEN 
           Q = -Q 
        ELSE 
           P = -P 
        ENDIF
        S = E 
        E = D 
        IF ( ( 2.0 * P < 3.0 * M * Q - ABS( TOL * Q ) ) .AND. ( P < ABS( 0.5 * S * Q ) ) ) THEN                             
           D = P / Q 
        ELSE 
           E = M 
           D = E 
        ENDIF
     ENDIF

     A  = B 
     FA = FB 
     IEXPA = IEXPB 

     IF ( ABS( D ) > TOL) THEN 
        B = B + D 
     ELSE 
        IF ( M > 0.0 ) THEN 
           B = B + TOL 
        ELSE 
           B = B - TOL 
        ENDIF
     ENDIF

     CALL FUNCT( B, FB, IEXPB ) 
     IF ( ( FB > 0.0 ) .EQV. ( FC > 0.0 ) ) GOTO 2000 
     GOTO 3000 
  ENDIF

  X = B 

  RETURN 
END SUBROUTINE ZBRENTX
