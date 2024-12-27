SUBROUTINE MERGEV( X, NX, Y, NY, Z, NZ )

  ! Merge two vectors X, Y into Z (duplicates are removed)

  IMPLICIT NONE
  INTEGER NX, NY, NZ, IX, IY, IZ
  REAL X(*), Y(*), Z(*)

  IX = 1
  IY = 1
  IZ = 1

  ! Loop to successively take values from vector X or Y

  DO WHILE ( IX <= NX .OR. IY <= NY )

     IF      ( IY > NY ) THEN             ! No more in vector Y, take from X 
        Z( IZ ) = X( IX )
        IX = IX + 1
        IZ = IZ + 1
     ELSE IF ( IX > NX ) THEN             ! No more in vector X, take from Y
        Z( IZ ) = Y( IY )
        IY = IY + 1
        IZ = IZ + 1
     ELSE IF ( X( IX ) <= Y( IY ) ) THEN  ! Vector X has the next largest
        Z( IZ ) = X( IX )
        IX = IX + 1
        IZ = IZ + 1
     ELSE                                 ! Vector Y has the next largest
        Z( IZ ) = Y( IY )
        IZ = IZ + 1
        IY = IY + 1
     END IF

     ! Check that we didn't just add a duplicate
     IF ( IZ > 2 ) THEN
        IF ( Z( IZ - 1 ) == Z( IZ - 2 ) ) THEN
           IZ = IZ - 1
        END IF
     END IF
  END DO

  NZ = IZ - 1

  RETURN
END SUBROUTINE MERGEV
