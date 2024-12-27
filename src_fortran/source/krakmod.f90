MODULE krakmod

   SAVE

   INTEGER,          PARAMETER :: MaxMedium = 500, NSets = 5
   REAL    (KIND=8), PARAMETER :: pi = 3.1415926535898D0
   COMPLEX (KIND=8), PARAMETER :: i = ( 0.0, 1.0 ) 

   INTEGER FirstAcoustic, LastAcoustic, NMedia, NV( NSets ), ISet, M, LRECL, ModeCount, Mode, IProf
   REAL    (KIND=8)              :: ET( NSets ), HV( NSets ), rhoT, rhoB, CMin, CLow, CHigh, &
                                    Freq, Omega2, RMax, BumDen, eta, xi
   REAL    (KIND=8), ALLOCATABLE :: EVMat( :, : ), Extrap( :, : ), VG( : )
   COMPLEX (KIND=8)              :: CPT, CST, CPB, CSB
   COMPLEX (KIND=8), ALLOCATABLE :: k( : )

   ! media properties
   INTEGER       ::  LOC(   MaxMedium ), NG( MaxMedium ), N(     MaxMedium )
   REAL (KIND=8) ::  Depth( MaxMedium ), H(  MaxMedium ), SIGMA( MaxMedium )
   CHARACTER     ::  Mater( MaxMedium )*8, TopOpt*8, BotOpt*8, Title*80
   REAL (KIND=8), ALLOCATABLE :: B1( : ), B1C( : ), B2( : ), B3( : ), B4( : ), RHO( : )

END MODULE krakmod
