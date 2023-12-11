FUNCTION KUPING( sigma, eta1SQ, rho1, eta2SQ, rho2, P, U ) 

  ! Evalutates the imaginary perturbation due to interfacial roughness using the Kuperman-Ingenito formulation                           
  ! P is the pressure at the interface                                
  ! U is P'/rho       "   "      "                                    

  REAL    (KIND=8  ) :: sigma, rho1, rho2
  COMPLEX (KIND=8)   :: KUPING, DEL, P, U, eta1, eta2, eta1SQ, eta2SQ, SCATRT, A11, A12, A21, A22
  COMPLEX, PARAMETER :: i = (0.0, 1.0)

  KUPING = 0.0 
  IF ( sigma == 0.0 ) RETURN 
  eta1 = SCATRT( eta1SQ ) 
  eta2 = SCATRT( eta2SQ ) 
  DEL  = rho1 * eta2 + rho2 * eta1 

  IF ( DEL /= 0.0 ) THEN
     A11 = 0.5 * ( eta1SQ - eta2SQ ) - ( rho2 * eta1SQ - rho1 * eta2SQ ) * ( eta1 + eta2 ) / DEL
     A12 =   i * ( rho2 - rho1 )**2 * eta1 * eta2 / DEL
     A21 =  -i * ( rho2 * eta1SQ - rho1 * eta2SQ )**2 / ( rho1 * rho2 * DEL )
     A22 = 0.5 * ( eta1SQ - eta2SQ ) + ( rho2 - rho1 ) * eta1 * eta2 * ( eta1 + eta2 ) / DEL

     KUPING = -sigma**2 * ( -A21 * P**2 + ( A11 - A22 ) * P * U + A12 * U**2 )
  ENDIF

  RETURN 
END FUNCTION KUPING

!**********************************************************************C

FUNCTION SCATRT( Z ) 

  ! Root for interfacial scatter                                      

  COMPLEX ( KIND=8 ) :: SCATRT, Z

  IF ( REAL( Z ) >= 0.0 ) THEN 
     SCATRT = SQRT( Z ) 
  ELSE 
     SCATRT = -( 0.0, 1.0 ) * SQRT( -Z ) 
  ENDIF

  RETURN 
END FUNCTION SCATRT
