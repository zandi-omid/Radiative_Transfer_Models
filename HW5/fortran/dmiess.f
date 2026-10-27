      SUBROUTINE DMIESS(  RO,      RFR,     RFI,     THETD,     JX,
     2                    QEXT,    QSCAT,   CTBRQS,  ELTRMX,    PI,
     3                    TAU,     CSTHT,   SI2THT,  ACAP, QBS, IT,
     4                    LL,      R,       RE2,     TMAG2,     WVNO  )
c
      implicit real*8  (a-h,o-z)
C
      COMMON / WARRAY / W(3,30000)
C
      COMPLEX*16    FNAP,      FNBP,      ACAP(LL),     W,
     2              FNA,       FNB,       RF,           RRF,
     3              RRFX,      WM1,       FN1,          FN2,
     4              TC1,       TC2,       WFN(2),       Z(4),
     5              K1,        K2,        K3,
     6              RC,        U(8),         DH1,
     7              DH2,       DH4,       P24H24,       P24H21,
     8              PSTORE,    HSTORE,    DUMMY,        DUMSQ
C
      real*8        T(5),      TA(4),     TB(2),        TC(2),
     2              TD(2),     TE(2),     PI( 3,IT ),   TAU( 3,IT ),
     3              CSTHT(IT), THETD(IT), SI2THT(IT),   ELTRMX( 4,IT,2 )
C
      EQUIVALENCE   (  FNA,TB(1) ),   ( FNB,TC(1) ),
     2              ( FNAP,TD(1) ),   (FNBP,TE(1) )
C
      IFLAG = 1
C     Equal core/shell material with equal radii is exactly homogeneous.
C     Avoid removable singularities in the unnecessary coated-sphere algebra.
      IF ( R .EQ. RO .AND. RE2 .EQ. RFR .AND.
     1     TMAG2 .EQ. RFI ) IFLAG = 2
      IF ( R/RO .LT. 1.0E-06 )   IFLAG = 2
C
      IF ( JX .LE. IT )   GO TO 20
         WRITE( 6,7 )
         WRITE( 6,6 )
         STOP 30
   20 RF =  DCMPLX( RFR,  -RFI )
      RC =  DCMPLX( RE2,-TMAG2 )
      X  =  RO * WVNO
C     Do not terminate at an isolated small coefficient before the tail.
      NSTOP = CEILING( X + 4.0D0 * X**(1.0D0/3.0D0) + 2.0D0 )
      K1 =  RC * WVNO
      K2 =  RF * WVNO
      K3 =  DCMPLX( WVNO, 0.0d+00 )
      Z(1) =  K2 * RO
      Z(2) =  K3 * RO
      Z(3) =  K1 * R
      Z(4) =  K2 * R
C
      X1   =  DREAL( Z(1) )
      Y1   =  DIMAG( Z(1) )
      X4   =  DREAL( Z(4) )
      Y4   =  DIMAG( Z(4) )
      RRF  =  1.0 / RF
      RX   =  1.0 / X
      RRFX =  RRF * RX
      T(1) =  ( X**2 ) * ( RFR**2 + RFI**2 )
      T(1) =  DSQRT( T(1) )
      NMX1 =  1.30 * T(1)
C
      IF ( NMX1 .LE. LL-1 )   GO TO 21
         WRITE(6,8)
         STOP 32
   21 NMX2 = 1.1 * T(1)
      IF ( NMX1 .GT.  200 )   GO TO 22
         NMX1 = 200
         NMX2 = 190
C
   22 ACAP( NMX1+1 )  =  ( 0.0,0.0 )
      IF ( IFLAG .EQ. 2 )   GO TO 26
         DO 29   N = 1,3
   29    W( N,NMX1+1 )  =  ( 0.0,0.0 )
   26 CONTINUE
      DO 23   N = 1,NMX1
         NN = NMX1 - N + 1
         ACAP(NN) = (NN+1) * RRFX - 1.0 / ( (NN+1) * RRFX + ACAP(NN+1) )
C
         IF ( IFLAG .EQ. 2 )   GO TO 23
            DO 31   M = 1,3
   31       W( M,NN ) = (NN+1) / Z(M+1)  -
     1                       1.0 / (  (NN+1) / Z(M+1)  +  W( M,NN+1 )  )
   23 CONTINUE
C
      DO 30   J = 1,JX
      IF ( THETD(J) .LT. 0.0 )  THETD(J) =  ABS( THETD(J) )
      IF ( THETD(J) .GT. 0.0 )  GO TO 24
      CSTHT(J)  = 1.0
      SI2THT(J) = 0.0
      GO TO 30
   24 IF ( THETD(J) .GE. 90.0 )  GO TO 25
      T(1)      =  ( 3.14159265359 * THETD(J) ) / 180.0
      CSTHT(J)  =  DCOS( T(1) )
      SI2THT(J) =  1.0 - CSTHT(J)**2
      GO TO 30
   25 IF ( THETD(J) .GT. 90.0 )  GO TO 28
      CSTHT(J)  =  0.0
      SI2THT(J) =  1.0
      GO TO 30
   28 WRITE( 6,5 )  THETD(J)
      WRITE( 6,6 )
      STOP 34
   30 CONTINUE
C
      DO 35  J = 1,JX
      PI(1,J)  =  0.0
      PI(2,J)  =  1.0
      TAU(1,J) =  0.0
      TAU(2,J) =  CSTHT(J)
   35 CONTINUE
C
C   INITIALIZATION OF HOMOGENEOUS SPHERE
      T(1)   =  DCOS(X)
      T(2)   =  DSIN(X)
      WM1    =  DCMPLX( T(1),-T(2) )
      WFN(1) =  DCMPLX( T(2), T(1) )
      TA(1)  =  T(2)
      TA(2)  =  T(1)
      WFN(2) =  RX * WFN(1) - WM1
      TA(3)  =   DREAL(WFN(2))
      TA(4)  =   DIMAG(WFN(2))
C
      IF ( IFLAG .EQ. 2 )   GO TO 560
C
      N = 1
C
C INITIALIZATION PROCEDURE FOR STRATIFIED SPHERE BEGINS HERE
C
      SINX1   =  DSIN( X1 )
      SINX4   =  DSIN( X4 )
      COSX1   =  DCOS( X1 )
      COSX4   =  DCOS( X4 )
      EY1     =  DEXP( Y1 )
      E2Y1    =  EY1 * EY1
      EY4     =  DEXP( Y4 )
      EY1MY4  =  DEXP( Y1 - Y4 )
      EY1PY4  =  EY1 * EY4
      EY1MY4  =  DEXP( Y1 - Y4 )
      AA  =  SINX4 * ( EY1PY4 + EY1MY4 )
      BB  =  COSX4 * ( EY1PY4 - EY1MY4 )
      CC  =  SINX1 * ( E2Y1 + 1.0 )
      DD  =  COSX1 * ( E2Y1 - 1.0 )
      DENOM   =  1.0  +  E2Y1 * ( 4.0 * SINX1 * SINX1 - 2.0 + E2Y1 )
      REALP   =  ( AA * CC  +  BB * DD ) / DENOM
      AMAGP   =  ( BB * CC  -  AA * DD ) / DENOM
      DUMMY   =  DCMPLX( REALP, AMAGP )
      AA  =  SINX4 * SINX4 - 0.5
      BB  =  COSX4 * SINX4
      P24H24  =  0.5 + DCMPLX( AA,BB ) * EY4 * EY4
      AA  =  SINX1 * SINX4  -  COSX1 * COSX4
      BB  =  SINX1 * COSX4  +  COSX1 * SINX4
      CC  =  SINX1 * SINX4  +  COSX1 * COSX4
      DD  = -SINX1 * COSX4  +  COSX1 * SINX4
      P24H21  =  0.5 * DCMPLX( AA,BB ) * EY1 * EY4  +
     2           0.5 * DCMPLX( CC,DD ) * EY1MY4
      DH4  =  Z(4) / ( 1.0 + ( 0.0,1.0 ) * Z(4) )  -  1.0 / Z(4)
      DH1  =  Z(1) / ( 1.0 + ( 0.0,1.0 ) * Z(1) )  -  1.0 / Z(1)
      DH2  =  Z(2) / ( 1.0 + ( 0.0,1.0 ) * Z(2) )  -  1.0 / Z(2)
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
C
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0,-1.0 )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
C
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) /
     2               ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) /
     2               ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
c 
c temporary print
c
c     print 742,  n
c     print 744,  u(1), u(2)
c     print 745,  u(3), u(4) 
c     print 746,  u(5), u(6) 
c     print 748,  u(7), u(8)
c     print 752,  fna, fnb

c
      GO TO 561
  560 TC1  =  ACAP(1) * RRF  +  RX
      TC2  =  ACAP(1) * RF   +  RX
      FNA  =  ( TC1 * TA(3)  -  TA(1) ) / ( TC1 * WFN(2)  -  WFN(1) )
      FNB  =  ( TC2 * TA(3)  -  TA(1) ) / ( TC2 * WFN(2)  -  WFN(1) )
C
  561 CONTINUE
      FNAP = FNA
      FNBP = FNB
      T(1) = 1.50
C
      TB(1) = T(1) * TB(1)
      TB(2) = T(1) * TB(2)
      TC(1) = T(1) * TC(1)
      TC(2) = T(1) * TC(2)
      DO 60 J = 1,JX
          ELTRMX(1,J,1) = TB(1) * PI(2,J) + TC(1) * TAU(2,J)
          ELTRMX(2,J,1) = TB(2) * PI(2,J) + TC(2) * TAU(2,J)
          ELTRMX(3,J,1) = TC(1) * PI(2,J) + TB(1) * TAU(2,J)
          ELTRMX(4,J,1) = TC(2) * PI(2,J) + TB(2) * TAU(2,J)
          ELTRMX(1,J,2) = TB(1) * PI(2,J) - TC(1) * TAU(2,J)
          ELTRMX(2,J,2) = TB(2) * PI(2,J) - TC(2) * TAU(2,J)
          ELTRMX(3,J,2) = TC(1) * PI(2,J) - TB(1) * TAU(2,J)
          ELTRMX(4,J,2) = TC(2) * PI(2,J) - TB(2) * TAU(2,J)
   60 CONTINUE
C
      QEXT   = 2.0 * ( TB(1) + TC(1))
      QSCAT  = ( TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2 ) / 0.75
C
      QBSR   = -2.0*(TC(1) - TB(1))
      QBSI   = -2.0*(TC(2) - TB(2))
      RMM    = -1.0
C
      CTBRQS = 0.0
      N = 2
   65 T(1) = 2*N - 1
      T(2) =   N - 1
      T(3) = 2*N + 1
      DO 70  J = 1,JX
          PI(3,J)  = ( T(1) * PI(2,J) * CSTHT(J) - N * PI(1,J) ) / T(2)
          TAU(3,J) = CSTHT(J) * ( PI(3,J) - PI(1,J) )  -
     1                          T(1) * SI2THT(J) * PI(2,J)  +  TAU(1,J)
   70 CONTINUE
C
C   HERE SET UP HOMOGENEOUS SPHERE
      WM1    =  WFN(1)
      WFN(1) =  WFN(2)
      TA(1)  =  DREAL(WFN(1))
      TA(2)  =  DIMAG(WFN(1))
      WFN(2) =  T(1) * RX * WFN(1)  -  WM1
      TA(3)  =  DREAL(WFN(2))
      TA(4)  =  DIMAG(WFN(2))
C
      IF ( IFLAG .EQ. 2 )   GO TO 1000
C
C   HERE SET UP STRATIFIED SPHERE
C
      DH2  =  - N / Z(2)  +  1.0 / ( N / Z(2) - DH2 )
      DH4  =  - N / Z(4)  +  1.0 / ( N / Z(4) - DH4 )
      DH1  =  - N / Z(1)  +  1.0 / ( N / Z(1) - DH1 )
      PSTORE  =  ( DH4 + N / Z(4) )  *  ( W(3,N) + N / Z(4) )
      P24H24  =  P24H24 / PSTORE
      HSTORE  =  ( DH1 + N / Z(1) )  *  ( W(3,N) + N / Z(4) )
      P24H21  =  P24H21 / HSTORE
      PSTORE  =  ( ACAP(N) + N / Z(1) )  /  ( W(3,N) + N / Z(4) )
      DUMMY   =  DUMMY * PSTORE
      DUMSQ   =  DUMMY * DUMMY
C
      U(1) =  K3 * ACAP(N)  -  K2 * W(1,N)
      U(2) =  K3 * ACAP(N)  -  K2 * DH2
      U(3) =  K2 * ACAP(N)  -  K3 * W(1,N)
      U(4) =  K2 * ACAP(N)  -  K3 * DH2
      U(5) =  K1 *  W(3,N)  -  K2 * W(2,N)
      U(6) =  K2 *  W(3,N)  -  K1 * W(2,N)
      U(7) =  ( 0.0,-1.0 )  *  ( DUMMY * P24H21 - P24H24 )
      U(8) =  TA(3) / WFN(2)
C
      FNA  =  U(8) * ( U(1)*U(5)*U(7)  +  K1*U(1)  -  DUMSQ*K3*U(5) ) /
     2               ( U(2)*U(5)*U(7)  +  K1*U(2)  -  DUMSQ*K3*U(5) )
      FNB  =  U(8) * ( U(3)*U(6)*U(7)  +  K2*U(3)  -  DUMSQ*K2*U(6) ) /
     2               ( U(4)*U(6)*U(7)  +  K2*U(4)  -  DUMSQ*K2*U(6) )
c
c temporary output to get u functions and a,b coefficients
c
c     print 742,  n
c     print 744,  u(1), u(2)
c     print 745,  u(3), u(4) 
c     print 746,  u(5), u(6) 
c     print 748,  u(7), u(8)
c     print 752,  fna, fnb

  742 format( /, 1x, 'Output for N = ',  i6 )
  744 format(    5x, 'U1,2:   ',  3( 1p2e14.6, 5x) )
  745 format(    5x, 'U3,4:   ',  3( 1p2e14.6, 5x) )
  746 format(    5x, 'U5,6:   ',  3( 1p2e14.6, 5x) )
  748 format(    5x, 'U7,8:   ',  3( 1p2e14.6, 5x) )
  752 format(    5x, 'An, Bn: ',  2( 1p2e14.6, 5x) )

C
 1000 CONTINUE
      TC1  =  ACAP(N) * RRF  +  N * RX
      TC2  =  ACAP(N) * RF   +  N * RX
      FN1  =  ( TC1 * TA(3)  -  TA(1) ) /  ( TC1 * WFN(2) - WFN(1) )
      FN2  =  ( TC2 * TA(3)  -  TA(1) ) /  ( TC2 * WFN(2) - WFN(1) )
      M    =  WVNO * R
      IF ( IFLAG .EQ. 2 )   GO TO 1001
      IF ( N .LT. M )   GO TO 1002
      IF ( CDABS(  ( FN1-FNA ) / FN1  )  .LT.  1.0E-09   .AND.
     1     CDABS(  ( FN2-FNB ) / FN2  )  .LT . 1.0E-09  )    IFLAG = 2
      IF ( IFLAG .EQ. 1 )   GO TO 1002
 1001 FNA  =  FN1
      FNB  =  FN2
C
 1002 CONTINUE
      T(5)  =  N
      T(4)  =  T(1) / ( T(5) * T(2) )
      T(2)  =  (  T(2) * ( T(5) + 1.0 )  ) / T(5)
C
      CTBRQS  =  CTBRQS  +  T(2) * ( TD(1) * TB(1)  +  TD(2) * TB(2)  +
     1                               TE(1) * TC(1)  +  TE(2) * TC(2) )
     2                   +  T(4) * ( TD(1) * TE(1)  +  TD(2) * TE(2) )
      QEXT    =   QEXT  +  T(3) * ( TB(1) + TC(1) )
      T(4)    =  TB(1)**2 + TB(2)**2 + TC(1)**2 + TC(2)**2
      QSCAT   =  QSCAT  +  T(3) * T(4)
C
      RMM     =  -RMM
      QBSR    =  QBSR + T(3)*RMM*(TC(1) - TB(1))
      QBSI    =  QBSI + T(3)*RMM*(TC(2) - TB(2))
C
      T(2)    =  N * (N+1)
      T(1)    =  T(3) / T(2)
      K = (N/2)*2
      DO 80 J = 1,JX
       ELTRMX(1,J,1) = ELTRMX(1,J,1)+T(1)*(TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,1) = ELTRMX(2,J,1)+T(1)*(TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,1) = ELTRMX(3,J,1)+T(1)*(TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,1) = ELTRMX(4,J,1)+T(1)*(TC(2)*PI(3,J)+TB(2)*TAU(3,J))
      IF ( K .EQ. N )   GO TO 75
       ELTRMX(1,J,2) = ELTRMX(1,J,2)+T(1)*(TB(1)*PI(3,J)-TC(1)*TAU(3,J))
       ELTRMX(2,J,2) = ELTRMX(2,J,2)+T(1)*(TB(2)*PI(3,J)-TC(2)*TAU(3,J))
       ELTRMX(3,J,2) = ELTRMX(3,J,2)+T(1)*(TC(1)*PI(3,J)-TB(1)*TAU(3,J))
       ELTRMX(4,J,2) = ELTRMX(4,J,2)+T(1)*(TC(2)*PI(3,J)-TB(2)*TAU(3,J))
      GO TO 80
   75  ELTRMX(1,J,2) =ELTRMX(1,J,2)+T(1)*(-TB(1)*PI(3,J)+TC(1)*TAU(3,J))
       ELTRMX(2,J,2) =ELTRMX(2,J,2)+T(1)*(-TB(2)*PI(3,J)+TC(2)*TAU(3,J))
       ELTRMX(3,J,2) =ELTRMX(3,J,2)+T(1)*(-TC(1)*PI(3,J)+TB(1)*TAU(3,J))
       ELTRMX(4,J,2) =ELTRMX(4,J,2)+T(1)*(-TC(2)*PI(3,J)+TB(2)*TAU(3,J))
   80 continue
C
      IF ( T(4) .LT. 1.0E-10 .AND. N .GE. NSTOP ) GO TO 100
      N = N + 1
      DO 90 J = 1,JX
         PI(1,J)   =   PI(2,J)
         PI(2,J)   =   PI(3,J)
         TAU(1,J)  =  TAU(2,J)
         TAU(2,J)  =  TAU(3,J)
   90 CONTINUE
      FNAP  =  FNA
      FNBP  =  FNB
      IF ( N .LE. NMX2 )   GO TO 65
         WRITE( 6,8 )

         st1 =  ( X**2 ) * ( RFR**2 + RFI**2 )
         write( 6,888 ) x, rfr, rfi, st1, n, nmx2
         write( 6,889 ) t(4)

         IF ( T(4) .LT. 1.0E-10 )   GO TO 100
         STOP 36
  100 DO 120 J = 1,JX
      DO 120 K = 1,2
         DO  115  I= 1,4
         T(I)  =  ELTRMX(I,J,K)
  115    CONTINUE
         ELTRMX(2,J,K)  =      T(1)**2  +  T(2)**2
         ELTRMX(1,J,K)  =      T(3)**2  +  T(4)**2
         ELTRMX(3,J,K)  =  T(1) * T(3)  +  T(2) * T(4)
         ELTRMX(4,J,K)  =  T(2) * T(3)  -  T(4) * T(1)
  120 CONTINUE
      T(1)    =    2.0 * RX**2
C
C QBS IS THE BACK SCATTER CROSS SECTION
C
      PIG   = ACOS(-1.0)
      RXP4  = RX*RX/(4.0*PIG)
      QBS   = RXP4*(QBSR**2 + QBSI**2)
C
      QEXT    =   QEXT * T(1)
      QSCAT   =  QSCAT * T(1)
      CTBRQS  =  2.0 * CTBRQS * T(1)
C
      RETURN
C
    5 FORMAT( 10X, 'THE VALUE OF THE SCATTERING ANGLE IS GREATER THAN
     1 90.0 DEGREES. IT IS ', E15.4 )
    6 FORMAT( // 10X, 'PLEASE READ COMMENTS.' // )
    7 FORMAT( // 10X, 'THE VALUE OF THE ARGUMENT JX IS GREATER THAN IT')
    8 FORMAT( // 10X, 'THE UPPER LIMIT FOR ACAP IS NOT ENOUGH. SUGGEST
     1 GET DETAILED OUTPUT AND MODIFY SUBROUTINE' // )
c
  888 format( 1h0, 'Output values of x, rfr, rfi, st1, n, nmx2:',  /,
     2         1x,  1p4e12.4, 3x, 2i10  )
  889 format( 1h0, 'Convergence criterion = ', 1pe14.5 )
c
      END



