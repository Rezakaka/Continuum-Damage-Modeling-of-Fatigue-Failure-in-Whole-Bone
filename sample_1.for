      SUBROUTINE UAMP(
     *     ampName, time, ampValueOld, dt, nProps, props, nSvars, 
     *     svars, lFlagsInfo,
     *     nSensor, sensorValues, sensorNames, jSensorLookUpTable, 
     *     AmpValueNew, 
     *     lFlagsDefine,
     *     AmpDerivative, AmpSecDerivative, AmpIncIntegral,
     *     AmpDoubleIntegral)
C
      INCLUDE 'ABA_PARAM.INC'

C     time indices
      parameter (iStepTime        = 1,
     *           iTotalTime       = 2,
     *           nTime            = 2)
C     flags passed in for information
      parameter (iInitialization   = 1,
     *           iRegularInc       = 2,
     *           iCuts             = 3,
     *           ikStep            = 4,
     *           nFlagsInfo        = 4)
C     optional flags to be defined
      parameter (iComputeDeriv       = 1,
     *           iComputeSecDeriv    = 2,
     *           iComputeInteg       = 3,
     *           iComputeDoubleInteg = 4,
     *           iStopAnalysis       = 5,
     *           iConcludeStep       = 6,
     *           nFlagsDefine        = 6)
      dimension time(nTime), lFlagsInfo(nFlagsInfo),
     *          lFlagsDefine(nFlagsDefine)
      dimension jSensorLookUpTable(*)
      dimension sensorValues(nSensor), svars(nSvars), props(nProps)
      character*80 sensorNames(nSensor)
      character*80 ampName

      !user coding to define AmpValueNew, and 
      !optionally lFlagsDefine, AmpDerivative, AmpSecDerivative, 
      !AmpIncIntegral, AmpDoubleIntegral
c     get sensor value
      AmpValueNew=0
      load=1260.0
      displc  = GetSensorValue('NODE_TRACKED_HISTORY',
     *                              jSensorLookUpTable,
     *                              sensorValues)
      IF (TIME(2) .eq. 0) THEN 
       threshold=1
      ELSE IF (TIME(2) .eq. 2) THEN 
        svars(1)=abs(load/displc)!initial stiffness value
        threshold=1
      ELSE
        svars(2)=abs(load/displc)! updated stiffness value
        threshold=svars(2)/svars(1)
      End IF

      if (threshold < 0.75 .AND. threshold > 0.0) then 
        lFlagsDefine(iConcludeStep)=1
      end if
C---------------------------------------------------------------------------------------
C End of UAMP code
C---------------------------------------------------------------------------------------
      RETURN
      END
C---------------------------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C

      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
      DIMENSION STATEV(NSTATV),
     1 DRPLDE(NTENS),
     2 TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3),
     4 JSTEP(4)
       
      DOUBLE PRECISION STRESS(NTENS),DDSDDE(NTENS,NTENS),
     1 DDSDDT(NTENS),STRAN(NTENS),DSTRAN(NTENS),strs(3,3),egv(3),
     2 PS(3), AN(3, 3)
      INTEGER el_no,deltN,no_cyc
      REAL Y_m,deltD,e_m,t_0_c
      INTEGER Eset(100, 2)
      REAL mat(100)
      DOUBLE PRECISION t_E,maxps,minps,maxsh,E
      DOUBLE PRECISION stff(6,6)


      ! Initialize the Eset matrix and materials array
      Eset(1,1)=1;Eset(1,2)=78
      Eset(2,1)=79;Eset(2,2)=121
      Eset(3,1)=122;Eset(3,2)=159
      Eset(4,1)=160;Eset(4,2)=198
      Eset(5,1)=199;Eset(5,2)=236
      Eset(6,1)=237;Eset(6,2)=279
      Eset(7,1)=280;Eset(7,2)=311
      Eset(8,1)=312;Eset(8,2)=358
      Eset(9,1)=359;Eset(9,2)=399
      Eset(10,1)=400;Eset(10,2)=444
      Eset(11,1)=445;Eset(11,2)=483
      Eset(12,1)=484;Eset(12,2)=537
      Eset(13,1)=538;Eset(13,2)=596
      Eset(14,1)=597;Eset(14,2)=650
      Eset(15,1)=651;Eset(15,2)=700
      Eset(16,1)=701;Eset(16,2)=782
      Eset(17,1)=783;Eset(17,2)=895
      Eset(18,1)=896;Eset(18,2)=1163
      Eset(19,1)=1164;Eset(19,2)=1642
      Eset(20,1)=1643;Eset(20,2)=2508
      Eset(21,1)=2509;Eset(21,2)=3937
      Eset(22,1)=3938;Eset(22,2)=5823
      Eset(23,1)=5824;Eset(23,2)=8171
      Eset(24,1)=8172;Eset(24,2)=10552
      Eset(25,1)=10553;Eset(25,2)=13587
      Eset(26,1)=13588;Eset(26,2)=15872
      Eset(27,1)=15873;Eset(27,2)=16549
      Eset(28,1)=16550;Eset(28,2)=17024
      Eset(29,1)=17025;Eset(29,2)=17522
      Eset(30,1)=17523;Eset(30,2)=18059
      Eset(31,1)=18060;Eset(31,2)=18541
      Eset(32,1)=18542;Eset(32,2)=19088
      Eset(33,1)=19089;Eset(33,2)=19666
      Eset(34,1)=19667;Eset(34,2)=20248
      Eset(35,1)=20249;Eset(35,2)=20926
      Eset(36,1)=20927;Eset(36,2)=21645
      Eset(37,1)=21646;Eset(37,2)=22370
      Eset(38,1)=22371;Eset(38,2)=23199
      Eset(39,1)=23200;Eset(39,2)=24168
      Eset(40,1)=24169;Eset(40,2)=25212
      Eset(41,1)=25213;Eset(41,2)=26308
      Eset(42,1)=26309;Eset(42,2)=27431
      Eset(43,1)=27432;Eset(43,2)=28621
      Eset(44,1)=28622;Eset(44,2)=29815
      Eset(45,1)=29816;Eset(45,2)=30993
      Eset(46,1)=30994;Eset(46,2)=32216
      Eset(47,1)=32217;Eset(47,2)=33450
      Eset(48,1)=33451;Eset(48,2)=34684
      Eset(49,1)=34685;Eset(49,2)=35947
      Eset(50,1)=35948;Eset(50,2)=37144
      Eset(51,1)=37145;Eset(51,2)=38320
      Eset(52,1)=38321;Eset(52,2)=39534
      Eset(53,1)=39535;Eset(53,2)=40715
      Eset(54,1)=40716;Eset(54,2)=41867
      Eset(55,1)=41868;Eset(55,2)=43054
      Eset(56,1)=43055;Eset(56,2)=44231
      Eset(57,1)=44232;Eset(57,2)=45375
      Eset(58,1)=45376;Eset(58,2)=46484
      Eset(59,1)=46485;Eset(59,2)=47543
      Eset(60,1)=47544;Eset(60,2)=48580
      Eset(61,1)=48581;Eset(61,2)=49682
      Eset(62,1)=49683;Eset(62,2)=50766
      Eset(63,1)=50767;Eset(63,2)=51761
      Eset(64,1)=51762;Eset(64,2)=52782
      Eset(65,1)=52783;Eset(65,2)=53829
      Eset(66,1)=53830;Eset(66,2)=54783
      Eset(67,1)=54784;Eset(67,2)=55771
      Eset(68,1)=55772;Eset(68,2)=56734
      Eset(69,1)=56735;Eset(69,2)=57708
      Eset(70,1)=57709;Eset(70,2)=58648
      Eset(71,1)=58649;Eset(71,2)=59454
      Eset(72,1)=59455;Eset(72,2)=60306
      Eset(73,1)=60307;Eset(73,2)=61069
      Eset(74,1)=61070;Eset(74,2)=61812
      Eset(75,1)=61813;Eset(75,2)=62476
      Eset(76,1)=62477;Eset(76,2)=63124
      Eset(77,1)=63125;Eset(77,2)=63791
      Eset(78,1)=63792;Eset(78,2)=64404
      Eset(79,1)=64405;Eset(79,2)=65022
      Eset(80,1)=65023;Eset(80,2)=65670
      Eset(81,1)=65671;Eset(81,2)=66249
      Eset(82,1)=66250;Eset(82,2)=66915
      Eset(83,1)=66916;Eset(83,2)=67678
      Eset(84,1)=67679;Eset(84,2)=68552
      Eset(85,1)=68553;Eset(85,2)=69663
      Eset(86,1)=69664;Eset(86,2)=70856
      Eset(87,1)=70857;Eset(87,2)=72021
      Eset(88,1)=72022;Eset(88,2)=73273
      Eset(89,1)=73274;Eset(89,2)=74511
      Eset(90,1)=74512;Eset(90,2)=75628
      Eset(91,1)=75629;Eset(91,2)=76717
      Eset(92,1)=76718;Eset(92,2)=77794
      Eset(93,1)=77795;Eset(93,2)=78638
      Eset(94,1)=78639;Eset(94,2)=79370
      Eset(95,1)=79371;Eset(95,2)=79924
      Eset(96,1)=79925;Eset(96,2)=80315
      Eset(97,1)=80316;Eset(97,2)=80522
      Eset(98,1)=80523;Eset(98,2)=80600
      Eset(99,1)=80601;Eset(99,2)=80626
      Eset(100,1)=80627;Eset(100,2)=80637

      mat(1)=2.5630
      mat(2)=2.5630
      mat(3)=2.5630
      mat(4)=2.5630
      mat(5)=2.5630
      mat(6)=2.5630
      mat(7)=2.5630
      mat(8)=2.5630
      mat(9)=2.5630
      mat(10)=2.5630
      mat(11)=2.5630
      mat(12)=2.5630
      mat(13)=2.5630
      mat(14)=2.5630
      mat(15)=2.5630
      mat(16)=2.5630
      mat(17)=2.5630
      mat(18)=2.5630
      mat(19)=2.5630
      mat(20)=2.5630
      mat(21)=2.5630
      mat(22)=2.5630
      mat(23)=2.5630
      mat(24)=2.5630
      mat(25)=2.5630
      mat(26)=2.5630
      mat(27)=2.5630
      mat(28)=2.5630
      mat(29)=2.5630
      mat(30)=4.4735
      mat(31)=35.0939
      mat(32)=88.0744
      mat(33)=160.4725
      mat(34)=250.6528
      mat(35)=357.4994
      mat(36)=480.1862
      mat(37)=618.0541
      mat(38)=770.5592
      mat(39)=937.2493
      mat(40)=1117.7335
      mat(41)=1311.6626
      mat(42)=1518.7360
      mat(43)=1738.6901
      mat(44)=1971.2641
      mat(45)=2216.2417
      mat(46)=2473.4047
      mat(47)=2742.5755
      mat(48)=3023.5784
      mat(49)=3316.2509
      mat(50)=3620.4131
      mat(51)=3935.9641
      mat(52)=4262.7449
      mat(53)=4600.6237
      mat(54)=4949.4976
      mat(55)=5309.2382
      mat(56)=5679.7466
      mat(57)=6060.8985
      mat(58)=6452.6350
      mat(59)=6854.8259
      mat(60)=7267.3961
      mat(61)=7690.2451
      mat(62)=8123.3053
      mat(63)=8566.4960
      mat(64)=9019.7220
      mat(65)=9482.9351
      mat(66)=9956.0614
      mat(67)=10439.0130
      mat(68)=10931.9378
      mat(69)=11434.3831
      mat(70)=11946.4684
      mat(71)=12468.0868
      mat(72)=12999.2574
      mat(73)=13540.0004
      mat(74)=14090.0768
      mat(75)=14649.5037
      mat(76)=15218.1675
      mat(77)=15796.2164
      mat(78)=16383.4041
      mat(79)=16979.8793
      mat(80)=17585.3919
      mat(81)=18199.9564
      mat(82)=18823.7230
      mat(83)=19456.3014
      mat(84)=20097.8402
      mat(85)=20748.2166
      mat(86)=21407.4437
      mat(87)=22075.5349
      mat(88)=22752.5040
      mat(89)=23438.0860
      mat(90)=24132.4321
      mat(91)=24835.4146
      mat(92)=25546.9041
      mat(93)=26267.1945
      mat(94)=26996.0151
      mat(95)=27733.5198
      mat(96)=28479.2910
      mat(97)=29233.6251
      mat(98)=29996.3894
      mat(99)=30767.4492
      mat(100)=31546.9594

C       Get the Young's modulus for a specific element
        el_no = noel
        CALL get_Y(el_no, Eset, mat, Y_m)
        
        IF(TIME(2) .eq. 0)THEN   !KSTEP: Step number      
          STATEV(1)=Y_m
          E=STATEV(1)
          STATEV(2)=0 !Total damage is zero during first step. (damage is updated at the end of each step)
        ELSE 
          E = STATEV(1) !use the updated Youngs modulus to calculate stress
        END IF


C       Von Mises strain incriment: e_m       
       e_m = SQRT(0.5 * ((STRAN(1)-STRAN(2))**2+
     1   (STRAN(2)-STRAN(3))**2+(STRAN(3)-STRAN(1))**2) +
     2   3.0 *(STRAN(4)**2 + STRAN(5)**2 + STRAN(6)**2))

       deltN=10000 ! time step
       IF (TIME(1).LT.200) THEN !after 100 cycles change deltN to 1000
            deltN=10000
       ELSE
            deltN=10000
       END IF
       
       STATEV(3)=e_m
       
       deltD=deltN*(10.950)*(e_m**3.2410)
       
       STATEV(4)=deltD
       STATEV(2)=deltD+STATEV(2) !update the total damage
       STATEV(1)=Y_m*(1-STATEV(2)) ! Y_m: undamaged Youngs modulus or E_0, STATEV(1): updated Youngs modulus
       E=STATEV(1)

C       Do not allow Youngs modulus of less than t_E (tested with 1, 0.01, 0, 1e-6: no sensitivity was seen)
       t_E=1
       IF(E .LE. t_E)THEN
        STATEV(1)=t_E
        E=t_E
       END IF
       STATEV(5)=1-E/Y_m !damage
       
C     Stiffness matrix       
C     Stiffness matrix       
       CALL ANISOM(E, stff)
       DDSDDE=stff

       DO j=1, 6
             STATEV(j+8)=DDSDDE(j,j)
       END DO

C      Calculation of STRESSES 
C
C      STRESS(NTENS)
C      This array is passed in as the stress tensor at the beginning of the increment and must be updated
C      in this routine to be the stress tensor at the end of the increment.
C      
C      DSTRAN(NTENS)
C      Array of strain increments.
      t_0_c=TIME(1)
      no_cyc=INT(MOD(t_0_c,2.0))
      t=12
      DO j=1, NTENS
            !STATEV(t)=DSTRAN(j)
            !STATEV(t+6)=STRESS(j)
            !STATEV(t+12)=STRAN(j)
            STRESS(j)=STRESS(j)+DDSDDE(j, 1)*DSTRAN(1)+
     1          DDSDDE(j, 2)*DSTRAN(2)+
     2          DDSDDE(j, 3)*DSTRAN(3)+DDSDDE(j, 4)*DSTRAN(4)+
     3          DDSDDE(j, 5)*DSTRAN(5)+DDSDDE(j, 6)*DSTRAN(6)
            t=t+1
      END DO

C---------------------------------------------------------------------------------------
C End of UMAT code
C---------------------------------------------------------------------------------------
      RETURN
      END

      SUBROUTINE get_Y(el_no, Eset, mat, Y_m)
      INTEGER el_no
      INTEGER Eset(100, 2)  ! Eset matrix (100x2)
      REAL mat(100)   ! Materials array (100x1)
      REAL Y_m   ! The Young's modulus to return
      INTEGER i

      ! Initialize young_modulus to 0 or some error value if not found
      Y_m = 0.0
      ! Search through each bin to find the range that contains the element
      DO 10 i = 1, 100
          IF (el_no .GE. Eset(i, 1) .AND. el_no .LE. Eset(i, 2)) THEN
              ! The element is in the current bin, return the corresponding Young's modulus
              Y_m = mat(i)
              RETURN
          END IF
  10  CONTINUE

      RETURN
      END SUBROUTINE get_Y

      SUBROUTINE ANISOM(E3, C)

      IMPLICIT NONE

      DOUBLE PRECISION E3
      DOUBLE PRECISION C(6,6)

      DOUBLE PRECISION E1, E2
      DOUBLE PRECISION G12, G23, G31
      DOUBLE PRECISION NU12, NU23, NU31
      DOUBLE PRECISION NU21, NU32, NU13
      DOUBLE PRECISION S(6,6)
      INTEGER I, J

      E1   = 0.574D0 * E3
      E2   = 0.577D0 * E3
      G12  = 0.195D0 * E3
      G23  = 0.265D0 * E3
      G31  = 0.216D0 * E3
      NU12 = 0.427D0
      NU23 = 0.234D0
      NU31 = 0.40D0

      NU21 = NU12 * (E2 / E1)
      NU32 = NU23 * (E3 / E2)
      NU13 = NU31 * (E1 / E3)

      DO I = 1,6
      DO J = 1,6
            S(I,J) = 0.0D0
      END DO
      END DO

      S(1,1) = 1.0D0 / E1
      S(2,2) = 1.0D0 / E2
      S(3,3) = 1.0D0 / E3

      S(1,2) = -NU12 / E1
      S(2,1) = -NU21 / E2
      S(2,3) = -NU23 / E2
      S(3,2) = -NU32 / E3
      S(3,1) = -NU31 / E3
      S(1,3) = -NU13 / E1
      S(4,4) = 1.0D0 / G23
      S(5,5) = 1.0D0 / G31
      S(6,6) = 1.0D0 / G12

      CALL INV6X6(S, C)

      RETURN
      END SUBROUTINE ANISOM


      SUBROUTINE INV6X6(A, AINV)

      IMPLICIT NONE

      DOUBLE PRECISION A(6,6), AINV(6,6)
      DOUBLE PRECISION WORK(6,6)
      DOUBLE PRECISION PIVOT, FACTOR
      INTEGER I, J, K, PIVOTR

      DO I = 1,6
      DO J = 1,6
            WORK(I,J) = A(I,J)
            AINV(I,J) = 0.0D0
      END DO
      AINV(I,I) = 1.0D0
      END DO

      DO K = 1, 6

      PIVOT = DABS(WORK(K,K))
      PIVOTR = K
      DO I = K+1, 6
            IF(DABS(WORK(I,K)) .GT. PIVOT) THEN
            PIVOT = DABS(WORK(I,K))
            PIVOTR = I
            END IF
      END DO

      IF(PIVOTR .NE. K) THEN
            DO J = 1,6
            FACTOR         = WORK(K,J)
            WORK(K,J)      = WORK(PIVOTR,J)
            WORK(PIVOTR,J) = FACTOR

            FACTOR         = AINV(K,J)
            AINV(K,J)      = AINV(PIVOTR,J)
            AINV(PIVOTR,J) = FACTOR
            END DO
      END IF

      PIVOT = WORK(K,K)
      IF(DABS(PIVOT) .LT. 1.0D-15) THEN
            WRITE(*,*) 'ERROR: Matrix singular or near singular in INV6X6'
            STOP
      END IF

      DO J = 1,6
            WORK(K,J) = WORK(K,J) / PIVOT
            AINV(K,J) = AINV(K,J) / PIVOT
      END DO

      DO I = 1,6
            IF(I .NE. K) THEN
            FACTOR = WORK(I,K)
            DO J = 1,6
                  WORK(I,J) = WORK(I,J) - FACTOR*WORK(K,J)
                  AINV(I,J) = AINV(I,J) - FACTOR*AINV(K,J)
            END DO
            END IF
      END DO

      END DO

      RETURN
      END SUBROUTINE INV6X6
