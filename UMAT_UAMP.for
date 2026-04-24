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
      
      OPEN(UNIT=7, FILE='D:\damage\tibia\Rabbit\Damage\inps\debug_output.txt', STATUS='UNKNOWN', ACTION='WRITE')
      WRITE(7,*) 'svars(1):', svars(1), ' svars(2):', svars(2), ' CURRENT U:', displc
      CLOSE(7)

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
      Eset(1,1)=1;Eset(1,2)=10
      Eset(2,1)=11;Eset(2,2)=16
      Eset(3,1)=17;Eset(3,2)=25
      Eset(4,1)=26;Eset(4,2)=35
      Eset(5,1)=36;Eset(5,2)=46
      Eset(6,1)=47;Eset(6,2)=57
      Eset(7,1)=58;Eset(7,2)=69
      Eset(8,1)=70;Eset(8,2)=84
      Eset(9,1)=85;Eset(9,2)=105
      Eset(10,1)=106;Eset(10,2)=127
      Eset(11,1)=128;Eset(11,2)=147
      Eset(12,1)=148;Eset(12,2)=173
      Eset(13,1)=174;Eset(13,2)=192
      Eset(14,1)=193;Eset(14,2)=221
      Eset(15,1)=222;Eset(15,2)=271
      Eset(16,1)=272;Eset(16,2)=380
      Eset(17,1)=381;Eset(17,2)=629
      Eset(18,1)=630;Eset(18,2)=1194
      Eset(19,1)=1195;Eset(19,2)=2190
      Eset(20,1)=2191;Eset(20,2)=3558
      Eset(21,1)=3559;Eset(21,2)=5343
      Eset(22,1)=5344;Eset(22,2)=7591
      Eset(23,1)=7592;Eset(23,2)=9904
      Eset(24,1)=9905;Eset(24,2)=12625
      Eset(25,1)=12626;Eset(25,2)=15689
      Eset(26,1)=15690;Eset(26,2)=16810
      Eset(27,1)=16811;Eset(27,2)=17359
      Eset(28,1)=17360;Eset(28,2)=17887
      Eset(29,1)=17888;Eset(29,2)=18408
      Eset(30,1)=18409;Eset(30,2)=18957
      Eset(31,1)=18958;Eset(31,2)=19559
      Eset(32,1)=19560;Eset(32,2)=20149
      Eset(33,1)=20150;Eset(33,2)=20788
      Eset(34,1)=20789;Eset(34,2)=21461
      Eset(35,1)=21462;Eset(35,2)=22216
      Eset(36,1)=22217;Eset(36,2)=23000
      Eset(37,1)=23001;Eset(37,2)=23817
      Eset(38,1)=23818;Eset(38,2)=24641
      Eset(39,1)=24642;Eset(39,2)=25568
      Eset(40,1)=25569;Eset(40,2)=26583
      Eset(41,1)=26584;Eset(41,2)=27628
      Eset(42,1)=27629;Eset(42,2)=28704
      Eset(43,1)=28705;Eset(43,2)=29785
      Eset(44,1)=29786;Eset(44,2)=30925
      Eset(45,1)=30926;Eset(45,2)=32063
      Eset(46,1)=32064;Eset(46,2)=33291
      Eset(47,1)=33292;Eset(47,2)=34441
      Eset(48,1)=34442;Eset(48,2)=35673
      Eset(49,1)=35674;Eset(49,2)=36945
      Eset(50,1)=36946;Eset(50,2)=38186
      Eset(51,1)=38187;Eset(51,2)=39420
      Eset(52,1)=39421;Eset(52,2)=40682
      Eset(53,1)=40683;Eset(53,2)=41912
      Eset(54,1)=41913;Eset(54,2)=43154
      Eset(55,1)=43155;Eset(55,2)=44435
      Eset(56,1)=44436;Eset(56,2)=45681
      Eset(57,1)=45682;Eset(57,2)=46847
      Eset(58,1)=46848;Eset(58,2)=48000
      Eset(59,1)=48001;Eset(59,2)=49202
      Eset(60,1)=49203;Eset(60,2)=50460
      Eset(61,1)=50461;Eset(61,2)=51659
      Eset(62,1)=51660;Eset(62,2)=52817
      Eset(63,1)=52818;Eset(63,2)=53928
      Eset(64,1)=53929;Eset(64,2)=54987
      Eset(65,1)=54988;Eset(65,2)=55997
      Eset(66,1)=55998;Eset(66,2)=56954
      Eset(67,1)=56955;Eset(67,2)=57864
      Eset(68,1)=57865;Eset(68,2)=58709
      Eset(69,1)=58710;Eset(69,2)=59532
      Eset(70,1)=59533;Eset(70,2)=60312
      Eset(71,1)=60313;Eset(71,2)=61090
      Eset(72,1)=61091;Eset(72,2)=61812
      Eset(73,1)=61813;Eset(73,2)=62530
      Eset(74,1)=62531;Eset(74,2)=63226
      Eset(75,1)=63227;Eset(75,2)=63911
      Eset(76,1)=63912;Eset(76,2)=64554
      Eset(77,1)=64555;Eset(77,2)=65243
      Eset(78,1)=65244;Eset(78,2)=65982
      Eset(79,1)=65983;Eset(79,2)=66670
      Eset(80,1)=66671;Eset(80,2)=67393
      Eset(81,1)=67394;Eset(81,2)=68112
      Eset(82,1)=68113;Eset(82,2)=68831
      Eset(83,1)=68832;Eset(83,2)=69565
      Eset(84,1)=69566;Eset(84,2)=70379
      Eset(85,1)=70380;Eset(85,2)=71197
      Eset(86,1)=71198;Eset(86,2)=72085
      Eset(87,1)=72086;Eset(87,2)=73088
      Eset(88,1)=73089;Eset(88,2)=74202
      Eset(89,1)=74203;Eset(89,2)=75456
      Eset(90,1)=75457;Eset(90,2)=76794
      Eset(91,1)=76795;Eset(91,2)=78260
      Eset(92,1)=78261;Eset(92,2)=79675
      Eset(93,1)=79676;Eset(93,2)=80834
      Eset(94,1)=80835;Eset(94,2)=81720
      Eset(95,1)=81721;Eset(95,2)=82308
      Eset(96,1)=82309;Eset(96,2)=82669
      Eset(97,1)=82670;Eset(97,2)=82856
      Eset(98,1)=82857;Eset(98,2)=82967
      Eset(99,1)=82968;Eset(99,2)=83017
      Eset(100,1)=83018;Eset(100,2)=83032

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
      mat(30)=25.2824
      mat(31)=70.9686
      mat(32)=135.2472
      mat(33)=216.4333
      mat(34)=313.4090
      mat(35)=425.3557
      mat(36)=551.6292
      mat(37)=691.7046
      mat(38)=845.1497
      mat(39)=1011.5705
      mat(40)=1190.6433
      mat(41)=1382.0875
      mat(42)=1585.6214
      mat(43)=1801.0320
      mat(44)=2028.0940
      mat(45)=2266.5993
      mat(46)=2516.3901
      mat(47)=2777.2804
      mat(48)=3049.1279
      mat(49)=3331.7949
      mat(50)=3625.1169
      mat(51)=3928.9875
      mat(52)=4243.2749
      mat(53)=4567.8647
      mat(54)=4902.6834
      mat(55)=5247.5765
      mat(56)=5602.4639
      mat(57)=5967.2749
      mat(58)=6341.8907
      mat(59)=6726.2589
      mat(60)=7120.2555
      mat(61)=7523.8360
      mat(62)=7936.9074
      mat(63)=8359.4465
      mat(64)=8791.2721
      mat(65)=9232.3817
      mat(66)=9682.7238
      mat(67)=10142.3822
      mat(68)=10610.9820
      mat(69)=11088.6910
      mat(70)=11575.2271
      mat(71)=12071.0508
      mat(72)=12575.4205
      mat(73)=13088.7338
      mat(74)=13610.7515
      mat(75)=14141.4894
      mat(76)=14681.0944
      mat(77)=15229.1915
      mat(78)=15785.7942
      mat(79)=16351.0488
      mat(80)=16924.9716
      mat(81)=17507.1778
      mat(82)=18097.6788
      mat(83)=18696.6217
      mat(84)=19304.0207
      mat(85)=19919.7540
      mat(86)=20543.6977
      mat(87)=21175.7261
      mat(88)=21815.9875
      mat(89)=22464.6331
      mat(90)=23121.2595
      mat(91)=23785.8762
      mat(92)=24458.6337
      mat(93)=25139.4029
      mat(94)=25828.3359
      mat(95)=26525.0189
      mat(96)=27229.6024
      mat(97)=27942.2393
      mat(98)=28662.7975
      mat(99)=29391.1431
      mat(100)=30127.1402
      
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

       deltN=100 ! time step
       IF (TIME(1).LT.200) THEN !after 100 cycles change deltN to 1000
            deltN=100
       ELSE
            deltN=100
       END IF
       
       STATEV(3)=e_m
       
       deltD=deltN*(33.70)*(e_m**3.2410)
       
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
