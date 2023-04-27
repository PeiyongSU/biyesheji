C
C       This program was designed to obtain EHL solutions of
C       elliptic contact problems by means of multi-grid method,
C       including MLMI to evaluate elastic deformation. 
C三维数组可看做是若干张二维表数据
C       This version is for 2005 EHL Workshop at Qingdao Tech U.
C
        PROGRAM PEHL2
        PARAMETER(M=3,NA=256,NB=256,KD=1,PI=3.1415926,NT1=100,NT2=100)
        DIMENSION X(M,0:NA),Y(M,0:NB),DX(M),DY(M),NX(M),NY(M),
     &            P(M,0:NA,0:NB),H(M,0:NA,0:NB),
     &            FP(M,0:NA,0:NB),FG(M),XS(0:NA),YS(0:NB),
     &            PS(0:NA,0:NB),HS(0:NA,0:NB),FPS(0:NA,0:NB),
     &            PB(M,0:NA,0:NB),POLD(0:NA,0:NB),MC(M),WW(M),
     &            RES(M,0:NA,0:NB),VISC(0:NA,0:NB),DENS(0:NA,0:NB),
     &            EPS(0:NA,0:NB),RH1(0:NA,0:NB),RESD(0:NA,0:NB),
     &            DB(M,0:NA,0:NB),DBS(0:NA,0:NB)
	  DIMENSION XX(0:NT2),YY(0:NT2),HMIN(0:NT2),HMINY(0:NT2),
     &            PMAX(0:NT2),TIMES(0:NT2),H0T(0:NT2),DEFOR(0:NA,0:NB)
	  
	  COMPLEX(8)::UDB(M,0:2*NA,0:2*NB)
	  COMPLEX(8)::YYD(0:2*NA)
C
        OPEN(15,FILE='PEHL2.DAT',STATUS='UNKNOWN')
	  INPUT=0
	  IF(INPUT.EQ.1) THEN
	  OPEN(20,FILE='DATA1.D',STATUS='OLD',FORM='UNFORMATTED')
	  READ(20) P,H00
	  CLOSE(20)
	  END IF
C
        OMIG1=0.3
        OMIG2=0.05
	  FSP=5.0
        MGAMA=2
        ISW=0
        MV1=2
        MV2=1
        MV0=19
	  RD=1.0
C
        R1=0.014
	  R2=0.01405
        RX=R1*R2/(R2-R1)
	  RY=RX*RD
	  R0=RX*RY/(RX+RY)
C	  PRINT *,RX,RY,R0
C	  PAUSE
	  RK=1.0
        XIN=-PI/2
        XOUT=PI/2
	  YIN=-PI/2
        YOUT=PI/2
C	 
	  Ue=1.0
C	  E1=2.1E11
C	  V1=0.3
	  E1=1.0E9
	  V1=0.4
	  E0=2.0/((1-V1**2)/E1)
        ALFA=2.19E-8
	  G0=ALFA*E0
        VIS0=0.01
        Z0=ALFA/(LOG(VIS0)+9.67)/5.1E-9
C
	  U0=VIS0*Ue/E0/RX
	  write(*,*) 'U0=',U0
	  W=2000.0
	  W0=W/E0/RX/RX
	  WRITE(*,*) 'W=',W
        F2=(1.0+(PI/2.0-1.0)*RD**(-1.0238))*(1.0+0.0486*RD**
     &      (-1.3358)*(LOG(RD))**1.0997)
	  B=(6.0*RK**2*F2*R0*W/PI/E0)**(1.0/3.0)
	  A=B/RK
	  PHERTZ=1.5*W/PI/A/B
	  AL=8.0*PI*RX*RK*U0/W0/A
	  AL1=2.0*PHERTZ*RX/PI/E0/A
C
        CALL MLGRID(X,Y,DX,DY,NX,NY,XIN,XOUT,YIN,YOUT,M,NA,NB)
        
        CALL COEFFI(DB,DX,DY,NX,NY,M,NA,NB,KD)
C
C	  IF(INPUT.EQ.0) THEN
C        CALL PMGIVE(P,X,Y,RK,M,NA,NB)
C	  END IF
        FFT_FORWARD=1
	  DO 80 K=1,M
	  DO 15 I=0,NX(K)
	  DO 10 J=0,NY(K)
	  UDB(K,I,J)=DB(K,I,J)
   10   CONTINUE
   15   CONTINUE
        ND=2*NX(K)
	  NE=2*NY(K)
	  DO 25 I=0,NX(K)
	  DO 20 J=NY(K)+1,NE-1
	  UDB(K,I,J)=UDB(K,I,NE-J)
   20   CONTINUE
   25   CONTINUE
        DO 35 J=0,NE-1
	  DO 30 I=NX(K)+1,ND-1
	  UDB(K,I,J)=UDB(K,ND-I,J)
   30   CONTINUE
   35   CONTINUE
        DO 50 I=0,ND-1
	  DO 40 J=0,NE-1
	  YYD(J)=UDB(K,I,J)
   40   CONTINUE
        CALL FCFFT(YYD,PI,NE,FFT_FORWARD)
	  DO 45 J=0,NE-1
	  UDB(K,I,J)=YYD(J)
   45   CONTINUE
   50   CONTINUE
        DO 65 J=0,NE-1
	  DO 55 I=0,ND-1
	  YYD(I)=UDB(K,I,J)
   55   CONTINUE
        CALL FCFFT(YYD,PI,ND,FFT_FORWARD)
	  DO 60 I=0,ND-1
	  UDB(K,I,J)=YYD(I)
   60   CONTINUE
   65   CONTINUE
   80   CONTINUE
	  H00=0.0

        DO 500 NT=0,NT2
	  PRINT *,'********************TIME STEP=',NT,'*****************'
C       有量纲时间T=1S
	  T=1.0
C       无量纲时间T0
	  T0=T*ue/a
	  PRINT *,'T0=',T0
C       有量纲时间间隔
        DT=T/NT1
	  TIMES(NT)=NT*DT
C       无量纲时间间隔
	  DT0=T0/NT1
	  IF(NT.EQ.0) THEN
	  DO 110 I=0,NA
	  DO 100 J=0,NB
	  H(M,I,J)=(1.0E-5)/A/A*RX+0.5*X(M,I)**2+0.5*Y(M,J)**2
  100   CONTINUE
  110   CONTINUE
	  GOTO 210
	  ENDIF
	  
	  DO 130 I=0,NX(M)
        DO 120 J=0,NY(M)
        FP(M,I,J)=-H(M,I,J)/DT0
  120   CONTINUE
  130   CONTINUE
        FG(M)=2.0*PI*RK/3.0
C
        DO 200 NUMBER=1,30
C
        DO 160 I=0,NX(M)
        DO 150 J=0,NY(M)
        POLD(I,J)=P(M,I,J)
  150   CONTINUE
  160   CONTINUE
        CALL CIRCLE(P,H,X,Y,FP,FG,DX,DY,NX,NY,MC,RES,WW,
     &       PB,PS,HS,XS,YS,FPS,VISC,DENS,EPS,RH1,
     &       RESD,H00,MV1,MV2,MV0,RD,AL,AL1,A1,A2,A3,A4,Z0,
     &       M,NA,NB,MGAMA,OMIG1,OMIG2,DB,DBS,
     &       KD,ISW,PI,DT0,UDB)
C
        ERRP=0.0
        SUMP=0.0
        DO 180 I=0,NX(M)
        DO 170 J=0,NY(M)
        ERRP=ERRP+ABS(P(M,I,J)-POLD(I,J))
        SUMP=SUMP+P(M,I,J)
  170   CONTINUE
  180   CONTINUE
        ERRP=ERRP/SUMP
        WR=0.0
        DO 182 I=1,NX(M)-1,2
        DO 181 J=1,NY(M)-1,2
        WR=WR+16.0*P(M,I,J)
     &    +4.0*(P(M,I-1,J)+P(M,I+1,J)+P(M,I,J-1)+P(M,I,J+1))
     &    +P(M,I-1,J-1)+P(M,I+1,J-1)+P(M,I-1,J+1)+P(M,I+1,J+1)
  181   CONTINUE
  182   CONTINUE
        WR=WR*DX(M)*DY(M)/9.0
        WT=2.0*PI*RK/3.0
        ERRW=(WR-WT)/WT
C
        IF(ERRP.GT.5.0E-3.AND.ISW.EQ.0) THEN
        ISW=0
        ELSE
        ISW=1
        END IF
C
        IF(ISW.EQ.1) H00=H00+FSP*OMIG2*ERRW
C
        WRITE(15,185) NUMBER,ERRP,ABS(ERRW)
  185   FORMAT(1X,' NUMBER=',I3,1X,'ERRP=',F10.6,1X,'ERRW=',F10.6)
        WRITE(*,*) 'ERRP=',ERRP,' ERRW=',ERRW, ' No.',NUMBER
C
        WRITE(*,*) '*********** OK ***********'
        IF((ERRP.LT.0.001).AND.(ABS(ERRW).LT.0.001)) GOTO 210
  200   CONTINUE
C
  210   XX=0.0
        YY=0.0
        XHY=0.0
	  ICY=NB/2
        HMIN(NT)=H(M,0,0)
        HMINY(NT)=H(M,0,ICY)
	  PMAX(NT)=P(M,1,1)
	  DO 260 I=0,NA
        IF(HMINY(NT).GT.H(M,I,ICY)) THEN
        HMINY(NT)=H(M,I,ICY)
        XHY=X(M,I)
        END IF
	  DO 250 J=0,NB
        IF(HMIN(NT).GT.H(M,I,J)) THEN
        HMIN(NT)=H(M,I,J)
        XX(NT)=X(M,I)
        YY(NT)=Y(M,J)
        END IF
	  if(PMAX(NT).lt.P(M,i,j)) then
	  PMAX(NT)=P(M,i,j)
	  end if
  250   CONTINUE
  260   CONTINUE
	  H0T(NT)=H00
C	  OPEN(25,FILE='DATA.D',STATUS='UNKNOWN',FORM='UNFORMATTED')
C	  WRITE(25) P,H00
C	  CLOSE(25)	 
	  DO 400 I=0,NA
	  DO 300 J=0,NB
	  DEFOR(I,J)=(H(M,I,J)-H00-0.5*(X(M,I)**2-Y(M,J)**2))*A*A/Rx
  300   CONTINUE
  400   CONTINUE 
	  IF(NT.EQ.1) THEN
        CALL RESULT0(X,Y,P,H,M,NA,NB,XIN,XOUT,YIN,YOUT,PH,RK,RX,A,ICY,
     &  PHERTZ,DEFOR,PI)
	  ELSEIF(NT.EQ.25) THEN
        CALL RESULT1(X,Y,P,H,M,NA,NB,XIN,XOUT,YIN,YOUT,PH,RK,RX,A,ICY,
     &  PHERTZ,DEFOR,PI)
	  ELSEIF(NT.EQ.50) THEN
        CALL RESULT2(X,Y,P,H,M,NA,NB,XIN,XOUT,YIN,YOUT,PH,RK,RX,A,ICY,
     &  PHERTZ,DEFOR,PI)
	  ELSEIF(NT.EQ.100) THEN
	  CALLRESULT3(X,Y,P,H,M,NA,NB,XIN,XOUT,YIN,YOUT,PH,RK,RX,A,ICY,
     &  PHERTZ,DEFOR,PI)
	  ENDIF
  500   CONTINUE
C
        OPEN(16,FILE='HPIA.DAT')
	  DO 610 K=0,NT2
	  WRITE(16,650) TIMES(K),HMIN(K)*A*A/RX*1.0E6,HMINY(K), !μm
     &                PMAX(K),H0T(K)      !MPa
c*A*A/RX   *PHERTZ*1.0E-6 *A*A/RX*1.0E6
  610   CONTINUE
  650   FORMAT(5(1X,E12.5))
	  CLOSE(15)
        END
C
C
        SUBROUTINE CIRCLE(P,H,X,Y,FP,FG,DX,DY,NX,NY,MC,RES,
     &             WW,PB,PS,HS,XS,YS,FPS,VISC,DENS,
     &             EPS,RH1,RESD,H00,MV1,MV2,MV0,RD,AL,AL1,
     &             A1,A2,A3,A4,Z0,M,NA,NB,MGAMA,OMIG1,OMIG2,
     &             DB,DBS,KD,ISW,PI,DT0,UDB)
        DIMENSION P(M,0:NA,0:NB),H(M,0:NA,0:NB),X(M,0:NA),Y(M,0:NB),
     &            FP(M,0:NA,0:NB),FG(M),DX(M),DY(M),
     &            NX(M),NY(M),MC(M),RES(M,0:NA,0:NB),
     &            WW(M),PB(M,0:NA,0:NB),PS(0:NA,0:NB),HS(0:NA,0:NB),
     &            XS(0:NA),YS(0:NB),FPS(0:NA,0:NB),DENS(0:NA,0:NB),
     &            EPS(0:NA,0:NB),RH1(0:NA,0:NB),RESD(0:NA,0:NB),
     &            VISC(0:NA,0:NB),DB(M,0:NA,0:NB),DBS(0:NA,0:NB)
	  COMPLEX(8)::UDB(M,0:2*NA,0:2*NB)
        DO 10 K=1,M
        MC(K)=0
   10   CONTINUE
        K=M
   20   NXS=NX(K)
        NYS=NY(K)
        KEY=1
        CALL SOLVER(P,H,X,Y,FP,FG,DX,DY,RES,WW,HS,PS,XS,YS,
     &       FPS,VISC,DENS,EPS,RH1,RESD,K,H00,NXS,
     &       NYS,MV1,RD,AL,AL1,A1,A2,A3,A4,Z0,OMIG1,OMIG2,M,
     &       NA,NB,KEY,DB,DBS,KD,ISW,PI,DT0,UDB)
        MC(K)=MC(K)+1
        CALL RESTR(P,M,NA,NB,K,NX,NY)
        CALL RESTR(RES,M,NA,NB,K,NX,NY)
        DO 40 I=0,NX(K-1)
        DO 30 J=0,NY(K-1)
        PB(K-1,I,J)=P(K-1,I,J)
   30   CONTINUE
   40   CONTINUE
        NXS1=NX(K-1)
        NYS1=NY(K-1)
        CALL FUNCAL(FP,FG,X,Y,DX,DY,P,H,RES,WW,XS,YS,PS,
     &       HS,VISC,DENS,EPS,RH1,FPS,K,NXS1,NYS1,H00,
     &       RD,AL,AL1,A1,A2,A3,A4,Z0,M,NA,NB,DB,DBS,
     &       KD,PI,DT0,UDB)
        K=K-1
        IF(K.GT.1) GOTO 20
        KEY=0
        NXS=NX(K)
        NYS=NY(K)
        CALL SOLVER(P,H,X,Y,FP,FG,DX,DY,RES,WW,HS,PS,XS,YS,
     &       FPS,VISC,DENS,EPS,RH1,RESD,K,H00,NXS,
     &       NYS,MV0,RD,AL,AL1,A1,A2,A3,A4,Z0,OMIG1,OMIG2,M,
     &       NA,NB,KEY,DB,DBS,KD,ISW,PI,DT0,UDB)
   50   CALL EXTEND(P,PB,NX,NY,RES,K,M,NA,NB)
        K=K+1
        NXS=NX(K)
        NYS=NY(K)
        CALL SOLVER(P,H,X,Y,FP,FG,DX,DY,RES,WW,HS,PS,XS,YS,
     &       FPS,VISC,DENS,EPS,RH1,RESD,K,H00,NXS,
     &       NYS,MV2,RD,AL,AL1,A1,A2,A3,A4,Z0,OMIG1,OMIG2,M,
     &       NA,NB,KEY,DB,DBS,KD,ISW,PI,DT0,UDB)
        IF(K.EQ.M) THEN
        WRITE(*,*) 'CYCLE ENDED'
        GOTO 60
        ELSE IF(MC(K).EQ.MGAMA) THEN
        MC(K)=0
        GOTO 50
        ELSE
        GOTO 20
        END IF
   60   RETURN
        END
C
C
        SUBROUTINE SOLVER(P,H,X,Y,FP,FG,DX,DY,RES,WW,HS,PS,
     &             XS,YS,FPS,VISC,DENS,EPS,RH1,RESD,
     &             K,H00,NXS,NYS,MV,RD,AL,AL1,A1,A2,A3,A4,Z0,
     &             OMIG1,OMIG2,M,NA,NB,KEY,DB,DBS,KD,ISW,PI,DT0,UDB)
        DIMENSION P(M,0:NA,0:NB),H(M,0:NA,0:NB),X(M,0:NA),Y(M,0:NB),
     &            FP(M,0:NA,0:NB),FG(M),DX(M),DY(M),
     &            RES(M,0:NA,0:NB),WW(M),
     &            HS(0:NA,0:NB),PS(0:NA,0:NB),XS(0:NA),YS(0:NB),
     &            FPS(0:NA,0:NB),VISC(0:NA,0:NB),
     &            DENS(0:NA,0:NB),EPS(0:NA,0:NB),RH1(0:NA,0:NB),
     &            RESD(0:NA,0:NB),DB(M,0:NA,0:NB),DBS(0:NA,0:NB)
	  COMPLEX(8)::UDB(M,0:2*NA,0:2*NB)
        DO 10 I=0,NXS
        XS(I)=X(K,I)
   10   CONTINUE
        DO 20 J=0,NYS
        YS(J)=Y(K,J)
   20   CONTINUE
        DO 50 I=0,NXS
        DO 30 J=0,NYS
        PS(I,J)=P(K,I,J)
        FPS(I,J)=FP(K,I,J)
   30   CONTINUE
   50   CONTINUE
        DO 55 I=0,NXS
        DO 54 J=0,NYS
        DBS(I,J)=DB(K,I,J)
   54   CONTINUE
   55   CONTINUE
        DXS=DX(K)
        DXS2=DXS**2
        DYS=DY(K)
        DYS2=DYS**2
        D0=DB(K,0,0)
        D1=DB(K,1,0)
C
        DO 150 MCC=1,MV
C
        CALL FILMTH(HS,PS,XS,YS,DBS,H00,RD,AL1,NA,NB,NXS,
     &       NYS,M,K,KD,PI,DB,UDB)
        DO 70 I=0,NXS
        DO 60 J=0,NYS
        EPS(I,J)=HS(I,J)**3/AL
        RH1(I,J)=HS(I,J)
   60   CONTINUE
   70   CONTINUE
        DO 80 I=0,NXS
	  PS(I,0)=0.0
        PS(I,NYS)=0.0
   80   CONTINUE
        DO 100 J=1,NYS-1
        PS(0,J)=0.0
        D00=D0
        D10=D1
        DO 90 I=1,NXS-1
        EPSXM=0.5*(EPS(I-1,J)+EPS(I,J))
        EPSXP=0.5*(EPS(I+1,J)+EPS(I,J))
        EPSXI=EPSXM+EPSXP
        EPSYM=0.5*(EPS(I,J-1)+EPS(I,J))
        EPSYP=0.5*(EPS(I,J+1)+EPS(I,J))
        EPSYI=EPSYM+EPSYP
        RESS=FPS(I,J)
     &       -(EPSXM*PS(I-1,J)-EPSXI*PS(I,J)+EPSXP*PS(I+1,J))/DXS2
     &       -(EPSYM*PS(I,J-1)-EPSYI*PS(I,J)+EPSYP*PS(I,J+1))/DYS2
     &       +RH1(I,J)/DT0
        DLDP=-EPSXI/DXS2-EPSYI/DYS2-AL1*DB(K,0,0)/DT0
C-AL1*D00/DXS
        DELT=RESS/DLDP
        PS(I,J)=PS(I,J)+OMIG1*DELT
        IF(PS(I,J).LT.0.0) PS(I,J)=0.0
   90   CONTINUE
        PS(NXS,J)=0.0
  100   CONTINUE
C
        IF(K.EQ.1.AND.MOD(MCC,5).EQ.0.AND.ISW.EQ.0) THEN
C
        
	  WK=0.0
        DO 130 I=1,NXS-1,2
        DO 120 J=1,NYS-1,2
        WK=WK+16.0*PS(I,J)
     &       +4.0*(PS(I-1,J)+PS(I+1,J)+PS(I,J-1)+PS(I,J+1))
     &       +PS(I-1,J-1)+PS(I+1,J-1)+PS(I-1,J+1)+PS(I+1,J+1)
  120   CONTINUE
  130   CONTINUE
        WK=WK*DXS*DYS/9.0
        ERRW=WK-FG(K)
        H00=H00+OMIG2*ERRW
        WRITE(*,*) 'ERRW=',ERRW,' H00=',H00
        END IF
  150   CONTINUE
        IF((K.EQ.1).OR.(KEY.EQ.0)) GOTO 200
        CALL RESCAL(RESD,PS,HS,XS,YS,FPS,VISC,DENS,EPS,
     &       RH1,H00,NXS,NYS,DXS,DYS,RD,AL,AL1,A1,A2,A3,
     &       A4,Z0,NA,NB,DBS,M,K,KD,PI,DB,DT0,UDB)
        DO 170 I=0,NXS
        DO 160 J=0,NYS
        RES(K,I,J)=RESD(I,J)
  160   CONTINUE
  170   CONTINUE
        WS=0.0
        DO 190 I=1,NXS-1,2
        DO 180 J=1,NYS-1,2
        WS=WS+16.0*PS(I,J)
     &       +4.0*(PS(I-1,J)+PS(I+1,J)+PS(I,J-1)+PS(I,J+1))
     &       +PS(I-1,J-1)+PS(I+1,J-1)+PS(I-1,J+1)+PS(I+1,J+1)
  180   CONTINUE
  190   CONTINUE
        WS=WS*DXS*DYS/9.0
        WW(K)=WS
  200   DO 220 I=0,NXS
        DO 210 J=0,NYS
        P(K,I,J)=PS(I,J)
        H(K,I,J)=HS(I,J)
  210   CONTINUE
  220   CONTINUE
        END
C
C
        SUBROUTINE RESCAL(RESD,PS,HS,XS,YS,FPS,VISC,DENS,
     &             EPS,RH1,H00,NXS,NYS,DXS,DYS,RD,AL,AL1,
     &             A1,A2,A3,A4,Z0,NA,NB,DBS,
     &             M,K,KD,PI,DB,DT0,UDB)
        DIMENSION RESD(0:NA,0:NB),PS(0:NA,0:NB),HS(0:NA,0:NB),
     &            XS(0:NA),YS(0:NB),FPS(0:NA,0:NB),
     &            VISC(0:NA,0:NB),DENS(0:NA,0:NB),
     &            EPS(0:NA,0:NB),RH1(0:NA,0:NB),DBS(0:NA,0:NB),
     &            DB(M,0:NA,0:NB)
	  COMPLEX(8)::UDB(M,0:2*NA,0:2*NB)
        CALL FILMTH(HS,PS,XS,YS,DBS,H00,RD,AL1,NA,NB,NXS,
     &       NYS,M,K,KD,PI,DB,UDB)
        DXS2=DXS**2
        DYS2=DYS**2
        DO 20 I=0,NXS
        DO 10 J=0,NYS
        EPS(I,J)=HS(I,J)**3/AL
        RH1(I,J)=HS(I,J)
   10   CONTINUE
   20   CONTINUE
        DO 30 I=0,NXS
        RESD(I,0)=0.0
   30   CONTINUE
        DO 60 J=1,NYS-1
        RESD(0,J)=0.0
        DO 50 I=1,NXS-1
        EPSXM=0.5*(EPS(I-1,J)+EPS(I,J))
        EPSXP=0.5*(EPS(I+1,J)+EPS(I,J))
        EPSXI=EPSXM+EPSXP
        EPSYM=0.5*(EPS(I,J-1)+EPS(I,J))
        EPSYP=0.5*(EPS(I,J+1)+EPS(I,J))
        EPSYI=EPSYM+EPSYP
        RESD(I,J)=FPS(I,J)
     &       -(EPSXM*PS(I-1,J)-EPSXI*PS(I,J)+EPSXP*PS(I+1,J))/DXS2
     &       -(EPSYM*PS(I,J-1)-EPSYI*PS(I,J)+EPSYP*PS(I,J+1))/DYS2
     &       +RH1(I,J)/DT0
   50   CONTINUE
        RESD(NXS,J)=0.0
   60   CONTINUE
        DO 70 I=0,NXS
        RESD(I,NYS)=0.0
   70   CONTINUE
        END
C
C
        SUBROUTINE FUNCAL(FP,FG,X,Y,DX,DY,P,H,RES,WW,XS,YS,PS,
     &             HS,VISC,DENS,EPS,RH1,FPS,K,NXS1,NYS1,H00,
     &             RD,AL,AL1,A1,A2,A3,A4,Z0,M,NA,NB,DB,DBS,
     &             KD,PI,DT0,UDB)
        DIMENSION FP(M,0:NA,0:NB),FG(M),X(M,0:NA),
     &            Y(M,0:NB),DX(M),DY(M),
     &            P(M,0:NA,0:NB),H(M,0:NA,0:NB),RES(M,0:NA,0:NB),
     &            WW(M),PS(0:NA,0:NB),XS(0:NA),
     &            YS(0:NB),HS(0:NA,0:NB),VISC(0:NA,0:NB),
     &            DENS(0:NA,0:NB),EPS(0:NA,0:NB),RH1(0:NA,0:NB),
     &            FPS(0:NA,0:NB),DB(M,0:NA,0:NB),
     &            DBS(0:NA,0:NB)
	  COMPLEX(8)::UDB(M,0:2*NA,0:2*NB)
        K1=K-1
        DO 10 I=0,NXS1
        XS(I)=X(K1,I)
   10   CONTINUE
        DO 20 J=0,NYS1
        YS(J)=Y(K1,J)
   20   CONTINUE
        DO 50 I=0,NXS1
        DO 30 J=0,NYS1
        PS(I,J)=P(K1,I,J)
   30   CONTINUE
   50   CONTINUE
        DO 55 I=0,NA
        DO 54 J=0,NB
        DBS(I,J)=DB(K1,I,J)
   54   CONTINUE
   55   CONTINUE
C
        DXS=DX(K1)
        DXS2=DXS**2
        DYS=DY(K1)
        DYS2=DYS**2
        CALL FILMTH(HS,PS,XS,YS,DBS,H00,RD,AL1,NA,NB,NXS1,
     &       NYS1,M,K1,KD,PI,DB,UDB)
C
        DO 80 I=0,NXS1
        DO 70 J=0,NYS1
        EPS(I,J)=HS(I,J)**3/AL
        RH1(I,J)=HS(I,J)
   70   CONTINUE
   80   CONTINUE
C
        DO 90 I=0,NXS1
        FPS(I,0)=0.0
   90   CONTINUE
        DO 110 J=1,NYS1-1
        FPS(0,J)=0.0
        DO 100 I=1,NXS1-1
        EPSXM=0.5*(EPS(I-1,J)+EPS(I,J))
        EPSXP=0.5*(EPS(I+1,J)+EPS(I,J))
        EPSXI=EPSXM+EPSXP
        EPSYM=0.5*(EPS(I,J-1)+EPS(I,J))
        EPSYP=0.5*(EPS(I,J+1)+EPS(I,J))
        EPSYI=EPSYM+EPSYP
        FPS(I,J)=RES(K1,I,J)
     &       +(EPSXM*PS(I-1,J)-EPSXI*PS(I,J)+EPSXP*PS(I+1,J))/DXS2
     &       +(EPSYM*PS(I,J-1)-EPSYI*PS(I,J)+EPSYP*PS(I,J+1))/DYS2
     &       -RH1(I,J)/DT0
  100   CONTINUE
        FPS(NXS1,J)=0.0
  110   CONTINUE
        DO 120 I=0,NXS1
        FPS(I,NYS1)=0.0
  120   CONTINUE
        WG=0.0
        DO 160 I=1,NXS1-1,2
        DO 150 J=1,NYS1-1,2
        WG=WG+16.0*PS(I,J)
     &       +4.0*(PS(I-1,J)+PS(I+1,J)+PS(I,J-1)+PS(I,J+1))
     &       +PS(I-1,J-1)+PS(I+1,J-1)+PS(I-1,J+1)+PS(I+1,J+1)
  150   CONTINUE
  160   CONTINUE
        WG=WG*DXS*DYS/9.0
        FG(K1)=WG+FG(K)-WW(K)
        DO 180 I=0,NXS1
        DO 170 J=0,NYS1
        FP(K1,I,J)=FPS(I,J)
  170   CONTINUE
  180   CONTINUE
        END
C
C
        SUBROUTINE MLGRID(X,Y,DX,DY,NX,NY,XIN,XOUT,YIN,YOUT,M,NA,NB)
        DIMENSION X(M,0:NA),Y(M,0:NB),DX(M),DY(M),NX(M),NY(M)
        NX(1)=NA/(2**(M-1))
        NY(1)=NB/(2**(M-1))
        DO 10 K=2,M
        NX(K)=NX(K-1)*2
        NY(K)=NY(K-1)*2
   10   CONTINUE
        DO 20 K=1,M
        DX(K)=(XOUT-XIN)/NX(K)
        DY(K)=(YOUT-YIN)/NY(K)
   20   CONTINUE
        DO 50 K=1,M
        DO 30 I=0,NX(K)
        X(K,I)=XIN+DX(K)*I
   30   CONTINUE
        DO 40 I=0,NY(K)
        Y(K,I)=YIN+DY(K)*I                 
   40   CONTINUE
   50   CONTINUE
        END
C
C
        SUBROUTINE COEFFI (DB,DX,DY,NX,NY,M,NA,NB,KD)
        DIMENSION DX(M),DY(M),NX(M),NY(M),DB(M,0:NA,0:NB)
        DOUBLE PRECISION A1,A2,B1,B2,D1,D2,D3,D4,D5,D6,D7,D8,
     &            E1,E2,E3,E4
        DO 210 LLA=M,1,-1
        DLTX=DX(LLA)
        DAA=DLTX/2.0
        DLTY=DY(LLA)
        DBB=DLTY/2.0
        DO 20 I=0,NX(LLA)
        A1=DBLE(I*DLTX+DAA)
        A2=DBLE(I*DLTX-DAA)
        DO 10 J=0,NY(LLA)
        B1=DBLE(J*DLTY+DBB)
        B2=DBLE(J*DLTY-DBB)
        D1=A1+SQRT(A1*A1+B1*B1)
        D2=A2+SQRT(A2*A2+B1*B1)
        D3=A2+SQRT(A2*A2+B2*B2)
        D4=A1+SQRT(A1*A1+B2*B2)
        D5=B1+SQRT(A1*A1+B1*B1)
        D6=B2+SQRT(A1*A1+B2*B2)
        D7=B2+SQRT(A2*A2+B2*B2)
        D8=B1+SQRT(A2*A2+B1*B1)
        E1=B1*(LOG(D1)-LOG(D2))
        E2=B2*(LOG(D3)-LOG(D4))
        E3=A1*(LOG(D5)-LOG(D6))
        E4=A2*(LOG(D7)-LOG(D8))
        DB(LLA,I,J)=REAL((E1+E2)+(E3+E4))
   10   CONTINUE
   20   CONTINUE
  210   CONTINUE
	  END
C
C
        SUBROUTINE PMGIVE(P,X,Y,RK,M,NA,NB)
        DIMENSION P(M,0:NA,0:NB),X(M,0:NA),Y(M,0:NB)
        EK=RK**2
        DO 20 I=0,NA
        DO 10 J=0,NB
        PJ=1.0-X(M,I)**2-Y(M,J)**2/EK
        IF(PJ.LE.1.0E-20) THEN
        P(M,I,J)=0.0
        ELSE
        P(M,I,J)=SQRT(PJ)
        END IF
   10   CONTINUE
   20   CONTINUE
        END
C
C
        SUBROUTINE FILMTH(HS,PS,XS,YS,DBS,H00,RD,AL1,NA,NB,NXS,
     &             NYS,M,ML,KD,PI,DB,UDB)
        DIMENSION HS(0:NA,0:NB),PS(0:NA,0:NB),XS(0:NA),YS(0:NB),
     &            DBS(0:NA,0:NB),DB(M,0:NA,0:NB)
	  COMPLEX(8)::UDB(M,0:2*NA,0:2*NB)
        COMPLEX(8)::UPS(0:2*NA,0:2*NB)
	  COMPLEX(8)::UDS(0:2*NA,0:2*NB)
	  COMPLEX(8)::YYP(0:2*NB)
	  COMPLEX(8)::YYD(0:2*NB)
	  COMPLEX(8)::YYDE(0:2*NB)
	  COMPLEX(8)::UDE(0:2*NA,0:2*NB)
C先将时域内的压力和变形系数赋到复数域
	  
	  FFT_FORWARD=1
	  FFT_INVERSE=-1
	  ND=2*NXS
	  NE=2*NYS
	  DO 20 I=0,NXS
	  DO 10 J=0,NYS
	  UPS(I,J)=PS(I,J)
   10   CONTINUE
   20   CONTINUE
C将频域内的压力扩展到2*NXS和2*NYS
        DO 40 I=0,NXS
	  DO 30 J=NYS+1,NE-1
	  UPS(I,J)=0
   30   CONTINUE
   40   CONTINUE
        DO 60 J=0,NE-1
	  DO 50 I=NXS+1,ND-1
	  UPS(I,J)=0
   50   CONTINUE
   60   CONTINUE
C对复数域内的变形系数和压力进行DFT
C先对列进行FFT
	  DO 130 I=0,ND-1
	  DO 110 J=0,NE-1
	  YYP(J)=UPS(I,J)
  110   CONTINUE
        CALL FCFFT(YYP,PI,NE,FFT_FORWARD)
	  DO 120 J=0,NE-1
	  UPS(I,J)=YYP(J)
  120   CONTINUE
  130   CONTINUE
C再对行进行FFT
        DO 160 J=0,NE-1
	  DO 140 I=0,ND-1
	  YYP(I)=UPS(I,J)
  140   CONTINUE
        CALL FCFFT(YYP,PI,ND,FFT_FORWARD)
	  DO 150 I=0,ND-1
	  UPS(I,J)=YYP(I)
  150   CONTINUE
  160   CONTINUE
C在频域内进行相乘计算卷积，即计算频域内的无量纲变形
        DO 180 I=0,ND-1
	  DO 170 J=0,NE-1
	  UDE(I,J)=UDB(ML,I,J)*UPS(I,J)
  170   CONTINUE
  180   CONTINUE
C对计算得到的UDE进行IFFT傅里叶逆变换
C先对每一列作IFFT
        DO 210 I=0,ND-1
	  DO 190 J=0,NE-1
	  YYDE(J)=UDE(I,J)
  190   CONTINUE
        CALL FCFFT(YYDE,PI,NE,FFT_INVERSE)
	  DO 200 J=0,NE-1
	  UDE(I,J)=YYDE(J)
  200   CONTINUE
  210   CONTINUE
C再对每一行作IFFT
        DO 240 J=0,NE-1
	  DO 220 I=0,ND-1
	  YYDE(I)=UDE(I,J)
  220   CONTINUE
        CALL FCFFT(YYDE,PI,NE,FFT_INVERSE)
	  DO 230 I=0,ND-1
	  UDE(I,J)=YYDE(I)
  230   CONTINUE
  240   CONTINUE
	  DO 260 I=0,NXS
	  DO 250 J=0,NYS
	  HS(I,J)=REAL(UDE(I,J)/ND/NE)
  250   CONTINUE
  260   CONTINUE
        DO 280 I=0,NXS
        DO 270 J=0,NYS
        HS(I,J)=H00+0.5*(XS(I)**2+YS(J)**2/RD)+AL1*HS(I,J)
  270   CONTINUE
  280   CONTINUE
        END
C
        SUBROUTINE INTER(U,M,NA,NB,K,NX,NY)
        DIMENSION U(M,0:NA,0:NB),NX(M),NY(M)
        K1=K+1
        DO 20 I=0,NX(K)
        DO 10 J=0,NY(K)
        U(K1,2*I,2*J)=U(K,I,J)
   10   CONTINUE
   20   CONTINUE
        DO 40 J=0,NY(K1),2
        DO 30 I=1,NX(K1)-1,2
        U(K1,I,J)=0.5*(U(K1,I-1,J)+U(K1,I+1,J))
   30   CONTINUE
   40   CONTINUE
        DO 60 I=0,NX(K1),2
        DO 50 J=1,NY(K1)-1,2
        U(K1,I,J)=0.5*(U(K1,I,J-1)+U(K1,I,J+1))
   50   CONTINUE
   60   CONTINUE
        DO 80 I=1,NX(K1)-1,2
        DO 70 J=1,NY(K1)-1,2
        U(K1,I,J)=0.5*(U(K1,I-1,J)+U(K1,I+1,J))
   70   CONTINUE
   80   CONTINUE
        END
C
C
        SUBROUTINE RESTR(U,M,NA,NB,K,NX,NY)
        DIMENSION U(M,0:NA,0:NB),NX(M),NY(M)
        K1=K-1
        DO 20 I=1,NX(K1)-1
        NI=2*I
        DO 10 J=1,NY(K1)-1
        NJ=2*J
        U(K1,I,J)=0.25*U(K,NI,NJ)+
     &           0.125*(U(K,NI+1,NJ)+U(K,NI-1,NJ)+
     &                  U(K,NI,NJ+1)+U(K,NI,NJ-1))+
     &          0.0625*(U(K,NI+1,NJ+1)+U(K,NI-1,NJ-1)+
     &                  U(K,NI+1,NJ-1)+U(K,NI-1,NJ+1))
   10   CONTINUE
   20   CONTINUE
        NU=NX(K)
        ND=NX(K1)
        DO 30 J=0,NY(K1)
        U(K1,0,J)=U(K,0,2*J)
        U(K1,ND,J)=U(K,NU,2*J)
   30   CONTINUE
        NU=NY(K)
	  ND=NY(K1)
        DO 40 I=1,NX(K1)-1
        U(K1,I,0)=U(K,2*I,0)
	  U(K1,I,ND)=U(K,2*I,NU)
   40   CONTINUE
        END
C
C
        SUBROUTINE EXTEND(P,PB,NX,NY,DP,K,M,NA,NB)
        DIMENSION P(M,0:NA,0:NB),PB(M,0:NA,0:NB),NX(M),NY(M),
     &            DP(M,0:NA,0:NB)
        DO 20 I=0,NX(K)
        DO 10 J=0,NY(K)
        DP(K,I,J)=P(K,I,J)-PB(K,I,J)
   10   CONTINUE
   20   CONTINUE
        CALL INTER(DP,M,NA,NB,K,NX,NY)
        DO 50 I=0,NX(K+1)
        DO 30 J=0,NY(K+1)
        P(K+1,I,J)=P(K+1,I,J)+DP(K+1,I,J)
   30   CONTINUE
   50   CONTINUE
        END
C
C
        SUBROUTINE VISDEN(VISC,DENS,PS,NXS,NYS,A1,A2,A3,A4,Z0,NA,NB)
        DIMENSION VISC(0:NA,0:NB),DENS(0:NA,0:NB),PS(0:NA,0:NB)
        DO 20 I=0,NXS
        DO 10 J=0,NYS
        VISC(I,J)=1.0
CEXP(A1*((1.0+A2*PS(I,J))**Z0-1.0))
        DENS(I,J)=1.0
C1.0+A3*PS(I,J)/(1.0+A4*PS(I,J))
   10   CONTINUE
   20   CONTINUE
        END
C
        
	  
C       下面为一维快速傅里叶变换
	  SUBROUTINE FCFFT(X,PI,N1,FORBACK)
	  COMPLEX(8)::X(0:N1-1)
	  COMPLEX(8)::CTMP
	  NCUR=SIZE(X)
   10	  NTEMP=NCUR
	  E=-2.0*PI/NCUR
	  NCUR=NCUR/2
	  DO J=0,NCUR-1
	  DO I=J,N1-1,NTEMP
	  ITMP=I+NCUR
	  CTMP=X(I)-X(ITMP)
	  X(I)=X(I)+X(ITMP)
	  X(ITMP)=CTMP*EXP(FORBACK*CMPLX(0.0,E*J))
	  ENDDO
	  ENDDO
	  IF(NCUR.GT.1) GOTO 10
   	  J=0
	  DO I=0,N1-2
	  IF(I.LT.J)THEN
	  CTMP=X(J)
	  X(J)=X(I)
	  X(I)=CTMP
	  ENDIF
	  K=N1/2
	  DO WHILE(K.LE.J)
	  J=J-K
	  K=K/2
	  ENDDO
	  J=J+K
	  ENDDO
	  END
C
        SUBROUTINE RESULT0(X,Y,P,H,M,NA,NB,XIN,XOUT,YIN,YOUT,PH,
     &             RK,RX,A,ICY,PHERTZ,DEFOR,PI)
        DIMENSION X(M,0:NA),Y(M,0:NB),P(M,0:NA,0:NB),H(M,0:NA,0:NB),
     &            DEFOR(0:NA,0:NB)
        OPEN(9,FILE='DEFOR0.DAT')
        OPEN(10,FILE='P0.DAT')
        OPEN(11,FILE='H0.DAT')
        OPEN(12,FILE='PHY0.DAT')
        MK=1
        DO 6 I=0,NA
        DO 5 J=0,NB
        H(M,I,J)=H(M,I,J)*A*A/RX*1.0E6  
	  P(M,I,J)=P(M,I,J)*PHERTZ*1.0E-6
    5   CONTINUE
    6   CONTINUE
C
        XL=-PI/2
        XR=PI/2
        YL=-PI/2
        YR=PI/2
C
        DLTX=(XOUT-XIN)/NA
        DLTY=(YOUT-YIN)/NB
        NXL=INT((XL-XIN)/DLTX+0.1)
        NXR=INT((XR-XIN)/DLTX+0.1)
        NYL=INT((YL-YIN)/DLTY+0.1)
        NYR=INT((YR-YIN)/DLTY+0.1)
C
        WRITE(9,100) 0,(Y(M,J),J=NYL,NYR,MK)
        WRITE(10,100) 0,(Y(M,J),J=NYL,NYR,MK)
        WRITE(11,100) 0,(Y(M,J),J=NYL,NYR,MK)
  100   FORMAT(I12,10000(1X,F12.5))
        DO 110 I=NXL,NXR,MK
        WRITE(9,120) X(M,I),(DEFOR(I,J),J=NYL,NYR,MK)
        WRITE(10,120) X(M,I),(P(M,I,J),J=NYL,NYR,MK)
        WRITE(11,120) X(M,I),(H(M,I,J),J=NYL,NYR,MK)
  110   CONTINUE
  120   FORMAT(F12.5,10000(1X,E12.5))
C
        DO 130 I=0,NA
        WRITE(12,140) X(M,I),P(M,I,ICY),H(M,I,ICY)
  130   CONTINUE
  140   FORMAT(1X,F15.7,E15.7,E15.7)
        CLOSE(9)
        CLOSE(10)
        CLOSE(11)
        CLOSE(12)
        END
C
        SUBROUTINE RESULT1(X,Y,P,H,M,NA,NB,XIN,XOUT,YIN,YOUT,PH,
     &             RK,RX,A,ICY,PHERTZ,DEFOR,PI)
        DIMENSION X(M,0:NA),Y(M,0:NB),P(M,0:NA,0:NB),H(M,0:NA,0:NB),
     &            DEFOR(0:NA,0:NB)
        OPEN(9,FILE='DEFOR1.DAT')
        OPEN(10,FILE='P1.DAT')
        OPEN(11,FILE='H1.DAT')
        OPEN(12,FILE='PHY1.DAT')
        MK=1
        DO 6 I=0,NA
        DO 5 J=0,NB
        H(M,I,J)=H(M,I,J)*A*A/RX*1.0E6  
	  P(M,I,J)=P(M,I,J)*PHERTZ*1.0E-6
C*A*A/RX*1.0E6  *PHERTZ*1.0E-6
    5   CONTINUE
    6   CONTINUE
C
        XL=-PI/2
        XR=PI/2
        YL=-PI/2
        YR=PI/2
C
        DLTX=(XOUT-XIN)/NA
        DLTY=(YOUT-YIN)/NB
        NXL=INT((XL-XIN)/DLTX+0.1)
        NXR=INT((XR-XIN)/DLTX+0.1)
        NYL=INT((YL-YIN)/DLTY+0.1)
        NYR=INT((YR-YIN)/DLTY+0.1)
C
        WRITE(9,100) 0,(Y(M,J),J=NYL,NYR,MK)
        WRITE(10,100) 0,(Y(M,J),J=NYL,NYR,MK)
        WRITE(11,100) 0,(Y(M,J),J=NYL,NYR,MK)
  100   FORMAT(I12,10000(1X,F12.5))
        DO 110 I=NXL,NXR,MK
        WRITE(9,120) X(M,I),(DEFOR(I,J),J=NYL,NYR,MK)
        WRITE(10,120) X(M,I),(P(M,I,J),J=NYL,NYR,MK)
        WRITE(11,120) X(M,I),(H(M,I,J),J=NYL,NYR,MK)
  110   CONTINUE
  120   FORMAT(F12.5,10000(1X,E12.5))
C
        DO 130 I=0,NA
        WRITE(12,140) X(M,I),P(M,I,ICY),H(M,I,ICY)
  130   CONTINUE
  140   FORMAT(1X,F15.7,E15.7,E15.7)
        CLOSE(9)
        CLOSE(10)
        CLOSE(11)
        CLOSE(12)
        END
C
        SUBROUTINE RESULT2(X,Y,P,H,M,NA,NB,XIN,XOUT,YIN,YOUT,PH,
     &             RK,RX,A,ICY,PHERTZ,DEFOR,PI)
        DIMENSION X(M,0:NA),Y(M,0:NB),P(M,0:NA,0:NB),H(M,0:NA,0:NB),
     &            DEFOR(0:NA,0:NB)
        OPEN(9,FILE='DEFOR2.DAT')
	  OPEN(10,FILE='P2.DAT')
        OPEN(11,FILE='H2.DAT')
        OPEN(12,FILE='PHY2.DAT')
        MK=1
        DO 6 I=0,NA
        DO 5 J=0,NB
        H(M,I,J)=H(M,I,J)*A*A/RX*1.0E6  
	  P(M,I,J)=P(M,I,J)*PHERTZ*1.0E-6
C*A*A/RX*1.0E6  *PHERTZ*1.0E-6
    5   CONTINUE
    6   CONTINUE
C
        XL=-PI/2
        XR=PI/2
        YL=-PI/2
        YR=PI/2
C
        DLTX=(XOUT-XIN)/NA
        DLTY=(YOUT-YIN)/NB
        NXL=INT((XL-XIN)/DLTX+0.1)
        NXR=INT((XR-XIN)/DLTX+0.1)
        NYL=INT((YL-YIN)/DLTY+0.1)
        NYR=INT((YR-YIN)/DLTY+0.1)
C
        WRITE(9,100) 0,(Y(M,J),J=NYL,NYR,MK)
	  WRITE(10,100) 0,(Y(M,J),J=NYL,NYR,MK)
        WRITE(11,100) 0,(Y(M,J),J=NYL,NYR,MK)
  100   FORMAT(I12,10000(1X,F12.5))
        DO 110 I=NXL,NXR,MK
        WRITE(9,120) X(M,I),(DEFOR(I,J),J=NYL,NYR,MK)
	  WRITE(10,120) X(M,I),(P(M,I,J),J=NYL,NYR,MK)
        WRITE(11,120) X(M,I),(H(M,I,J),J=NYL,NYR,MK)
  110   CONTINUE
  120   FORMAT(F12.5,10000(1X,E12.5))
C
        DO 130 I=0,NA
        WRITE(12,140) X(M,I),P(M,I,ICY),H(M,I,ICY)
  130   CONTINUE
  140   FORMAT(1X,F15.7,E15.7,E15.7)
        CLOSE(9)
	  CLOSE(10)
        CLOSE(11)
        CLOSE(12)
        END
C	  	  
        SUBROUTINE RESULT3(X,Y,P,H,M,NA,NB,XIN,XOUT,YIN,YOUT,PH,
     &             RK,RX,A,ICY,PHERTZ,DEFOR,PI)
        DIMENSION X(M,0:NA),Y(M,0:NB),P(M,0:NA,0:NB),H(M,0:NA,0:NB),
     &            DEFOR(0:NA,0:NB)
        OPEN(9,FILE='DEFOR3.DAT')
	  OPEN(10,FILE='P3.DAT')
        OPEN(11,FILE='H3.DAT')
        OPEN(12,FILE='PHY3.DAT')
        MK=1
        DO 6 I=0,NA
        DO 5 J=0,NB
        H(M,I,J)=H(M,I,J)*A*A/RX*1.0E6  
	  P(M,I,J)=P(M,I,J)*PHERTZ*1.0E-6
C*A*A/RX*1.0E6  *PHERTZ*1.0E-6
    5   CONTINUE
    6   CONTINUE
C
        XL=-PI/2
        XR=PI/2
        YL=-PI/2
        YR=PI/2
C
        DLTX=(XOUT-XIN)/NA
        DLTY=(YOUT-YIN)/NB
        NXL=INT((XL-XIN)/DLTX+0.1)
        NXR=INT((XR-XIN)/DLTX+0.1)
        NYL=INT((YL-YIN)/DLTY+0.1)
        NYR=INT((YR-YIN)/DLTY+0.1)
C
        WRITE(9,100) 0,(Y(M,J),J=NYL,NYR,MK)
	  WRITE(10,100) 0,(Y(M,J),J=NYL,NYR,MK)
        WRITE(11,100) 0,(Y(M,J),J=NYL,NYR,MK)
  100   FORMAT(I12,10000(1X,F12.5))
        DO 110 I=NXL,NXR,MK
        WRITE(9,120) X(M,I),(DEFOR(I,J),J=NYL,NYR,MK)
	  WRITE(10,120) X(M,I),(P(M,I,J),J=NYL,NYR,MK)
        WRITE(11,120) X(M,I),(H(M,I,J),J=NYL,NYR,MK)
  110   CONTINUE
  120   FORMAT(F12.5,10000(1X,E12.5))
C
        DO 130 I=0,NA
        WRITE(12,140) X(M,I),P(M,I,ICY),H(M,I,ICY)
  130   CONTINUE
  140   FORMAT(1X,F15.7,E15.7,E15.7)
        CLOSE(9)
	  CLOSE(10)
        CLOSE(11)
        CLOSE(12)
        END
        