	
	!In this subroutine, different surface type can be obtained
	SUBROUTINE SURFACE_TYPE(ISHAP,IQ,AMPT,M,N,SURFACE,XS,YS,DX,DY,R0,RX,RY,A0,B0,H0)
	IMPLICIT NONE
	INTEGER ISHAP,IQ,M,N,I,J,NX,NY,IR,JR,K,L
	REAL*8 XS,YS,DX,DY,R0,RX,RY,A0,B0,H0
	REAL*8 XEND1,YEND1,DXR,DYR,XX,YY,A1,A2,A3,A4,XD,YD
	REAL*8 SURFACE(M,N),AMPT,RASP,RR0,RR1,X,Y,ZMIN
	REAL*8 RRA1,RRQ1,RRT1,SUM,PI,Y1,Y2,R,RL
	REAL*8,ALLOCATABLE::Z0(:,:),ZH(:,:),HEIGHT(:,:)
	CHARACTER Int2Str*4
	PI  = 4.0*ATAN(1.0D0)
	ALLOCATE(Z0(M,N),ZH(M,N))
	!-------------------------------------------------------------
	IF(ISHAP == 0) THEN			!FOR SMOOTH SPHERE SURFACE
		DO I=1,M
			DO J=1,N
				Z0(I,J)=0.0
			ENDDO
		ENDDO
	END IF
	!-------------------------------------------------------------
	IF(ISHAP == 1) THEN			!SINGLE SEMI-SPHERE
		RASP = 0.2				!RADIUS OF ASPERITY,[A0]		
		RR0  = RASP*RASP
		DO I=1,M
			X=XS+DX*(I-1)
			DO J=1,N
				Y=YS+DY*(J-1)
				RR1=X*X+Y*Y
				IF(RR1.GE.RR0) RR1=RR0
				Z0(I,J)=-(AMPT*R0/A0/B0)*SQRT(RR0-RR1)/RASP
			ENDDO
		ENDDO
	ENDIF
	!-------------------------------------------------------------
	IF(ISHAP == 2) THEN			!TRANVERSE CONSINE
		DO I=1,M
			X=XS+DX*(I-1)
			DO J=1,N
				Y=YS+DY*(J-1)
		!		Z0(I,J) = (AMPT*R0/A0/B0)*COS(2*3.1415926*X/(0.32))	!TRANSVERSE ROUGHNESS
				Z0(I,J) = (AMPT*R0/A0/B0)*COS(2*3.1415926*Y/(0.50))	!LONGITUDE ROUGHNESS
			ENDDO
		ENDDO
	ENDIF
	!--------------------------------------------------------------
	IF(ISHAP == 3) THEN			!a transversal ridge from Venner's paper
		RASP=0.7
		DO I=1,M
			X=XS+DX*(I-1)
			DO J=1,N
				Y=YS+DY*(J-1)
				Z0(I,J)=(AMPT*R0/A0/B0)*EXP(-10.*X**2/RASP**2)&
							           *COS(2*PI*X/RASP)
			ENDDO
		ENDDO
	ENDIF
	!-------------------------------------------------------------
	IF(ISHAP == 4) THEN			!ISOTROPIC SINUSOIDAL SURFACE
		DO I=1,M
			X=XS+DX*(I-1)
			DO J=1,N
				Y=YS+DY*(J-1)
				Z0(I,J) = -(AMPT)*COS(2*PI*X/(1.0))&
										 *COS(2*PI*Y/(1.0))
			ENDDO
		ENDDO
	ENDIF
	!----------------------------------------------------------------
	IF(ISHAP == 5) THEN			!single asperity
		RASP=0.4
		DO I=1,M
			X=XS+DX*(I-1)
			DO J=1,N
				Y=YS+DY*(J-1)
				Z0(I,J)=(AMPT*R0/A0/B0)*EXP(-10.*(X**2+Y**2)/RASP**2)&
							           *COS(2*PI*sqrt(X**2+Y**2)/RASP)
			ENDDO
		ENDDO
	ENDIF
	!----------------------------------------------------------------
	IF(ISHAP < 0) THEN			!REAL ENGINEERING SURFACE
		OPEN(1,FILE='SurfaceGeometry.dat')
			READ(1,*) NX,NY,XEND1,YEND1
			ALLOCATE(HEIGHT(0:NX-1,0:NY-1))
			DO I=0,NX-1
				READ(1,*) (HEIGHT(I,J),J=0,NY-1)
			ENDDO
!			READ(1,*) RRA1,RRQ1,RRT1
		CLOSE(1)

300		FORMAT(1X,5E16.8)
400		FORMAT(2X,2I4,4X,4E16.8)  
		DXR=XEND1/A0/(NX-1)
		DYR=YEND1/A0/(NY-1)
		DO I=1,M
			X=DX*(I-1)
			DO J=1,N
				Y=DY*(J-1)
				IF(Y>YEND1/A0) THEN
					Y=Y-YEND1/A0
				ENDIF
				XX=-DXR
				DO K=0,NX-1
					XX=XX+DXR
					IR=K
					IF (XX.GT.X+1.0E-10) GO TO 20
				ENDDO
20				YY=-DYR
				DO L=0,NY-1
					YY=YY+DYR
					JR=L
					IF (YY.GE.Y+1.0E-10) GO TO 40
				ENDDO
40				A1=HEIGHT(IR-1,JR-1)
				A2=HEIGHT(IR  ,JR-1)-A1
				A3=HEIGHT(IR-1,JR  )-A1
				A4=HEIGHT(IR  ,JR  )-A1-A2-A3
				XD=(X-XX+DXR)/DXR
				YD=(Y-YY+DYR)/DYR
				Z0(I,J)=A1+A2*XD+A3*YD+A4*XD*YD
			ENDDO
		ENDDO
		SUM=0.
		DO I=1,M
			DO J=1,N
				SUM=SUM+Z0(I,J)*Z0(I,J)
			ENDDO
		ENDDO
		SUM = SUM/M/N
		SUM = SQRT(SUM)
		DO I=1,M
			DO J=1,N
				Z0(I,J)=AMPT*Z0(I,J)/A0
			ENDDO
		ENDDO
		DEALLOCATE(HEIGHT)
	ENDIF
	!----------------------------------------------------------------
	!OPEN(20,FILE='Roughness.DAT')
	!	DO J=1,N
	!		WRITE(20,'('//Int2Str(M,4)//'F16.7)')(Z0(I,J),I=1,M)  	!nondimensional roughness(/A0)
	!	ENDDO
	!CLOSE(20)
	!----------------------------------------------------------------
	IF(IQ == 1) THEN	!IQ=1 FOR POINT CONTACTS
		DO I=1,M
			X=XS+DX*(I-1)
			DO J=1,N
				Y=(YS+DY*(J-1))
				RR1=X*X/RX*A0+Y*Y/RY*A0
				RR0=H0+RR1/2.
				ZH(I,J)=RR0    !Nondimensional
			ENDDO
		ENDDO
	ENDIF
	IF(IQ == 2) THEN	!IQ=2 for cylinder of finite length
		RL=1.0*a0	!Unit: mm
		R=24.		!Unit: mm
		!RL=1.0
		!R =30
		DO I=1,N
			X=XS+DX*(I-1)
			DO J=1,M
				Y=YS+DY*(J-1)
				Y1=Dabs(DABS(YS)-DABS(Y))*B0
				IF(Y1<=RL) THEN
					Y2=	R-DSQRT(R**2-(RL-Y1)**2)
				ELSE
					Y2=0.0
				ENDIF
			!	write(*,*) y2
				RR1=R0*(X*X/(RX-Y2))
				RR0=H0+RR1/2.+Y2*R0/a0/b0
				ZH(I,J)=RR0
			ENDDO
		!	pause
		ENDDO
	ENDIF
	IF(IQ == 3) THEN	!IQ=2 FOR FLAT CONTACT
		DO I=1,N
			DO J=1,M
				ZH(I,J)=0
			ENDDO
		ENDDO
	ENDIF
	IF(IQ == 4) THEN	!IQ=4 for cylinder of finite length
		RL=1.0*a0	!Unit: mm
		R=24.		!Unit: mm
		DO I=1,N
			X=XS+DX*(I-1)
			DO J=1,M
				Y=YS+DY*(J-1)
			!	y2=dabs(x)-0.2/a0
				if(dabs(y)>=0.2/a0) then
					RR1=(dabs(y)-0.2/a0)/(1.8/a0)
				else
					rr1=0.0
				endif
				RR0=H0+RR1*R0/a0/b0
				ZH(I,J)=RR0
			ENDDO
		ENDDO
	ENDIF

	!-- SUPERPOSE ROUGHNESS TO THE SPHERICAL SURFACE --
	ZMIN=0.0
	DO I=1,M
		DO J=1,N
			Z0(I,J)=ZH(I,J)+Z0(I,J)
			IF(Z0(I,J) < ZMIN) ZMIN=Z0(I,J)
		ENDDO
	ENDDO
	DO I=1,M
		DO J=1,N
			SURFACE(I,J)=Z0(I,J)-ZMIN
		ENDDO
	ENDDO

!	OPEN(20,FILE='GEOMETRY.DAT')
!		DO J=1,N
!			WRITE(20,'('//Int2Str(M,4)//'F16.7)')(SURFACE(I,J),I=1,M) !nondimensional geometry(/A0)
!		ENDDO
!	CLOSE(20)
	DEALLOCATE(Z0,ZH)
	END
