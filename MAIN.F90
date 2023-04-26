!************************************************************************************************************************
!	The source code is developed to calculate contact problem for a layer-substrate system 
!   with an imperfectly bonded interface: spring-like condition.
!   The sliding contact problem of two elastic solids with each of them covered by a coating can be solved efficiently.
!   The contact surface can be smooth or rough.
!   The code was developed by Zhanjiang Wang and Hao Yu in 2017, China.
!   The formulas and derivations can be found in the paper listed below,
!   1. Wang, Z., Yu, H., Wang, Q., 2017, "Layer-substrate System with an Imperfectly Bonded Interface: Spring-like Condition," International Journal of Mechanical Sciences.
!   Some related references are also listed as follows,
!   2. Wang, Z., Yu, H., Wang, Q., 2017, "Layer-substrate System with an Imperfectly Bonded Interface: Coupled Dislocation-Like and Force-Like Conditions," International Journal of Solids and Structures, 122-123, pp. 91-109.
!   3. Yu, C., Wang, Z., and Wang, Q., 2014, "Analytical Frequency Response Functions for Contact of Multilayered Materials," Mechanics of Materials, 76, pp. 102-120.
!   4. Wang, Z., Yu, C., and Wang, Q., 2015, "An Efficient Method for Solving Three-dimensional Fretting Contact Problems Involving Multilayered or Functionally Graded Materials," International Journal of Solids and Structures, 66, pp. 46-61.
!   5. Liu, S., Wang, Q., 2002, "Studying Contact Stress Fields Caused by Surface Tractions with a Discrete Convolution and Fast Fourier Transform Algorithm," ASME Journal of Tribology 124(1), 36-45.
!   6. O¡¯Sullivian, T.C., King, R.B., 1988, "Sliding Contact Stress-field Due to a Spherical Indenter on a Layered Elastic Half-space," ASME Journal of Tribology 110 (2), 235-240.
!   The code can be used for free for non-commercial and academic use.
!************************************************************************************************************************
    PROGRAM DRY_CONTACT_FFT2
	IMPLICIT NONE
	INTEGER M,N,M2,N2,L,L1,L2,ITERATION
	PARAMETER(M=129,N=M,L=100,M2=(M-1)*2,N2=(N-1)*2,L1=L/2,L2=L/2)  !M*N*L1: grids in coating; M*N*L2: grids in substrate
	INTEGER ISHAP,IQ,I,J,K,ITERATION_NUM
	REAL*8 WLD,R11,R12,R21,R22,EM1,EM2,E1,E2,E0,E00,EV,EVV,A0,B0,P0,R0,AKE,RX,RY
	REAL*8 PI,XS,XE,YS,YE,DX,DY,FRC
	REAL*8 H0,AMPT,RASP,W,ERROR_P,ERROR_L
	REAL*8 P(M,N),GAP(M,N),Z0(M,N)
	DOUBLE COMPLEX ZZFFT1(M2,N2),ZZFFT2(M2,N2),ZZFFT(M2,N2)	!This is for influence coefficient
	CHARACTER FileIn*20,Int2Str*4,Formatt*11
	INTEGER NGaussQ,alias,exp
	REAL*8 EM1_s,EM2_s,E1_s,E2_s,h_c1,h_c2
    REAL*8 t(3)  !0<t<=+oo interfacial boundary parameters
	REAL*8 Z1,Z2
    !Stress_XYZ(1,M,N,L),...,(6,M,N,L)--Stress components 11,22,33,12,13,23
	REAL*8 Stress_XYZ(6,M,N,L),VonMises_Stress(M,N,L)
    
	FORMATT='('//Int2Str(M,4)//'F16.7)'
	PI  = 4.0*ATAN(1.0D0)

	!---------determine some parameters-------------------------------------------
	!Working conditions, iteration constant, and the geomerty of contacting body
	WLD  = 14.652		!Applied load [N]
	FRC  = 0.0          !Friction coefficient
	IQ   = 1			!IQ=1 for point contact;=2 for flat contact
	ISHAP= 0 !5			!0->smooth sphere;1->semisphere;2->sin surface
							!<0->measured surface
	AMPT = 1.0          !Amplitude of artificial surface; rms for real surface
	ERROR_P=0.1E-6		!Convergence criterion of pressure
	ERROR_L=0.5E-4		!Convergence criterion for load balance
	ITERATION_NUM=2000	!Maximum iteration number
	R11 = 10.		    !Radius of the ball 1, [mm]
	R12 = 10.		    !Radius of the ball 1, [mm]
	R21 = 1.0E35        !Radius of the ball 2, [mm]
	R22 = 1.0E35        !Radius of the ball 2, [mm]
    
    !Dimensionless compauting domain (XS,XE)*(YS,YE)
	XS  = -2.0
	XE  =  2.0
	YS  = -2.0
	YE  =  2.0

    !Material parameters
	EM1 = 0.3		    !Poisson's ratio of coating in body 1
    EM1_s=0.3           !Poisson's ratio of substrate in body 1
	EM2 = 0.3		    !Poisson's ratio of coating in body 2
    EM2_s=0.3           !Poisson's ratio of substrate in body 2
	E1  = 2.1E35		!Young's modulus of coating in body 1, [MPa]
	E1_s =2.1E35        !Young's modulus of substrate in body 1, [MPa]
	E2  = 0.5E5	        !Young's modulus of coating in body 2, [MPa]
	E2_s =1.0E5         !Young's modulus of substrate in body 2, [MPa]
	h_c1=0.2            !Dimensionless coating thickness in body 1, [h_c1/a0]
	h_c2=1.0            !Dimensionless coating thickness in body 2, [h_c2/a0]
    
	NGaussQ=64  !Number of quadrature points for the Gaussian integral when calculate FRF in the original point. 
	alias=0     !AL in Eq.(27) in  Wang et al. IJSS 2017
	EXP=8       !Gamma in Eq.(27) in  Wang et al. IJSS 2017
    
    !Dimensionless interfacial jumping coefficients T(i) (i=1,2,3)
    !The tangential spring coefficients in x and y directions are defined to be equal, i.e., T(1) = T(2)
    T(1)=10.0 !1.0E35 !100.0 !28.58849 !1.0E35 !T(1)=T(2)
    T(2)=T(1) !1.0E35
    T(3)=10.0 !1.0E35 !28.58849 !-3000 !1.0E35

	!Determine hertz contact parameters
	CALL HERTZ_PARAMETER(WLD,R11,R12,R21,R22,EM1_s,EM2_s,E1_s,E2_s,A0,B0,RX,RY,R0,P0,AKE,E0,E00,EV,EVV)
    
	!Determin mesh size
	DX  = (XE-XS)/(M-1)		!Mesh size in X, nondimensional
	DY  = (YE-YS)/(N-1)		!Mesh size in Y, nondimensional
	
	WRITE(*,*)
	WRITE(*,*) '          ****** NUMERICAL ANALYSIS OF ELASTIC CONTACT ******'
	WRITE(*,250) WLD,RX,RY,A0,B0,AKE,P0,FRC
	
	!Obtain initial gap
	!Note the initial minimum gap should be zero 
	IF(ISHAP==0) AMPT=0.0
	H0  = 0.0	!Gap at the contact center

    CALL SURFACE_TYPE(ISHAP,IQ,AMPT,M,N,Z0,XS,YS,DX,DY,R0,RX,RY,A0,B0,H0)
	
	!Determine initial pressure the subroutine is in initial press.f90
	CALL PINITIAL(M,N,P,XS,YS,DX,DY)
	
    !Dimensionless Young's modulus	
    E1  = E1/P0
	E1_s =E1_s/P0
	E2  = E2/P0
	E2_s =E2_s/P0
    
	!Calculate influence coefficients relating pressure to surface displacement
    !ZZFFT1--the discrete Fourier transform of the ICs relating pressure to surface displacement u3
    CALL Multilayer_Displacement_IC(M-1,N-1,M2,N2,ZZFFT1,NGaussQ,alias,exp,XS,XE,YS,YE&
	                               ,E1_s,E1,EM1_s,EM1,h_c1,t,9)   ! For body 1
    CALL Multilayer_Displacement_IC(M-1,N-1,M2,N2,ZZFFT2,NGaussQ,alias,exp,XS,XE,YS,YE&
	                               ,E2_s,E2,EM2_s,EM2,h_c2,t,9)   ! For body 2
    
    !Body 1 + Body 2
    DO I=1,M2
	  DO J=1,N2
	    ZZFFT(I,J)=ZZFFT1(I,J)+ZZFFT2(I,J)
	  ENDDO
	ENDDO

	WRITE(*,*)
	WRITE(*,*)'  **** DETERMINE PRESSURE AND DEFORMATION BY CONJUGATE GRADIENT METHOD ****'
	WRITE(*,*) '       ITERATION             ERR       '

	W=2.0/3.0*PI
    CALL CGM(M,N,P,GAP,Z0,W,DX,DY,ERROR_P,ERROR_L,ITERATION_NUM,M2,N2,ZZFFT)
	OPEN(100,FILE='P.DAT')
    DO J=1,N
		WRITE(100,FORMATT) (P(I,J), I=1,M)  !Dimensionless Pressure
	ENDDO
	CLOSE(100)

	WRITE(*,*)
	WRITE(*,*)'  **** EVALUATE SUB VON MISES STRESS ****'
    Z1=h_c2   !Computational domain in z direction for the coating
	Z2=h_c2   !Computational domain in z direction for the substrate
	NGaussQ=64
	alias=0
	EXP=4
    
    !Stress calculation
	CALL Multilayer_Shear_Stress(P,FRC,Stress_XYZ,M,N,M2,N2,L,L1,L2,NGaussQ,alias,exp,&
	                             XS,XE,YS,YE,Z1,Z2,E2_s,E2,EM2_s,EM2,h_c2,t)  !For body 2

	CALL VonMises_Stress_All(Stress_XYZ,VonMises_Stress,m,n,l)

	OPEN(6,FILE='VonMises_Stress.DAT')
	OPEN(7,FILE='StressXX.DAT')
	OPEN(8,FILE='StressYY.DAT')
	OPEN(9,FILE='StressZZ.DAT')
	OPEN(10,FILE='StressXY.DAT')
	OPEN(11,FILE='StressXZ.DAT')
 	OPEN(12,FILE='StressYZ.DAT')
		DO K=L,1,-1
			WRITE(6,FORMATT)(VonMises_Stress(I,(N+1)/2,K),I=1,M) !Dimensionless Stress in xoz plane (j=(N+1)/2)
			WRITE(7,FORMATT)(STRESS_XYZ(1,I,(N+1)/2,K),I=1,M)
			WRITE(8,FORMATT)(STRESS_XYZ(2,I,(N+1)/2,K),I=1,M)
			WRITE(9,FORMATT)(STRESS_XYZ(3,I,(N+1)/2,K),I=1,M)
			WRITE(10,FORMATT)(STRESS_XYZ(4,I,(N+1)/2,K),I=1,M)
 			WRITE(11,FORMATT)(STRESS_XYZ(5,I,(N+1)/2,K),I=1,M)
			WRITE(12,FORMATT)(STRESS_XYZ(6,I,(N+1)/2,K),I=1,M)         
		ENDDO
	CLOSE(6)
	CLOSE(7)
	CLOSE(8)
	CLOSE(9)
	CLOSE(10)
	CLOSE(11)
	CLOSE(12)

	!OPEN(7,FILE='Surface_Stress_XX.DAT')
	!OPEN(8,FILE='Interface_Stress_XZ.DAT')
	!	DO J=1,N
    !  		WRITE(7,FORMATT)(Stress_XYZ(1,I,J,1),I=1,M)
	!		WRITE(8,FORMATT)(Stress_XYZ(5,I,J,L1),I=1,M)
	!	ENDDO
    !CLOSE(7)
	!CLOSE(8)

	!--------------------FORMAT PART FOR OUTPUT AND DISPLAY------------------------------------
250 FORMAT(/,20X,'LOAD=',F10.5,' N',/,20X,'RX  =',F10.5,' MM',5X,&
			'RY =',F10.5,' MM',/,20X,'A0  =',F10.5,' MM',5X,'B0 =',F10.5,&
			' MM', /,20X,'AKE =',F10.5,8X,'P0 =',F10.5,' N/MM2',/,20X,&
			'FRC =',F10.5,/,20X) 
	
	END