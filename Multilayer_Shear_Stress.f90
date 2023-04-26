   ! Calculate stresses caused by pressure P and shear traction (FRC*P)
	Subroutine Multilayer_Shear_Stress(P,FRC,Stress_XYZ,M,N,M2,N2,L,L1,L2,NGaussQ,al,exp,&
                                       XS,XE,YS,YE,Z1,Z2,E_s,E_c,Po_s,Po_c,h_c,t)
	implicit none
    integer NGaussQ,al,exp,M,N,M2,N2,L1,L2,L,layer,strN,SectionN
	real*8  P(M,N),Stress_XYZ(6,M,N,L)
	real*8	XS,XE,YS,YE,dx,dy,Z1,Z2,zPosition
	integer i,j,k
    real*8 FRC,E_s,E_c,Po_s,Po_c,h_c,t(3)
	double complex,allocatable::P_exp(:,:),P_fft(:,:),IC_S(:,:)
	double complex,allocatable:: stressFFT(:,:),stress(:,:)
	allocate (P_exp(M2,N2),P_fft(M2,N2),stressFFT(M2,N2),stress(M2,N2),IC_S(M2,N2))

	DO I=1,M2
		DO J=1,N2
			IF(I.LE.M .AND. J.LE.N) THEN
				P_exp(I,J)=CMPLX(P(I,J),0.0)
			ELSE
				P_exp(I,J)=CMPLX(0.0, 0.0)
			ENDIF
        ENDDO
	ENDDO
	CALL DFFT2D (M2,N2,P_exp,M2,P_FFT,M2)

    layer=1
	DO SectionN=1,L1
	   DO StrN=1,6
	   zPosition=(SectionN-1)*Z1/(L1-1.0)
       Call Multilayer_Shear_Stress_IC(FRC,M-1,N-1,M2,N2,IC_s,NGaussQ,al,exp,XS,XE,YS,YE,E_s,E_c,Po_s,Po_c,h_c,t,layer,zPosition,strN)
       DO i=1,M2
		  DO j=1,N2
			StressFFT(i,j)=IC_s(i,j)*p_fft(i,j)
		  ENDDO
	   ENDDO
	   CALL DFFT2B(M2,N2,stressFFT,M2,stress,M2)
	   DO i=1,M
		  DO j=1,N
		    Stress_XYZ(StrN,i,j,SectionN)=real(stress(i,j))/m2/n2
		   ENDDO
	   ENDDO
	   ENDDO
    ENDDO

    layer=2
	DO SectionN=1,L2
	   DO StrN=1,6
	   zPosition=SectionN*Z2/L2
       Call Multilayer_Shear_Stress_IC(FRC,M-1,N-1,M2,N2,IC_s,NGaussQ,al,exp,XS,XE,YS,YE,E_s,E_c,Po_s,Po_c,h_c,t,layer,zPosition,strN)
       DO i=1,M2
		  DO j=1,N2
			StressFFT(i,j)=IC_s(i,j)*p_fft(i,j)
		  ENDDO
	   ENDDO
	   CALL DFFT2B(M2,N2,stressFFT,M2,stress,M2)
	   DO i=1,M
		  DO j=1,N
		    k=L1+SectionN
		    Stress_XYZ(StrN,i,j,k)=real(stress(i,j))/m2/n2
		   ENDDO
	   ENDDO
	   ENDDO
    ENDDO
	
	DEALLOCATE (P_exp,P_fft,stress,stressFFT,IC_S)
	END

	!Calculate the discrete Fourier transform of the ICs (ICU) relating pressure to stresses
	SUBROUTINE Multilayer_Shear_Stress_IC(FRC,M,N,M2,N2,ICU,NGaussQ,al,exp,XS,XE,YS,YE,E_s,E_c,Po_s,Po_c,h_c,t,layer,zPosition,strN)
	IMPLICIT NONE
    INTEGER NGaussQ,al,exp,M,N,M2,N2,layer,strN
	REAL*8	XS,XE,YS,YE,dx,dy,zPosition
	Double COMPLEX ICU(M2,N2)
	INTEGER NX,NY,Me,Ne,lx,ly,ND,i,j,k,l
	REAL*8	pi,area,wnyquist(2)
	REAL*8  wx0,wy0,wx,wy,dwx,dwy,Gpx,Gpy,FTPulse,erfc_r
    real*8 FRC,E_s,E_c,Po_s,Po_c,h_c,t(3)
	Double COMPLEX FRF_Shear_Stress

	Double COMPLEX	sumY,ctemp,frf,c,temp
	Double COMPLEX,ALLOCATABLE:: v_fAll(:,:)
	REAL*8,ALLOCATABLE::ICS(:,:),QW(:),QX(:),QXFIX(:)
	ALLOCATE (v_fAll(M*exp,N*exp),ICS(M2,N2),QW(NGaussQ), QX(NGaussQ), QXFIX(2))

    PI  = 4.0*ATAN(1.0D0)

	dx=(XE-XS)/M
	dy=(YE-YS)/N 
	area=dx*dy

	ND=2
    Me=M*exp
	Ne=N*exp	
	NX=2*M
	NY=2*N
	wnyquist(1)=2.*pi/dx
	wnyquist(2)=2.*pi/dy
	dwx=wnyquist(1)/(Me-1)
	dwy=wnyquist(2)/(Ne-1)
    CALL gauleg(-1.0d0,1.0d0,QX,QW,NGaussQ)  !Calculate the points and weights for the Gauss formulas
!	CALL DGQRUL (NGaussQ, 1, DBLE(0.0), DBLE(0.0), 0, QXFIX, QX, QW)

    DO i=1,Me	
		DO j=1,Ne
 			IF(i.gt.Me/2+1) THEN
				wx0=-dwx*(Me-i+1) 
			ELSE
				wx0= dwx*(i-1)
			ENDIF
			IF(j.gt.Ne/2+1) THEN
				wy0=-dwy*(Ne-j+1)
			ELSE
				wy0= dwy*(j-1)
			ENDIF
			!Calculate FRF
			frf=0.0
			DO lx=-al,al
				DO ly=-al,al
					wx=wx0-lx*wnyquist(1)
					wy=wy0-ly*wnyquist(2)
					IF(DABS(wx).le.dwx.and.DABS(wy).le.dwy) THEN !Calculate FRF when wx=wy=0
						ctemp=0.0
						DO k=1,NGaussQ
							sumY=0.0
							DO l=1,NGaussQ
								GPx=wx+dwx/2.*QX(k)
								GPy=wy+dwy/2.*QX(l)
								c=FRF_Shear_Stress(FRC,GPx,GPy,E_s,E_c,Po_s,Po_c,h_c,t,layer,zPosition,strN)
								sumY=sumY+QW(l)*c&
										*FTPulse(GPx,dx)*FTPulse(GPy,dy)/dx/dy
							ENDDO
							ctemp=ctemp+sumY*QW(k)
						ENDDO
						ctemp=ctemp/4.0
                    ELSE
					    c=FRF_Shear_Stress(FRC,wx,wy,E_s,E_c,Po_s,Po_c,h_c,t,layer,zPosition,strN)
						ctemp=c*FTPulse(wx,dx)*FTPulse(wy,dy)/dx/dy
					ENDIF
					frf=frf+ctemp
				ENDDO
			ENDDO
			v_fAll(i,j)=frf
		ENDDO
	ENDDO

	CALL DFFT2B (ME,NE,v_fALL,ME,v_fALL,ME)

	DO i=1,Me
		DO j=1,Ne
			v_fALL(i,j)=v_fALL(i,j)/Me/Ne
		ENDDO
	ENDDO
	!Anti-wrap around order
	DO i=1,Me/2
		DO j=1,Ne
			temp=v_fALL(i,j)
			v_fALL(i,j)=v_fALL(i+Me/2,j)
			v_fALL(i+Me/2,j)=temp	  
		ENDDO
	ENDDO
	DO i=1,Me
		DO j=1,Ne/2
			temp=v_fALL(i,j)
			v_fALL(i,j)=v_fALL(i,j+Ne/2)
			v_fALL(i,j+Ne/2)=temp
		ENDDO
	ENDDO
	DO i=1,M*ND
		DO j=1,N*ND
			ICs(i,j)=real(v_fALL(i+Me/2-M,j+Ne/2-N))
		ENDDO
    ENDDO

    !Wrap around order
	DO i=1,M*ND/2
		DO j=1,N*ND
			temp=ICs(i,j)
			ICs(i,j)=ICs(i+M*ND/2,j)
			ICs(i+M*ND/2,j)=temp	  
		ENDDO
	ENDDO
	DO i=1,M*ND
		DO j=1,N*ND/2
			temp=ICs(i,j)
			ICs(i,j)=ICs(i,j+N*ND/2)
			ICs(i,j+N*ND/2)=temp
		ENDDO
	ENDDO
	DO i=1,M2
		DO j=1,N2
			ICU(i,j)=ICs(i,j)
		ENDDO
	ENDDO
	CALL DFFT2D (M2,N2,ICU,M2,ICU,M2)

    END
    
!******Frequency response functions for Stresses***********************************************
!The formulas can be found in the paper (Wang et al., IJMS 2017).
	DOUBLE COMPLEX FUNCTION FRF_Shear_Stress(FRC,wx,wy,E_s,E_c,Po_s,Po_c,h_c,t,layer,z,strN)  !E_s--Young's modulus of substrate; E_c--Young's modulus of coating
	implicit none
	integer layer,strN
	real*8  z
	real*8 FRC,wx,wy,E_s,E_c,Po_s,Po_c,h_c,Po,t(3) 	!Po_s--Poisson's ratio of substrate; Po_c--Poisson's ratio of coating; h_c--coating thickness.
	real*8 alpha,thet
    real*8 t1,t2,t3
	real*8 G_C,G_S,G_CS,G  !G_C--Shear modulus of coating; G_S--Shear modulus of substrate
	double complex B1,BB1,B1x,BB1x,B2,BB2,B2x,BB2x   !B1x--the derivative of B1 with respect to wx; BB1x--the derivative of BB1 with respect to wx;
    double complex BB1_T,BB1x_T,B2_T,B2x_T !BB1_T=exp(2*alpha*h_c)*BB1; BB1x_T=exp(2*alpha*h_c)*BB1x;B2_T=exp(alpha*h_c)*B2; B2x_T=exp(alpha*h_c)*B2x 
    double complex B,BB,Bx,BBx
	double complex A1,AA1,C1,CC1,AA1_T,CC1_T,A2,A2_T,AA2,C2,C2_T,CC2,A,AA,C,CC !AA1_T=exp(2*alpha*h_c)*AA1; CC1_T=exp(2*alpha*h_c)*CC1; A2_T=exp(alpha*h_c)*A2; C2_T=exp(alpha*h_c)*C2
    double complex S1,S2,S3,S4,S5,S6,S3_T,S4_T,S5_T,S6_T !S2_T=exp(2.0*alpha*h_c)*S2; S3_T=exp(alpha*h_c)*S3; S4_T=exp(alpha*h_c)*S4; S5_T=exp(alpha*h_c)*S5; S6_T=exp(alpha*h_c)*S6
    real*8 k1,k2,k3,k4,k5,k6,k7,k8,k1_T,k3_T,k5_T,k7_T,k10_T,k13_T,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18
    double complex Sa,Sb,Sa_T,Sb_T,Sc,Sd !Sb_T=exp(2*alpha*h_c)*Sb
    double complex i
    i=(0,1.0)
    
    t1=t(1)
    t2=t(2)
    t3=t(3)

	alpha=dsqrt(wx*wx+wy*wy)
	thet=dexp(-2.0*alpha*h_c)
    G_C=E_c/2.0/(1.0+Po_c)
	G_S=E_s/2.0/(1.0+Po_s)    
    G_CS=E_c*(1.0+Po_s)/E_s/(1.0+Po_c)

    BB1_T=FRC*(t1-alpha*G_C-t1*G_CS)/(2.0*alpha*(1.0-Po_c)*((t1-alpha*G_C-t1*G_CS)*thet+t1+alpha*G_C+t1*G_CS))
	BB1=BB1_T*thet
	B1=BB1-FRC/(2.0*alpha*(1.0-Po_c))

    BB1x_T=FRC*wx/(2.0*alpha**2*(1.0-Po_c))
	BB1x_T=BB1x_T*(-2.0*t1*G_C-(2.0*h_c+1.0/alpha)*(t1-alpha*G_C-t1*G_CS)*(t1+alpha*G_C+t1*G_CS)-thet/alpha*(t1-alpha*G_C-t1*G_CS)**2)   
    BB1x_T=BB1x_T/((t1-alpha*G_C-t1*G_CS)*thet+t1+alpha*G_C+t1*G_CS)**2
	BB1x=BB1x_T*thet
    B1x=BB1x+FRC*wx/(2.0*(1.0-Po_c)*alpha**3)

	B2_T=(1.0-Po_c)*(thet-1.0)/(1.0-Po_s)*BB1_T-FRC/2.0/alpha/(1.0-Po_s)
	B2=B2_T*dexp(-alpha*h_c)
    B2x_T=-wx*h_c*(1.0-Po_c)*(thet+1.0)/alpha/(1.0-Po_s)*BB1_T &
          +(1.0-Po_c)*(thet-1.0)/(1.0-Po_s)*BB1x_T &
          +FRC*(1.0+alpha*h_c)*wx/2.0/alpha/alpha/alpha/(1.0-Po_s)
	B2x=B2x_T*dexp(-alpha*h_c)
    BB2x=0.0 
    
    S1=-1.0-2.0*i*wx*(1.0-Po_c)*(B1+BB1)-i*alpha**2*(B1x+BB1x)
	S1=S1/alpha
	S2=-wx*(B1-BB1)-alpha*alpha*(B1x-BB1x)
    S2=S2*i/alpha
	S3_T=-2.0*wx*h_c*alpha*(B1+BB1_T)+wx*(2.0+t1*h_c/G_C)*(B1-BB1_T)-t1*alpha/G_C*(B1x+BB1x_T)  &
         +2.0*alpha**2*(B1x-BB1x_T)+t1*alpha/G_S*B2x_T
    S3_T=S3_T*i/alpha
    S3=S3_T*dexp(-alpha*h_c)
    S4_T=-wx*alpha*(t3*h_c/G_C+4.0*(1.0-Po_c))*(B1+BB1_T)+wx*(t3/G_C+2.0*alpha**2*h_c)*(B1-BB1_T)&
         +t3*alpha**2/G_C*(B1x-BB1x_T)-2.0*alpha**3*(B1x+BB1x_T)-t3*wx/G_S*B2_T-t3*alpha**2/G_S*B2x_T
    S4_T=S4_T*i/alpha
    S4=S4_T*dexp(-alpha*h_c)
	S5_T=(1.0-alpha*h_c)*wx*B1-(1.0+alpha*h_c)*wx*BB1_T+alpha**2*(B1x-BB1x_T)-wx*B2_T-alpha**2*B2x_T
	S5_T=-S5_T*i/alpha
    S5=S5_T*dexp(-alpha*h_c)
	S6_T=wx*alpha*h_c*(B1-BB1_T)-2.0*wx*(1.0-po_c)*(B1+BB1_T)-alpha*alpha*(B1x+BB1x_T)+2.0*wx*(1.0-Po_s)*B2_T+alpha**2*B2x_T
    S6_T=S6_T*i/alpha
    S6=S6_T*dexp(-alpha*h_c)

    k1_T=(t3*G_CS*alpha-t3*alpha+2.0*G_C*alpha**2)/(2.0-2.0*Po_s)
    k1=k1_T*thet
    k2=(t3*alpha+2.0*G_C*alpha**2+(3.0-4.0*Po_s)*t3*G_CS*alpha)/(2.0-2.0*Po_s)
    k3_T=(t3*G_CS*(3.0-2.0*Po_c-2.0*Po_s+alpha*h_c)-t3*(3.0-4.0*Po_c)+4.0*G_C*alpha*(1.0-Po_c)+alpha*h_c*(2.0*G_C*alpha-t3))/(2.0-2.0*Po_s)   
    k3=k3_T*thet
    k4=(alpha*h_c*(2.0*G_C*alpha+t3)-t3*(3.0-4.0*Po_c)-4.0*G_C*alpha*(1.0-Po_c)-t3*G_CS*(5.0-6.0*Po_s-(3.0-4.0*Po_s)*(2.0*Po_c+alpha*h_c)))/(2.0-2.0*Po_s)
    k5_T=t1*alpha-2.0*G_C*alpha**2-t1*G_CS*alpha
    k5=k5_T*thet
    k6=t1*alpha+2.0*G_C*alpha**2+t1*G_CS*alpha*(3.0-4.0*Po_s)
    k7_T=alpha*h_c*(t1-2.0*G_C*alpha)-2.0*G_C*alpha*(1.0-2.0*Po_c)+t1*G_CS*(2.0*Po_c-2.0*Po_s-alpha*h_c)    
    k7=k7_T*thet
    k8=alpha*h_c*(t1+2.0*G_C*alpha)-2.0*G_C*alpha*(1.0-2.0*Po_c)-t1*G_CS*(2.0*(1.0-2.0*Po_s)*(1.0-2.0*Po_c-alpha*h_c)+(2.0-2.0*Po_c-2.0*Po_s-alpha*h_c))   
    k9=k1+k2
    k10_T=k3_T-k1_T/alpha*(1.0-2.0*Po_c)
    k10=k10_T*thet
    k11=k4-k1/alpha*(1.0-2.0*Po_c)
    k12=k5+k6
    k13_T=k7_T-k5_T/alpha*(1.0-2.0*Po_c)
    k13=k13_T*thet
    k14=k8-k5/alpha*(1.0-2.0*Po_c)
    k15=k10-k9/2.0/alpha
    k16=k11-k9/2.0/alpha*(4.0*Po_c-3.0)
    k17=k13-k12/2.0/alpha
    k18=k14-k12/2.0/alpha*(4.0*Po_c-3.0)
    
    Sa_T=G_C/(2.0-2.0*Po_s)*S4_T-t3*G_CS*(1.0-2.0*Po_s)/(2.0-2.0*Po_s)*S5_T+t3*G_CS*S6_T
    Sa=Sa_T*thet
    Sb_T=alpha*G_C*S3_T-t1*G_CS*(2.0-2.0*Po_s)*S5_T+t1*G_CS*(1.0-2.0*Po_s)*S6_T
    Sb=Sb_T*thet
    Sc=Sa-k1/alpha*S2-k9/2.0/alpha*(S1-S2)
    Sd=Sb-k5/alpha*S2-k12/2.0/alpha*(S1-S2)
    
    if(layer==1) then
	  CC1_T=(k10-k9/2.0/alpha)*(Sb_T-k5_T/alpha*S2)-(k13-k12/2.0/alpha)*(Sa_T-k1_T/alpha*S2)&
            +(k9*k13_T-k10_T*k12)*(S1-S2)/2.0/alpha
      CC1_T=CC1_T/(k15*k18-k16*k17)
      CC1=CC1_T*thet
      C1=(Sc-k16*CC1)/k15
      
      AA1_T=(k10_T*S1-Sa_T)/k15+(k16/k15-(4.0*Po_c-3.0))*CC1_T
      AA1_T=AA1_T/2.0/alpha
      AA1=AA1_T*thet
      A1=(S2+alpha*AA1-(1.0-2.0*po_c)*C1-CC1*(1.0-2.0*po_c))/alpha
      
      B=B1*dexp(-alpha*z)
!	  BB=BB1*dexp(alpha*z)
	  BB=BB1_T*dexp(-2.0*alpha*h_c+alpha*z)    !Avoid exponential blow up
	  Bx=B1x*dexp(-alpha*z)
!	  BBx=BB1x*dexp(alpha*z)
	  BBx=BB1x_T*dexp(-2.0*alpha*h_c+alpha*z)   !Avoid exponential blow up

  	  A=A1*dexp(-alpha*z)
!	  AA=AA1*dexp(alpha*z)
      AA=AA1_T*dexp(-2.0*alpha*h_c+alpha*z)     !Avoid exponential blow up

	  C=C1*dexp(-alpha*z)
!	  CC=CC1*dexp(alpha*z)
	  CC=CC1_T*dexp(-2.0*alpha*h_c+alpha*z)     !Avoid exponential blow up

	  Po=Po_c
      G=G_C
      
    else if(layer==2) then
      CC1_T=(k10-k9/2.0/alpha)*(Sb_T-k5_T/alpha*S2)-(k13-k12/2.0/alpha)*(Sa_T-k1_T/alpha*S2)&
            +(k9*k13_T-k10_T*k12)*(S1-S2)/2.0/alpha
      CC1_T=CC1_T/(k15*k18-k16*k17)
      CC1=CC1_T*thet
      C1=(Sc-k16*CC1)/k15
      
      AA1_T=(k10_T*S1-Sa_T)/k15+(k16/k15-(4.0*Po_c-3.0))*CC1_T
      AA1_T=AA1_T/2.0/alpha
      AA1=AA1_T*thet
      A1=(S2+alpha*AA1-(1.0-2.0*po_c)*C1-CC1*(1.0-2.0*po_c))/alpha
     
	  A2_T=(t1*alpha-2.0*G_C*alpha**2)*A1+(t1*alpha+2.0*G_C*alpha**2)*AA1_T+((t1-2.0*G_C*alpha)*alpha*h_c-2.0*G_C*alpha*(1.0-2.0*Po_c))*C1 &
          +(alpha*h_c*(t1+2.0*G_C*alpha)-2.0*G_C*alpha*(1.0-2.0*Po_c))*CC1_T-alpha*G_C*S3_T
      A2_T=A2_T/t1/G_CS/alpha
      A2=A2_T*dexp(-alpha*h_c)
      C2_T=alpha*A1+alpha*AA1_T+(2.0*(1.0-Po_c)+alpha*h_c)*C1-(2.0*(1.0-Po_c)-alpha*h_c)*CC1_T-alpha*A2_T-S6_T
      C2_T=C2_T/(2.0-2.0*Po_s)
      C2=C2_T*dexp(-alpha*h_c)

      B=B2*dexp(-alpha*z)
	  BB=0.0 !BB2*dexp(alpha*z)
	  Bx=B2x*dexp(-alpha*z)
	  BBx=0.0 !BB2x*dexp(alpha*z)
	  A=A2*dexp(-alpha*z)
	  AA=0.0
	  C=C2*dexp(-alpha*z)
	  CC=0.0

	  Po=Po_s
      G=G_S
	endif

    GOTO (1,2,3,4,5,6), strN
1   FRF_Shear_Stress=2.0*i*(Po-2.0)*wx*(B+BB)+i*wx**3*z/alpha*(B-BB)-i*wx**2*(Bx+BBx)&
                     -wx*wx*(A+AA)+2.0*alpha*Po*(C-CC)-z*wx*wx*(C+CC)
    RETURN
2   FRF_Shear_Stress=-2.0*i*Po*wx*(B+BB)+i*z*wy**2*wx/alpha*(B-BB)-i*wy**2*(Bx+BBx)&
                     -wy*wy*(A+AA)+2.0*alpha*Po*(C-CC)-z*wy*wy*(C+CC)
	RETURN
3   FRF_Shear_Stress=2.0*i*(1.0-Po)*wx*(B+BB)-i*z*alpha*wx*(B-BB)+i*alpha**2*(Bx+BBx)&
                     +alpha*alpha*(A+AA)+2.0*(1.0-Po)*alpha*(C-CC)+z*alpha*alpha*(C+CC)
    RETURN               
4   FRF_Shear_Stress=-2.0*i*wy*(1.0-Po)*(B+BB)+i*wx**2*wy*z/alpha*(B-BB)-i*wx*wy*(Bx+BBx)&
                     -wx*wy*(A+AA)-z*wx*wy*(C+CC)
	RETURN
5   FRF_Shear_Stress=-wx**2*z*(B+BB)+(2.0*alpha*(1.0-Po)+wx**2/alpha)*(B-BB)+alpha*wx*(Bx-BBx)&
                     -i*(wx*alpha*(A-AA)+wx*(1.0-2.0*Po)*(C+CC)+z*wx*alpha*(C-CC))
	RETURN
6   FRF_Shear_Stress=(1.0/alpha-z)*wx*wy*B-(1.0/alpha+z)*wx*wy*BB+alpha*wy*(Bx-BBx)&
                     -i*(wy*alpha*(A-AA)+(1.0-2.0*Po)*wy*(C+CC)+z*wy*alpha*(C-CC))
	RETURN

	END FUNCTION
!********************************************************************************************!
