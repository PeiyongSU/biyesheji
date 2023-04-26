

	!Calculate the discrete Fourier transform of the ICs (ICU) relating pressure to displacements
	SUBROUTINE Multilayer_Displacement_IC(M,N,M2,N2,ICU,NGaussQ,al,exp,XS,XE,YS,YE,E_s,E_c,Po_s,Po_c,h_c,t,Number)
	IMPLICIT NONE
    INTEGER NGaussQ,al,exp,M,N,M2,N2,Number
	REAL*8	XS,XE,YS,YE,dx,dy
	Double COMPLEX ICU(M2,N2)
	INTEGER NX,NY,Me,Ne,lx,ly,ND,i,j,k,l
	REAL*8	pi,area,wnyquist(2)
	REAL*8  wx0,wy0,wx,wy,dwx,dwy,Gpx,Gpy,FTPulse,erfc_r
    real*8 E_s,E_c,Po_s,Po_c,h_c,t(3)
	Double COMPLEX Multilayer_FRF_Surf

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

    CALL gauleg(-1.0d0,1.0d0,QX,QW,NGaussQ) !Calculate the points and weights for the Gauss formulas
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
								c=Multilayer_FRF_Surf(GPx,GPy,E_s,E_c,Po_s,Po_c,h_c,t,Number)
								sumY=sumY+QW(l)*c&
										*FTPulse(GPx,dx)*FTPulse(GPy,dy)/dx/dy
							ENDDO
							ctemp=ctemp+sumY*QW(k)
						ENDDO
						ctemp=ctemp/4.0
					ELSE
					    c=Multilayer_FRF_Surf(wx,wy,E_s,E_c,Po_s,Po_c,h_c,t,Number)
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
!	Anti-wrap around order
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
	DEALLOCATE (v_fAll,ICs)
	END

!*********FTPulse Eq.(26) (one dimension) in Wang et al. IJSS 2017 ******************************************************* 
	REAL*8 FUNCTION FTPulse(wX,Duration)
	Implicit none
	REAL*8 wX,Duration
	IF(wX.EQ.0.) THEN
		FTPulse=Duration
	ELSE
		FTPulse=2.*DSIN(wX*Duration/2.)/wX
	ENDIF
	END FUNCTION

!******Frequency response functions for surface displacements ***********************************************
!The formulas can be found in the paper (Wang et al., IJMS 2017).
	DOUBLE COMPLEX FUNCTION Multilayer_FRF_Surf(wx,wy,E_s,E_c,Po_s,Po_c,h_c,t,Number)   !E_s--Young's modulus of substrate; E_c--Young's modulus of coating
	implicit none                                               
	integer Number
	real*8 FRC,wx,wy,E_s,E_c,Po_s,Po_c,h_c	!Po_s--Poisson's ratio of substrate; Po_c--Poisson's ratio of coating; h_c--coating thickness.
	real*8 alpha,thet
    real*8 t(3),t1,t2,t3
    real*8 G_C,G_S,G_CS  !G_C--Shear modulus of coating; G_S--Shear modulus of substrate
	double complex B1,BB1,B1x,BB1x,B2,BB2,B2x,BB2x   !B1x--the derivative of B1 with respect to wx; BB1x--the derivative of BB1 with respect to wx;
    double complex BB1_T,BB1x_T,B2_T,B2x_T !BB1_T=exp(2*alpha*h_c)*BB1; BB1x_T=exp(2*alpha*h_c)*BB1x;B2_T=exp(alpha*h_c)*B2; B2x_T=exp(alpha*h_c)*B2x 
    double complex B,BB,Bx,BBx
	double complex S1,S2,S3,S4,S5,S6,S3_T,S4_T,S5_T,S6_T !S2_T=exp(2.0*alpha*h_c)*S2,S3_T=exp(alpha*h_c)*S3,S4_T=exp(alpha*h_c)*S4; S5_T=exp(alpha*h_c)*S5; S6_T=exp(alpha*h_c)*S6
    real*8 k1,k2,k3,k4,k5,k6,k7,k8,k1_T,k3_T,k5_T,k7_T,k10_T,k13_T,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18
    double complex A1,AA1,C1,CC1,AA1_T,CC1_T
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
 
    if(Number==3.OR.Number==9) then
        FRC=0.0
    elseif(Number==1.OR.Number==4) then
        FRC=1.0
    endif

    !When FRC=0, all B terms equal to 0    
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
   
    if(Number==3.OR.Number==9) then
        S1=-1.0-2.0*i*wx*(1.0-Po_c)*(B1+BB1)-i*alpha**2*(B1x+BB1x)
    	S1=S1/alpha
    elseif(Number==1.OR.Number==4) then
        S1=-2.0*i*wx*(1.0-Po_c)*(B1+BB1)-i*alpha**2*(B1x+BB1x)
    	S1=S1/alpha   
    endif

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
    
    CC1_T=(k10-k9/2.0/alpha)*(Sb_T-k5_T/alpha*S2)-(k13-k12/2.0/alpha)*(Sa_T-k1_T/alpha*S2)&
          +(k9*k13_T-k10_T*k12)*(S1-S2)/2.0/alpha
    CC1_T=CC1_T/(k15*k18-k16*k17)
    CC1=CC1_T*thet
    C1=(Sc-k16*CC1)/k15
      
    AA1_T=(k10_T*S1-Sa_T)/k15+(k16/k15-(4.0*Po_c-3.0))*CC1_T
    AA1_T=AA1_T/2.0/alpha
    AA1=AA1_T*thet
    A1=(S2+alpha*AA1-(1.0-2.0*po_c)*C1-CC1*(1.0-2.0*po_c))/alpha
    
    if(Number==9) then
    Multilayer_FRF_Surf=-alpha*(A1-AA1)-(3.0-4.0*Po_c)*(C1+CC1)
   	Multilayer_FRF_Surf=Multilayer_FRF_Surf/2.0/G_C
	return
	endif

	END FUNCTION
!********************************************************************************************