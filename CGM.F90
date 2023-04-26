	
	!This subroutine is to use cgm to obtain the pressure distribution based on the
	!description in Polonsky and Keer, Wear 231(1999)
	SUBROUTINE CGM(M,N,P,GAP,S,W,DX,DY,ERRP,ERRL,ITERATION_NUM,M2,N2,ZZFFT)
	IMPLICIT NONE
	INTEGER M,N,M2,N2,ITERATION,OK,I,J,ITERATION_NUM,NE
	REAL*8  PI,W,AREA,PMIN,DELTA,G_OLD,LOAD,ERRP,ERRL,DX,DY
	REAL*8  APPROACH,ERR,AVEP,G_NEW,BETA,TOU1,TOU2,TOU,SE,SEM,CORR,GBAR,RBAR
	REAL*8  P(M,N),GAP(M,N),S(M,N)
	REAL*8,ALLOCATABLE::T(:,:),U(:,:),R(:,:),POLD(:,:)
	DOUBLE COMPLEX ZZFFT(M2,N2)

	PI  = 4.0*ATAN(1.0D0)
	PMIN  = 2.E-7		 
	DELTA = 0.0		 ! 'DELTA' in Polonsky and Keer, Wear 1999
	G_OLD = 1.0		 ! 'G' in Polonsky and Keer, Wear 1999
	ERR   = 1.0
	AREA  = DX*DY
	AVEP  = W/AREA
	LOAD  = 0.0
	
	ALLOCATE(T(M,N),U(M,N),R(M,N),POLD(M,N))
    
	!Initialize conjugate gradient direction as zero
	DO I=1,M
		DO J=1,N
			T(I,J)=0.0
		ENDDO
	ENDDO

	!Begin CGM iteration
	ITERATION=1
	DO WHILE (ERR.GT.ERRP.OR.ABS(LOAD-W)/W.GT.ERRL)
		IF(ITERATION > ITERATION_NUM) THEN
			WRITE(*,*) ' '
			WRITE(*,*) 'MAXIMUM ITERATION NUMBER IS REACHED!'
			WRITE(*,*) 'SORRY, IT IS DIFFCULT TO CONVERGENCE'
			RETURN
		ENDIF
	!Gap calculation
	!1> Deformation calculation base on fft technique
	CALL DEFORMATIONFFT(M,N,P,DX,DY,U,M2,N2,ZZFFT)
	
	!2> Obtain gap
		DO I=1,M
			DO J=1,N
				GAP(I,J)=S(I,J)+U(I,J)
			ENDDO
		ENDDO
		
	!Calculation average value in contact area
		CALL AVERAGEINCONTACTAREA(M,N,GAP,GBAR,P)
	
	!Obtain new gap by minus average value gbar
		DO I=1,M
			DO J=1,N
				GAP(I,J)=GAP(I,J)-GBAR
			ENDDO
		ENDDO
		APPROACH=GBAR
	
	!Apply the inequality constrains
		DO I=1,M
			DO J=1,N
				IF (P(I,J) ==0.0 .AND. GAP(I,J) < 0.0) THEN
					P(I,J)=PMIN
					DELTA=0.0
				ENDIF
			ENDDO
		ENDDO
		IF(DELTA == 0.0) THEN
			CALL AVERAGEINCONTACTAREA(M,N,GAP,GBAR,P)
			DO I=1,M
				DO J=1,N
					GAP(I,J)=GAP(I,J)-GBAR
				ENDDO
			ENDDO
			APPROACH=APPROACH+GBAR
		ENDIF

		G_NEW=0.0
		DO I=1,M
			DO J=1,N
				IF(P(I,J) > 0.0) THEN
					G_NEW=G_NEW+GAP(I,J)**2          
	   			ENDIF
			ENDDO
		ENDDO

		!Renew the conjugate gradient direction	tij
		BETA=DELTA*G_NEW/G_OLD
		G_OLD=G_NEW
		DELTA=1.0
		DO I=1,M
			DO J=1,N
				IF(P(I,J) > 0.0) THEN
					T(I,J)=GAP(I,J)+T(I,J)*BETA
				ELSE 
					T(I,J)=0.0
				END IF
			ENDDO
		ENDDO

		!Calculation the convolution between tij and influence coefficients
		CALL DEFORMATIONFFT(M,N,T,DX,DY,U,M2,N2,ZZFFT)
		
		!Obtain rij
		DO I=1,M
			DO J=1,N
  				R(I,J)=U(I,J)
			ENDDO
		ENDDO
		
		!Calculation the average value of rij
		CALL AVERAGEINCONTACTAREA(M,N,R,RBAR,P)
		!Renew the matrix rij
		DO I=1,M
			DO J=1,N
				R(I,J)=R(I,J)-RBAR
			ENDDO
		ENDDO

		TOU1=0.
		TOU2=0.
		DO I=1,M
			DO J=1,N
				IF(P(I,J) > 0.0) THEN
					TOU1=TOU1+GAP(I,J)*T(I,J)
					TOU2=TOU2+R(I,J)*T(I,J)
				ENDIF
			ENDDO
		ENDDO
		IF(TOU2.GT.0.0) THEN
			TOU=TOU1/TOU2
		ELSE
			TOU=0.5
		ENDIF
	
	!Store the current pij as pold for subsequent error estimation
		DO I=1,M
			DO J=1,N
				POLD(I,J)=P(i,J)
			ENDDO
		ENDDO

	
        !OK loop
		OK=0
		DO WHILE (OK == 0)	
			DO I=1,M
				DO J=1,N
					IF(P(I,J) > 0.0) THEN
						P(I,J)=P(I,J)-TOU*T(I,J)
					ENDIF
				ENDDO
			ENDDO

			DO I=1,M
				DO J=1,N
					IF(P(I,J) < 0.0) P(I,J)=0.0
				ENDDO
			ENDDO

			SE=0.0
			DO I=1,M
				DO J=1,N
					IF(P(I,J) > 0.0) THEN
						SE=SE+P(I,J)
					ENDIF
				ENDDO
			ENDDO
			CORR=AVEP/SE
			DO I=1,M
				DO J=1,N
					IF(P(I,J) > 0.0) THEN
						P(I,J)=P(I,J)*CORR
					ENDIF
				ENDDO
			ENDDO
			OK=1
		ENDDO  !END OK loop

		LOAD=0.0
		DO I=1,M
			DO J=1,N
				LOAD=LOAD+P(I,J)*AREA
			ENDDO
		ENDDO

		ERR=0.0
		DO I=1,M
			DO J=1,N
				ERR=ERR+ABS(P(I,J)-POLD(I,J))
			ENDDO
		ENDDO
		ERR=ERR*AREA/W
		IF(ITERATION == ITERATION/4*4) THEN
			WRITE(*,*) '',ITERATION,'       ',ERR
		ENDIF
		ITERATION=ITERATION+1
    ENDDO  
			
	WRITE(*,*) ' '
	WRITE(*,*) 'TOTAL ITERATION:', ITERATION

	DEALLOCATE(T,U,R,POLD)

    END

    
    SUBROUTINE AVERAGEINCONTACTAREA(M,N,A,AVE,P)
	IMPLICIT NONE
    INTEGER M,N,COUNT,I,J
    REAL*8 A(M,N),P(M,N)
	REAL*8 AVE

	AVE=0.0
	COUNT=0
    DO I=1,M
		DO J=1,N
			IF(P(I,J).GT.0.0) THEN
				AVE=AVE+A(I,J)          
				COUNT=COUNT+1
	   		END IF
		ENDDO
	ENDDO
	IF(COUNT.GT.0) THEN
		AVE=AVE/COUNT
	ELSE
		AVE=1.0E+10
	ENDIF
    
	END