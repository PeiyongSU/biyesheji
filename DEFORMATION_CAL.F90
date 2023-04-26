	!Calculate displacement V caused by pressure P based on fft algorithm
    !ZZFFT is the discrete Fourier transform of the IC relating pressure to surface displacement
	SUBROUTINE DEFORMATIONFFT(M,N,P,DX,DY,V,M2,N2,ZZFFT)
	IMPLICIT NONE
    INTEGER M,N,M2,N2	
	REAL*8 P(M,N),V(M,N)
	DOUBLE COMPLEX ZZFFT(M2, N2)
	REAL*8, ALLOCATABLE::DEF(:,:)
	DOUBLE COMPLEX, ALLOCATABLE:: PPFFT(:,:), DFFT(:,:)
	REAL*8 DX,DY,PI,C0
	PI  = 4.0*ATAN(1.0D0)

	!M--the number of mesh in x in extended computation domain
	!N--the number of mesh in y in extended computation domain
	
	!C0 is a constant for fft
   	C0=1.0/M2/N2
	!Allocate temporary matrix for fft
	ALLOCATE (PPFFT(M2,N2), DFFT(M2,N2))
	
	!Pressure zero-padding and perform fft to pressure
	CALL PRESSUREFFT(M,N,P,M2,N2,PPFFT)
	
	!Perform piecewise multiplicaton in frequence domain (convolution theorem)
	CALL DEFORMFFT(M2,N2,ZZFFT,PPFFT,DFFT)
	
	!Perform inverse fft to obtain the deformation in computation 
	!Domain by discarding the term beyond computation domain
	CALL DEFORMATION(C0,M2,N2,DFFT,M,N,V)

	DEALLOCATE (PPFFT, DFFT)

	END

	SUBROUTINE PRESSUREFFT(M,N,P,M2,N2,PPFFT)
	REAL*8 P(M,N)
	DOUBLE COMPLEX PPFFT(M2,N2)
	DOUBLE COMPLEX,ALLOCATABLE:: P_RECT(:,:)
	ALLOCATE (P_RECT(M2,N2))

	!Zero-pad to pressure in extended domain and convert into complex type
	DO I=1,M2
		DO J=1,N2
			IF(I.LE.M .AND. J.LE.N) THEN
				P_RECT(I,J)=CMPLX(P(I,J),0.0)
			ELSE
				P_RECT(I,J)=CMPLX(0.0, 0.0)
			ENDIF
          ENDDO
	ENDDO
	
	!Perfom fft the obtained pressure and store in complex array ppfft 
	CALL DFFT2D (M2,N2,P_RECT,M2, PPFFT,M2)
	DEALLOCATE (P_RECT)

	END

	SUBROUTINE DEFORMFFT(M2,N2,ZZFFT,PFFT,DFFT)
	DOUBLE COMPLEX ZZFFT(M2,N2)
	DOUBLE COMPLEX PFFT(M2,N2)
	DOUBLE COMPLEX DFFT(M2,N2)
	
	!Convolution theorem (perform piece-wise multiplication)
	!To obtain the deformations in frequence domain. 
	DO I=1,M2
		DO J=1,N2
            DFFT(I,J)=PFFT(I,J)*ZZFFT(I,J)
		ENDDO
	ENDDO

	END

	SUBROUTINE DEFORMATION(C0,M2,N2,DFFT,M,N,DEF)
	IMPLICIT NONE
	INTEGER M,N,M2,N2,I,J
	REAL*8 C0
	REAL*8 DEF(M,N)
	DOUBLE COMPLEX DFFT(M2,N2)
	DOUBLE COMPLEX,ALLOCATABLE:: DEFORM(:,:)
	ALLOCATE (DEFORM(M2,N2))
	
	!Perform inverse fft
	CALL DFFT2B(M2,N2,DFFT,M2,DEFORM,M2)
	DO I=1,M
		DO J=1,N
			DEF(I,J)=REAL(DEFORM(I,J))*C0
		ENDDO
	ENDDO
	DEALLOCATE (DEFORM)

	END

