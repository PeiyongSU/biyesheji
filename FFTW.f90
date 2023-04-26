	!Two demensional fft---complex to complex 
	subroutine DFFT2D (M2,N2,in,notuse1,out,notuse2)  
 !   use, intrinsic :: iso_c_binding
     use MKL_DFTI
     implicit none
     integer M2,N2,notuse1,notuse2
     Integer :: Status
     type(DFTI_DESCRIPTOR), POINTER :: hand
     DOUBLE COMPLEX in(M2,N2),out(M2,N2)
     status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,[M2,N2])
     status = DftiCommitDescriptor(hand)
     status = DftiComputeForward(hand, in(:,1))
     out=in
     status = DftiFreeDescriptor(hand)
    end
    
	!Two dimensional inverse fft---complex to complex 
	subroutine DFFT2B (M2,N2,in,notuse1,out,notuse2)  
 !   use, intrinsic :: iso_c_binding
     use MKL_DFTI
     implicit none
     integer M2,N2,notuse1,notuse2
     Integer :: Status
     type(DFTI_DESCRIPTOR), POINTER :: hand
     DOUBLE COMPLEX in(M2,N2),out(M2,N2)
     status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,[M2,N2])
     status = DftiCommitDescriptor(hand)
     status = DftiComputeBackward(hand, in(:,1))
     out=in
     status = DftiFreeDescriptor(hand)
    end
    
	!Three demensional fft---complex to complex 
	subroutine DFFT3F (M2,N2,L2,in,notuse1,notuse2,out,notuse3,notuse4)  
 !   use, intrinsic :: iso_c_binding
     use MKL_DFTI
     implicit none
     integer M2,N2,L2,notuse1,notuse2,notuse3,notuse4
     Integer :: Status
     type(DFTI_DESCRIPTOR), POINTER :: hand
     DOUBLE COMPLEX in(M2,N2,L2),out(M2,N2,L2)
     status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,[M2,N2,L2])
     status = DftiCommitDescriptor(hand)
     status = DftiComputeForward(hand, in(:,1,1))
     out=in
     status = DftiFreeDescriptor(hand)
    end
    
    !Three demensional fft---complex to complex 
	subroutine DFFT3B (M2,N2,L2,in,notuse1,notuse2,out,notuse3,notuse4)  
 !   use, intrinsic :: iso_c_binding
     use MKL_DFTI
     implicit none
     integer M2,N2,L2,notuse1,notuse2,notuse3,notuse4
     Integer :: Status
     type(DFTI_DESCRIPTOR), POINTER :: hand
     DOUBLE COMPLEX in(M2,N2,L2),out(M2,N2,L2)
     status = DftiCreateDescriptor(hand, DFTI_DOUBLE, DFTI_COMPLEX, 2,[M2,N2,L2])
     status = DftiCommitDescriptor(hand)
     status = DftiComputeBackward(hand, in(:,1,1))
     out=in
     status = DftiFreeDescriptor(hand)
    end
    