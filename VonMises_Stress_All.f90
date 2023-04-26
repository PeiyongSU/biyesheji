!all subsurface VonMise stress
    SUBROUTINE 	VonMises_Stress_All(stress,VonMises_Stress,m,n,l)
	implicit none
    integer m,n,l,i,j,k
	real*8 stress(6,m,n,l),VonMises_Stress(m,n,l)
	do i=1,m
	  do j=1,n
	    do k=1,l
    	VonMises_Stress(i,j,k)=sqrt(((stress(1,i,j,k)-stress(2,i,j,k))**2+(stress(1,i,j,k)-stress(3,i,j,k))**2+&
	            (stress(2,i,j,k)-stress(3,i,j,k))**2+6*(stress(4,i,j,k)**2+stress(5,i,j,k)**2+stress(6,i,j,k)**2))/6.0)
		enddo
      enddo
	enddo

    END
    
    SUBROUTINE gauleg(x1,x2,x,w,n)
    INTEGER n
    REAL*8 x1,x2,x(n),w(n)
    DOUBLE PRECISION EPS
    PARAMETER (EPS=3.d-14)
    INTEGER i,j,m
    REAL*8 p1,p2,p3,pp,xl,xm,z,z1
    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    do 12 i=1,m
    z=cos(3.14159265358979323846d0*(i-.25d0)/(n+.5d0))
1       continue
        p1=1.d0
        p2=0.d0
        do 11 j=1,n
        p3=p2
        p2=p1
        p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11        continue
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
    if(abs(z-z1).gt.EPS)goto 1
    x(i)=xm-xl*z
    x(n+1-i)=xm+xl*z
    w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
    w(n+1-i)=w(i)
12    continue
    return
    END