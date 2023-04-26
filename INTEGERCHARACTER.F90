	!Transform integer to character
	function Int2Str(nNum,nLen)
    IMPLICIT REAL*8 (A-H,O-Z)
	CHARACTER*(4) Int2Str
		Int2Str=''
		nTemp=nNum
		IChZero=ichar('0')
		do 10 i=1,nLen
			j=10**(nLen-i)
			IBit=nTemp/j
			Int2Str(i:i)=char(IChZero+IBit)
			ntemp=ntemp-IBit*j
10		continue
	end function

	function Leng(nNum,nLen)
	IMPLICIT REAL*8 (A-H,O-Z)
		leng=0
		do i=1,nLen
			j=10**(nLen-i)
			IBit=nNum/j
			IF(IBit.NE.0) leng=leng+1
		enddo
	end function
