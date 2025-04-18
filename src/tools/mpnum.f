	subroutine mpnum(seq, n, marginp, numseq)
	integer n, numseq(n), i
	character*1, seq(n)
	double precision marginp(4)
	
	marginp(1)=0.0d0
	marginp(2)=0.0d0
	marginp(3)=0.0d0
	marginp(4)=0.0d0
	
	do i=1, n
	   if (seq(i). eq. 'A'.or. seq(i). eq. 'a') then
	     marginp(1)=marginp(1)+1.d0
	      numseq(i)=1
	   end if	
	   if (seq(i). eq. 'C'.or. seq(i). eq. 'c') then
	     marginp(2)=marginp(2)+1.d0
	      numseq(i)=2
	   end if
	   	
	   if (seq(i). eq. 'G'.or. seq(i). eq. 'g') then
	     marginp(3)=marginp(3)+1.d0
	      numseq(i)=3
	   end if
	   
	   if (seq(i). eq. 'U'.or. seq(i). eq. 'u'
     &		.or.seq(i). eq. 'T'.or. seq(i). eq. 't' ) then
	      marginp(4)=marginp(4)+1.d0
	      numseq(i)=4
	   end if
	   
	end do
	
	marginp(1)=marginp(1)/dble(n)
	marginp(2)=marginp(2)/dble(n)
	marginp(3)=marginp(3)/dble(n)
	marginp(4)=marginp(4)/dble(n)
	
	return
	end 					   		
	     
	      	
