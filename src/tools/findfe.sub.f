***** This is the F77 program for computing the free energy
***** of a given structure. The gcg connect file is chosen
***** to be the input file format because we need in this
***** program the sequence and the structure. Both bp and
***** helix formats present only the structure, thus we
***** will still need an additional input file for the seq.
***********************************************************

	subroutine f77findfe(n, numseq, fbp, h5, h3, hnum, seq,
     &		eall, pardir, pardirlen)

	integer n, i

	integer numseq(n), fbp(n),
     &		h5(n), h3(n), hnum(n),
     &		hindex

        double precision marginp(4),
     &		eall, estack, eext, ehp, ebg, eint, emloop

        character*1 seq(n)
c        character*80 a(300), b

        integer pardirlen
	character*(*) pardir

      	external mpnum, data, energy

c        print *, pardir(1:pardirlen)

	do i=1, n
	  if(seq(i). eq. 'T'. or. seq(i). eq .'t' 
     &		.or.seq(i). eq .'U'.or.seq(i). eq .'u' ) then
	    seq(i)='U'
c 2010-03-17 (Adam) Data integrity check
	    if (h3(i) .gt. n .or. h3(1) .lt. 0) then
            write(0, *) "invalid h3(" , i, ")=", h3(i)
        end if
	    if (h5(i) .gt. n .or. h5(i) .lt. 0) then
            write(0, *) "invalid h5(" , i, ")=", h5(i)
        end if
	    if (fbp(i) .gt. n .or. fbp(i) .lt. 0) then
            write(0, *) "invalid fbp(" , i, ")=", fbp(i)
        end if
	  endif

c    2004/12/20 (Clarence):
c       Replace 'N' by 'A' and its complement by 'T'. See
c       additional notes in ../notes/ directory for the reason
c       in choosing A to replace N.
c
	  if (seq(i). eq. 'A'. or. seq(i). eq .'a'
     &		.or. seq(i). eq. 'N'. or. seq(i). eq .'n') then
	     seq(i) ='A'
	  end if
	  
	  if (seq(i). eq. 'C'. or. seq(i). eq .'c') then
	     seq(i)='C'
	  end if
	  
	  if (seq(i). eq. 'G'. or. seq(i). eq .'g') then
	     seq(i)='G'
	  end if
	end do

c        print first 10 and last 10 nt and see input data if n >=10;
c               print thr wholes sequence if n < 10

c        write(0, *) ""
c        if (n.ge.10) then

c        do i=1, 10
c          write(0, *)  i, seq(i)
c        end do
c        do i=n-9, n
c           write(0, *) i, seq(i)
c        end do

c        else
c        do i=1, n
c          write(0, *)  i, seq(i)
c        end do

c        end if


	call mpnum(seq, n, marginp, numseq)
c
c-- Ivan Auger - 2/25/10 - do not call data if "pardirlen" is 0 or less
	if (pardirlen .gt. 0) then
	    call data(pardir, pardirlen)
	end if

c		    structure in Zuker structure input format and total free energy
c			reset initial values
	   hindex=0
	   do i=1, n
	      h5(i)=0
	      h3(i)=0
	      hnum(i)=0
	   end do
	   
	   do i=1, n

c		3 cases, fbp(i)> i, < i or =0 (neq 0 in first 2 cases)	   
c		for the helix with i > fbp(i)  it is a redundant copy
	    
	   if (i.gt. fbp(i). or. fbp(i).eq. 0) then
	      go to 99
	   end if 
c					new helix for i=1	     	   
	   if (i.eq. 1. and. fbp(i).ne.0) then	         
	      hindex=1
	      h5(1)=1
	      h3(1)=fbp(1)
	      hnum(1)=1
	      go to 99	         
	   end if

c					new helix for any i

           if (i.ge.2) then

	     if (fbp(i-1).eq.0. and. fbp(i).ne.0) then
	       hindex=hindex+1
	       h5(hindex)=i
	       h3(hindex)=fbp(i)
	       hnum(hindex)=1
	       go to 99
	     end if

c		can be continuing helix or a new helix if fbp(i-1) is the end
c		of a "reduntant" helix or it is a bulge	   
	   
	     if (fbp(i-1).ne.0. and. fbp(i).ne.0) then
	     
	        if ((i-1).lt. fbp(i-1).and.
     & 			(fbp(i-1).eq. (fbp(i)+1))) then
	          hnum(hindex)=hnum(hindex)+1
	        else
	          hindex=hindex+1
	          h5(hindex)=i
	          h3(hindex)=fbp(i)
	          hnum(hindex)=1
	        end if
	        
	      end if
	    end if

99	   continue
	   
	   end do
	   

c			compute free energy of the sampled structure	   

	call energy(n, 	numseq, h5, h3, hnum,
     &		eall, estack, eext, ehp, ebg, eint, emloop)

c       2010-02-15 (Adam) Perform some operation on eall. print eall or divide
c          eall by 1. This tells the f90 -fast compiler that we are indeed
c          interested in this variable, don't optimize it away.
c          I dare you to remove it.
        eall = eall / 1.0
c	fe(isample)=eall

c        write(0, *) ""
cc        write(0, *) "Sequence file = ", infile
c        write(0, 4) "Sequence length = ", n
c        write(0, 5) "Free energy     = ", eall
c        write(0, *) ""

c4	format(A, i9)
c5	format(A, f9.2)

	return
c	stop
	end


	include 'mpnum.f'
	include 'data.f'
	include 'energy.f'
