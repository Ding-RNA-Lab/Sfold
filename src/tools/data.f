c

        
	subroutine data(pardir, pardirlen)

        character*(*) pardir
        integer pardirlen
	character*20, a1
       	character*6, a
	character*5, a2
	
	
	character*80, b, c*120
	double precision bonus, eh3(4,4,4,4,4), eh4(4,4,4,4,4,4),
     &		e(16, 16), es(4, 4, 4, 4), ed3(4,4,4), ed5(4,4,4),
     &		ehs(4, 4, 4, 4), eis(4, 4, 4, 4),
     &		eh(30), eb(30), ei(30),  e21(96, 24), e22(576,16),
     &		e11(24, 24), ei11(4,4,4,4,4,4), 
     &		ei22a(4,4,4,4,4,4,4), ei22c(4,4,4,4,4,4,4),
     &		ei22g(4,4,4,4,4,4,4), ei22u(4,4,4,4,4,4,4),     
     &		e21s(4,4,4,4,4,4,4),  e21d(4,4,4,4,4,4,4)
     
	double precision  eh3ex(4,4,4,4,4), eh4ex(4,4,4,4,4,4),
     &		 esex(4, 4, 4, 4), ed3ex(4,4,4), ed5ex(4,4,4),
     &		ehsex(4, 4, 4, 4), eisex(4, 4, 4, 4),
     &		ehex(30), ebex(30), eiex(30),  
     &		 ei11ex(4,4,4,4,4,4), 
     &		ei22aex(4,4,4,4,4,4,4), ei22cex(4,4,4,4,4,4,4),
     &		ei22gex(4,4,4,4,4,4,4), ei22uex(4,4,4,4,4,4,4),     
     &		e21sex(4,4,4,4,4,4,4),  e21dex(4,4,4,4,4,4,4),
     &		scale, am, bm, cm, slope, 
     &		tphb, hgu, h3c, ha, hb, iasy, 
     &		aex, bex, cex, slopeex, 
     &		tphbex, hguex, h3cex, haex, hbex, iasyex, ex3     
          
	integer baseid, i1,i2,i3,i4,i5,i6, i7, i8
	integer i,iblock,j,m,l,k,irb,icb,ilb,isb
	
	common /stack/ es /dangle/ ed3, ed5 /HBI/ eh, eb, ei
     &		/Hstack/ ehs /Istack/ eis /int11/ ei11
     &		/int21/ e21s, e21d/ int22
     &		/ ei22a,  ei22c, ei22g, ei22u 
     &		/tetra/ eh4 / tri /eh3	
     
	common /stackex/ esex /dangleex/ ed3ex, ed5ex 
     &		/HBIex/ ehex, ebex, eiex
     &		/Hstackex/ ehsex /Istackex/ eisex /int11ex/ ei11ex
     &		/int21ex/ e21sex, e21dex/ int22ex
     &		/ ei22aex,  ei22cex, ei22gex, ei22uex 
     &		/tetraex/ eh4ex / triex /eh3ex
     &		/other/aex, bex, cex, slopeex, 
     &		tphbex, hguex, h3cex, haex, hbex, iasyex, ex3
          	     
	external baseid
c
        integer buflen, maxfnlen
        parameter(buflen=1000)
        parameter(maxfnlen=200)
        character*(buflen) bufstr

        if (pardirlen+maxfnlen .gt. buflen) then
          write(0, *)  " Error: length of parameter directory too long"
          stop
        end if
	
        bufstr = pardir(1:pardirlen) // '/stack.dat'
        open(unit=11, file=bufstr(1:len_trim(bufstr)),
     &    action='read',status='old')
        bufstr = pardir(1:pardirlen) // '/dangle.dat'
        open(unit=12, file=bufstr(1:len_trim(bufstr)),
     &    action='read',status='old')
        bufstr = pardir(1:pardirlen) // '/loop.dat'
        open(unit=13, file=bufstr(1:len_trim(bufstr)),
     &    action='read',status='old')
        bufstr = pardir(1:pardirlen) // '/tstackh.dat'
        open(unit=14, file=bufstr(1:len_trim(bufstr)),
     &    action='read',status='old')	
        bufstr = pardir(1:pardirlen) // '/tstacki.dat'
        open(unit=15, file=bufstr(1:len_trim(bufstr)),
     &    action='read',status='old')
        bufstr = pardir(1:pardirlen) // '/int11.dat'
        open(unit=16, file=bufstr(1:len_trim(bufstr)),
     &    action='read',status='old')	
        bufstr = pardir(1:pardirlen) // '/int21.dat'
        open(unit=17, file=bufstr(1:len_trim(bufstr)),
     &    action='read',status='old')	
        bufstr = pardir(1:pardirlen) // '/int22.dat'
        open(unit=18, file=bufstr(1:len_trim(bufstr)),
     &    action='read',status='old')
        bufstr = pardir(1:pardirlen) // '/tloop.dat'
        open(unit=19, file=bufstr(1:len_trim(bufstr)),
     &    action='read',status='old')
        bufstr = pardir(1:pardirlen) // '/triloop981.dat'
        open(unit=20, file=bufstr(1:len_trim(bufstr)),
     &    action='read',status='old')
	data
     &		slope, am, bm, cm,
     &		tphb, hgu, h3c, ha, hb,
     &		iasy /
     &		1.079d0, 3.4d0, 0.0d0, 0.4d0,
     &		0.5d0, -2.2d0, 1.4d0, 0.3d0, 1.6d0,
     &		0.5d0 /
     
	scale=1000.d0/(1.9872d0*310.15d0)
	slopeex=dexp(-slope*scale)
	aex=dexp(-am*scale)
	bex=dexp(-bm*scale)
	cex=dexp(-cm*scale)
	tphbex=dexp(-tphb*scale)
	hguex=dexp(-hgu*scale)	
	h3cex=dexp(-h3c*scale)	
	haex=dexp(-ha*scale)
	hbex=dexp(-hb*scale)
	iasyex=dexp(-iasy*scale)		
	ex3=dexp(-3.d0*scale)
			
	do i=1, 15
	  read (11, 11) b
	end do
	
11	format(A80)
12	format(A120)	
	do iblock=1, 4
	  do i=1, 11
	    read (11, 11) b
	  end do
	  do i=1, 4 
	    read (11, 11) b
	    do j=1, 16
	    
	      if (b((j-1)*5+1:j*5).eq. ' .   ' .or.
     &	          b((j-1)*5+1:j*5).eq. '  .  ' .or. 
     &	          b((j-1)*5+1:j*5).eq. '   . ' .or.      	    
     &	          b((j-1)*5+1:j*5).eq. '    .' ) then	   	    
	        b((j-1)*5+1:j*5)='999.0'
	      end if
	      read (b((j-1)*5+1:j*5), '(f5.1)') e((iblock-1)*4+i, j)
	    end do
	  end do
	end do
	
	
	
 	do i=1, 4
 	   do j=1, 4
 	      do m=1, 4
 	         do l=1, 4
 	         es(i, m, l, j)=e(4*(i-1)+m, 4*(j-1)+l)
 	         
 	         if(es(i, m, l, j).gt.998.d0) then
 	            esex(i, m, l, j)=0.d0
 	            else
 	            esex(i, m, l, j)=dexp(-es(i, m, l, j)*scale)
 	         end if
 	         
 	         end do
 	      end do
 	   end do
 	end do
	do iblock=1, 8
	  do i=1, 11
	    read (12, 11) b
	  end do
	    do j=1, 16
	      if (b((j-1)*5+1:j*5).eq. ' .   ') then
	        b((j-1)*5+1:j*5)='000.0'	        
	      end if
	      read (b((j-1)*5+1:j*5), '(f5.1)') e(iblock, j)
	      
	    end do
	end do
	
	
	do j=1, 4
	   do i=1, 4
	     do k=1, 4
	     	ed3(i,j,k)=e(j, (i-1)*4+k)
	     	ed5(i,j,k)=e(j+4, (i-1)*4+k)
	
 	        ed3ex(i,j,k)=dexp(-ed3(i,j,k)*scale)
 		ed5ex(i,j,k)=dexp(-ed5(i,j,k)*scale)
    	 	         	         	     	
	     end do
	   end do
	end do
	
	do i=1, 4
	  read (13, 11) b
	end do
	
	do i=1, 30
	  read (13, 11) b
	  
	  if (b(16:18).eq.' . ') then
	    b(14:18)='999.0'
	  end if
	  if (b(51:53).eq.' . ') then
	    b(49:53)='999.0'
	  end if
	  read (b(14:18), '(f5.1)') ei(i)
	  read (b(34:36), '(f5.1)') eb(i)	  	  	  
	  read (b(49:53), '(f5.1)') eh(i)
	  if(ei(i).gt.998.d0) then
 	       eiex(i)=0.d0
 	       else
 	       eiex(i)=dexp(-ei(i)*scale)
 	  end if
 	  	     
	  if(eb(i).gt.998.d0) then
 	       ebex(i)=0.d0
 	       else
 	       ebex(i)=dexp(-eb(i)*scale)
 	  end if
 	  
	  if(eh(i).gt.998.d0) then
 	       ehex(i)=0.d0
 	       else
 	       ehex(i)=dexp(-eh(i)*scale)
 	  end if
 	   	  
	end do
	do i=1, 15
	  read (14, 11) b
	end do
		
	do iblock=1, 4
	  do i=1, 11
	    read (14, 11) b
	  end do
	  do i=1, 4 
	    read (14, 11) b
	    do j=1, 16
	    
	      if (b((j-1)*5+1:j*5).eq. ' .   ' .or.
     &	          b((j-1)*5+1:j*5).eq. '  .  ' .or. 
     &	          b((j-1)*5+1:j*5).eq. '   . ' .or.      	    
     &	          b((j-1)*5+1:j*5).eq. '    .' ) then
     	    
	        b((j-1)*5+1:j*5)='999.0'
	      end if
	      read (b((j-1)*5+1:j*5), '(f5.1)') e((iblock-1)*4+i, j)
	    end do
	  end do
	end do 
	
 	do i=1, 4
 	   do j=1, 4
 	      do m=1, 4
 	         do l=1, 4
 	         ehs(i, m, l, j)=e(4*(i-1)+m, 4*(j-1)+l)
 	         
 	         if(ehs(i, m, l, j).gt.998.d0) then
 	            ehsex(i, m, l, j)=0.d0
 	            else
 	            ehsex(i, m, l, j)=dexp(-ehs(i, m, l, j)*scale)
 	         end if	
 	                  
 	         end do
 	      end do
 	   end do
 	end do
 	
	  
	do i=1, 15
	  read (15, 11) b
	end do
		
	do iblock=1, 4
	  do i=1, 11
	    read (15, 11) b
	  end do
	  do i=1, 4 
	    read (15, 11) b
	    do j=1, 16
	    
	      if (b((j-1)*5+1:j*5).eq. ' .   ' .or.
     &	          b((j-1)*5+1:j*5).eq. '  .  ' .or. 
     &	          b((j-1)*5+1:j*5).eq. '   . ' .or.      	    
     &	          b((j-1)*5+1:j*5).eq. '    .' ) then
     
	        b((j-1)*5+1:j*5)='999.0'
	      end if
	      read (b((j-1)*5+1:j*5), '(f5.1)') e((iblock-1)*4+i, j)
	    end do
	  end do
	end do 
	
 	do i=1, 4
 	   do j=1, 4
 	      do m=1, 4
 	         do l=1, 4
 	         eis(i, m, l, j)=e(4*(i-1)+m, 4*(j-1)+l)
	         
 	         if(eis(i, m, l, j).gt.998.d0) then
 	            eisex(i, m, l, j)=0.d0
 	            else
 	            eisex(i, m, l, j)=dexp(-eis(i, m, l, j)*scale)
 	         end if 
 	         	         
 	         end do
 	      end do
 	   end do
 	end do	    	  
	
	
	do i=1, 21
	  read (16, 12) c
	end do
	
	do iblock=1, 6
	  do i=1, 11
	    read (16, 12) c
	  end do
	  
	  do i=1, 4
	    read (16,'(24f5.1)') (e11((iblock-1)*4+i,j),  j=1, 24) 
	  end do 
	end do 
	
	do i1=1, 4
	 do i2=1, 4
	  do i3=1, 4
	   do i4=1, 4
	    do i5=1, 4
	     do i6=1, 4
	      ei11(i1,i2,i3,i4,i5,i6)=999.d0
	     end do
	    end do
	   end do
	  end do
	 end do
	end do
	do irb=1, 6
	 if(irb.eq.1) then
	   i1=1
	   i2=4
	 end if
	 if(irb.eq.2) then
	   i1=2
	   i2=3
	 end if	 
	 if(irb.eq.3) then
	   i1=3
	   i2=2
	 end if
	 if(irb.eq.4) then
	   i1=4
	   i2=1
	 end if
	 if(irb.eq.5) then
	   i1=3
	   i2=4
	 end if
	 if(irb.eq.6) then
	   i1=4
	   i2=3
	 end if	
	 
	 do i5=1, 4
	  do i6=1, 4
	   
	    do icb=1, 6
	     if (icb.eq.1) then
	       i3=1
	       i4=4
	     end if
	     if (icb.eq.2) then
	       i3=2
	       i4=3
	     end if	     
	     if (icb.eq.3) then
	       i3=3
	       i4=2
	     end if	     
	     if (icb.eq.4) then
	       i3=4
	       i4=1
	     end if
	     if (icb.eq.5) then
	       i3=3
	       i4=4
	     end if
	     if (icb.eq.6) then
	       i3=4
	       i4=3
	     end if
	     
		ei11(i1,i2,i3,i4,i5,i6)=
     &		e11((irb-1)*4+i5, (icb-1)*4+i6)
     	         
 	         if(ei11(i1,i2,i3,i4,i5,i6).gt.998.d0) then
 	            ei11ex(i1,i2,i3,i4,i5,i6)=0.d0
 	            else
 	   	    ei11ex(i1,i2,i3,i4,i5,i6)=
     &			dexp(-ei11(i1,i2,i3,i4,i5,i6)*scale)
 	         end if
 	         
	    end do
	  end do
	 end do
	end do       
 		    	
	do i=1, 20
	  read (17, 12) c
	end do
	
	do iblock=1, 24
	  do i=1, 10
	    read (17, 12) c
	  end do
	  
	  do i=1, 4
	    read (17,'(24f5.1)') (e21((iblock-1)*4+i,j),  j=1, 24) 
	  end do 
	end do 
	
	do i1=1, 4
	 do i2=1, 4
	  do i3=1, 4
	   do i4=1, 4
	    do i5=1, 4
	     do i6=1, 4
	      do i7=1, 4
	      e21s(i1,i2,i3,i4,i5,i6,i7)=999.d0
	      e21d(i1,i2,i3,i4,i5,i6,i7)=999.d0
	      end do
	     end do
	    end do
	   end do
	  end do
	 end do
	end do
	do ilb=1, 6
	 if(ilb.eq.1) then
	   i1=1
	   i2=4
	 end if
	 if(ilb.eq.2) then
	   i1=2
	   i2=3
	 end if	 
	 if(ilb.eq.3) then
	   i1=3
	   i2=2
	 end if
	 if(ilb.eq.4) then
	   i1=4
	   i2=1
	 end if
	 if(ilb.eq.5) then
	   i1=3
	   i2=4
	 end if
	 if(ilb.eq.6) then
	   i1=4
	   i2=3
	 end if	
	 
	 do i6=1, 4
	  do i5=1, 4
	   
	    do icb=1, 6
	     if (icb.eq.1) then
	       i3=1
	       i4=4
	     end if
	     if (icb.eq.2) then
	       i3=2
	       i4=3
	     end if	     
	     if (icb.eq.3) then
	       i3=3
	       i4=2
	     end if	     
	     if (icb.eq.4) then
	       i3=4
	       i4=1
	     end if
	     if (icb.eq.5) then
	       i3=3
	       i4=4
	     end if
	     if (icb.eq.6) then
	       i3=4
	       i4=3
	     end if
	     
	      do i7=1, 4
		e21s(i1,i2,i3,i4,i5,i6,i7)=
     &		e21((ilb-1)*16+(i6-1)*4+i5, (icb-1)*4+i7)
		e21d(i4,i3,i2,i1,i6,i7,i5)=e21s(i1,i2,i3,i4,i5,i6,i7)
 	         
 	         if(e21s(i1,i2,i3,i4,i5,i6,i7).gt.998.d0) then
 	            e21sex(i1,i2,i3,i4,i5,i6,i7) =0.d0
 	            else
 	            e21sex(i1,i2,i3,i4,i5,i6,i7)=
     &			dexp(-e21s(i1,i2,i3,i4,i5,i6,i7)*scale)
 	         end if
	         
 	         if(e21d(i4,i3,i2,i1,i6,i7,i5).gt.998.d0) then
 	            e21dex(i4,i3,i2,i1,i6,i7,i5) =0.d0
 	            else
 	            e21dex(i4,i3,i2,i1,i6,i7,i5)=
     &			dexp(-e21d(i4,i3,i2,i1,i6,i7,i5)*scale)
 	         end if
 	          	         			
	      end do
	    end do
	  end do
	 end do
	end do       
		           	     	     		    	   	   	  	   	 	 	
	do i=1, 31
	  read (18, 11) b
	end do
	
	do iblock=1, 36
	  do i=1, 10
	    read (18, 11) b
	  end do
	  
	  do i=1, 16
	    read (18,'(16f5.1)') (e22((iblock-1)*16+i,j),  j=1, 16) 
	  end do 
	end do 
	do i1=1, 4
	 do i2=1, 4
	  do i3=1, 4
	   do i4=1, 4
	    do i5=1, 4
	     do i6=1, 4
	      do i7=1, 4
	       
	       ei22a(i1,i2,i3,i4,i5,i6,i7)=999.d0
	       ei22c(i1,i2,i3,i4,i5,i6,i7)=999.d0
	       ei22g(i1,i2,i3,i4,i5,i6,i7)=999.d0	       
	       ei22u(i1,i2,i3,i4,i5,i6,i7)=999.d0
	       	       	       
	      end do
	     end do
	    end do
	   end do
	  end do
	 end do
	end do
	do ilb=1, 6
	 if(ilb.eq.1) then
	   i1=1
	   i2=4
	 end if
	 if(ilb.eq.2) then
	   i1=2
	   i2=3
	 end if	 
	 if(ilb.eq.3) then
	   i1=3
	   i2=2
	 end if
	 if(ilb.eq.4) then
	   i1=4
	   i2=1
	 end if
	 if(ilb.eq.5) then
	   i1=3
	   i2=4
	 end if
	 if(ilb.eq.6) then
	   i1=4
	   i2=3
	 end if	
	 do isb=1, 6
	  if(isb.eq.1) then
	    i3=1
	    i4=4
	  end if
	  if(isb.eq.2) then
	    i3=2
	    i4=3
	  end if	 
	  if(isb.eq.3) then
	    i3=3
	    i4=2
	  end if
	  if(isb.eq.4) then
	    i3=4
	    i4=1
	  end if
	  if(isb.eq.5) then
	    i3=3
	    i4=4
	  end if
	  if(isb.eq.6) then
	    i3=4
	    i4=3
	  end if
	  
	   do i5=1, 4
	    do i8=1, 4
	     do i6=1,4
	      do i7=1, 4
	        if (i8.eq.1) then
	        ei22a(i1,i2,i3,i4,i5,i6,i7)=
     &		e22((ilb-1)*96+(isb-1)*16+(i5-1)*4+i8, (i6-1)*4+i7)
     		end if
     		
	        if (i8.eq.2) then
	        ei22c(i1,i2,i3,i4,i5,i6,i7)=
     &		e22((ilb-1)*96+(isb-1)*16+(i5-1)*4+i8, (i6-1)*4+i7)
     		end if
	        if (i8.eq.3) then
	        ei22g(i1,i2,i3,i4,i5,i6,i7)=
     &		e22((ilb-1)*96+(isb-1)*16+(i5-1)*4+i8, (i6-1)*4+i7)
     		end if
     		
	        if (i8.eq.4) then
	        ei22u(i1,i2,i3,i4,i5,i6,i7)=
     &		e22((ilb-1)*96+(isb-1)*16+(i5-1)*4+i8, (i6-1)*4+i7)
     		end if     		     		
 	         
 	         if(ei22a(i1,i2,i3,i4,i5,i6,i7).gt.998.d0) then
 	            ei22aex(i1,i2,i3,i4,i5,i6,i7)=0.d0
 	            else
 	            ei22aex(i1,i2,i3,i4,i5,i6,i7)=
     &			dexp(-ei22a(i1,i2,i3,i4,i5,i6,i7)*scale)
 	         end if
  	         
 	         if(ei22c(i1,i2,i3,i4,i5,i6,i7).gt.998.d0) then
 	            ei22cex(i1,i2,i3,i4,i5,i6,i7)=0.d0
 	            else
 	            ei22cex(i1,i2,i3,i4,i5,i6,i7)=
     &			dexp(-ei22c(i1,i2,i3,i4,i5,i6,i7)*scale)
 	         end if
 	         
 	         if(ei22g(i1,i2,i3,i4,i5,i6,i7).gt.998.d0) then
 	            ei22gex(i1,i2,i3,i4,i5,i6,i7)=0.d0
 	            else
 	            ei22gex(i1,i2,i3,i4,i5,i6,i7)=
     &			dexp(-ei22g(i1,i2,i3,i4,i5,i6,i7)*scale)
 	         end if
  	         
 	         if(ei22u(i1,i2,i3,i4,i5,i6,i7).gt.998.d0) then
 	            ei22uex(i1,i2,i3,i4,i5,i6,i7)=0.d0
 	            else
 	            ei22uex(i1,i2,i3,i4,i5,i6,i7)=
     &			dexp(-ei22u(i1,i2,i3,i4,i5,i6,i7)*scale)
 	         end if
 	              		              	 	              		              		     		
      	      end do
      	     end do
      	    end do
      	   end do
      	 end do
      	end do
  	
	      		
	do i=1, 2
	  read (19, 191) a1
	end do
191	format(A20)
		
	do i1=1, 4
	 do i2=1, 4
	  do i3=1, 4
	   do i4=1, 4
	    do i5=1, 4	
	     do i6=1, 4
	    	eh4(i1,i2,i3,i4,i5,i6)=eh(4)  	         
 	        eh4ex(i1,i2,i3,i4,i5,i6)=
     &			dexp(-eh4(i1,i2,i3,i4,i5,i6)*scale)	        	         	    	
	     end do
	    end do
	   end do
	  end do
	 end do
	end do
	
	do i=1, 30
	
	  read (19, 192) a, bonus
	  i1=baseid(a(1:1))
	  i2=baseid(a(2:2))	
	  i3=baseid(a(3:3))
	  i4=baseid(a(4:4))	
	  i5=baseid(a(5:5))
	  i6=baseid(a(6:6))
	
	  eh4(i1,i2,i3,i4,i5,i6)=eh4(i1,i2,i3,i4,i5,i6)+bonus 	         	         
 	  eh4ex(i1,i2,i3,i4,i5,i6)=
     &			dexp(-eh4(i1,i2,i3,i4,i5,i6)*scale) 	  
 	         
	end do
192	format(A6,1x, f4.1)
	do i=1, 2
	  read (20, 191) a1	
	end do
	
	do i1=1, 4
	 do i2=1, 4
	  do i3=1, 4
	   do i4=1, 4
	    do i5=1, 4	
	    	eh3(i1,i2,i3,i4,i5)=eh(3) 
 	        eh3ex(i1,i2,i3,i4,i5)=dexp(-eh3(i1,i2,i3,i4,i5)*scale) 	         	    	
	    end do
	   end do
	  end do
	 end do
	end do	
	do i=1, 20
	
	  read (20, 202) a2, bonus
			
	  i1=baseid(a2(1:1))
	  i2=baseid(a2(2:2))	
	  i3=baseid(a2(3:3))
	  i4=baseid(a2(4:4))	
	  i5=baseid(a2(5:5))
	  eh3(i1,i2,i3,i4,i5)=eh3(i1,i2,i3,i4,i5)+bonus
	         
 	  eh3ex(i1,i2,i3,i4,i5)=dexp(-eh3(i1,i2,i3,i4,i5)*scale)
 	         
	end do
		
202	format(1x, A5, 3x, f4.1)
	
	close(11)
	close(12)
	close(13)		
	close(14)
	close(15)
	close(16)
	close(17)
	close(18)		
	close(19)
	close(20)
	
	return
	end
	
	function baseid(base)
	character*1, base
	integer baseid
	
	if (base. eq. 'A'.or.base . eq. 'a') then
	  baseid=1
	end if
	
	if (base. eq. 'C'.or.base . eq. 'c') then 
	  baseid=2
	end if	
	
	if (base. eq. 'G'.or.base . eq. 'g') then 
	  baseid=3
	end if	
	if (base. eq. 'U'.or. base. eq. 'u') then 
	  baseid=4
	end if	
	
	if (base. eq. 'T'.or. base. eq. 't') then 
	  baseid=4
	end if
	
	return
	end 	
