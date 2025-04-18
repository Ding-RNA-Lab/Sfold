C!!!	do not input non-canonical BPS, efn of Mfold treat them as unpaired
C		i.e., Zuker' Madison format
	subroutine  energy(n, numseq, h5, h3, hnum,
     &		eall, estack, eext, ehp, ebg, eint, emloop)

	integer n, numseq(n), h5(n), h3(n), hnum(n),
     &		intbase(n), index, ind5, ind3, index2,
     & 		hp5, hp3, bext5, bext3, bint5, bint3, 
     &		iext5,iext3, iint5, iint3, nh, m5(n), m3(n),
     &		i, j
     
	double precision a, b, c, slope,
     &		tphb, hgu, h3c, ha, hb, iasy, 
     &		es(4, 4, 4, 4), ed3(4,4,4), ed5(4,4,4),
     &		ehs(4, 4, 4, 4), eis(4, 4, 4, 4),
     &		e21s(4,4,4,4,4,4,4),  e21d(4,4,4,4,4,4,4),
     &		ei22a(4,4,4,4,4,4,4), ei22c(4,4,4,4,4,4,4),
     &		ei22g(4,4,4,4,4,4,4), ei22u(4,4,4,4,4,4,4),      
     &		ei11(4,4,4,4,4,4),
     &		eh4(4,4,4,4,4,4), 
     &		eh(30), eb(30), ei(30),
     &		eall, estack, eext, ehp, ebg, eint, emloop
            
	double precision hpe, be, inte, mloope
	
	common /stack/ es /dangle/ ed3, ed5 /HBI/ eh, eb, ei
     &		/Hstack/ ehs /Istack/ eis /int11/ ei11
     &		/int21/ e21s, e21d/ int22
     &		/ ei22a,  ei22c, ei22g, ei22u
     &		/tetra/ eh4
	data
     &		slope, a, b, c,
     &		tphb, hgu, h3c, ha, hb,
     &		iasy /
     &		1.079d0, 3.4d0, 0.0d0, 0.4d0,
     &		0.5d0, -2.2d0, 1.4d0, 0.3d0, 1.6d0,
     &		0.5d0 /
     
ccc					1  total stacking E
	estack=0.d0
	
	do i=1, n
	   if (h5(i).ne.0. and. hnum(i). gt. 1) then
	   
	        do j=1, hnum(i)-1
	          estack=estack+
     &		es( numseq(h5(i)+j-1), numseq(h5(i)+j), 
     &		numseq(h3(i)-j), numseq(h3(i)-j+1) )
     	 	end do
    	   end if
    	end do
   
ccc			external: part 1, 2  dangling energies at 5'&3' ends
	eext=0.d0
	if (h5(1).gt.1) then 	
	   eext=eext+ed5(numseq(h5(1)),numseq(h3(1)),numseq(h5(1)-1))
        end if
     
	index=1
	do i=2, n
	   if (h3(i).eq.0 ) then
	      go to 10
	   else
	      if (h3(i).gt. h3(index)) then
	          index=i
	      end if
	    end if
	end do
	   
10 	continue
	if (h3(index).lt.n) then 	
	   eext=eext+
     &	   ed3(numseq(h5(index)),numseq(h3(index)),numseq(h3(index)+1))
        end if	
 
        
		
	do i=1, n
	   intbase(i)=1
	end do
	
	do i=1, n
	   if (h5(i). ne.0) then
	      do j=1, hnum(i)
	      	 intbase(h5(i)+j-1)=0
	      	 intbase(h3(i)-j+1)=0
	      end do
	      
	      do j=h5(i)+hnum(i), h3(i)-hnum(i)
	         intbase(j)=0
	      end do
	   end if
	end do
   
C				add dangling energies
	do i=h3(1)+1, h5(index)-1
	   if (intbase(i-1).eq.0. and. intbase(i).eq.1. 
     &			and.intbase(i+1).eq.0) then
	     do j=1, n
	        if(h3(j).eq. (i-1)) then
	          ind5=j
	        end if
	        
	        if(h5(j).eq. (i+1)) then
	          ind3=j
	        end if
	     end do
	     	      
	     eext=eext+min(ed3(numseq(h5(ind5)),numseq(h3(ind5)),numseq(i)), 
     &		ed5(numseq(h5(ind3)), numseq(h3(ind3)), numseq(i)) )
     	     
	   end if
	   
	   if (intbase(i-1).eq.0. and. intbase(i).eq.1. 
     &			and.intbase(i+1).eq.1) then
	     do j=1, n
	        if(h3(j).eq. (i-1)) then
	          ind5=j
	        end if
	     end do
     		
	     eext=eext+ed3(numseq(h5(ind5)),numseq(h3(ind5)),numseq(i))
	   end if
 	   if (intbase(i-1).eq.1. and. intbase(i).eq.1. 
     &			and.intbase(i+1).eq.0) then    
	     do j=1, n
	        if(h5(j).eq. (i+1)) then
	          ind3=j
	        end if
	     end do
	     eext=eext+ed5(numseq(h5(ind3)), numseq(h3(ind3)), numseq(i)) 
	        
	   end if
	   
	end do
	
	index=0
	
	do i=1, n-1
	   if (h5(i).ne.0. and. h5(i+1).eq.0 ) then
	      index=i
	      go to 25
	   end if
	end do
25	continue	
	   if ( ((numseq(h5(1)).eq.1).and.(numseq(h3(1)).eq.4)).or.
     &	     ((numseq(h5(1)).eq.4).and.(numseq(h3(1)).eq.1)).or.
     &	     ((numseq(h5(1)).eq.3).and.(numseq(h3(1)).eq.4)).or.
     &	     ((numseq(h5(1)).eq.4).and.(numseq(h3(1)).eq.3)) ) then
	        eext=eext+tphb
	   end if 	   
	do i=2, index
	    do j=1, i-1
	      if ( h5(j).lt. h5(i). and. 
     &		   h3(i). lt. h3(j)) then
     		 go to 27
     	      end if
     	    end do	   	   
	
	   if ( ((numseq(h5(i)).eq.1).and.(numseq(h3(i)).eq.4)).or.
     &	     ((numseq(h5(i)).eq.4).and.(numseq(h3(i)).eq.1)).or.
     &	     ((numseq(h5(i)).eq.3).and.(numseq(h3(i)).eq.4)).or.
     &	     ((numseq(h5(i)).eq.4).and.(numseq(h3(i)).eq.3)) ) then
	        eext=eext+tphb
	   end if 	   
27	continue	   
	end do
ccc				compute for H, B, I, M loop cases
	ehp=0.d0
	ebg=0.d0
	eint=0.d0
	emloop=0.d0
	
	do i= 1, index
	   
	   do j=i+1, index
	      if ((h5(i)+hnum(i)-1).lt. h5(j). and. 
     &		   h3(j). lt. (h3(i)-hnum(i)+1)	) then
     		 index2=j
     		 go to 30
     	      end if	   	   
	   end do
	   hp5=h5(i)+hnum(i)-1
	   hp3=h3(i)-hnum(i)+1
	   ehp=ehp+
     &	  	hpe(hp5, hp3, n, numseq,eh, eh4, ehs, 
     &				slope, hgu, h3c, ha, hb, tphb)     
	   go to 100
		
30	   nh=2
	   m5(1)=h5(i)+hnum(i)-1
	   m3(1)=h3(i)-hnum(i)+1
	   m5(2)=h5(index2)
	   m3(2)=h3(index2)	
	   do j=index2+1, index
	      if (h3(index2).lt. h5(j). and. 
     &		   h3(j). lt. (h3(i)-hnum(i)+1) ) then
		 nh=nh+1
		 m5(nh)=h5(j)    
		 m3(nh)=h3(j)
		 index2=j
	      end if
	   end do
	   if (nh.eq.2) then
	       if (m5(2).eq. (m5(1)+1). or. m3(2).eq.(m3(1)-1) ) then
		   bext5=m5(1)
     		   bext3=m3(1)
     		   bint5=m5(2)
     	           bint3=m3(2)
		   ebg=ebg+
     &		be(bext5, bext3, bint5, bint3, 
     &			n, numseq, eb, es, slope, tphb)
               else 
     	 	  iext5=m5(1)
     		  iext3=m3(1)
     		  iint5=m5(2)
     		  iint3=m3(2)
		  eint=eint+
     &		inte(iext5,iext3, iint5, iint3, n, numseq, slope, 
     &	 iasy, ei11,ei22a, ei22c,ei22g,ei22u, e21s, e21d, eis, ei)
     	       end if
     	   end if
     	     
	  if (nh.gt.2) then
	   emloop=emloop+mloope(n, numseq,nh, m5, m3, a, b,c, tphb, ed3, ed5)
	  end if
100	continue
	  
	end do
	eall=estack+eext+ehp+ebg+eint+emloop
	
	
	return
	end  
	
ccc		function to evaluate hairpin enery with closing hp5-hp3
	function hpe(hp5, hp3, n, numseq,eh, eh4, ehs, 
     &	slope, hgu, h3c, ha, hb, tphb)
	
	integer n, numseq(n), hp5, hp3, k
	double precision hpe, eh(30), eh4(4,4,4,4,4,4), ehs(4,4,4,4),
     &			slope, hgu, h3c, ha, hb, tphb
	
	if ((hp3-hp5-1). gt. 30) then
	    hpe=eh(30)+slope*dlog(dble(hp3-hp5-1)/30.d0)
	else
	    hpe=eh(hp3-hp5-1)
	end if
	if ((hp3-hp5-1).eq.3) then
	
	   if ( ((numseq(hp5).eq.1).and.(numseq(hp3).eq.4)).or.
     &		((numseq(hp5).eq.4).and.(numseq(hp3).eq.1)).or.
     &		((numseq(hp5).eq.3).and.(numseq(hp3).eq.4)).or.
     &	       ((numseq(hp5).eq.4).and.(numseq(hp3).eq.3)) ) then
     		  hpe=hpe+tphb   
           end if 	
	   
	end if
	if ((hp3-hp5-1).eq.4) then
	   hpe=eh4( numseq(hp5), numseq(hp5+1),numseq(hp5+2),
     &		    numseq(hp5+3),numseq(hp5+4), numseq(hp5+5) )
        end if	   
	      if ((hp3-hp5-1).ge.4) then
	         hpe=hpe+
     &	         ehs(numseq(hp5), numseq(hp5+1), 
     &	         numseq(hp3-1), numseq(hp3))
     	      end if
     	      
	
C			DHM et al say 2 GG in BP, this is not so
C			by mfold, 
     	      if( (hp5-2).ge.1.and.(hp3+2).le.n.and.
     &		numseq(hp5).eq. 3. and. numseq(hp3).eq.4. and.     	      
     &		numseq(hp5-1).eq. 3. and. numseq(hp5-2).eq.3 
     &			) then
     		hpe=hpe+hgu
     	      end if 
     	      
  	          	            	      
	      if ((hp3-hp5-1).eq.3.and. numseq(hp5+1).eq.2. 
     &			and. numseq(hp5+2).eq.2.	      
     &			and. numseq(hp5+3).eq.2) then
     		hpe=hpe+h3c
     		go to 5
     	      end if  
	      if ((hp3-hp5-1).gt.3) then
	        do k=hp5+1, hp3-1
	          if(numseq(k).ne.2) then
	            go to 5
	          end if
	        end do
	        hpe=hpe+ha*dble(hp3-hp5-1)+hb
	      end if
	      
5	      continue
	return
	end 
	
ccc				function to evaluate bulge energ with  
ccc		 		ext BP	bext5-bext3, int BP: bint5, bint3
	function be(bext5, bext3, bint5, bint3, 
     &			n, numseq, eb, es, slope, tphb)
    	integer n, numseq(n), bext5, bext3, bint5, bint3 
    	double precision be, eb(30), es(4,4,4,4), slope, tphb
    	
 
    	if ((bint5-bext5-1+bext3-bint3-1). gt. 30) then 
    	   be=eb(30)+slope*dlog(dble(bint5-bext5-1+bext3-bint3-1)/30.d0)
    	else
    	   be=eb(bint5-bext5-1+bext3-bint3-1)
    	end if
	
	if ((bint5-bext5-1+bext3-bint3-1).eq.1) then 
	   be=be+es(numseq(bext5), numseq(bint5), 
     &				numseq(bint3), numseq(bext3))	 
        end if
   
        
ccc		for size >1,	add AU, UA, GU , UG penalty to both ends  
	if   ((bint5-bext5-1+bext3-bint3-1). gt. 1) then
	   if ( ((numseq(bext5).eq.1).and.(numseq(bext3).eq.4)).or.
     &		((numseq(bext5).eq.4).and.(numseq(bext3).eq.1)).or.
     &		((numseq(bext5).eq.3).and.(numseq(bext3).eq.4)).or.
     &	((numseq(bext5).eq.4).and.(numseq(bext3).eq.3)) ) then
     		     be=be+tphb
           end if 
            
 	   if ( ((numseq(bint5).eq.1).and.(numseq(bint3).eq.4)).or.
     &		((numseq(bint5).eq.4).and.(numseq(bint3).eq.1)).or.
     &		((numseq(bint5).eq.3).and.(numseq(bint3).eq.4)).or.
     &	((numseq(bint5).eq.4).and.(numseq(bint3).eq.3)) ) then
     		     be=be+tphb
     		     
           end if
        end if
     
	return
	end 
ccc				function to evaluate interio loop energy with  
ccc		 		ext BP	iext5-iext3, int BP: iint5, iint3  
	function inte(iext5,iext3, iint5, iint3, n, numseq, slope, 
     &	 iasy, ei11,ei22a, ei22c,ei22g,ei22u, e21s, e21d, eis, ei)	
     	integer iext5,iext3, iint5, iint3, n, numseq(n)
     	double precision inte, slope,iasy, 
     &		e21s(4,4,4,4,4,4,4),  e21d(4,4,4,4,4,4,4),
     &		ei22a(4,4,4,4,4,4,4), ei22c(4,4,4,4,4,4,4),
     &		ei22g(4,4,4,4,4,4,4), ei22u(4,4,4,4,4,4,4),      
     &		 ei11(4,4,4,4,4,4),
     &		ei(30), eis(4, 4, 4, 4)
		    if(iint5.eq.(iext5+2). and. iint3.eq.(iext3-2)) then
		       inte=
     &			ei11(numseq(iext5), numseq(iext3),
     &			    numseq(iext5+2), numseq(iext3-2),	         
     &			    numseq(iext5+1), numseq(iext3-1) )
     			go to 6		      		      
		    end if 
		    				    	
		    if (iint5.eq.(iext5+2).and. iint3.eq.(iext3-3) ) then
		      	inte=
     &			e21s(numseq(iext5), numseq(iext3),
     &			numseq(iext5+2), numseq(iext3-3),	         
     &			numseq(iext5+1), 	
     &			numseq(iext3-2), numseq(iext3-1))
     			go to 6
     		    end if  
     		    
		    if (iint5.eq.(iext5+3).and.iint3.eq.(iext3-2) ) then
		      	inte=
     &			e21d( numseq(iext5), numseq(iext3),
     &			numseq(iext5+3), numseq(iext3-2),	         
     &			numseq(iext5+1), numseq(iext5+2),	
     &			numseq(iext3-1) )
     			go to 6
     		    end if 
		    if (iint5.eq.(iext5+3).and.iint3.eq.(iext3-3) ) then
		    
		    	if (numseq(iext3-1) .eq. 1) then
		      	inte=
     &			ei22a(numseq(iext5), numseq(iext3),
     &			numseq(iext5+3), numseq(iext3-3),	         
     &			numseq(iext5+1), numseq(iext5+2),	
     &			numseq(iext3-2) )
     			end if
     			
		    	if (numseq(iext3-1) .eq. 2) then
		      	inte=
     &			ei22c(numseq(iext5), numseq(iext3),
     &			numseq(iext5+3), numseq(iext3-3),	         
     &			numseq(iext5+1), numseq(iext5+2),	
     &			numseq(iext3-2) )
     			end if
     			
		    	if (numseq(iext3-1) .eq. 3) then
		      	inte=
     &			ei22g(numseq(iext5), numseq(iext3),
     &			numseq(iext5+3), numseq(iext3-3),	         
     &			numseq(iext5+1), numseq(iext5+2),	
     &			numseq(iext3-2) )
     			end if 
     			 
		    	if (numseq(iext3-1) .eq. 4) then
		      	inte=
     &			ei22u(numseq(iext5), numseq(iext3),
     &			numseq(iext5+3), numseq(iext3-3),	         
     &			numseq(iext5+1), numseq(iext5+2),	
     &			numseq(iext3-2) )
     			end if     			
     			   			     			
     			go to 6
     		    end if     		        		       	
     		    
     		    if ( ((iint5-iext5-1+iext3-iint3-1).gt.4 ).or. 
     &     	        (iint5.eq.(iext5+2). and.iint3. eq. (iext3-4)) .or.
     &		     (iint5 .eq.(iext5+4). and.iint3. eq. (iext3-2)) ) then
		       if ((iint5-iext5-1+iext3-iint3-1).gt.30 )  then
		           inte=ei(30)+
     &	            slope*dlog(dble(iint5-iext5-1+iext3-iint3-1)/30.d0)
     		       else
		            inte=ei(iint5-iext5-1+iext3-iint3-1)
		       end if 
			if (iint5.eq.(iext5+2).or.iint3 . eq. (iext3-2)) then
			   inte=inte
     &		          +eis(numseq(iext5), 1, 1, numseq(iext3))
     &		+min(3.d0, iasy*abs(dble((iint5-iext5-1)-(iext3-iint3-1))) )
     &			  +eis(numseq(iint3), 1, 1, numseq(iint5))
      		        else
     		      inte=inte
     &		          +eis(numseq(iext5), numseq(iext5+1), 
     &				numseq(iext3-1), numseq(iext3))
     &		+min(3.d0,iasy*abs(dble((iint5-iext5-1)-(iext3-iint3-1))) )
     &			  +eis(numseq(iint3), numseq(iint3+1), 
     &			      numseq(iint5-1), numseq(iint5))
                       end if
     	            end if        
     	            
6	continue
    	return
     	end 
     	    	
ccc					M-loop, nh: # of helices, 
	function mloope(n, numseq,  nh, m5, m3, a, b, c, tphb, ed3, ed5)
	integer n, numseq(n),nh, m5(n), m3(n), i
	double precision mloope, a, b,c, tphb, ed3(4,4,4), ed5(4,4,4)
	
	mloope=a+nh*c
	do i=1, nh
 	      if ( ((numseq(m5(i)).eq.1).and.(numseq(m3(i)).eq.4)).or.
     &		((numseq(m5(i)).eq.4).and.(numseq(m3(i)).eq.1)).or.
     &		((numseq(m5(i)).eq.3).and.(numseq(m3(i)).eq.4)).or.
     &		((numseq(m5(i)).eq.4).and.(numseq(m3(i)).eq.3)) ) then
     		  mloope=mloope+tphb 
              end if	
     
              	
	   if (i.eq. 1)  then
	      mloope=mloope+(m5(2)-m5(1)-1)*b
              if (m5(2).eq.(m5(1)+2)) then 
                 mloope=mloope+
     &		  min( ed3(numseq(m3(1)), numseq(m5(1)), numseq(m5(1)+1)), 
     &		       ed5(numseq(m5(2)), numseq(m3(2)), numseq(m5(2)-1)) ) 
              end if
              
              if (m5(2).gt.(m5(1)+2)) then 
                 mloope=mloope+
     &			ed3(numseq(m3(1)), numseq(m5(1)), numseq(m5(1)+1))+ 
     &			ed5(numseq(m5(2)), numseq(m3(2)), numseq(m5(2)-1)) 
              end if 
           end if
           
	   if (i.gt. 1. and. i. lt. nh)  then
	      mloope=mloope+(m5(i)-m3(i-1)-1)*b
              if (m3(i).eq.(m5(i+1)-2)) then 
                 mloope=mloope+
     &	      min( ed3(numseq(m5(i)), numseq(m3(i)), numseq(m3(i)+1)), 
     &		   ed5(numseq(m5(i+1)), numseq(m3(i+1)), numseq(m5(i+1)-1)) ) 
              end if
              
              if (m3(i).lt.(m5(i+1)-2)) then 
                 mloope=mloope+
     &		   ed3(numseq(m5(i)),numseq(m3(i)), numseq(m3(i)+1))+ 
     &		   ed5(numseq(m5(i+1)), numseq(m3(i+1)), numseq(m5(i+1)-1))
              end if 
           end if           
	   if ( i. eq. nh)  then
	      mloope=mloope+(m3(1)-m3(nh)-1)*b
	                    
              if (m3(nh).eq.(m3(1)-2)) then 
                 mloope=mloope+
     &		min( ed3(numseq(m5(nh)), numseq(m3(nh)), numseq(m3(nh)+1)), 
     &		     ed5(numseq(m3(1)), numseq(m5(1)), numseq(m3(1)-1)) ) 
              end if
              
              if (m3(nh).lt.(m3(1)-2)) then 
                 mloope=mloope+
     &		   ed3(numseq(m5(nh)), numseq(m3(nh)), numseq(m3(nh)+1))+ 
     &		   ed5(numseq(m3(1)), numseq(m5(1)), numseq(m3(1)-1))
              end if 
           end if 
           
   	end do
   	       
        return
        end                                     
