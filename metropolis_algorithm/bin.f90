

			!Date: 27/08/2015
			! Binning
			! Bin size=s
			! Smallest number in the input RN data= start
			! Biggest number in the input RN data= enda
			! Total bins=a=(enda-start)/s
			! Total number of RNs in the input file to be binned= n
		 
			PROGRAM shaktiman
			IMPLICIT NONE
			REAL::x(5000000),bin(500000,3)=0,s,start,enda
			INTEGER::i,a,m,n,temp,p,j

			s=0.1
			start=-10                      	!INPUTS
			enda=10
			n=10000

			a=(enda-start)/s
			OPEN(UNIT=10,FILE='metroout.dat')
			OPEN(UNIT=20,FILE='binout.dat')
			READ(10,*),(x(i),i=1,n)			
			DO i=1,a
				bin(i,1)=start
				start=start+s
				bin(i,3)=NINT(start*(10))
			ENDDO			
			DO i=1,n		
				p=x(i)/s   
				DO j=1,a
					IF(bin(j,3)==p)THEN
					temp=j
					ENDIF
				ENDDO
				bin(temp,2)=bin(temp,2)+1 				
			ENDDO		  
			DO i=1,a
				WRITE(20,*),bin(i,1),bin(i,2),bin(i,3)
			ENDDO
			END PROGRAM		
				

			
