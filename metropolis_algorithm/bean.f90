


			!Date: 27/08/2015
			! Binning random numbers for plotting
			! All RNs should be positive
			! Enter bin size (s) and number of RNs to be binned(n)
			
			PROGRAM shaktiman
			IMPLICIT NONE
			REAL::x(5000000),bin(5000000)=0,s=2
			INTEGER::i,p,n=10000	
			OPEN(UNIT=10,FILE='metroout.dat')
			OPEN(UNIT=20,FILE='beanout.dat')			
			READ(10,*),(x(i),i=1,n)			
			DO i=1,n		
				p=x(i)/s   
				bin(p)=bin(p)+1 				
			ENDDO			  
			DO i=1,n
				WRITE(20,*),i,bin(i)
			ENDDO
			END PROGRAM		
				

			
