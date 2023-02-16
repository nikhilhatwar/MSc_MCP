


			!Program to find average quantities from ising model output
			!Enter the value of n = total no. of data points in the file 
			!is, js starting values after equilibriation
			! k is dummy variable	
			!kb = boltzman const

			PROGRAM avenger
			IMPLICIT NONE
			REAL,ALLOCATABLE::E(:),M(:)
			INTEGER::i,j,n,k,g,is,js
			REAL::mE,mM,msE,msM,X,C,t,kb
	
			n=100000				!Enter the value of n
			is=50000			!Enter i starting value
			js=50000			!Enter j starting value
			kb=1				!
			mE=0				!Mean Energy
			mM=0				!mean magnetization
			msE=0				!mean square energy
			msM=0				!mean square magnetization


			OPEN(UNIT=10,FILE='energy.dat')
			OPEN(UNIT=20,FILE='mag.dat')
			ALLOCATE(E(n))
			ALLOCATE(M(n))
			
			DO i=1,n
				READ(10,*),k,E(i),t
				READ(20,*),k,M(i)
			ENDDO
		
			DO i=is,n			
				mE= mE + E(i)
				msE= msE+ E(i)**2
			ENDDO

			mE=mE/REAL(n-is)	!mean energy per spin
			msE=msE/REAL(n-is)	!mean sq energy per spin

			C=(msE-mE**2)/(kb*t*t)		! Specific heat
			
			DO j=js,n			
				mM= mM + M(j)
				msM= msM + m(j)**2
			ENDDO

			mM=mM/REAL(n-js)	!mean magnetization per spin
			msM=msM/REAL(n-js)	!mean sq magnetization per spin

			X=(msM-mM**2)/(kb*t*t)		! suceptibility
			
			PRINT*,t,mE,mM,C,X
		
			END PROGRAM














