
			! PRN :- 2013P027	
			! 6/10/2015
			! Objective:- To plot temp vs magnetization, temp vs suceptibility 
			! and temp vs specific heat for 2D Ising Model

			! p1,p2 = temporary integer varibles 
			! L =length of the 2D array of spins s(n:n)
			! omcc = 1 monte carlo cycles
			! tmcc = total number of monte carlo cycles
			! Je = strength of spin-spin interaction
			! t = temperature in reduced units
			! TE,TM = Total Energy and Total Magnetization
			! E1,E2.. energy of present spin due to nearest neighbours
			! r1,r2.. Random Numbers

			! Enter the file name for UNIT 40 is 'size of lattice'
			! Enter the value of 'is' and 'js' in avengers subroutine
	!_____________________________________________________________________________

			MODULE goku
			INTEGER,ALLOCATABLE::s(:,:)
			REAL,ALLOCATABLE::E(:),M(:)
			INTEGER::L,omcc,tmcc,p1,p2,n
			REAL::Je,t,mE,mM,C,X
			END MODULE
	!_____________________________________________________________________________

			PROGRAM dising
			USE goku
			IMPLICIT NONE
			INTEGER::i,j,h
			REAL::r1,r2,r3,w,delE,dels,TE,TM,E1,E2,tem1,tem2
			
			OPEN(UNIT=10,FILE='sing.dat')
			OPEN(UNIT=20,FILE='energy.dat')
			OPEN(UNIT=30,FILE='mag.dat')
			OPEN(UNIT=40,FILE='80.dat')
			READ(10,*),L,Je,t,omcc,tmcc
			ALLOCATE(s(0:L+1,0:L+1),E(tmcc),M(tmcc))	
			
			!Initialization
			DO i=1,L
				DO j=1,L
				s(i,j)=1
				CALL pbc(i,j)
				ENDDO
			ENDDO

			DO h=1,200

				!Total energy calculation
				tem1=0
				DO i=1,L
					DO j=1,L
			tem1=tem1+s(i,j)*s(i-1,j)+s(i,j)*s(i,j-1)+s(i,j)*s(i+1,j)+s(i,j)*s(i,j+1)
					ENDDO
				ENDDO
				TE=-Je*tem1/2
				
			!Total magnetization calculation
			!	TM=SUM(s(1:L,1:L))
			!	PRINT*,TE,ABS(TM)
			!	WRITE(20,*),0,TE
			!	WRITE(30,*),0,TM
	
				DO n=1,tmcc
					tem2=0
					DO i=1,omcc
			
						! call a random spin
						CALL RANDOM_NUMBER(r1)
						CALL RANDOM_NUMBER(r2)
						p1=1+r1*L
						p2=1+r2*L
						
						! flip the spin and find delE
						CALL energy(E1)					
						s(p1,p2)=-s(p1,p2)
						CALL pbc(p1,p2)
						CALL energy(E2)
					
						delE=E2-E1
						
						!Metropolis
						IF(delE<=0)GOTO 10
						IF(delE>0)THEN
							CALL RANDOM_NUMBER(r3)
							dels=delE/Je
							w=EXP(-dels/t)
		
							IF(r3<w)THEN 
								GOTO 10
							ELSE
								s(p1,p2)=-s(p1,p2)
								CALL pbc(p1,p2)
								GOTO 20
							ENDIF
						ENDIF		
		
			10			tem2=tem2+delE
			20			tem2=tem2
					ENDDO
	
				TE=TE+tem2
				TM=SUM(s(1:L,1:L))
	
			!	WRITE(20,*),n,TE
			!	WRITE(30,*),n,ABS(TM)
				E(n)=TE
				M(n)=ABS(TM)
				ENDDO
			
			CALL avengers()

			PRINT*,t,mE,ABS(mM),C,X

			WRITE(40,*),t,mE,ABS(mM),C,X

			t=t+0.2	

			IF(t>6)EXIT	
			
			ENDDO

			END PROGRAM
	!_________________________________________________________________________

		!to calculate energy due to only neighbouring interaction
			SUBROUTINE energy(En)
			USE goku
			REAL, INTENT(INOUT)::En

			En=     s(p1,p2)*s(p1+1,p2) + s(p1,p2)*s(p1,p2+1)
			En= En + s(p1,p2)*s(p1-1,p2) + s(p1,p2)*s(p1,p2-1)
			En=-Je*En

			END SUBROUTINE energy
	!_________________________________________________________________________________
		
		!subroutine to apply periodic boundary condition
			SUBROUTINE pbc(p,q)
			USE goku
			INTEGER,INTENT(IN)::p,q
		
			IF(p==1) s(L+1,q)=s(p,q)
			IF(q==1) s(p,L+1)=s(p,q)
			IF(p==L) s(0,q)=s(p,q)
			IF(q==L) s(p,0)=s(p,q)

			END SUBROUTINE pbc	
	!__________________________________________________________________________________

			! is,js = data point number after which equilibriation starts
			! k is dummy variable	
			! kb = boltzman const

			SUBROUTINE avengers()
			USE goku
			IMPLICIT NONE
		
			INTEGER::i,j,k,g,is,js
			REAL::msE,msM,kb
			
			is=n/2				!Enter i starting value
			js=n/2				!Enter j starting value
			kb=1				!Reduced unit
			mE=0				!Mean Energy
			mM=0				!mean magnetization
			msE=0				!mean square energy
			msM=0				!mean square magnetization

			DO i=is,n			
				mE= mE + E(i)
				msE= msE + E(i)**2
			ENDDO

			mE=mE/REAL(n-is)		!mean energy
			msE=msE/REAL(n-is)		!mean sq energy 

			C=(msE-mE**2)/(kb*t*t)		! Specific heat

			C=C/(L*L)			!specific heat per spin

			mE= mE/(L*L)			!mean energy per spin
			
			DO j=js,n			
				mM= mM + M(j)
				msM= msM + M(j)**2
			ENDDO

			mM=mM/REAL(n-js)		!mean magnetization 
			msM=msM/REAL(n-js)		!mean sq magnetization

			X=(msM-mM**2)/(kb*t*t)		!suceptibility
			
			X=X/(L*L)			!suceptibility per spin 

			mM = mM/(L*L)			!mean magnetization per spin

			END SUBROUTINE

	   !___________________________________________________________________________________________			
			


			
			
			
			
