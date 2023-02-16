
			! 6/10/2015
			! 2D Ising Model
			! p1,p2 = temporary varibles
			! n = number of spins in the linear array
			! omcc = 1 monte carlo cycles
			! tmcc = total number of monte carlo cycles
			! TE = Total Energy
			! TM = Total Magnetization
!	__________________________________________________________________________

			MODULE scooby
			INTEGER::n,omcc,tmcc
			INTEGER,ALLOCATABLE,DIMENSION(:,:)::s
			REAL::t,Je
			END MODULE
!	__________________________________________________________________________

			PROGRAM ultimatrix
			USE scooby
			IMPLICIT NONE
			INTEGER::l,i,j,p1,p2,dels
			REAL::r1,r2,w,oldE,newE,delE,a,TE,TM
			
			OPEN(UNIT=10,FILE='sing.dat')
			OPEN(UNIT=20,FILE='energy.dat')
			OPEN(UNIT=30,FILE='mag.dat')
			READ(10,*),n,Je,t,omcc,tmcc
			ALLOCATE(s(n,n))	

			DO i=1,n															
				DO j=1,n
					s(i,j)=1
				ENDDO
			ENDDO
			CALL energy(TE)
			CALL magnet(TM)
			print*,TE,TM
			DO l=1,tmcc
				DO i=1,omcc
					CALL RANDOM_NUMBER(r1)
					CALL RANDOM_NUMBER(r2)
					p1=r1*n
					p2=r2*n
					IF(p1==0 .OR. p2==0)THEN
						p1=1
						p2=1
					ENDIF
					oldE=0
					oldE=s(p1,p2)*s(p1,p2+1)+oldE
					oldE=s(p1,p2)*s(p1,p2-1)+oldE
					oldE=s(p1,p2)*s(p1+1,p2)+oldE
					oldE=s(p1,p2)*s(p1-1,p2)+oldE
					oldE=-J*oldE

					s(p1,p2)=-s(p1,p2)
					
					newE=0
					newE=s(p1,p2)*s(p1,p2+1)+newE
					newE=s(p1,p2)*s(p1,p2-1)+newE
					newE=s(p1,p2)*s(p1+1,p2)+newE
					newE=s(p1,p2)*s(p1-1,p2)+newE
					newE=-J*newE
				
					delE=newE-oldE
					
					IF(delE<=0)THEN 
						GOTO 10
					ELSE
						CALL RANDOM_NUMBER(a)
						dels=delE/Je
						w=exp(-dels/t)	
						IF(a<=w)THEN 
							GOTO 10
						ELSE
						s(p1,p2)=-s(p1,p2)
						ENDIF
					ENDIF		
		10		p1=p1
				ENDDO
				CALL energy(TE)
				CALL magnet(TM)
				
				WRITE(20,*),l,TE,t
				WRITE(30,*),l,TM
			ENDDO
			print*,TE,TM
			END PROGRAM
!	__________________________________________________________________________

			SUBROUTINE energy(TE)
			USE scooby
			REAL, INTENT(INOUT)::TE
			INTEGER::i,j,temp
			TE=0
			temp=0
			DO i=2,n-1
				DO j=2,n-1
					temp=temp+s(i,j)*s(i,j+1)
					temp=temp+s(i,j)*s(i,j-1)
					temp=temp+s(i,j)*s(i+1,j)
					temp=temp+s(i,j)*s(i-1,j)
					TE=-temp*Je
				ENDDO
			ENDDO
			TE=TE/REAL(n*n)
			END SUBROUTINE energy
!	__________________________________________________________________________

			SUBROUTINE magnet(TM)
			USE scooby
			REAL, INTENT(INOUT)::TM
			INTEGER::i,j,temp
			TM=0
			TM=SUM(s)
			TM=TM/REAL(n*n)
			END SUBROUTINE magnet				
!	__________________________________________________________________________			
			
			
			
