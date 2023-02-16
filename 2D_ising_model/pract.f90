


			MODULE groot
			INTEGER,ALLOCATABLE::s(:,:)
			INTEGER::n,omcc,tmcc
			REAL::t,B,Je
			END MODULE

	!____________________________________________________________________

			PROGRAM ising
			USE groot
			IMPLICIT NONE
			INTEGER::l,i,j,p1,p2,dels
			REAL::r1,r2,w,oldE,newE,delE,a,TE,TM
			
			
			n=3
			ALLOCATE(s(n,n))	

			DO i=1,3,0.5
			

			print*,i

			ENDDO

			END PROGRAM
