
			!date: 17-09-2015
			! Importance sampling integration only with
			! input of RNs of required distribution
			! put appropriate normalization costant if required

			PROGRAM kamehameha
			IMPLICIT NONE
			REAL::r(500000),x(500000),fn,f,p,nr,dr,ratio=0,y
			INTEGER::i,n
			OPEN(UNIT=10,FILE='impin.dat')
			OPEN(UNIT=20,FILE='impout.dat')
			n=10000
			READ(10,*),(r(i),i=1,n)
			DO i=1,n
				y=x(i)
				nr=f(y)
				dr=p(y)
				ratio=(nr/dr)+ratio
			ENDDO		
			fn=(1/n)*ratio
			PRINT*,'Value of Integration is',fn							
			END PROGRAM

			REAL FUNCTION f(y)
				f=exp(-y*y)
			END FUNCTION f

			REAL FUNCTION p(y)
				p=exp(-y)
			END FUNCTION p

		
