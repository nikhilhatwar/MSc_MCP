


			!PRN: 2013P027
			!DATE: 24/11/2015
			! Objective : To find first excited state energy for 1D finite 
			! potential well.		
		
		
			! V = Potential:- V(x)=-V0 for -a<x<a ;V(x) = 0 elsewhere
			! 13/11/2015
			! V0 = depth of finite potential well
			! a = length of potential (-a to a)
			! h = increment in 'x' (mesh parameter)
			! xm = initial value of 'x' to start with
			! N = State of system (gs=0, 1st excited=1..)
			! xc = matching point
			! espt = trial eigen value of given state to start with
			! cl,cr = the 2nd value of y from left and right side
			! y() = is the wavefunction array
			! x() = x variable array
			! V() = potential array
			! er() = array to store error
			! esp(:) array to store energy values
			! ei = increment in trail energy value in iteration
			! en = trail energy increment counter variable
			! e1,e2,e3 = temporary variables for energy
			! er1,er2,er3 = temporary variable for error
	
			! we are using reduced units as (mass)/(planckconst)**2 =1
		
		!_____________________________________________________________________________________
	
			MODULE omnitrix
			REAL,ALLOCATABLE::y(:),x(:),V(:),esp(:)
			INTEGER::n,Ns,xcn,en,i,j
			REAL::V0,a,h,xm,del,xc,cl,cr,yf,yb,er(100000)
			REAL::e1,e2,e3,er1,er2,er3
	
			REAL::gn,gnp1,gnp2,gnm1,gnm2
			REAL*8::fn,fnp1,fnp2,fnm1,fnm2
			
			END MODULE
		!_____________________________________________________________________________________
	
			PROGRAM Shenron
			USE omnitrix
			IMPLICIT NONE
			
			REAL::ei,espt
			INTEGER::k
			
			OPEN(UNIT=10,FILE='in.dat')
			OPEN(UNIT=20,FILE='v.dat')
			OPEN(UNIT=30,FILE='y.dat')
			OPEN(UNIT=40,FILE='e.dat')
			OPEN(UNIT=50,FILE='g.dat')
			OPEN(UNIT=60,FILE='ye.dat')
	
			READ(10,*),V0,a,h,xm,del,Ns,xc,espt,cl,cr,ei
	
			n=xm/h		
			ALLOCATE(x(-n:n),y(-n:n),V(-n:n),esp(n))
				
			DO i=-n,n
				y(i)=0.0
				x(i)=i*h-xm+n*h
				V(i)=0
				IF(x(i)>=-a .AND. x(i)<=a)V(i)=V0
				WRITE(20,*),x(i),V(i)
			ENDDO
	
			xcn=(xc)/h
			y(-n)=0.0	! for forward relation
			y(-n+1)=cl
			y(n)=0.0	! for backward relation
			y(n-1)=cr
			e1=V0
			e2=0
			en=ABS(V0)/ei
			
			DO j=1,en	
		
				CALL hulk(espt)	
				DO i=-n,n
					WRITE(30,*),x(i),y(i)
				ENDDO
				er(j)=yf-yb
				esp(j)=espt
				WRITE(40,*),esp(j),er(j)	
		10		espt=espt+ei
				
			ENDDO
	
			!Bisection method:-
			
			DO i=1,en
				IF((er(i)>0 .AND. er(i+1)<0) .OR. (er(i)<0 .AND. er(i+1)>0))THEN
				e1=esp(i)
				e2=esp(i+1)		
				EXIT
				ENDIF
			ENDDO
	
			! eigen value lies between e1 and e2
	
			DO j=1,10000
				CALL hulk(e1)
				er1=yf-yb
				CALL hulk(e2)
				er2=yf-yb
				e3=(e1+e2)/2
				CALL hulk(e3)
				er3=yf-yb
				IF(er3>0)e1=e3
				IF(er3<0)e2=e3			
				IF(((ABS(e1)-ABS(e2))/ABS(e1))<del)EXIT
			ENDDO	
		
			CALL hulk(e3)
			DO i=-n,n
				WRITE(60,*),x(i),y(i)
			ENDDO
	
			PRINT*,'Eigenvalue for energy state =',Ns,'is',e3
			END PROGRAM
	    	!____________________________________________________________________________________
	
			SUBROUTINE hulk(espt)
			USE omnitrix
	
				!Forward relation
			
				DO i=-n,xcn,+1
				gn=  2*(espt-V(i))			!gn = g(n)
				gnp1=2*(espt-V(i+1))			!gnp1=f(n+1)
				gnp2=2*(espt-V(i+2))			!gnp2=f(n+2)
				fn=  1+((gn*(h*h))/12.0)		!fn=f(n)
				fnp1=1+((gnp1*(h*h))/12.0)		!fnp1 = f(n+1) 
				fnp2=1+((gnp2*(h*h))/12.0)		!fnp2 = f(n+2)
				
				y(i+2)= ((12-10*fnp1)*y(i+1) - fn*y(i))/fnp2
	
				IF(i==xcn)yf=(y(i+2)-y(i+1))/(h*y(i+2))
	
				ENDDO
	
				!Backward relation		
	
				DO i=n,xcn,-1
				gn = 2*(espt-V(i))
				gnm1=2*(espt-V(i-1))			!gnm1 = g(n-1)
				gnm2=2*(espt-V(i-2))
				fn=  1+((gn*(h*h))/12.0)		
				fnm1=1+((gnm1*(h*h))/12.0)		!fnm1 = f(n-1)
				fnm2=1+((gnm2*(h*h))/12.0)
				y(i-2) = ((12-10*fnm1)*y(i-1) - fn*y(i))/fnm2
				IF(i==xcn)yb=(y(i-1)-y(i-2))/(h*y(i-1))
				ENDDO
			
			END SUBROUTINE
	
		!________________________________________________________________________________________
			
	
		





	
		
