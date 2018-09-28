program b3_fft
implicit none
real*8::T,Tr,rmax,b3,dr,dk,pi,fr
real*8:: fk(1024,1024)
integer::i,j,n
T=100
open(unit=6,file='b3.dat')
do while (T<=500)
Tr=T/117.7
pi=4*atan(1.0)
rmax=50.0
n=2**10
dr=rmax/float(n)
dk=pi/(dr*n)
do j=1,n
    do i=1,n
       fk(i,j)=((2/pi)**0.5)*(1./(j*dk))*(dr*i)*fr(i*dr,Tr)*sin(i*dr*j*dk)
     enddo
     
enddo
b3=0
do j=1,n 
   do i=1,n
        b3=b3-((2*pi)**1.5)*(1./3.)*((j*dk)**2)*(fk(i,j))**3
    enddo    
enddo
write(6,100)T,Tr,10*b3/float(n)
100 format(f9.2,f15.4,f15.4)
!print*,T,Tr,b3/float(n) 
T=T+5
enddo 
close(6) 
end program b3_fft

function fr(r,Tr)
implicit none
real*8::r,Tr,sigma,epsil,epsilon4,s2,s6,u,x,fr
sigma=1.0
epsil=1.0
epsilon4=epsil*4
s2=(sigma*sigma)/(r*r)
s6=s2*s2*s2
u=epsilon4*s6*(s6-1.0)
x=-u/Tr
if (abs(x)<0.01) then
    fr=x+x*x/2+(1./6.)*x**3+(1./24.)*(x**4)+(1./120.)*x**5
else
    fr=exp(x)-1.0
end if
end function fr        

          
