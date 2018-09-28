program B3_MC
implicit none
real*8::x2,x3,y2,y3,z2,z3,r12,r13,r23,f,s,error
real*8::T,Tr,constante1,constante2,rc,pi,B3,suma
integer::i,N
print*,'Ingrese el numero de corridas, N:'
read*,N
print*,'Ingrese la temperatura del gas:'
read*,T
print*,'Ingrese la constante del gas:'
read*,constante1
open(unit=6,file='virial3.dat')
pi=4*atan(1.0)
rc=1.92

constante2=-(1.0/3.0)*(16.0*pi*pi)/float(N)
do while (T<=450)
suma=0.0
s=0.0
 do i=1,N
   Tr=T/constante1
   x2=rc*rand()
   y2=rc*rand()
   z2=rc*rand()
   x3=rc*rand()
   y3=rc*rand()
   z3=rc*rand()
   r12=sqrt(x2*x2+y2*y2+z2*z2)
   r13=sqrt(x3*x3+y3*y3+z3*z3)
   r23=sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
   suma=suma+rc*rc*f(r12,Tr)*f(r13,Tr)*f(r23,Tr)*r12**2*r13**2
   s=s+suma**2   
 enddo
   B3=constante2*suma
   error=-constante2*sqrt((s/N-(suma/N)**2)/N)
   
   
   write(6,100)T,Tr,B3,error
   100 format(f9.4,f15.4,f15.4,f15.4)
   !print*,Tr,B3
   T=T+2.0
enddo   
close(6)   
end program B3_MC     

function f(r,Temp)
implicit none
real*8::r,f,Temp,s2,s6,x,u
s2=1.0/(r*r)
s6=s2*s2*s2
u=4*s6*(s6-1.0)
x=-u/Temp
if (abs(x)<0.01) then
    f=x+(x*x/2.0)+(x**3/6.0)+(x**4/24.0)+(x**5/120.0)
else
   f=exp(x)-1.0
endif
end function f       
