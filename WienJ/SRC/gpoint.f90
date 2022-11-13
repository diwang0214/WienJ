subroutine gpoint(spt,weight)
implicit none
integer :: n
real*8 :: spt(3,378)
real*8 :: weight(378)
real*8 :: pi,flake
real*8 :: z,rxy,phi
real*8 :: xx(14)
real*8 :: sum
integer :: i,j,index

pi=4.0d0*atan(1.0d0)

n=0
call grule(14,xx,weight)
do i=8,14
j=i-7
weight(i)=weight(j)
xx(i)=-xx(j)
end do
flake=exp(-4.0d0)
do j=1,14
z=xx(j)
rxy=sqrt(1.0d0-z*z)
do i=0,26
n=n+1
phi=pi*(dble(2*i)/dble(27)+flake)
index=14*i+j
weight(index)=weight(j)
spt(1,index)=cos(phi)*rxy
spt(2,index)=sin(phi)*rxy
spt(3,index)=z
end do
end do
!if(n!=378)write(*,*)n
sum=0.0d0
do j=1,n
sum=sum+weight(j)
end do
sum=4.0d0*pi/sum
weight=weight*sum
return
end

subroutine grule(n,x,w)
implicit none
integer :: n
real*8 :: x(n),w(n)
real*8 :: pi
integer :: m
integer :: i,it,k
real*8 :: e1,t,x0,pkm1,pk,t1,pkp1,den,d1,dpn,d2pn,d3pn,d4pn,u,v,h,p,dp,fx
pi=4.0d0*atan(1.0d0)
m=(n+1)/2
e1=n*(n+1)
do i=1,m
t=(4*i-1)*pi/(4*n+2)
x0=(1.0d0-(1.0d0-1.0d0/n)/(8.0d0*n*n))*cos(t)
do it=1,3
pkm1=1.0d0
pk=x0
do k=2,n
t1=x0*pk
pkp1=t1-pkm1-(t1-pkm1)/k+t1
pkm1=pk
pk=pkp1
end do
den=1.0d0-x0*x0
d1=n*(pkm1-x0*pk)
dpn=d1/den
d2pn=(2.0d0*x0*dpn-e1*pk)/den
d3pn=(4.0d0*x0*d2pn+(2.0d0-e1)*dpn)/den
d4pn=(6.0d0*x0*d3pn+(6.0d0-e1)*d2pn)/den
u=pk/dpn
v=d2pn/dpn
h=-u*(1.0d0+0.5d0*u*(v+u*(v*v-u*d3pn/(3.0d0*dpn))))
p=pk+h*(dpn+0.5d0*h*(d2pn+h/3.0d0*(d3pn+0.25d0*h*d4pn)))
dp=dpn+h*(d2pn+0.5d0*h*(d3pn+h*d4pn/3.0d0))
h=h-p/dp
x0=x0+h
end do
x(i)=x0
fx=d1-h*e1*(pk+0.5d0*h*(dpn+h/3.0d0*(d2pn+0.25d0*h*(d3pn+0.2d0*h*d4pn))))
w(i)=2.0d0*(1.0d0-x(i)*x(i))/(fx*fx)
end do
if(m+m.gt.n)x(m)=0.0d0
return
end
