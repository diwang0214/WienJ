subroutine readalm(itape,filename,f,e,almblm)
use cl
use bl1
use bl4
use bl9
implicit none

integer :: itape
character*22 :: filename
real*8 :: f(nmax,kmax)
real*8 :: e(nmax,kmax)
complex*16 :: almblm(5,(lmax+1)*(lmax+1),nmax,kmax,jspin)
complex*16 :: almblmtmp(5,(lmax+1)*(lmax+1),nmax,kmax)
real*8 :: eferm
integer :: nlo(0:lmax,nat)
integer :: kn(kmax)
real*8 :: weight
integer :: i,j,k,l,m,n
integer :: na,mu,index,nm,na2
character*80 :: kong

nm=0
f=0.0d0
e=0.0d0
almblm=(0.0d0,0.0d0)
open(unit=itape,file=filename)
read(itape,*)kong
read(itape,*)kong
read(itape,*)kong
read(itape,*)eferm
do na=1,nat
do l=0,lmax
read(itape,*)kong
read(itape,*)nlo(l,na)
do j=1,nlo(l,na)
read(itape,*)kong
end do
end do
do k=1,kmax
read(itape,*)kong
read(itape,*)kong
read(itape,*)l,m,kn(k)
!write(*,*)kn(k,na)
do n=1,kn(k)
read(itape,*)weight,e(n,k)
!f(n,k)=weight/kweight(k)*dfloat(kmax2)
if(e(n,k).gt.eferm)then
f(n,k)=0.0d0  
else
f(n,k)=1.0d0
end if
end do
do mu=1,mult(na)  !
!nm=mu
!do n=1,na-1
!nm=nm+mult(n)
!end do
!write(*,*)nm
read(itape,*)kong
do n=1,kn(k)
index=0
do l=0,lmax
do m=-l,l
index=index+1
read(itape,*)almblmtmp(1,index,n,k),almblmtmp(2,index,n,k)
do j=1,nlo(l,na)
read(itape,*)almblmtmp(j+2,index,n,k)
end do
end do
end do
end do
do na2=1,jspin
if(na.eq.atom(na2).and.mu.eq.amu(na2))then
almblm(:,:,:,k,na2)=almblmtmp(:,:,:,k)
if(k.eq.1)then
write(6,*)na
nm=nm+1
end if
end if
end do
end do
end do
if(nm.eq.jspin)exit
end do
close(itape)
return
end
