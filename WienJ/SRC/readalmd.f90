subroutine readalmd(itape,filename,k1,k2)
use cl
use bl1
use bl4
use bl5
use bl9
implicit none

integer :: itape
character*22 :: filename
integer :: k1,k2
real*8 :: f(nmax)
real*8 :: e(nmax)
!complex*16 :: almblm(5,(lmax+1)*(lmax+1),nmax*ndif,kmax,jspin)
complex*16 :: almblmtmp(5,(lmax+1)*(lmax+1),nmax)
real*8 :: eferm
!integer :: nlo(0:lmax,nat)
!integer :: kn(kmax)
real*8 :: weight
integer :: i,j,k,l,m,n
integer :: spin,na,mu,index,nm,nk,na2
character*80 :: kong

spin=itape-10
nk=2*jspin
nm=0
!f=0.0d0
!e=0.0d0
!almblm=(0.0d0,0.0d0)
open(unit=itape,file=filename)
read(itape,*)
read(itape,*)
read(itape,*)
read(itape,*)eferm
do na=1,nat
do l=0,lmax
read(itape,*)
read(itape,*)
do j=1,nlo(l,na,spin)
read(itape,*)
end do
end do
do k=1,kmax
read(itape,*)
read(itape,*)
read(itape,*)
!write(*,*)kn(k,na)
do n=1,kn2(k,spin)
read(itape,*)weight,e(n)
!f(n)=weight/kweight(k)*dfloat(kmax2)
if(e(n).gt.eferm)then
f(n)=0.0d0  
else
f(n)=1.0d0
end if
end do
do mu=1,mult(na)  !
!nm=mu
!do n=1,na-1
!nm=nm+mult(n)
!end do
!write(*,*)nm
read(itape,*)kong
do n=1,kn2(k,spin)
index=0
do l=0,lmax
do m=-l,l
index=index+1
read(itape,*)almblmtmp(1,index,n),almblmtmp(2,index,n)
do j=1,nlo(l,na,spin)
read(itape,*)almblmtmp(j+2,index,n)
end do
end do
end do
end do
do na2=1,jspin
if(na.eq.atom(na2).and.mu.eq.amu(na2))then
if(k.eq.k1)then
fermi2(:,1,spin)=f
energy2(:,1,spin)=e
almblm2(:,:,:,1,na2,spin)=almblmtmp
nm=nm+1
end if
if(k.eq.k2)then
fermi2(:,2,spin)=f
energy2(:,2,spin)=e
almblm2(:,:,:,2,na2,spin)=almblmtmp
nm=nm+1
end if
end if
end do
if(nm.eq.nk)goto 5
end do
end do
!if(nm.eq.jspin)exit
end do
5 continue
close(itape)
return
end
