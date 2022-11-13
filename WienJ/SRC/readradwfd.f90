subroutine readradwfd(itape,filename,radwf,nlo,r)
use cl
use bl1
use bl9
implicit none
integer :: itape
character*22 :: filename
real*8 :: radwf(2,5,0:lmax,781,jspin)
real*8 :: radwftmp(2,5,0:lmax,781)
integer :: nlo(0:lmax,nat)
real*8 :: r(781,jspin)
real*8 :: rtmp(781)
integer :: i,j,k,l,na,kk,na2,nm
character*80 :: kong

nm=0
radwf=0.0d0
open(unit=itape,file=filename)
do na=1,nat
read(itape,*)kong
do l=0,lmax
do i=1,nlo(l,na)
read(itape,*)kong
do j=1,781
read(itape,*)rtmp(j),radwftmp(1,i+2,l,j),radwftmp(2,i+2,l,j)
end do
end do
read(itape,*)kk
!write(6,*)kk
do j=1,781
read(itape,*)rtmp(j),radwftmp(1,1,l,j),radwftmp(1,2,l,j)
end do
read(itape,*)kong
do j=1,781
read(itape,*)rtmp(j),radwftmp(2,1,l,j),radwftmp(2,2,l,j)
end do
end do
do na2=1,jspin
if(na.eq.atom(na2))then
r(:,na2)=rtmp
radwf(:,:,:,:,na2)=radwftmp
nm=nm+1
end if
end do
if(nm.eq.jspin)exit
end do
close(itape)
return
end
