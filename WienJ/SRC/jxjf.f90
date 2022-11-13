subroutine jxjf
use cl
use bl1
use bl5
use bl6
use bl7
use bl8
use bl9
implicit none
real*8,parameter :: cin=1.0d0/137.0359895d0**2
integer :: spin1,spin2,na
integer :: l1,l2,i1,i2,lm3
integer :: i,j,k,l,m,n
logical :: alive
integer :: istat
real*8 :: uu

spin1=1
spin2=2
uuu=0.0d0
do na=1,jspin
do l1=0,lmax
do i1=1,2+nlo2(l1,na,spin1)
do l2=0,lmax
do i2=1,2+nlo2(l2,na,spin2)
do lm3=1,lmmax2(na)
uu=(radwf2(1,i1,l1,1,na,spin1)*radwf2(1,i2,l2,1,na,spin2)+cin*radwf2(2,i1,l1,1,na,spin1)*radwf2(2,i2,l2,1,na,spin2))*rr2(1,na)*(dx(atom(na))+1.0d0)/3.0d0*vlm2(1,lm3,na)
j=2
do
uu=uu+(radwf2(1,i1,l1,j,na,spin1)*radwf2(1,i2,l2,j,na,spin2)+cin*radwf2(2,i1,l1,j,na,spin1)*radwf2(2,i2,l2,j,na,spin2))*dx(atom(na))*rr2(j,na)*4.0d0/3.0d0*vlm2(j,lm3,na)
j=j+1
if(j.eq.781)exit
uu=uu+(radwf2(1,i1,l1,j,na,spin1)*radwf2(1,i2,l2,j,na,spin2)+cin*radwf2(2,i1,l1,j,na,spin1)*radwf2(2,i2,l2,j,na,spin2))*dx(atom(na))*rr2(j,na)*2.0d0/3.0d0*vlm2(j,lm3,na)
j=j+1
end do
uu=uu+(radwf2(1,i1,l1,781,na,spin1)*radwf2(1,i2,l2,781,na,spin2)+cin*radwf2(2,i1,l1,781,na,spin1)*radwf2(2,i2,l2,781,na,spin2))*dx(atom(na))*rr2(781,na)/3.0d0*vlm2(781,lm3,na)

uuu(lm3,i2,l2,i1,l1,na)=uu
!if(na.eq.2.and.abs(uu).gt.1.0d-8)write(6,*)na,l1,i1,l2,i2,lm3
end do
end do
end do
end do
end do
end do

inquire(file=foldername(1:flen)//'dat.vorbup',exist=alive)
nlorb=0
if(alive)then

open(unit=7,file=foldername(1:flen)//'dat.vorbup')
open(unit=8,file=foldername(1:flen)//'dat.vorbdn')

read(7,*,iostat=istat)
if(istat/=0)then
close(7,status='delete')
close(8,status='delete')
else
close(7)
open(unit=7,file=foldername(1:flen)//'dat.vorbup')

read(7,*)nmod,nsp,natorb
read(8,*)nmod,nsp,natorb
do k=1,natorb
read(7,*)iatom,m
read(8,*)iatom,m
do n=1,m
read(7,*)l
read(8,*)l
do i=-l,l
do j=-l,l
read(7,*)rval,cval
vorbup=cmplx(rval,cval)
read(8,*)rval,cval
vorbdn=cmplx(rval,cval)
do na=1,jspin
if(iatom.eq.atom(na))then
nlorb(na)=m
lorb(n,na)=l
vorb(j,i,n,na)=(vorbup-vorbdn)/2.0d0
end if
end do
end do
end do
end do
end do
close(7)
close(8)

orbuu=0.0d0
do na=1,jspin
do n=1,nlorb(na)
l=lorb(n,na)
do i1=1,2+nlo2(l,na,spin1)
do i2=1,2+nlo2(l,na,spin2)
uu=(radwf2(1,i1,l,1,na,spin1)*radwf2(1,i2,l,1,na,spin2)+cin*radwf2(2,i1,l,1,na,spin1)*radwf2(2,i2,l,1,na,spin2))*rr2(1,na)*(dx(atom(na))+1.0d0)/3.0d0
j=2
do
uu=uu+(radwf2(1,i1,l,j,na,spin1)*radwf2(1,i2,l,j,na,spin2)+cin*radwf2(2,i1,l,j,na,spin1)*radwf2(2,i2,l,j,na,spin2))*dx(atom(na))*rr2(j,na)*4.0d0/3.0d0
j=j+1
if(j.eq.781)exit
uu=uu+(radwf2(1,i1,l,j,na,spin1)*radwf2(1,i2,l,j,na,spin2)+cin*radwf2(2,i1,l,j,na,spin1)*radwf2(2,i2,l,j,na,spin2))*dx(atom(na))*rr2(j,na)*2.0d0/3.0d0
j=j+1
end do
uu=uu+(radwf2(1,i1,l,781,na,spin1)*radwf2(1,i2,l,781,na,spin2)+cin*radwf2(2,i1,l,781,na,spin1)*radwf2(2,i2,l,781,na,spin2))*dx(atom(na))*rr2(781,na)/3.0d0

orbuu(i2,i1,n,na)=uu
end do
end do
end do
end do
write(6,*)'orb!'
end if

else
write(6,*)'no orb'
end if

return
end
