subroutine caljqso(r1)
use cl
use bl1
use bl2
use bl4
use bl5
use bl6
use bl7
use bl8
use bl9
use wmpi
implicit none
real*8 :: r1(3)
real*8 :: sk3(3)
integer :: k1,k2,n1,n2,spin1,spin2,na,q
integer :: l1,m1,l2,m2,lm1,lm2,lm3,i1,i2
integer :: i,j,k,l,n
integer*8 :: nd,nx,m
real*8,parameter :: cin=1.0d0/137.0359895d0**2
real*8 :: uu
complex*16 :: wf1,wf2
real*8 :: j1tmp,j1,j3tmp
real*8 :: sktmp1(3),sktmp2(3),sktmp3(3)
complex*16 :: j2tmp,j2(2,2),j3
complex*16 :: sumtmp(nxy),sumj
character*22 :: filename(2)
character*4  :: myStringq           

filename(1)=foldername(1:flen)//'dat.almblmup'
filename(2)=foldername(1:flen)//'dat.almblmdn'


spin1=1
spin2=2

allocate(fermi2(nmax,kmax2,2),energy2(nmax,kmax2,2))
allocate(almblm2(5,(lmax+1)*(lmax+1),nmax,kmax2,jspin,2))

call readalm(11,filename(1),fermi2(:,:,1),energy2(:,:,1),almblm2(:,:,:,:,:,1))
call readalm(12,filename(2),fermi2(:,:,2),energy2(:,:,2),almblm2(:,:,:,:,:,2))
write(6,*)'almblm read'


m=0
nd=0
nx=0
sumj=0.0d0 ! kn1=kn2,energy1=energy2,fermi1=fermi2
!open(unit=22,file='jq.txt')

do q=qmin+cpuid, qmax, num_cpu
write(myStringq,'(i4.4)') q
open(unit=22+cpuid,file='./data/jq_'//myStringq)
!write(22+cpuid,*)'now open it, q and cpuid =',q,cpuid

write(6,*)'sztest:q,cpuid',q,cpuid,stdout
!if(cpuid.eq.0) write(stdout,*)'sztest:q,cpuid',q,cpuid

if(q.eq.1)then
write(22+cpuid,*)kmax2,so,nxy
write(22+cpuid,*)r1
end if
sumtmp=0.0d0
do k1=1,kmax2
sktmp1=sk2(:,k1)+sk2(:,q)
do k2=1,kmax2
do i=1,3
sktmp2(i)=sktmp1(i)-sk2(i,k2)
end do
call dmatrixprod(br_dir,sktmp2,sktmp3,3)
if(abs(exp(2.0d0*pi*xs*sktmp3(1))-1.0d0).lt.1.d-8)then ! 2pi
if(abs(exp(2.0d0*pi*xs*sktmp3(2))-1.0d0).lt.1.d-8)then
if(abs(exp(2.0d0*pi*xs*sktmp3(3))-1.0d0).lt.1.d-8)then

do n1=1,kn2(k1,1)
do n2=1,kn2(k2,2)
j1tmp=fermi2(n1,k1,1)-fermi2(n2,k2,2)
if(abs(j1tmp).lt.1.d-8)cycle
m=m+1
j1=j1tmp/(energy2(n1,k1,1)-energy2(n2,k2,2)) !*kweight(k1)/nsym
if(abs(j1).gt.1.d3)then
nd=nd+1
!write(6,*)k1,k2,n1,n2
!write(6,*)j1
!write(6,*)fermi2(n1,k1,1),fermi2(n2,k2,2),energy2(n1,k1,1),energy2(n2,k2,2)
!stop
!cycle
else
nx=nx+1
end if

j2=0.0d0
na=1
do spin1=1,2
spin2=3-spin1
lm1=0
do l1=0,lmax
do m1=-l1,l1
lm1=lm1+1
lm2=0
do l2=0,lmax
do m2=-l2,l2
lm2=lm2+1

do lm3=1,lmmax2(na)
if(abs(real(yyy(lm3,lm2,lm1,na))).lt.1.0d-8.and.abs(aimag(yyy(lm3,lm2,lm1,na))).lt.1.0d-8)cycle
j2tmp=0.0d0
do i1=1,2+nlo2(l1,na,spin1)
do i2=1,2+nlo2(l2,na,spin2)
if(spin1.eq.1)then
j2tmp=j2tmp+conjg(almblm2(i1,lm1,n1,k1,na,spin1))*almblm2(i2,lm2,n2,k2,na,spin2)*uuu(lm3,i2,l2,i1,l1,na)
else
j2tmp=j2tmp+conjg(almblm2(i1,lm1,n1,k1,na,spin1))*almblm2(i2,lm2,n2,k2,na,spin2)*uuu(lm3,i1,l1,i2,l2,na)
end if
end do
end do
j2(spin1,1)=j2(spin1,1)+j2tmp*yyy(lm3,lm2,lm1,na) !*(-1.d0)**spin2
end do
end do
end do
end do
end do
do n=1,nlorb(na)
l=lorb(n,na)
lm1=l*l
do m1=-l,l
lm1=lm1+1
lm2=l*l
do m2=-l,l
lm2=lm2+1
j2tmp=0.0d0
do i1=1,2+nlo2(l,na,spin1)
do i2=1,2+nlo2(l,na,spin2)
if(spin1.eq.1)then
j2tmp=j2tmp+conjg(almblm2(i1,lm1,n1,k1,na,spin1))*almblm2(i2,lm2,n2,k2,na,spin2)*orbuu(i2,i1,n,na)
else
j2tmp=j2tmp+conjg(almblm2(i1,lm1,n1,k1,na,spin1))*almblm2(i2,lm2,n2,k2,na,spin2)*orbuu(i1,i2,n,na)
end if
end do
end do
j2(spin1,1)=j2(spin1,1)+j2tmp*vorb(m2,m1,n,na) !*(-1.0d0)**spin2
end do
end do
end do
end do

na=jspin
do spin1=1,2
spin2=3-spin1
lm1=0
do l1=0,lmax
do m1=-l1,l1
lm1=lm1+1
lm2=0
do l2=0,lmax
do m2=-l2,l2
lm2=lm2+1

do lm3=1,lmmax2(na)
if(abs(real(yyy(lm3,lm2,lm1,na))).lt.1.0d-8.and.abs(aimag(yyy(lm3,lm2,lm1,na))).lt.1.0d-8)cycle
j2tmp=0.0d0
do i1=1,2+nlo2(l1,na,spin1)
do i2=1,2+nlo2(l2,na,spin2)
if(spin1.eq.1)then
j2tmp=j2tmp+conjg(almblm2(i1,lm1,n2,k2,na,spin1))*almblm2(i2,lm2,n1,k1,na,spin2)*uuu(lm3,i2,l2,i1,l1,na)
else
j2tmp=j2tmp+conjg(almblm2(i1,lm1,n2,k2,na,spin1))*almblm2(i2,lm2,n1,k1,na,spin2)*uuu(lm3,i1,l1,i2,l2,na)
end if
end do
end do
j2(spin1,2)=j2(spin1,2)+j2tmp*yyy(lm3,lm2,lm1,na) !*(-1.d0)**spin2
end do
end do
end do
end do
end do
do n=1,nlorb(na)
l=lorb(n,na)
lm1=l*l
do m1=-l,l
lm1=lm1+1
lm2=l*l
do m2=-l,l
lm2=lm2+1
j2tmp=0.0d0
do i1=1,2+nlo2(l,na,spin1)
do i2=1,2+nlo2(l,na,spin2)
if(spin1.eq.1)then
j2tmp=j2tmp+conjg(almblm2(i1,lm1,n2,k2,na,spin1))*almblm2(i2,lm2,n1,k1,na,spin2)*orbuu(i2,i1,n,na)
else
j2tmp=j2tmp+conjg(almblm2(i1,lm1,n2,k2,na,spin1))*almblm2(i2,lm2,n1,k1,na,spin2)*orbuu(i1,i2,n,na)
end if
end do
end do
j2(spin1,2)=j2(spin1,2)+j2tmp*vorb(m2,m1,n,na) !*(-1.0d0)**spin2
end do
end do
end do
end do


do i=1,nxy
if(xy1(i).eq.'x'.and.xy2(i).eq.'x')then
sumtmp(i)=sumtmp(i)+j1*(j2(1,1)+j2(2,1))*(j2(1,2)+j2(2,2)) !xx
else if(xy1(i).eq.'y'.and.xy2(i).eq.'y')then
sumtmp(i)=sumtmp(i)-j1*(j2(1,1)-j2(2,1))*(j2(1,2)-j2(2,2)) !yy
else if(xy1(i).eq.'x'.and.xy2(i).eq.'y')then
sumtmp(i)=sumtmp(i)+xs*j1*(j2(1,1)+j2(2,1))*(j2(1,2)-j2(2,2)) !xy
else if(xy1(i).eq.'y'.and.xy2(i).eq.'x')then
sumtmp(i)=sumtmp(i)+xs*j1*(j2(1,1)-j2(2,1))*(j2(1,2)+j2(2,2))
else 
write(6,*)'x,y'
stop
end if
end do
!write(6,*)k1,k2,n1,n2
end do
end do
end if
end if
end if
end do
end do
!write(6,*)sk2(:,q)
write(22+cpuid,*)sk2(:,q)
do i=1,nxy
sumtmp(i)=sumtmp(i)/kmax2
write(22+cpuid,*)sumtmp(i)
end do



!write(22+cpuid,*)'now close it, q and cpuid =',q,cpuid
close(22+cpuid)
end do
!close(22)

write(6,*)m,nd,nx

return
end


