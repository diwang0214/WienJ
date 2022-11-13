subroutine readr2v
use cl
use bl1
use bl7
use bl9
implicit none
real*8 :: vlmup(781),vlmdn(781)
real*8,allocatable :: vlmtmp(:,:) ! max(lmmax(na))
real*8 :: c_kub(0:10,0:10)
integer :: lm1,i,j,na,na2,nm

c_kub=0.0d0
      c_kub(0,0)=1.d0
      c_kub(3,2)=1.d0
      c_kub(4,0)=.5d0*SQRT(7.d0/3.d0)
      c_kub(4,4)=.5d0*SQRT(5.d0/3.d0)
      c_kub(6,0)=.5d0*SQRT(.5d0)
      c_kub(6,2)=.25d0*SQRT(11.d0)
      c_kub(6,4)=-.5d0*SQRT(7.d0/2.d0)
      c_kub(6,6)=-.25d0*SQRT(5.d0)
      c_kub(7,2)=.5d0*SQRT(13.d0/6.d0)
      c_kub(7,6)=.5d0*SQRT(11.d0/6.d0)
      c_kub(8,0)=.125d0*SQRT(33.d0)
      c_kub(8,4)=.25d0*SQRT(7.d0/3.d0)
      c_kub(8,8)=.125d0*SQRT(65.d0/3.d0)
      c_kub(9,2)=.25d0*SQRT(3.d0)
      c_kub(9,4)=.5d0*SQRT(17.d0/6.d0)
      c_kub(9,6)=-.25d0*SQRT(13.d0)
      c_kub(9,8)=-.5d0*SQRT(7.d0/6.d0)
      c_kub(10,0)=.125d0*SQRT(65.D0/6.D0)
      c_kub(10,2)=.125d0*SQRT(247.D0/6.D0)
      c_kub(10,4)=-.25d0*SQRT(11.D0/2.D0)
      c_kub(10,6)=0.0625d0*SQRT(19.D0/3.D0)
      c_kub(10,8)=-.125d0*SQRT(187.D0/6.D0)
      c_kub(10,10)=-.0625d0*SQRT(85.d0)

nm=0
allocate(lmmax2(jspin))
open(unit=15,file=foldername(1:flen)//"dat.r2v")
open(unit=16,file=foldername(1:flen)//"dat.r2vdn")
read(15,'(/,/)')
read(16,'(/,/)')
do na=1,nat
read(15,'(15X,I3)')jatom
read(15,'(15X,I3,/,/)')lmmax
!write(*,*)lmmax(na)
read(16,'(15X,I3)')jatom
read(16,'(15X,I3,/,/)')lmmax
!write(*,*)lmmax(na)
if(na.eq.1)then
allocate(lm(2,lmmax+100))
allocate(lmm2(2,lmmax+100,jspin))
allocate(vlmtmp(781,lmmax+100))
allocate(vlm2(781,lmmax+100,jspin))
end if
vlmtmp=0.0d0
do lm1=1,lmmax
read(15,'(15X,I3,5X,I2,/)')lm(1,lm1),lm(2,lm1)
!write(*,*)lm(:,lm1,na)
read(15,'(3X,4ES19.12)')(vlmup(j),j=1,781)
read(15,'(/)')
read(16,'(15X,I3,5X,I2,/)')lm(1,lm1),lm(2,lm1)
!write(*,*)lm(:,lm1,na)
read(16,'(3X,4ES19.12)')(vlmdn(j),j=1,781)
read(16,'(/)')
do j=1,781
vlmtmp(j,lm1)=(vlmup(j)-vlmdn(j))/2.d0  ! /(r(j)**2)
end do
end do
read(15,'(/,/,/)')
read(16,'(/,/,/)')

do na2=1,jspin
if(na.eq.atom(na2))then
lmmax2(na2)=lmmax
lmm2(:,:,na2)=lm
if(iatnr(na).gt.0)then
i=1
do
if(i.gt.lmmax2(na2))exit
if(lmm2(1,i,na2).eq.0.and.lmm2(2,i,na2).eq.0)then
vlm2(:,i,na2)=vlmtmp(:,i)
i=i+1
else if(lmm2(1,i,na2).eq.-3.and.lmm2(2,i,na2).eq.2)then
vlm2(:,i,na2)=vlmtmp(:,i)
i=i+1
else if(lmm2(1,i,na2).eq.4.or.lmm2(1,i,na2).eq.6.or.lmm2(1,i,na2).eq.-7)then
vlm2(:,i,na2)=c_kub(abs(lmm2(1,i,na2)),lmm2(2,i,na2))*c_kub(abs(lmm2(1,i,na2)),lmm2(2,i,na2))*vlmtmp(:,i)
vlm2(:,i+1,na2)=c_kub(abs(lmm2(1,i,na2)),lmm2(2,i,na2)+4)*c_kub(abs(lmm2(1,i,na2)),lmm2(2,i,na2))*vlmtmp(:,i)
i=i+2
end if
end do
else
do i=1,lmmax2(na2)
do j=1,781
vlm2(j,i,na2)=vlmtmp(j,i)
!if(vlm2(j,i,na2).gt.1.0d-8)write(6,*)na2,i,j,vlm2(j,i,na2)
end do
end do
end if
nm=nm+1
end if
end do
if(nm.eq.jspin)exit
end do
!close(15)
!close(16)

return
end

