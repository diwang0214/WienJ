subroutine readklist
use cl
use bl4
implicit none
integer :: isx,isy,isz,idv
integer :: i,j
real*8 :: br(3,3)
real*8 :: tmp(3)
character*10 :: kname

!br=0.0d0
!do i=1,3
!do j=1,3
!br(j,i)=br_rec(i,j)
!end do
!end do

open(unit=10,file=foldername(1:flen)//"dat.klist")
i=1
read(10,"(A10,4I10,3F5.2,I10)")kname,isx,isy,isz,idv,tmp(1),e1,e2,kmax2
allocate(sk(3,kmax2+1),kweight(kmax2+1))
sk=0.0d0
kweight(1)=tmp(1)
do
tmp(1)=dble(isx)/dble(idv)
tmp(2)=dble(isy)/dble(idv)
tmp(3)=dble(isz)/dble(idv)
!call dmatrixprod(br_rec,tmp,sk(:,i),3)
sk(:,i)=tmp
i=i+1 
read(10,"(A10,4I10,F5.2)")kname,isx,isy,isz,idv,kweight(i)
if(kname.eq."END       ")exit
end do
kmax=i-1
close(10)
return
end
