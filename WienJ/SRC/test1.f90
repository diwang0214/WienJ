program test1
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
real*8 :: r1(3),r2(3),r12(3,2)  ! input
!character :: xy1,xy2   ! input
integer :: itape
character*22 :: filename(2)
character*80 :: kong
integer :: kongi
real*8 :: kongr
complex*16 :: kongc
integer :: i,j,k,l,m,n,index,index2
integer :: lm1,lm2,lm3
!real*8 :: mat(3,3,48)
complex*16 :: sumj,sum
logical :: alive
integer :: ierr                          
!integer,parameter :: stdout= 8
character*4  :: myStringq
integer :: rrmax
real*8,allocatable :: szrr(:,:),szrr2(:,:),szrr3(:,:),szrr4(:,:)
real*8,allocatable :: dr(:),dr2(:)
real*8,allocatable :: szsk(:,:)
complex*16,allocatable :: jq(:,:)
real*8 :: r(3),tmp,Rmax,rtau(3),sumjr(9)
real*8 :: jdmvalue(4)



!szadd
     ierr = 0
     cpuid= 0
     num_cpu= 1
     !> initial the environment of mpi
#if defined (MPI)
     call mpi_init(ierr)
     call mpi_comm_rank(mpi_cmw,cpuid,ierr)
     call mpi_comm_size(mpi_cmw,num_cpu,ierr)
#endif

!     if (cpuid==0) open(unit=stdout, file='szout')
     open(unit=stdout, file='szout')

     !> if mpi initial wrong, alarm
     if (cpuid==0.and.ierr.ne.0)then
        write(stdout,*)'mpi initialize wrong'
        stop
     endif

     !> print information for mpi
     if (cpuid==0) then
        write(stdout, '(1x, a, i5, a)')'You are using ', num_cpu, ' CPU cores'
        write(stdout, *)' '
     endif
!szaddend

pi=4.0d0*atan(1.0d0)

br_rec=0.0d0
br_dir=0.0d0
do i=1,3
br_rec(i,i)=1.0d0
br_dir(i,i)=1.0d0
end do

open(unit=6,file='dat.output')
write(6,*)'20220413-VERSION'
write(6,*)'FOR WIEN2k_18.2'

!inquire(file='../sz/dat.struct',exist=alive)
!if(alive)then
!foldername='../sz/'
!flen=6
!else
!inquire(file='../../sz/dat.struct',exist=alive)
!if(alive)then
!foldername='../../sz/'
!flen=9
!else
!write(6,*)'Can not find the input file'
!stop
!end if
!end if
!write(6,*)foldername


foldername='./'
flen=2

call readstruct

if(lattic(1:3).eq.'CXY')then
br_rec(1,1)=1.0d0
br_rec(1,2)=1.0d0
br_rec(1,3)=0.0d0
br_rec(2,1)=-1.0d0
br_rec(2,2)=1.0d0
br_rec(2,3)=0.0d0
br_rec(3,1)=0.0d0
br_rec(3,2)=0.0d0
br_rec(3,3)=1.0d0
br_dir(1,1)=0.5d0
br_dir(1,2)=-0.5d0
br_dir(1,3)=0.0d0
br_dir(2,1)=0.5d0
br_dir(2,2)=0.5d0
br_dir(2,3)=0.0d0
br_dir(3,1)=0.0d0
br_dir(3,2)=0.0d0
br_dir(3,3)=1.0d0
else if(lattic(1:3).eq.'CXZ')then
br_rec(1,1)=1.0d0
br_rec(1,2)=0.0d0
br_rec(1,3)=1.0d0
br_rec(2,1)=0.0d0
br_rec(2,2)=1.0d0
br_rec(2,3)=0.0d0
br_rec(3,1)=-1.0d0
br_rec(3,2)=0.0d0
br_rec(3,3)=1.0d0
br_dir(1,1)=0.5d0
br_dir(1,2)=0.0d0
br_dir(1,3)=-0.5d0
br_dir(2,1)=0.0d0
br_dir(2,2)=1.0d0
br_dir(2,3)=0.0d0
br_dir(3,1)=0.5d0
br_dir(3,2)=0.0d0
br_dir(3,3)=0.5d0
else if(lattic(1:3).eq.'CYZ')then
br_rec(1,1)=1.0d0
br_rec(1,2)=0.0d0
br_rec(1,3)=0.0d0
br_rec(2,1)=0.0d0
br_rec(2,2)=1.0d0
br_rec(2,3)=1.0d0
br_rec(3,1)=0.0d0
br_rec(3,2)=-1.0d0
br_rec(3,3)=1.0d0
br_dir(1,1)=1.0d0
br_dir(1,2)=0.0d0
br_dir(1,3)=0.0d0
br_dir(2,1)=0.0d0
br_dir(2,2)=0.5d0
br_dir(2,3)=-0.5d0
br_dir(3,1)=0.0d0
br_dir(3,2)=0.5d0
br_dir(3,3)=0.5d0
else if(lattic(1:1).eq.'F')then
br_rec=1.0d0
br_dir=0.5d0
do i=1,3
br_rec(i,i)=-1.0d0
br_dir(i,i)=0.0d0
end do
else if(lattic(1:1).eq.'B')then
br_rec=1.0d0
br_dir=0.5d0
do i=1,3
br_rec(i,i)=0.0d0
br_dir(i,i)=-0.5d0
end do
end if
write(6,*)'br_rec:'
write(6,*)br_rec
write(6,*)'br_dir'
write(6,*)br_dir

open(unit=97,file="dat.in")
read(97,*)atom(1),amu(1)
read(97,*)atom(2),amu(2)
!read(97,*)qk
qk=0
!read(97,*)qmin,qmax
qmin=0;qmax=1000
read(97,*)szfmafm
read(97,*)so
read(97,*)szjdm
nxy=4
if(so.eq.1)then
!read(97,*)nxy
nxy=4
allocate(xy1(nxy),xy2(nxy))
xy1(1)='x';xy2(1)='x'
xy1(2)='y';xy2(2)='y'
xy1(3)='x';xy2(3)='y'
xy1(4)='y';xy2(4)='x'
!do i=1,nxy
!read(97,*)xy1(i),xy2(i)
!write(6,*)xy1(i),xy2(i)
!end do
end if
close(97)

if(atom(1).eq.atom(2).and.amu(1).eq.amu(2))then
jspin=1
r1=0.0d0
else
jspin=2
do i=1,jspin
index=1
index2=amu(i)
do j=1,atom(i)-1
index=index+mult(j)
index2=index2+mult(j)
end do
do n=1,iord
m=1+iord-n
call dmatrixprod(dble(iz(:,:,m)),pos(:,index),r2,3)
r12(:,i)=r2+tau(:,m)-pos(:,index2)
call dmatrixprod(br_rec,r12(:,i),r2,3)
if(abs(exp(2.0d0*pi*xs*r2(1))-1.0d0).lt.1.d-8)then ! 2pi
if(abs(exp(2.0d0*pi*xs*r2(2))-1.0d0).lt.1.d-8)then
if(abs(exp(2.0d0*pi*xs*r2(3))-1.0d0).lt.1.d-8)then
write(6,*)m,r12(:,i)
exit
end if
end if
end if
end do
end do
r1=r12(:,2)-r12(:,1)
end if
write(6,*)jspin

call readklist
if(qmin.eq.0)then
qmin=1
qmax=kmax
else if(qmax.gt.kmax)then
write(6,*)'qmax>kmax'
stop
end if
write(6,*)kmax,kmax2
kmax2=kmax

allocate(sk2(3,kmax2))
allocate(nlo(0:lmax,nat,2))
!allocate(kn(kmax,nat,2))
!allocate(fermi(nmax,kmax,nat,2),energy(nmax,kmax,nat,2),almblm(5,(lmax+1)*(lmax+1),nmax,kmax,ndif,2))
allocate(nlo2(0:lmax,jspin,2))
allocate(kn2(kmax2,2))

nmax=0
filename(1)=foldername(1:flen)//'dat.almblmup'
filename(2)=foldername(1:flen)//'dat.almblmdn'
do lm3=1,2
itape=10+lm3
open(unit=itape,file=filename(lm3))
do i=1,4
read(itape,*)
end do
do lm1=1,nat
do l=0,lmax
read(itape,*)
read(itape,*)nlo(l,lm1,lm3)
do j=1,nlo(l,lm1,lm3)
read(itape,*)
end do
end do
do k=1,kmax
read(itape,*)
read(itape,*)
read(itape,*)l,m,kn2(k,lm3)
if(kn2(k,lm3).gt.nmax)nmax=kn2(k,lm3)
do n=1,kn2(k,lm3)
read(itape,*)
end do
do lm2=1,mult(lm1)
read(itape,*)
read(itape,*)
do n=1,kn2(k,lm3)
!index=0
do l=0,lmax
do m=-l,l
!index=index+1
read(itape,*)kongc,kongc
do j=1,nlo(l,lm1,lm3)
read(itape,*)
end do
end do
end do
end do
end do
end do
end do
close(itape)
end do
write(6,*)'nmax=',nmax

!allocate(fermi2(nmax,2,2),energy2(nmax,2,2),almblm2(5,(lmax+1)*(lmax+1),nmax,2,jspin,2))
!allocate(r(781,nat))
!allocate(radwf(2,5,0:lmax,781,nat,2))
allocate(rr2(781,jspin))
allocate(radwf2(2,5,0:lmax,781,jspin,2))
!call readalmd(11,"dat.almblmup",fermi2(:,:,1),energy2(:,:,1),almblm2(:,:,:,:,:,1),nlo(:,:,1),kn2(:,1))
!call readalmd(12,"dat.almblmdn",fermi2(:,:,2),energy2(:,:,2),almblm2(:,:,:,:,:,2),nlo(:,:,2),kn2(:,2))

do i=1,2
do n=1,jspin
nlo2(:,n,i)=nlo(:,atom(n),i)
end do
end do

!write(6,*)'almblm read'
filename(1)=foldername(1:flen)//"dat.radwfup"
filename(2)=foldername(1:flen)//"dat.radwfdn"
call readradwfd(13,filename(1),radwf2(:,:,:,:,:,1),nlo(:,:,1),rr2)
call readradwfd(14,filename(2),radwf2(:,:,:,:,:,2),nlo(:,:,2),rr2)
write(6,*)'radwf read'
call readr2v
write(6,*)'r2v read'

!lmmax2(1)=lmmax(atom(1))
!lmmax2(2)=lmmax(atom(2))
index=0
if(lmmax2(1).gt.lmmax2(jspin))then
m=lmmax2(1)
else
m=lmmax2(jspin)
end if
l=(lmax+1)*(lmax+1)
allocate(yy(m,378,jspin))
allocate(yyy(m,l,l,jspin))
allocate(uuu(m,5,0:lmax,5,0:lmax,jspin))
!allocate(vlm2(781,m,2))
call gpoint(spt,pweight)
do i=1,378
call ylm(spt(:,i),yl(:,i),8)
end do
do n=1,jspin
do i=1,378
do lm1=1,lmmax2(n)
call suml(yl(:,i),yy(lm1,i,n),lmm2(:,lm1,n))
end do
end do
do lm1=1,l
do lm2=1,l
do lm3=1,lmmax2(n)
sum=(0.0d0,0.0d0)
do i=1,378
sum=sum+conjg(yl(lm1,i))*yl(lm2,i)*yy(lm3,i,n)*pweight(i)
end do
if(abs(aimag(sum)).gt.1.0d-8)index=index+1
yyy(lm3,lm2,lm1,n)=sum
!if(abs(yyy(lm3,lm2,lm1,n)).gt.1.d-8)write(6,*)n,lm1,lm2,lm3
end do
end do
end do
end do
write(6,*)'yyy aimag=',index

!call k2k2(atom,amu)
do k=1,kmax
sk2(:,k)=sk(:,k)
end do

do n=1,jspin
do lm1=1,lmmax2(n)
do j=1,781
vlm2(j,lm1,n)=vlm2(j,lm1,n)/rr2(j,n)**2
end do
end do
end do
write(6,*)'cc'

call jxjf

!if(qk.eq.0)then
if(so.eq.1)then
call caljqso(r1)
else
call caljq(r1)
end if
!else
!if(so.eq.1)then
!call caljqso_bx(r1)
!else
!call caljq_bx(r1)
!end if
!end if
!write(6,*)'dd'


!note: need it, or the code will close when 1 core close
call mpi_finalize(ierr)

!close(6)
!stop

!new add  combine

pi=4.0d0*atan(1.0d0)

br_rec=0.0d0
br_dir=0.0d0
do i=1,3
br_rec(i,i)=1.0d0
br_dir(i,i)=1.0d0
end do

open(unit=6,file='j.output')

!inquire(file='../sz/dat.struct',exist=alive)
!if(alive)then
!foldername='../sz/'
!flen=6
!else
!inquire(file='../../sz/dat.struct',exist=alive)
!if(alive)then
!foldername='../../sz/'
!flen=9
!else
!write(*,*)'Can not find the struct file'
!stop
!end if
!end if

foldername='./'
flen=2

!call readstruct
call latgen
do i=1,3
write(6,'(3F11.8)')br1(:,i)
end do
br2=br1
if(lattic(1:1).eq.'R')then
do i=1,3
do j=1,3
br2(i,j)=br1(j,i)
end do
end do
end if

if(lattic(1:3).eq.'CXY')then
br_rec(1,1)=1.0d0
br_rec(1,2)=1.0d0
br_rec(1,3)=0.0d0
br_rec(2,1)=-1.0d0
br_rec(2,2)=1.0d0
br_rec(2,3)=0.0d0
br_rec(3,1)=0.0d0
br_rec(3,2)=0.0d0
br_rec(3,3)=1.0d0
br_dir(1,1)=0.5d0
br_dir(1,2)=-0.5d0
br_dir(1,3)=0.0d0
br_dir(2,1)=0.5d0
br_dir(2,2)=0.5d0
br_dir(2,3)=0.0d0
br_dir(3,1)=0.0d0
br_dir(3,2)=0.0d0
br_dir(3,3)=1.0d0
else if(lattic(1:3).eq.'CXZ')then
br_rec(1,1)=1.0d0
br_rec(1,2)=0.0d0
br_rec(1,3)=1.0d0
br_rec(2,1)=0.0d0
br_rec(2,2)=1.0d0
br_rec(2,3)=0.0d0
br_rec(3,1)=-1.0d0
br_rec(3,2)=0.0d0
br_rec(3,3)=1.0d0
br_dir(1,1)=0.5d0
br_dir(1,2)=0.0d0
br_dir(1,3)=-0.5d0
br_dir(2,1)=0.0d0
br_dir(2,2)=1.0d0
br_dir(2,3)=0.0d0
br_dir(3,1)=0.5d0
br_dir(3,2)=0.0d0
br_dir(3,3)=0.5d0
else if(lattic(1:3).eq.'CYZ')then
br_rec(1,1)=1.0d0
br_rec(1,2)=0.0d0
br_rec(1,3)=0.0d0
br_rec(2,1)=0.0d0
br_rec(2,2)=1.0d0
br_rec(2,3)=1.0d0
br_rec(3,1)=0.0d0
br_rec(3,2)=-1.0d0
br_rec(3,3)=1.0d0
br_dir(1,1)=1.0d0
br_dir(1,2)=0.0d0
br_dir(1,3)=0.0d0
br_dir(2,1)=0.0d0
br_dir(2,2)=0.5d0
br_dir(2,3)=-0.5d0
br_dir(3,1)=0.0d0
br_dir(3,2)=0.5d0
br_dir(3,3)=0.5d0
else if(lattic(1:1).eq.'F')then
br_rec=1.0d0
br_dir=0.5d0
do i=1,3
br_rec(i,i)=-1.0d0
br_dir(i,i)=0.0d0
end do
else if(lattic(1:1).eq.'B')then
br_rec=1.0d0
br_dir=0.5d0
do i=1,3
br_rec(i,i)=0.0d0
br_dir(i,i)=-0.5d0
end do
end if



open(unit=97,file="dat.in")
read(97,*)atom(1),amu(1)
read(97,*)atom(2),amu(2)
!read(97,*)qk
qk=0
!read(97,*)qmin,qmax
qmin=0;qmax=1000
read(97,*)szfmafm
read(97,*)so
read(97,*)szjdm
nxy=4
if(so.eq.1)then
!read(97,*)nxy
nxy=4
!allocate(xy1(nxy),xy2(nxy))
xy1(1)='x';xy2(1)='x'
xy1(2)='y';xy2(2)='y'
xy1(3)='x';xy2(3)='y'
xy1(4)='y';xy2(4)='x'
!do i=1,nxy
!read(97,*)xy1(i),xy2(i)
!write(6,*)xy1(i),xy2(i)
!end do
end if
read(97,*)Rmax
Rmax=Rmax/0.5291772083d0
close(97)



!open(unit=97,file="dat.in")
!read(97,*)atom(1),amu(1)
!read(97,*)atom(2),amu(2)
!read(97,*)qk
!read(97,*)qmin,qmax
!read(97,*)so
!nxy=1
!if(so.eq.1)then
!read(97,*)nxy
!allocate(xy1(nxy),xy2(nxy))
!do i=1,nxy
!read(97,*)xy1(i),xy2(i)
!!write(6,*)xy1(i),xy2(i)
!end do
!end if
!close(97)

i=0
do j=1,nat
do k=1,mult(j)
i=i+1
if(j.eq.atom(1).and.k.eq.amu(1))index=i
if(j.eq.atom(2).and.k.eq.amu(2))index2=i
end do
end do
do i=1,3
rtau(i)=pos(i,index2)-pos(i,index)
end do


open(unit=10,file='./data/jq_0001')
read(10,*)kmax,so,nxy
allocate(szsk(3,kmax),jq(kmax,nxy))
read(10,*)r1
close(10)

do k=1, kmax
write(myStringq,'(i4.4)') k
open(unit=10,file='./data/jq_'//myStringq)
if(k.eq.1)then
read(10,*)kmax,so,nxy
read(10,*)r1
end if
read(10,*)szsk(:,k)
do i=1,nxy
read(10,*)jq(k,i)
end do
end do

!open(unit=10,file='jq.txt')
!read(10,*)kmax,so,nxy
!allocate(sk(3,kmax),jq(kmax,nxy))
!read(10,*)r1
!do k=1,kmax
!read(10,*)sk(:,k)
!do i=1,nxy
!read(10,*)jq(k,i)
!end do
!end do
!close(10)

!write(*,*)'20171212'
!write(*,*)"please give Rmax (A): "
!read(*,*)Rmax
!Rmax=Rmax/0.5291772083d0

if(aa.lt.bb.and.aa.lt.cc)then
rrmax=INT(Rmax/aa*2.0d0)
else if(bb.lt.cc)then
rrmax=INT(Rmax/bb*2.0d0)
else
rrmax=INT(Rmax/cc*2.0d0)
end if

index=(2*rrmax+1)**3
allocate(szrr(3,index),szrr2(3,index),szrr3(3,index),szrr4(3,index))
allocate(dr(index),dr2(index))

index=0
do i=-rrmax,rrmax
do j=-rrmax,rrmax
do k=-rrmax,rrmax
index=index+1
r(1)=dble(i)
r(2)=dble(j)
r(3)=dble(k)
call dmatrixprod(br_dir,r,szrr(:,index),3)
szrr2(:,index)=szrr(:,index)+rtau
call dmatrixprod(br2,szrr2(:,index),r,3)
szrr3(:,index)=r
dr(index)=sqrt(r(1)**2+r(2)**2+r(3)**2)
end do
end do
end do

index2=0
do i=1,index
if(dr(i).lt.Rmax)then
index2=index2+1
szrr2(:,index2)=szrr(:,i)
szrr4(:,index2)=szrr3(:,i)
dr2(index2)=dr(i)
do j=index2-1,1,-1
if(dr2(j+1).lt.dr2(j))then
r=szrr2(:,j+1)
szrr2(:,j+1)=szrr2(:,j)
szrr2(:,j)=r
r=szrr4(:,j+1)
szrr4(:,j+1)=szrr4(:,j)
szrr4(:,j)=r
tmp=dr2(j+1)
dr2(j+1)=dr2(j)
dr2(j)=tmp
end if
end do
end if
end do



sumjr=0.0d0
open(unit=7,file='j.txt')
if(szjdm.eq.1)then
write(7,*)"    R               distance    J"
end if
if(szjdm.eq.2)then
write(7,*)"    R               distance    DM"
end if
if(szjdm.eq.3)then
write(7,*)"    R               distance    J          DM"
end if
do index=1,index2

r=szrr2(:,index)-r1
do j=1,nxy
sum=0.0d0
do k=1,kmax
tmp=0.0d0
do i=1,3
tmp=tmp+r(i)*szsk(i,k)
end do
sum=sum+jq(k,j)*exp(-2.0d0*pi*xs*tmp)
end do
sumjr(j)=sumjr(j)+real(sum/kmax)
jdmvalue(j)=real(sum/kmax*13606.0d0)
if(szfmafm.eq.0)then
jdmvalue(j)=jdmvalue(j)*(-1.0d0)
end if
end do
if(szjdm.eq.1 .and. dr2(index).gt.0.01d0)then
write(7,'(3f5.1,2f11.6)')szrr2(:,index),dr2(index)*0.5291772083d0,(jdmvalue(1)+jdmvalue(2))/2
end if
if(szjdm.eq.2 .and. dr2(index).gt.0.01d0)then
write(7,'(3f5.1,2f11.6)')szrr2(:,index),dr2(index)*0.5291772083d0,(jdmvalue(3)-jdmvalue(4))/2
end if
if(szjdm.eq.3 .and. dr2(index).gt.0.01d0)then
write(7,'(3f5.1,3f11.6)')szrr2(:,index),dr2(index)*0.5291772083d0,(jdmvalue(1)+jdmvalue(2))/2,(jdmvalue(3)-jdmvalue(4))/2
end if

end do
close(7)

do j=1,nxy
write(6,*)'sumjr=',sumjr(j)
write(6,*)'j(q=0)=',jq(1,j)
end do
close(6)



end
