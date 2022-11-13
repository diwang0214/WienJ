module cl
implicit none
integer,parameter :: lmax=3
integer,parameter :: l21=7
integer :: nmax
complex*16,parameter :: xs=(0.0d0,1.0d0)
real*8 :: pi
character*20 :: foldername
integer :: flen
end module cl

module bl1  ! readstruct
implicit none
logical :: ortho,inversion
integer :: nat,nsym,ndif
character*4 :: lattic,irel,cform
character*80 :: title
character*10,allocatable :: aname(:)
real*8 :: br1(3,3),aa,bb,cc,br2(3,3),avec(3,3)
real*8 :: a(3),alpha(3),pia(3),vol
integer,allocatable :: iatnr(:),mult(:),jri(:),isplit(:) !NATOM
real*8,allocatable :: r0(:),dx(:)                        !NATOM
real*8,allocatable :: rotloc(:,:,:),rmt(:),zz(:),vsph(:) !NATOM
real*8,allocatable :: rotij(:,:,:),tauij(:,:)            !NDIF
real*8,pointer :: pos(:,:)
end module bl1

module bl2  ! readstruct
integer :: iord
integer,allocatable :: iz(:,:,:),inum(:),kkk(:,:) !NSYM
integer,allocatable :: kzz(:,:),inst(:)             !NVAW
real*8,allocatable :: tau(:,:)                      !NSYM
end module bl2

!module bl3
!implicit none
!complex*16,allocatable :: matrix(:,:,:)
!real*8 :: da,db,dc
!integer :: dd
!end module bl3

module bl4  ! readklist
implicit none
integer :: kmax,kmax2
real*8,allocatable :: sk(:,:),kweight(:)
real*8 :: e1,e2
real*8,allocatable :: sk2(:,:)
real*8 :: br_rec(3,3),br_dir(3,3)
end module bl4

module bl5  ! readalmd,readradwfd
use cl
implicit none
integer,allocatable :: nlo(:,:,:) !,nlodn(:,:)
!integer,allocatable :: kn(:,:,:) !,kndn(:,:)  ! kmax
!real*8,allocatable :: fermi(:,:,:,:),energy(:,:,:,:) !,fdn(:,:,:),edn(:,:,:)
!complex*16,allocatable :: almblm(:,:,:,:,:,:) !,almblmdn(:,:,:,:,:)
integer,allocatable :: nlo2(:,:,:)
integer,allocatable :: kn2(:,:) !,kn2dn(:,:)
real*8,allocatable :: fermi2(:,:,:),energy2(:,:,:)
complex*16,allocatable :: almblm2(:,:,:,:,:,:)
!real*8,allocatable :: r(:,:)
real*8,allocatable :: rr2(:,:)
!real*8,allocatable :: radwf(:,:,:,:,:,:) !,radwfdn(:,:,:,:,:)
real*8,allocatable :: radwf2(:,:,:,:,:,:)
end module bl5

module bl6
use cl
implicit none
real*8 :: spt(3,378),pweight(378)
complex*16 :: yl(81,378)
complex*16,allocatable :: yy(:,:,:)
complex*8,allocatable :: yyy(:,:,:,:)
real*8,allocatable :: uuu(:,:,:,:,:,:)
end module bl6

module bl7  ! readr2v
implicit none
integer :: jatom,lmmax
integer,allocatable :: lmmax2(:)
integer,allocatable :: lm(:,:)
integer,allocatable :: lmm2(:,:,:)
!real*8,allocatable :: vlm(:,:,:)
real*8,allocatable :: vlm2(:,:,:)
end module bl7

module bl8  ! +U
implicit none
integer :: nmod,nsp,natorb,iatom
integer :: nlorb(2),lorb(3,2)
real*8 :: rval,cval
complex*16 :: vorbup,vorbdn
complex*16 :: vorb(-3:3,-3:3,3,2)
real*8 :: orbuu(5,5,3,2)
end module bl8

module bl9 ! input
implicit none
integer :: atom(2),amu(2)
integer :: jspin
character,allocatable :: xy1(:),xy2(:)
integer :: so,nxy,qk
integer :: qmin,qmax
integer :: szfmafm,szjdm
end module bl9


!mpi

!include 'mkl_dss.f90'
  module prec
     !>> A module controls the precision. 
     !> when the nnzmax is larger than 2,147,483,647 then li=8,
     !> otherwise, li=4. 
     !> warning: li=4 was tested, li=8 is not tested yet.
     integer,parameter :: li=4 ! long integer
     integer,parameter :: Dp=kind(1.0d0) ! double precision  
     integer,parameter :: stdout=237 ! long integer
  end module prec

  module wmpi
     use prec

#if defined (MPI)
     include 'mpif.h'
#endif

     integer :: cpuid  ! CPU id for mpi
     integer :: num_cpu  ! Number of processors for mpi

#if defined (MPI)
     integer, parameter :: mpi_in= mpi_integer
     integer, parameter :: mpi_dp= mpi_double_precision
     integer, parameter :: mpi_dc= mpi_double_complex
     integer, parameter :: mpi_cmw= mpi_comm_world
#endif 

     !> Define a structure containing information for doing communication
     type WTParCSRComm
  
        !> mpi communicator
        integer :: comm
  
        !> how many cpus that we need to send data on
        integer :: NumSends
  
        !> which cpus that we need to send data on
        integer, pointer :: SendCPUs(:)
  
        !> when before we send the vector data to other cpus, we need to get the 
        !> data which should be sent, then put these data to a array called 
        !> x_buf_data(:). The array SendMapElements(:) gives the position of the 
        !> data in vector that should be sent.
        integer, pointer :: SendMapStarts(:)
  
        !> with this array, we can select the vector data that should be sent
        integer(li), pointer :: SendMapElements(:)
  
        !> How many cpus that we need to recieve data from
        integer :: NumRecvs
  
        !> Which cpus that we need to recieve data from
        integer, pointer :: RecvCPUs(:)
  
        !> When recieved data from other cpus, we need to arrange those data into 
        !> an array. The length of this 
        integer, pointer :: RecvVecStarts(:)
     end type WTParCSRComm
  
     !> Define a structure containing information for doing communication
     type WTParVecComm
  
        !> mpi communicator
        integer :: comm
  
        !> how many cpus that we need to send data on
        integer :: NumSends
  
        !> which cpus that we need to send data on
        integer, pointer :: SendCPUs(:)
  
        !> when before we send the vector data to other cpus, we need to get the 
        !> data which should be sent, then put these data to a array called 
        !> x_buf_data(:). The array SendMapElements(:) gives the position of the 
        !> data in vector that should be sent.
        integer, pointer :: RecvMapStarts(:)
  
        !> with this array, we can select the vector data that should be sent
        integer(li), pointer :: RecvMapElements(:)
  
        !> How many cpus that we need to recieve data from
        integer :: NumRecvs
  
        !> Which cpus that we need to recieve data from
        integer, pointer :: RecvCPUs(:)
  
        !> When recieved data from other cpus, we need to arrange those data into 
        !> an array. The length of this 
        integer, pointer :: SendVecStarts(:)
  
        integer :: NumRowsDiag
        integer :: NumRowsOffd
  
        integer(li), pointer :: RowMapOffd(:)
        integer(li), pointer :: LocalIndexOffd(:)
        integer(li), pointer :: LocalIndexDiag(:)
     end type WTParVecComm
  
  
  
     !> define a handle for comm, that can be created and destroyed
     type WTCommHandle
  
        type(WTParCSRComm), pointer :: sendrecv
  
        integer :: numrequest
        integer, pointer :: mpirequest(:)
        
        complex(dp), pointer :: senddata(:)
        complex(dp), pointer :: recvdata(:)
  
     end type WTCommHandle
  
     integer               :: BasisStart
     integer               :: BasisEnd
  
     contains
  
     !> generate partition for any vector with a given length 
  
     !> generate local partition for any vector with a given length 
  
     !> given the send data, recieve data and sendrecv list, we can use
     !> this subroutine to send and recieve data from other cpus.
     !> when finished this subroutine calls, we need call WTCommHandleDestroy
     !> to check whether the send recv operation is finished
  
     !> define my mpi_allreduce for complex array

  end module wmpi

