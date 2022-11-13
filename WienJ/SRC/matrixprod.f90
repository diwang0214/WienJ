subroutine matrixprod(matrix,ylm_old,ylm_new,m)  ! 
    implicit none
    integer :: m
    complex*16 :: matrix(m,m)
    complex*16 :: ylm_old(m)
    complex*16 :: ylm_new(m)
    integer :: m1,m2

    ylm_new=(0.0d0,0.0d0)
    do m1=1,m
        do m2=1,m
            ylm_new(m1)=ylm_new(m1)+matrix(m2,m1)*ylm_old(m2)
        end do
    end do
    
    return
end
subroutine dmatrixprod(matrix,ylm_old,ylm_new,m)  ! 
    implicit none
    integer :: m
    real*8 :: matrix(m,m)
    real*8 :: ylm_old(m)
    real*8 :: ylm_new(m)
    integer :: m1,m2

    ylm_new=0d0
    do m1=1,m
        do m2=1,m
            ylm_new(m1)=ylm_new(m1)+matrix(m2,m1)*ylm_old(m2)
        end do
    end do
    
    return
end
