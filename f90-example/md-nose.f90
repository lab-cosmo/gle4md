! very simple module for Nosé-Hoover chains.
! tersely commented, as it's here only for comparison with GLE.
module md_nose
  implicit none
  integer nc
  real*8, allocatable ::  pi(:,:), xi(:,:)
  real*8 :: nhq, onq, nhcham
contains
  subroutine nhc_init(nchain, wopt, irnd)
    use md_tools, only: rang, ndim, kt
    implicit none
    integer, intent(inout):: irnd
    integer, intent(in) :: nchain
    real*8,  intent(in) :: wopt
    integer i, j
    nc=nchain
    allocate(pi(nc,ndim))
    allocate(xi(nc,ndim))

    xi=0.
    nhq=kt/(wopt*wopt)
    onq=1.d0/nhq

    do i=1,nc
    do j=1,ndim
      pi(i,j)=rang(irnd)*sqrt(kt)
    enddo
    enddo
    nhcham=0.d0  ! sets to zero nose bit of the conserved quantity
  end subroutine

  subroutine nhc_step(p, fulldt, mts) 
    use md_tools, only:  ndim, kt
    implicit none
    real*8, intent(inout)     :: p(:)
    real*8, intent(in)        :: fulldt
    integer, intent(in)       ::mts

    real*8 dt, dt2
    real*8 a, b, c
    integer k, j, i
    dt=fulldt/mts
    dt2=0.5d0*dt

    do j=1,ndim
    do k=1,mts
      ! evolve p
      b=onq*pi(1,j)
      a=b*dt2
      xi(1,j)=xi(1,j)+a
      p(j)=p(j)*dexp(-a)   

      ! half time-step on odd xi's and even pi's
      c=b*pi(1,j)-kt
      do i=3,nc,2
        b=onq*pi(i,j)
        a=b*dt2
        xi(i,j)=xi(i,j)+a
        if (b.eq.0.d0) then
          pi(i-1,j)=pi(i-1,j)+c*dt2
        else
          c=c/b
          pi(i-1,j)=c+(pi(i-1,j)-c)*dexp(-a)
        endif
        c=b*pi(i,j)-kt
      enddo
      if (mod(nc,2).eq.0) then
        pi(nc,j)=pi(nc,j)+c*dt2
      endif

      ! full time step on odd pi's and even xi's
      c=p(j)*p(j)-kt  ! mass is one! 
      if (nc.eq.1) then
        pi(1,j)=pi(1,j)+c*dt
      else
        do i=2,nc,2
          b=onq*pi(i,j)
          a=b*dt
          xi(i,j)=xi(i,j)+a
          if (b.eq.0.d0) then
            pi(i-1,j)=pi(i-1,j)+c*dt
          else
            c=c/b
            pi(i-1,j)=(pi(i-1,j)-c)*dexp(-a)+c
          endif
          c=b*pi(i,j)-kt
        enddo
        if (mod(nc,2).ne.0) then
          pi(nc,j)=pi(nc,j)+c*dt
        endif
      endif

      ! evolve p
      b=onq*pi(1,j)
      a=b*dt2
      xi(1,j)=xi(1,j)+a
      p(j)=p(j)*dexp(-a)   

      ! half time-step on odd xi's and even pi's
      c=b*pi(1,j)-kt
      do i=3,nc,2
        b=onq*pi(i,j) 
        a=b*dt2
        xi(i,j)=xi(i,j)+a
        if (b.eq.0.d0) then
          pi(i-1,j)=pi(i-1,j)+c*dt2
        else
          c=c/b
          pi(i-1,j)=c+(pi(i-1,j)-c)*dexp(-a)
        endif
        c=b*pi(i,j)-kt
      enddo
      if (mod(nc,2).eq.0) then
        pi(nc,j)=pi(nc,j)+c*dt2
      endif
    enddo
    enddo
  end subroutine

  subroutine nhc_cons()
    use md_tools, only:  ndim, kt
    implicit none
    integer i,j

    nhcham=0.d0
    do j=1,ndim
    do i=1,nc
      nhcham=nhcham+pi(i,j)*pi(i,j)*0.5*onq+kt*xi(i,j)
    end do
    end do
  end subroutine
end module md_nose
