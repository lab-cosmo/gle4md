! a couple of functions needed for dynamics, and shared variables
! obviously for any serious work better prng should be used

module md_tools
  implicit none
  integer ndim
  real*8 kt
  real*8, allocatable :: wtw(:,:)
contains

real*8 function rang(idum)
  implicit none
  integer, intent(inout) :: idum
  integer iset
  real(8) gset,v1,v2,fac,rsq
  save iset,gset
  data iset/0/
  if (iset.eq.0) then
1    v1=2.d0*ran2(idum)-1.d0
     v2=2.d0*ran2(idum)-1.d0
     rsq=v1**2+v2**2  
     if(rsq.ge.1.d0.or.rsq.eq.0.d0)goto 1 
     fac=sqrt(-2.d0*log(rsq)/rsq) 
     gset=v1 * fac
     rang=v2 * fac
     iset=1
  else 
     rang=gset
     iset=0 
  endif
  return
end function rang


real*8 function ran2(idum)
  implicit none
  real *8 x
  integer, intent(inout) :: idum
  integer iseed,i
  integer, allocatable :: seed(:)
    
  if(idum.le.0) then 
     idum=-idum
     call random_seed(size=iseed) 
     allocate(seed(iseed))
     do i=1,iseed  !ugly. once again, this is just a stub. you should use a GOOD prng!
       seed(i)=idum+i
     end do
     call random_seed(put=seed)
  endif
  call random_number(harvest=x)
  ran2=x
  return
end function ran2

subroutine force(q, f)
  real*8, intent(in) :: q(:)
  real*8, intent(out) :: f(:)
  integer i,j
  f=0.d0
  do i=1,ndim
    do j=1,ndim
      f(i)=f(i)-wtw(i,j)*q(j)
    end do
  end do
end subroutine

subroutine pot(q, v)
  real*8, intent(in) :: q(:)
  real*8, intent(out) :: v
  integer i,j
  v=0.d0
  do i=1,ndim
    do j=1,ndim
      v=v+q(i)*wtw(i,j)*q(j)
    end do
  end do
  
  v=v*0.5
end subroutine

subroutine kin(p, k)
  real*8, intent(in) :: p(:)
  real*8, intent(out) :: k
  integer i
  k=0.d0
  do i=1,ndim
    k=k+p(i)*p(i)
  enddo
  k=k*0.5d0
end subroutine

end module md_tools
