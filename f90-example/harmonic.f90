! *********************************************************
! * An awkwardly simple code for the dynamics of an       *
! * an n-dimensional harmonic oscillator, to demonstrate  *
! * the use of colored-noise Langevin equation in         * 
! * thermostatting.                                       *   
! *                                                       *
! * Code has been kept as modular as possible, and it     *
! * should be relatively straightforward to adapt it for  *
! * other serial MD codes. In parallel codes using domain *
! * decomposition, the additional degrees of freedom      *
! * should "follow" the corresponding atoms.              *
! *                                                       *
! * Feel free to copy, modify, lend, borrow, nuke this    *
! * code, which is licensed under GPLv3 [www.gnu.org]     *
! * e-mail me at michele dot ceriotti at gmail dot com    *
! *********************************************************

program harmonic
  use md_tools
  use md_nose
  use md_gle

  implicit none

  ! * options to be given in input file *
  ! seed              initial seed for PRNG
  ! temp              target temperature
  ! dt                timestep
  ! nstep             number of steps to be performed
  ! stride            output every stride frames
  ! ndof              dimensionality of the oscillator
  ! wfile             name of the file to read the hessian 
  !                   decomposition from (see below)
  ! traj              logical, whether to output p,q trajectory
  ! * thermostatting options * 
  ! wnw               "optimal" frequency for langevin (0. 0-> WN off) 
  ! glew              range shift freq. for GLE (0. 0-> GLE off) 
  ! nchains           number of NH chains to be used (0 -> NHC off)
  ! nhmts             timestep subdivision for NHC integration
  ! nhw               "optimal" frequency for NHC (Q=kt/nhw**2)
  ! * a note on thermostats: all the thermostats can be used at once,
  !   even if that makes not much sense. again, this is just a small
  !   test code
  ! * a note on units: mass is set to one, and k_B as well, so that
  !   at equilibrium <p^2>=omega^2 <q^2>=temp
  
  real*8 dum, dt, temp, nhw, wnw, glew
  integer argc, seed, nstep, stride, ndof, nchains, nhmts, tstride
  character *256 fname, wfile, prefix
  namelist /inp/ seed, wfile, dt, temp, nstep, stride, tstride, ndof, nchains, nhmts, nhw, wnw, glew

  real*8, allocatable :: q(:), p(:), f(:), wm(:,:)
  real*8 v, k, h ! yes, it's all we need!
  real*8 dt2
#ifdef USELIBS
  integer, allocatable :: ipiv(:)  ! needed to invert matrix
#endif
  integer irnd, istep, i, j

  ! reads command line
  argc = iargc()
  if ( argc .ne. 1 ) then
    write(6,*) '* Call me as: harmonic <input>'
    stop
  endif  

  ! reads input file
  tstride=0
  call getarg (1,fname)
  open(101,file=fname)
  read(101,inp)
  close (unit=101)
  irnd=-seed
  dt2=dt*0.5
  kt=temp
  ndim=ndof
  allocate(q(ndof))
  allocate(p(ndof))
  allocate(f(ndof))

  allocate(wtw(ndof,ndof))
  allocate(wm(ndof,ndof))
#ifdef USELIBS
  allocate(ipiv(ndof))
#endif 

  ! init random seed
  dum=ran2(irnd)
  
  ! reads square root of the hessian (v(q)=1/2 q^T W^T W q)
  open(101,file=wfile)
  read(101,*) wm
  close (unit=101)
  wtw=matmul(transpose(wm),wm)

  ! init momenta
  do i=1,ndof
    p(i)=rang(irnd)*sqrt(temp)
  enddo
  ! init coordinates (requires inverting wm ONLY WORKS IF WM IS SQRT(HESSIAN), NOT IF DONE WITH CHOLESKY!!!)
#ifdef USELIBS
  do i=1,ndof
    q(i)=rang(irnd)*sqrt(temp)
  enddo 
  call dgesv(ndof, 1, wm, ndof, ipiv, q, ndof, i)
#else
  ! define your own inversion if you wish. 
  ! I am lazy and just init to zero.
  q=0.
#endif  
  call force(q,f)  

  !initializes thermostats
  if(nchains .gt. 0) call nhc_init(nchains, nhw, irnd) 
  if(glew .gt. 0.d0) call gle_init(dt2,glew,irnd)
  if(wnw .gt. 0.d0) call wn_init(dt2,wnw)

  ! opens file for output
  if (tstride>0) open(102,file='traj-p.out') 
  if (tstride>0) open(103,file='traj-q.out') 
  open(104,file='statis.out') 

  ! we are already at the main dynamics loop!
  do istep=1,nstep

    !calls the active thermostats
    if(nchains .gt. 0)    call nhc_step(p,dt2,nhmts)
    if (wnw .gt. 0.d0 .or. glew .gt. 0.d0) then 
      call kin(p, k)
      langham=langham+k
      if (wnw .gt. 0.d0)  call wn_step(p,irnd)
      if (glew .gt. 0.d0) call gle_step(p,irnd)
      call kin(p, k)
      langham=langham-k
    endif

    ! hamiltonian step for dynamics
    p=p+f*dt2
    q=q+p*dt
    call force(q,f)
    p=p+f*dt2

    ! thermostats, second bit
    if (wnw .gt. 0.d0 .or. glew .gt. 0.d0) then 
      call kin(p, k)
      langham=langham+k
      if (wnw .gt. 0.d0)  call wn_step(p,irnd)
      if (glew .gt. 0.d0) call gle_step(p,irnd)
      call kin(p, k)
      langham=langham-k
    endif
    if(nchains .gt. 0)  call nhc_step(p,dt2,nhmts)

    ! computes properties & outputs
    if (mod(istep,stride).eq.0) then
      call pot(q, v)
      call kin(p, k)
      h=k+v
      if (nchains.gt.0) then
        call nhc_cons()
        h=h+nhcham
      end if
      if (wnw .gt. 0.d0 .or. glew .gt. 0.d0) then
        h=h+langham
      end if
      if (tstride.gt.0 .and. mod(istep,tstride).eq.0) then
        write(102,'(1e17.8 )', advance='NO') istep*dt
        write(103,'(1e17.8 )', advance='NO') istep*dt
        do i=1,ndof
          write(102,'( 1e17.8 )', advance='NO') p(i)
          write(103,'( 1e17.8 )', advance='NO') q(i)
        enddo
        write(102,*) ""
        write(103,*) ""
      endif
      write(104,'(4e17.8)') istep*dt, v/ndof, k/ndof, h/ndof
    endif
  enddo

  if (tstride>0) close(102)
  if (tstride>0) close(103)
  close(104)
end program harmonic
