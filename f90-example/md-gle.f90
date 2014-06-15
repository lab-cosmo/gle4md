! *********************************************************
! * here we define the functions to be called for         *
! * colored-noise thermostatting. incidentally, we also   *
! * write here a couple of functions for white-noise.     *
! *                                                       *
! * code is licensed under GPLv3 [www.gnu.org]            *
! * please consider citing the relevant papers (listed    *
! * below) if you use GLE in your simulations.            *
! *                                                       *
! * e-mail me at michele dot ceriotti at gmail dot com    *
! *********************************************************

module md_gle
  implicit none
  real*8, allocatable, save ::  gS(:,:), gT(:,:), gp(:,:), ngp(:,:)
  real*8 wnt, wns, langham
  integer ns
contains

  ! initialize white-noise thermostat. 
  ! the init and the propagator here are written with the same phylosophy 
  ! used for the full-fledged colored-noise stuff.
  subroutine wn_init(dt,wopt)
    use md_tools, only: kt
    implicit none
    real*8, intent(in) :: dt, wopt
    real*8 g
    g=2.*wopt
    wnt=exp(-dt*g)
    wns=sqrt(kt*(1.-wnt*wnt))
    langham=0.d0  ! sets to zero accumulator for langevin 'conserved' quantity
  end subroutine

  ! white-noise propagator. time-step has been set in wn_init
  subroutine wn_step(p, irnd)
    use md_tools, only: kt, rang, ndim
    implicit none
    real*8, intent(inout)  :: p(:)
    integer, intent(inout) :: irnd
    integer i
    do i=1,ndim
      p(i)=wnt*p(i)+wns*rang(irnd)
    end do
  end subroutine
  
  ! initialize gle_init
  subroutine gle_init(dt,wopt,irnd)
    use md_tools, only: kt, rang, ndim
    implicit none
    real*8, intent(in)  :: dt, wopt
    integer, intent(inout) :: irnd
    real *8, allocatable :: gA(:,:), gC(:,:), gr(:)
    integer i, j, k, h, cns, ios
    

    write(6,*) "# Initialization of GLE thermostat.                           "
    write(6,*) "# Please cite the relevant works among:                       "
    write(6,*) "#                                                             "
    write(6,*) "# M. Ceriotti, G. Bussi and M. Parrinello                     "
    write(6,*) "# Phys. Rev. Lett. 102, 020601 (2009)                         "
    write(6,*) "#                                                             "
    write(6,*) "# M. Ceriotti, G. Bussi and M. Parrinello                     "
    write(6,*) "# Phys. Rev. Lett. 103, 030603 (2009)                         "
    write(6,*) "#                                                             "
    write(6,*) "# M. Ceriotti, G. Bussi and M. Parrinello                     "
    write(6,*) "# J. Chem. Theory Comput. 6, 1170 (2010)                      "
    write(6,*) "#                                                             "

    !reads in matrices
    !reads A (in units of the "optimal" frequency of the fitting)
    open(121,file='GLE-A',status='OLD',iostat=ios)
    read(121,*) ns

    !allocate everything we need
    allocate(gA(ns+1,ns+1))
    allocate(gC(ns+1,ns+1))
    allocate(gS(ns+1,ns+1))
    allocate(gT(ns+1,ns+1))
    allocate(gp(ndim,ns+1))   
    allocate(ngp(ndim,ns+1))   
    allocate(gr(ns+1))
    
    if (ios.ne.0) write(0,*) "Error: could not read GLE-A file!"
    do i=1,ns+1
       read(121,*) gA(i,:)
    enddo
    close(121)

    ! gamma for a WN langevin will be 1/tau, which would make it optimal for w=1/(2tau) angular frequency. 
    ! we scale gA (which is expected to be fitted such that maximum fitted frequency is one) accordingly 
    gA=gA*wopt

    ! reads C (in K), or init to kT
    open(121,file='GLE-C',status='OLD',iostat=ios)
    if (ios.ne.0) then            
       write(6,*) "# Using canonical-sampling, Cp=kT"
       gC=0.
       do i=1,ns+1
          gC(i,i)=kt
       enddo
    else    
       write(6,*) "# Reading specialized Cp matrix"
       read(121,*) cns
       if (cns.ne.ns) write(0,*) " Error: size mismatch between given GLE-A and GLE-C!"
       do i=1,ns+1
          read(121,*) gC(i,:)
       enddo
    end if
    close(121)

    ! the deterministic part of the propagator is obtained in a second
    call matrix_exp(-dt*gA, ns+1,15,15,gT)
    
    ! the stochastic part is just as easy. we use gA as a temporary array
    gA=gC-matmul(gT,matmul(gC,transpose(gT)))
    call cholesky(gA, gS, ns+1)


    ! then, we must initialize the auxiliary vectors. we keep general - as we might be 
    ! using non-diagonal C to break detailed balance - and we use cholesky decomposition
    ! of C. again, since one would like to initialize correctly the velocities in 
    ! case of generic C, we use an extra slot for gp for the physical momentum, as we 
    ! could then use it to initialize the momentum in the calling code
    gA=gC   
    call cholesky(gA, gC, ns+1)
    
    do j=1,ndim
      do i=1,ns+1
        gr(i)=rang(irnd)
      enddo
      gp(j,:)=matmul(gC,gr)
    end do
    
    deallocate(gA)
    deallocate(gC)
    deallocate(gr)
    langham=0.d0  ! sets to zero accumulator for langevin 'conserved' quantity
  end subroutine gle_init

  ! the GLE propagator. 
  ! gT contains the deterministic (friction) part, and gS the stochastic (diffusion) part.
  ! gp(j,1) must be filled with the mass-scaled actual j-th momentum, and contains in
  ! gp(j,2:ns+1) the current values of additional momenta. 
  ! the matrix multiplies are performed on the whole array at once, and the new momentum
  ! is passed back to the caller, while the new s's are kept stored in gp.
  ! please note that one can avoid the double conversion between mass-scaled and actual
  ! momentum/velocity (as used by the calling code) by scaling with the mass the white-noise
  ! random numbers. 
  subroutine gle_step(p,irnd)
    use md_tools
    implicit none
    integer,intent(inout)  :: irnd
    real *8, intent(inout) :: p(:)
    integer i, j
    real*8 mfac, totm

    do j=1,ndim
       gp(j,1)=p(j)   !<-- if m!= 1, here a proper scaling must be performed
    enddo

#ifdef USELIBS
    call dgemm('n','t',ndim,ns+1,ns+1,1.0d0,gp,ndim,gT,ns+1,0.0d0,ngp,ndim)
#else
    ngp=transpose(matmul(gT,transpose(gp)))
#endif

    !now, must compute random part. 
    !first, fill up gp of random n
    do j=1,ndim
      do i=1,ns+1
        gp(j,i)=rang(irnd)     !<-- if m!= 1, alternatively one could perform the scaling here (check also init!)
      end do
    end do

#ifdef USELIBS    
    call dgemm('n','t',ndim,ns+1,ns+1,1.0d0,gp,ndim,gS,ns+1,1.0d0,ngp,ndim)
    gp=ngp
#else
    gp=ngp+transpose(matmul(gS,transpose(gp)))
#endif

    do j=1,ndim
      p(j)=gp(j,1)     !<-- if m!= 1, here a proper inverse scaling must be performed
    end do
  end subroutine gle_step

  ! matrix exponential by scale & square.
  ! one can diagonalize with lapack, but it's not worth it, as 
  ! we call this routine only once!      
  subroutine matrix_exp(M, n, j, k, EM)
    integer, intent(in)  :: n, j, k
    real*8, intent(in)   :: M(n,n)
    real*8, intent(out)   :: EM(n,n)
    
    real *8 :: tc(j+1), tmp(n,n), SM(n,n)
    integer p, i
    tc(1)=1
    do i=1,j
       tc(i+1)=tc(i)/dble(i)
    enddo
    
    !scale
    SM=M*(1./2.**k)
    EM=0.
    do i=1,n
       EM(i,i)=tc(j+1)
    enddo
    
    !taylor exp of scaled matrix
    do p=j,1,-1
       EM=matmul(SM,EM);
       do i=1,n
          EM(i,i)=EM(i,i)+tc(p)
       enddo
    enddo
    
    !square
    do p=1,k
       EM=matmul(EM,EM)
    enddo
  end subroutine matrix_exp
  
  ! brute-force "stabilized" cholesky decomposition.
  ! in practice, we compute LDL^T decomposition, and force
  ! to zero negative eigenvalues.
  subroutine cholesky(SST, S, n)
    integer, intent(in)  :: n
    real*8, intent(in)   :: SST(n,n)
    real*8, intent(out)   :: S(n,n)
    real *8 :: L(n,n), D(n,n) 
    integer i,j,k
    S=0.
    L=0.
    D=0.
    do i=1,n
       L(i,i)=1.0
       do j=1,i-1
          L(i,j)=SST(i,j);
          do k=1,j-1
             L(i,j)=L(i,j)-L(i,k)*L(j,k)*D(k,k)
          enddo
          if (D(j,j).ne. 0.0) then
            L(i,j)=L(i,j)/D(j,j) 
          else
            write(0,*) "Warning: zero eigenvalue in LDL^T decomposition."
            L(i,j)=0.
          end if
       enddo
       D(i,i)=SST(i,i)
       do k=1,i-1
          D(i,i)=D(i,i)-L(i,k)**2*D(k,k)
       end do
    enddo
    do i=1,n
       if ( D(i,i).ge. 0.0d0 ) then
         D(i,i)=sqrt(D(i,i))
       else
         write(0,*) "Warning: negative eigenvalue (",D(i,i),")in LDL^T decomposition."
         D(i,i)=0.0
       end if
    end do
    S=matmul(L,D)
  end subroutine cholesky
end module
