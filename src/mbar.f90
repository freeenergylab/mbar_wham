module mbar
use parameter
use extract
implicit none
private
public :: mbar2
contains
  subroutine mbar2
  implicit none
  integer(kind=4)                          ::i,j,k,l
  integer(kind=4)                          ::stat,system,numbw
  integer(kind=4),allocatable              ::Nk(:)
  real(kind=8)                             ::sum1,sum2
  real(kind=8),allocatable                 ::u(:,:,:)
  real(kind=8),allocatable                 ::ushift(:,:,:)
  real(kind=8),allocatable                 ::f(:),f_old(:),f_old_kcal(:)
  real(kind=8),allocatable                 ::deltaf(:)
  character(len=80)                        ::char_i
 
  allocate(Nk(nwindows))
  open(14,file='fort.13',status='old')
  do i=1,nwindows
    read(14,*) Nk(i)
  end do
  close(14)
  
  stat=system('rm fort.13')
   
  allocate(u(nwindows,nwindows,maxval(Nk)))
  allocate(ushift(nwindows,nwindows,maxval(Nk)))
  u=0.d0
  ushift=0.d0
  do i=1,nwindows
    do l=1,Nk(i)
      do k=1,nwindows
        u(i,k,l)=ub(i,k,l)
        ushift(i,i,l)=u(i,i,l)
      end do
    end do
  end do

!  print*,maxval(u),minval(u)

  do i=1,nwindows
    do l=1,Nk(i)
      do k=1,nwindows
        u(i,k,l)=u(i,k,l)-ushift(i,i,l)
      end do
    end do
  end do
  
!  print*,maxval(u),minval(u)

  allocate(f(nwindows))
  allocate(f_old(nwindows))
  allocate(f_old_kcal(nwindows))
  allocate(deltaf(nwindows))
  f=0.d0
!  f_old=0.d0 !first kind of initial guess (too long time needed!!!)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!initial guess from WHAM-fi in order to accelerate the convergence of MBAR!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(21,file='WHAM-fi.dat')  
  read(21,*) 
  do i=1,nwindows
    read(21,'(I2,2f13.4)') numbw,f_old_kcal(i),f_old(i)
  end do
  close(21)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  deltaf=0.d0
  open(15,file='MBAR-pmf.dat')
  do
    do i=1,nwindows 
      sum2=0.d0
      do j=1,nwindows
        do k=1,Nk(j)
          sum1=0.d0
          do l=1,nwindows
            sum1=sum1+Nk(l)*exp(f_old(l)-u(j,l,k)) !here u(j,l,k) represents all the energy value in the same simulation!                  
          end do
          !print*,'sum1=',sum1
          sum2=sum2+exp(-u(j,i,k))/sum1 !here u(j,i,k) represents the state i Hamiltonian energy of all frame structures in all the simulations!
        end do
        !print*,'sum2=',sum2
      end do  
      f(i)=-log(sum2)
    end do
    f=f-f(1)
    deltaf(1)=0.d0
    do i=2,nwindows 
      deltaf(i)=abs(f(i)-f_old(i))/abs(f_old(i))
    end do
    print*,"MBAR's convergence tolerance:",abs(maxval(deltaf))

    if(abs(maxval(deltaf)).lt.tolerance)then
      goto 111
    else
      f_old=f
    end if
  end do

111 do i=1,nwindows
      write(15,'(I2,a,f13.4,a)') i,"th window's free energy (relative to  1th window)= ", J_PER_CAL*(f(i))/beta,' kJ/mol'
    end do
  write(15,*) '------------------------------------------------------------------------'
  do i=1,nwindows
    write(15,'(I2,a,f13.4,a)') i,"th window's free energy (relative to  1th window)= ", f(i)/beta,' kcal/mol'
  end do
  write(15,*) '------------------------------------------------------------------------'
  close(15)

  open(16,file='MBAR-fi.dat')
  write(16,'(a)') "#windows,     in kcal/mol,    in reduced units"
  do i=1,nwindows
    write(16,'(I2,2f13.4)') i,(f(i)-f(1))/beta,(f(i)-f(1))
  end do
  close(16)

  if(allocated(Nk)) deallocate(Nk)
  if(allocated(u)) deallocate(u)
  if(allocated(f)) deallocate(f)
  if(allocated(f_old)) deallocate(f_old)
  if(allocated(deltaf)) deallocate(deltaf)
  end subroutine
end module
