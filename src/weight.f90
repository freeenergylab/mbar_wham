module weight
use parameter
use extract
implicit none
contains
  subroutine weight2
  implicit none
  integer(kind=4)                          ::i,j,k,l
  integer(kind=4)                          ::stat,system,numbw
  integer(kind=4),allocatable              ::Nk(:)
  real(kind=8)                             ::sum1
  real(kind=8),allocatable                 ::u(:,:,:)
  real(kind=8),allocatable                 ::ushift(:,:,:)
  real(kind=8),allocatable                 ::f(:),f_old(:),f_old_kcal(:)
  real(kind=8),allocatable                 ::weights(:,:,:)
  real(kind=8),allocatable                 ::one(:),two(:,:)
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
  f=0.d0
  open(21,file='WHAM-fi.dat')  
  read(21,*) 
  do i=1,nwindows
    read(21,'(I2,2f13.4)') numbw,f_old_kcal(i),f_old(i)
  end do
  close(21)

  allocate(weights(nwindows,nwindows,maxval(Nk)))
  do i=1,nwindows 
    do j=1,nwindows
      do k=1,Nk(j)
        sum1=0.d0
        do l=1,nwindows
          sum1=sum1+Nk(l)*exp(f_old(l)-u(j,l,k)) !here u(j,l,k) represents all the energy value in the same simulation!                  
        end do
        weights(j,i,k)=exp(f_old(i)-u(j,i,k))/sum1 !here u(j,i,k) represents the state i Hamiltonian energy of all frame structures in all the simulations!
      end do
    end do
  end do

  allocate(one(nwindows))
  do i=1,nwindows 
    one(i)=0.d0
    do j=1,nwindows
      do k=1,Nk(j)
        one(i)=one(i)+weights(j,i,k) 
      end do
    end do
    !print*,one(i)
  end do

  allocate(two(nwindows,maxval(Nk)))
  do j=1,nwindows
    do k=1,Nk(j)
      two(j,k)=0.d0
      do i=1,nwindows 
        two(j,k)=two(j,k)+Nk(i)*weights(j,i,k) 
      end do
      !print*,two(j,k)
    end do
  end do

  do i=1,nwindows 
    write(char_i,'(I5)') i
    open(15,file='./weights/weights'//trim(adjustl(char_i))//'.dat')
    do k=1,Nk(i) 
      !write(15,'(I5,31f20.10)') k,(weights(i,j,k),j=1,nwindows)
      write(15,'(I5,f20.10)') k,weights(i,i,k)
    end do
    close(15)  
  end do
  end subroutine
end module
