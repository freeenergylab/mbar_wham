module extract
use parameter
implicit none
real(kind=8),allocatable,save              ::ub(:,:,:)
contains
  subroutine extr
  implicit none
  integer(kind=4)                          ::i,k,l
  integer(kind=4)                          ::numb
  real(kind=8),allocatable                 ::x(:,:)
  real(kind=8),allocatable                 ::x_local(:),kw(:)
  character(len=80)                        ::char_i
  character(len=20),allocatable            ::filename(:)

  allocate(filename(nwindows))
  allocate(x_local(nwindows))
  allocate(kw(nwindows))
  open(11,file='./data/meta.dat',status='old')
  do i=1,nwindows
    read(11,'(a20,f8.1,f8.0)') filename(i),x_local(i),kw(i)
    !write(*,'(a20,f8.1,f8.0)') trim(adjustl(filename(i))),x_local(i),kw(i)
    !print*,trim(adjustl(filename(i)))
  end do
  close(11)
  
  allocate(x(nwindows,nframes))
  do i=1,nwindows
    open(12,file=trim(adjustl(filename(i))),status='old')
      do l=1,nframes
        read(12,'(I8,f10.4)') numb,x(i,l)    
        !write(*,'(I8,f10.4)') numb,x(i,l)    
      end do
    close(12)
  end do 
  
  allocate(ub(nwindows,nwindows,nframes))
  do i=1,nwindows
    do l=1,nframes
      do k=1,nwindows
        ub(i,k,l)=1.d0/2.d0*kw(k)*(x(i,l)-x_local(k))**2.d0*beta !in reduced potential
        !print*,ub(i,k,l)
      end do
    end do
  end do
  
  do i=1,nwindows
    write(char_i,'(I0)') i
    open(15,file="./energy/energy."//trim(adjustl(char_i))//".dat")
    do l=1,nframes
      write(15,'(31f13.4)') (ub(i,k,l),k=1,nwindows)
    end do
    close(15)
  end do

  open(13,file='fort.13')
  do i=1,nwindows
    write(13,*) nframes
  end do
  close(13)

!  print*,maxval(ub),minval(ub)

  end subroutine
end module
