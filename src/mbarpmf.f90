module mbarpmf
use parameter
implicit none
contains
  subroutine mbarpmf2
  implicit none
  integer(kind=4)                          ::i,k,l
  integer(kind=4)                          ::npoint,bindex
  integer(kind=4)                          ::numb,numbw
  integer(kind=4),allocatable              ::Ni(:),bins(:,:),sum_bins(:)
  real(kind=8)                             ::bin_width
  real(kind=8),allocatable                 ::x(:,:)
  real(kind=8),allocatable                 ::x_local(:),kw(:)
  real(kind=8),allocatable                 ::rc(:)
  real(kind=8),allocatable                 ::weight_factor(:,:),bias_factor(:,:),sum_bias_factor(:)
  real(kind=8),allocatable                 ::f_i(:),f_i_old(:),f_i_i(:),deltaf_i(:),p_l(:),pmf(:)
  real(kind=8),allocatable                 ::f_i_old_kcal(:)
  character(len=20),allocatable            ::filename(:)

  bin_width=(rcmax-rcmin)/real(binsnum)
  !print*,bin_width

  allocate(filename(nwindows))
  allocate(x_local(nwindows))
  allocate(kw(nwindows))
  open(11,file='./data/meta.dat',status='old')
  do i=1,nwindows
    read(11,'(a20,f8.1,f8.0)') filename(i),x_local(i),kw(i)
    !print*, filename(i),x_local(i),kw(i)
  end do
  close(11)
  
  allocate(x(nwindows,nframes))
  allocate(Ni(nwindows))
  allocate(bins(nwindows,binsnum))
  bins=0.d0
  do i=1,nwindows
    open(12,file=trim(adjustl(filename(i))),status='old')
    npoint=0
    do l=1,nframes
      read(12,'(I8,f10.4)') numb,x(i,l)    
      if (x(i,l)>rcmin.and.x(i,l)<=rcmax) then
        npoint=npoint+1
        bindex=int(abs(x(i,l)-rcmax)/bin_width)+1 !different system has different bindex.
        !print*,bindex
        bins(i,bindex)=bins(i,bindex)+1.d0
      end if
    end do
    !print*,bins(i,:)
    Ni(i)=npoint
    !print*,Ni(i)
    close(12)
  end do 

  allocate(rc(binsnum))
  allocate(weight_factor(nwindows,binsnum))
  allocate(bias_factor(nwindows,binsnum))
  weight_factor=0.d0
  bias_factor=0.d0
  do i=1,nwindows
    !print*,x_local(i)
    do k=1,binsnum
      rc(k)=(rcmax-bin_width/2.d0)-(k-1)*bin_width
      !print*,rc(k)
      weight_factor(i,k)=1.d0/2.d0*kw(i)*(rc(k)-x_local(i))**2
      bias_factor(i,k)=exp(-beta*weight_factor(i,k))
    end do
    !print*,weight_factor(i,:)
    !print*,bias_factor(i,:)
  end do

  allocate(f_i_old(nwindows))
  allocate(f_i_old_kcal(nwindows))
  allocate(sum_bins(binsnum))
  allocate(sum_bias_factor(binsnum))
  allocate(p_l(binsnum))
!  f_i_old=1.d0 !initial guess
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!f_i_old have been obtained by MBAR                                       !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(21,file='MBAR-fi.dat')  
  read(21,*) 
  do i=1,nwindows
    read(21,'(I2,2f13.4)') numbw,f_i_old(i),f_i_old_kcal(i)
    f_i_old(i)=exp(beta*f_i_old(i)) !mbar: e^{-beta*f_{i}} = wham: f_{i}^{-1}
  end do
  close(21)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do k=1,binsnum
    sum_bins(k)=0.d0
    sum_bias_factor(k)=0.d0
    do i=1,nwindows
      sum_bins(k)=bins(i,k)+sum_bins(k)
      sum_bias_factor(k)=Ni(i)*f_i_old(i)*bias_factor(i,k)+sum_bias_factor(k)
    end do
    p_l(k)=real(sum_bins(k))/sum_bias_factor(k)
  end do
  
  open(20,file='MBAR-pmf.dat')
  allocate(pmf(binsnum))
  do k=1,binsnum
    pmf(k)=(-1.d0/beta)*log(p_l(k)/(p_l(binsnum)))
    write(20,'(2f13.4)') rc(k),pmf(k)
  end do
  close(20)

  end subroutine
end module
