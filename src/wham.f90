module wham
use parameter
implicit none
contains
  subroutine wham2
  implicit none
  integer(kind=4)                          ::i,k,l
  integer(kind=4)                          ::npoint,bindex
  integer(kind=4)                          ::numb
  integer(kind=4),allocatable              ::Ni(:),bins(:,:),sum_bins(:)
  real(kind=8)                             ::bin_width
  real(kind=8),allocatable                 ::x(:,:)
  real(kind=8),allocatable                 ::x_local(:),kw(:)
  real(kind=8),allocatable                 ::rc(:)
  real(kind=8),allocatable                 ::weight_factor(:,:),bias_factor(:,:),sum_bias_factor(:)
  real(kind=8),allocatable                 ::f_i(:),f_i_old(:),f_i_i(:),deltaf_i(:),p_l(:),pmf(:)
  character(len=20),allocatable            ::filename(:)
  character(len=80)                        ::char_i

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
    write(char_i,'(I5)') i
    open(15,file='./bindexs/bindexs'//trim(adjustl(char_i))//'.dat')
    do l=1,nframes
      read(12,'(I8,f10.4)') numb,x(i,l)    
      if (x(i,l)>rcmin.and.x(i,l)<=rcmax) then
        npoint=npoint+1
        bindex=int(abs(x(i,l)-rcmax)/bin_width)+1 !different system has different bindex.
        !print*,bindex
        write(15,'(2I5)') l,bindex
        bins(i,bindex)=bins(i,bindex)+1.d0
      elseif(x(i,l)<=rcmin.or.x(i,l)>rcmax) then
        bindex=0
        write(15,'(2I5)') l,bindex  
      end if
    end do
    !print*,bins(i,:)
    Ni(i)=npoint
    !print*,Ni(i)
    close(12)
    close(15)
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

  allocate(f_i(nwindows))
  allocate(f_i_old(nwindows))
  allocate(deltaf_i(nwindows))
  allocate(f_i_i(nwindows))
  allocate(sum_bins(binsnum))
  allocate(sum_bias_factor(binsnum))
  allocate(p_l(binsnum))
  f_i=0.d0
  f_i_old=1.d0 !initial guess
  deltaf_i=0.d0
  do while(.true.) 
    do k=1,binsnum
      sum_bins(k)=0.d0
      sum_bias_factor(k)=0.d0
      do i=1,nwindows
        sum_bins(k)=bins(i,k)+sum_bins(k)
        sum_bias_factor(k)=Ni(i)*f_i_old(i)*bias_factor(i,k)+sum_bias_factor(k)
      end do
      p_l(k)=real(sum_bins(k))/sum_bias_factor(k)
      !print*,real(sum_bins(k)),sum_bias_factor(k),p_l(k)
    end do
    do i=1,nwindows
      f_i_i(i)=0.d0
      do k=1,binsnum
        f_i_i(i)=p_l(k)*bias_factor(i,k)+f_i_i(i)
      end do
      f_i(i)=1.d0/(f_i_i(i))
    end do
    do i=1,nwindows
      deltaf_i(i)=abs(f_i(i)-f_i_old(i))/abs(f_i_old(i))
    end do
    print*,"wham's convergence tolerance:",abs(maxval(deltaf_i))
   
    if (abs(maxval(deltaf_i))<=tolerance) then
      exit
    else
      f_i_old=f_i
    end if
  end do
  
  open(19,file='WHAM-fi.dat')
  write(19,'(a)') "#windows,     in kcal/mol,    in reduced units"
  do i=1,nwindows
    write(19,'(I2,2f13.4)') i,(1.d0/beta)*log(f_i(i)/f_i(1)),(1.d0/beta)*log(f_i(i)/f_i(1))*beta
  end do
  close(19)

  open(20,file='WHAM-pmf.dat')
  allocate(pmf(binsnum))
  do k=1,binsnum
    pmf(k)=(-1.d0/beta)*log(p_l(k)/(p_l(binsnum)))
    write(20,'(2f13.4)') rc(k),pmf(k)
  end do
  close(20)

  end subroutine
end module
