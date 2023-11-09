module parameter
implicit none
integer(kind=4),parameter          ::nwindows=31
integer(kind=4),parameter          ::nframes=10000
real(kind=8),parameter             ::kb=0.0019872041d0 ! Boltzmann's constant (kcal/mol/K)
real(kind=8),parameter             ::T=300.d0
real(kind=8),parameter             ::beta=1.d0/(kb*T)
real(kind=8),parameter             ::J_PER_CAL=4.184d0
real(kind=8),parameter             ::tolerance=1.0e-7 !Termination criterion for iteration calculations
real(kind=8),parameter             ::rcmin=-1.55d0
real(kind=8),parameter             ::rcmax=1.55d0
integer(kind=4),parameter          ::binsnum=60
end module
