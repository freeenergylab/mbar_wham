program main
use wham
use extract
use mbar
use mbarpmf
use weight
implicit none

call system('mkdir -p bindexs')
call wham2
call system('mkdir -p energy')
call extr
call mbar2
call mbarpmf2
call system('mkdir -p weights')
call weight2

end
