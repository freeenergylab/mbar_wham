#makefile

FC:=gfortran

OBJ:=main.o parameter.o extract.o mbar.o wham.o mbarpmf.o weight.o

%.o:%.f90
	$(FC) -c $<

all:main.x

main.x:$(OBJ)
	$(FC) -o $@ $^

main.o:extract.o mbar.o wham.o mbarpmf.o weight.o
extract.o:parameter.o
wham.o:parameter.o
mbar.o:parameter.o
mbarpmf.o:parameter.o
weight.o:parameter.o

clean:
	rm *.o *.mod 
