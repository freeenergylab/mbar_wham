#makefile

FC:=gfortran

OBJ:=main.o parameter.o extract.o wham.o

%.o:%.f90
	$(FC) -c $<

all:main.x

main.x:$(OBJ)
	$(FC) -o $@ $^

main.o:extract.o wham.o
extract.o:parameter.o
wham.o:parameter.o

clean:
	rm *.o *.mod 
