FORT    =  gfortran
#NICE    =  -fdefault-real-16 -fdefault-double-8
NICE    =  -fdefault-real-8 -fdefault-double-8 
LIB     =   /users/spencerbryngelson/downloads/lapack/liblapack.a /users/spencerbryngelson/downloads/lapack/librefblas.a
OPTS    =  -O3 
#DEBUG   = -O2 -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace
PROF    = 

OBJECTS =  cheb.o prms.o main.o data.o out.o ev.o

film: $(OBJECTS) makefile 
	$(FORT) $(LIB) -o film $(PROF) $(OBJECTS) $(LIB)

main.o: main.f90 cheb.o prms.o out.o ev.o makefile
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  main.f90 

ev.o: ev.f90 data.o prms.o cheb.o makefile
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  ev.f90 

out.o: out.f90 data.o prms.o cheb.o makefile
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  out.f90 

data.o: data.f90 prms.o cheb.o makefile
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  data.f90 

cheb.o: cheb.f90 prms.o makefile
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  cheb.f90

prms.o: prms.f90 makefile
	$(FORT) -c $(NICE) $(OPTS) $(PROF) $(DEBUG)  prms.f90

clean:
	rm -f *.o
#	rm le -f

clear:
	rm *.out* -f


