# Makefile

FC = gfortran

OBJS = main.o initialise.o extrap_and_print.o firstord_radextrap.o utils.o

radextrap: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o radextrap

