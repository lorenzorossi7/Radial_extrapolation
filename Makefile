# Makefile

#FC = gfortran
FC = mpif90

OBJS = main.o initialise.o extrap_and_print.o firstord_radextrap.o utils.o

EXEC = radextrap

$(EXEC): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

%.o: %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJS) $(EXEC)
