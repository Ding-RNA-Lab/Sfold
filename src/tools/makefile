CC	= cc
FC	= f90
CFLAGS	= -fast -m64
CLIBS = -lmopt
FFLAGS	= -e
#Project has been forked from main Sfold. Use local .f files.
FDIR	= .

disruptEn: disruptEn.o findfe.sub.o makefile
	$(FC) -o disruptEn disruptEn.o  findfe.sub.o $(CLIBS) $(CFLAGS) 

findfe.sub.o: $(FDIR)/findfe.sub.f 	mpnum.o data.o energy.o makefile
	$(FC) -c $(FDIR)/findfe.sub.f $(CFLAGS) $(FFLAGS)
mpnum.o: $(FDIR)/mpnum.f makefile
	$(FC) -c $(FDIR)/mpnum.f $(CFLAGS) $(FFLAGS)
data.o: $(FDIR)/data.f makefile
	$(FC) -c $(FDIR)/data.f $(CFLAGS) $(FFLAGS)
energy.o: $(FDIR)/energy.f makefile
	$(FC) -c $(FDIR)/energy.f $(CFLAGS) $(FFLAGS)
	
disruptEn.o: disruptEn.c makefile
	$(CC) -c disruptEn.c $(CFLAGS)

clean:
	rm -f *.o disruptEn

