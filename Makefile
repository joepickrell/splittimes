HOMEL=$(PWD)
TOP=bin
## binary for install  
BIN=$(HOMEL)/bin
DEBUG_OPTIONS=

NLIB=$(HOMEL)/lib/nicklib.a 
IDIR=$(HOMEL)/include/
INC=-I/opt/local/include
LIB=-L/opt/local/lib

CFLAGS= -c -g -I$(IDIR) -Wimplicit #
TWTAB=\"$(HOMEL)/tables/twtable\"
VPATH=.:nicksrc

PROGS= kimtrans gaurav test 
OBJ=strsubs.o sortit.o vsubs.o statsubs.o linsubs.o getpars.o xsearch.o gauss.o	gds.o


N1=kimtrans 
N1O=kimtrans.o  kimsubs.o wynn.o

N2=gaurav 
N2O=gaurav.o  kimsubs.o wynn.o

N3=test
N30=test.o Kimura.o kimsubs.o wynn.o

N4=splittimes
N40=SplitTimes.o gzstream.o CountData.o CmdLine.o Kimura.o kimsubs.o wynn.o

statsubs.o:     nicksrc/statsubs.c
	$(CC)  $(CFLAGS) -DTWTAB=$(TWTAB) -o statsubs.o nicksrc/statsubs.c

$(N1): nicklib $(N1O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(N1) $(N1O) $(NLIB) 
 
$(N2): nicklib $(N2O)
	gcc -I$(IDIR) $(DEBUG_OPTIONS) -lm  -o $(N2) $(N2O) $(NLIB) 

$(N3): nicklib $(N30)
	g++ -I$(IDIR) $(INC) $(LIB) -lm -lgsl  -o $(N3) $(N30) $(NLIB)

$(N4): nicklib $(N40)
	g++ -I$(IDIR) $(INC) $(LIB) -lz -lm -lgsl  -o $(N4) $(N40) $(NLIB)

libnick.a:	dirs tables  $(OBJ)
	ar -r libnick.a $(OBJ)
	ranlib libnick.a

nicklib:	libnick.a 
	cp libnick.a  $(NLIB)

tables:    
	echo "tables made"  > tables
	cp twtable  $(HOMEL)/tables
	
dirs:	
	mkdir -p  $(HOMEL)/lib
	mkdir -p  $(BIN)
	cp  *.h  $(IDIR).
	cp  nicksrc/*.h  $(IDIR)

%.o: %.cpp
	g++ -c $< -o $@ $(INC)

clean: 
	rm -f *.o 
	rm -f core

clobber: clean 
	rm -f $(PROGS)

