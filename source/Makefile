SRCS=	random.c main.c draw.c bmd.c samplemd.c sampletis.c movepath.c init.c

OBJS=	random.o main.o draw.o bmd.o samplemd.o sampletis.o movepath.o init.o


HFILES = 

#

CF =   cc
FFLAGS =   
LIBS =   -lm -framework OpenGL -framework Glut -framework Cocoa
OUT =	../run/bmd.run



opt:	$(SRCS:.c=.o)
	$(CF) $(FFLAGS) -o $(OUT) $(SRCS:.c=.o) $(LIBS)

prof:	$(SRCS:.c=.o)
	$(CF) $(FFLAGS) -pg -o $(OUT) $(SRCS) $(LIBS)

opt3:	$(SRCS)
	$(CF)  -O3 -o $(OUT) $(SRCS) $(LIBS)

deb:	$(SRCS)
	$(CF) -g -O2 -o rods.deb $(SRCS) $(LIBS)

$(SRCS):	$(HFILES)
$(OBJS):	$(HFILES)	







