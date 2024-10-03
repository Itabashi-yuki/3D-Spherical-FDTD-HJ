PREFIX = $(HOME)
LIB_DIR = $(PREFIX)/lib
INC_DIR = $(PREFIX)/include

LIB_SUFFIX = _20	

OBJS = fdtd3d.o update_E.o update_H.o current_source.o cal_obs_n0.o allocate.o output_E.o output_pal.o make_dir.o \
		update_E_PML.o update_D_PML.o update_H_PML.o update_J.o initialize_PML.o initialize_Plasma.o Geomag_IGRF.o
HEADERS = fdtd3d.h

OPTS = -O3 -Wall -std=gnu++17

all: main

main: $(OBJS)
	g++ -o $@ $(OBJS) -fopenmp -Wl,-R$(LIB_DIR) -L$(LIB_DIR) -ligrf$(LIB_SUFFIX)

%.o: %.cpp $(HEADERS)
	g++ -c $< $(OPTS) -fopenmp -I$(INC_DIR)

clean:
	rm -rf main *.o

.PHONY: all clean