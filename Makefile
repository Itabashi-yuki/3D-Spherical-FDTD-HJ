OBJS = fdtd3d.o update_E.o update_H.o current_source.o cal_obs_n0.o allocate.o output_E.o output_pal.o make_dir.o \
		update_E_PML.o update_H_PML.o initialize_PML.o
HEADERS = fdtd3d.h

OPTS = -O3 -Wall -std=gnu++17

all: main

main: $(OBJS)
	g++ -o $@ $(OBJS)

%.o: %.cpp $(HEADERS)
	g++ -c $< $(OPTS)

clean:
	rm -rf main *.o

.PHONY: all clean