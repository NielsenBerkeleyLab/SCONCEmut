# compiler
CC=gcc
# compiler flags
CFLAGS=-lm


all: SCONCEmut_sim

SCONCEmut_sim: SCONCEmut_sim.c
	$(CC) SCONCEmut_sim.c ranlib.c rnglib.c -o SCONCEmut_sim $(CFLAGS)

.PHONY: clean
clean:
	rm -f SCONCEmut_sim true_*healthy_cell* true_*cancer_cell* simu_*healthy_cell* simu_*cancer_cell*

