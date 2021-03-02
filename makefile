TARGET = kashikoma

SRC_DIR = src
SOURCES = main.c pre_process.c post_process.c solver.c

OBJECTS := $(addprefix $(SRC_DIR)/,$(SOURCES:.c=.o))

CC = gcc
#CC = icc
#CC = mpicc

CFLAGS = -g
#CFLAGS = -O3 -fopenmp

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) $^ -o $(TARGET)
	@echo "make completed!! Please execute!!"

.SUFFIXES : .o .c
.c.o:
	${CC} $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET) $(OBJECTS)

# header dependencies
src/main.o : src/func_header.h
src/pre_process.o : src/func_header.h
src/post_process.o : src/func_header.h
src/solver.o : src/func_header.h

