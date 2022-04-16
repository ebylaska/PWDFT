OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
        OPTIONS += -framework OpenCL
else
        OPTIONS += -l OpenCL -lm
endif

main: test3.c
	gcc -Wall -g test3.c -o test3 $(OPTIONS)

clean:
	rm -rf test3
