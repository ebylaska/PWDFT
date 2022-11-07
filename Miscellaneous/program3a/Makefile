OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
        OPTIONS += -framework OpenCL
else
        OPTIONS += -l OpenCL -lm
endif

main: test3a.c
	gcc -Wall -g test3a.c -o test3a $(OPTIONS)

clean:
	rm -rf test3a
