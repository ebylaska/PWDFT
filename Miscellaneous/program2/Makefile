OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
        OPTIONS += -framework OpenCL
else
        OPTIONS += -l OpenCL -lm
endif

main: test2.c
	gcc -Wall -g test2.c -o test2 $(OPTIONS)

clean:
	rm -rf test1
