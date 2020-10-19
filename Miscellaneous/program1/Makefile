OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
        OPTIONS += -framework OpenCL
else
        OPTIONS += -l OpenCL
endif

main: test1.c
	gcc -Wall -g test1.c -o test1 $(OPTIONS)

clean:
	rm -rf test1
