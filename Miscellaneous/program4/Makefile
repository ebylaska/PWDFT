OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
        OPTIONS += -framework OpenCL
else
        OPTIONS += -l OpenCL
endif

main: test4.c
	gcc -Wall -g test4.c -o test4 $(OPTIONS)

clean:
	rm -rf test4
