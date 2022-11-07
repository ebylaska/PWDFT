OS := $(shell uname)
OPTIONS:= 

ifeq ($(OS),Darwin)
        OPTIONS += -framework OpenCL
else
        OPTIONS += -l OpenCL -lm
endif

main: test5.cpp
	g++ -Wall -g test5.cpp -o test5 $(OPTIONS)

clean:
	rm -rf test5
