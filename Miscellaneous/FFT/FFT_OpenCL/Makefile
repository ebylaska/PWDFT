# Directories to compile
#DIR_1D = $(pwd)/1D
#DIR_1D = $(pwd)/2D
#DIR_1D = $(pwd)/3D

DIR_1D = $(dir $(mkfile_dir))1D
DIR_2D = $(dir $(mkfile_dir))2D
DIR_3D = $(dir $(mkfile_dir))3D

#DIR_2D := $(dir $(realpath $(firstword $(MAKEFILE_LIST))))/2D/
#DIR_3D := $(dir $(realpath $(firstword $(MAKEFILE_LIST))))/3D/

# Phony target
.PHONY: all

all:
	@(cd $(DIR_1D) ; $(MAKE))
	@(cd $(DIR_2D) ; $(MAKE))
	@(cd $(DIR_3D) ; $(MAKE))

# Clean target
clean:
	@(cd $(DIR_1D) ; $(MAKE) $@)
	@(cd $(DIR_2D) ; $(MAKE) $@)
	@(cd $(DIR_3D) ; $(MAKE) $@)
