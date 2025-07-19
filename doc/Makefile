# Makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Custom targets for PWDFT documentation

# Generate Doxygen XML
doxygen:
	doxygen Doxyfile

# Build complete documentation (Doxygen + Sphinx)
all: doxygen html

# Build HTML documentation
html: doxygen
	$(SPHINXBUILD) -b html $(SPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/html

# Build PDF documentation
pdf: doxygen
	$(SPHINXBUILD) -b latex $(SPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/latex
	$(MAKE) -C $(BUILDDIR)/latex all-pdf

# Clean all build artifacts
clean:
	rm -rf $(BUILDDIR)
	rm -rf ../build/doxygen

# Install required Python packages
install-deps:
	pip install sphinx sphinx-rtd-theme breathe

# Check if required tools are available
check-deps:
	@echo "Checking documentation dependencies..."
	@which sphinx-build > /dev/null || (echo "sphinx-build not found. Install with: pip install sphinx" && exit 1)
	@which doxygen > /dev/null || (echo "doxygen not found. Install with: brew install doxygen (macOS) or apt-get install doxygen (Ubuntu)" && exit 1)
	@python -c "import breathe" 2>/dev/null || (echo "breathe not found. Install with: pip install breathe" && exit 1)
	@echo "All dependencies found!"

# Development server (requires sphinx-autobuild)
serve:
	sphinx-autobuild $(SOURCEDIR) $(BUILDDIR)/html $(SPHINXOPTS) 