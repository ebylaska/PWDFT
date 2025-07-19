# PWDFT Documentation

This directory contains the complete documentation for PWDFT, built using Sphinx + Doxygen + Breathe.

## Quick Start

### Prerequisites

1. **Python packages:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Doxygen:**
   - macOS: `brew install doxygen`
   - Ubuntu: `sudo apt-get install doxygen`
   - CentOS: `sudo yum install doxygen`

### Building Documentation

1. **Check dependencies:**
   ```bash
   make check-deps
   ```

2. **Build HTML documentation:**
   ```bash
   make html
   ```

3. **Build PDF documentation:**
   ```bash
   make pdf
   ```

4. **View the documentation:**
   - HTML: Open `_build/html/index.html` in your browser
   - PDF: Find `_build/latex/PWDFT.pdf`

## Documentation Structure

```
doc/
├── conf.py              # Sphinx configuration
├── Doxyfile            # Doxygen configuration
├── index.rst           # Main documentation index
├── Makefile            # Build automation
├── requirements.txt    # Python dependencies
├── introduction/       # Project overview and theory
├── installation/       # Installation guides
├── tutorials/          # Step-by-step tutorials
├── theory/            # Theoretical background
├── api/               # Auto-generated API reference
├── faq/               # Frequently asked questions
└── developers/        # Developer documentation
```

## Contributing to Documentation

### Adding New Content

1. **Create new RST files** in the appropriate directory
2. **Add to toctree** in the parent `index.rst`
3. **Build and test** your changes

### Improving API Documentation

1. **Add Doxygen comments** to source code functions:
   ```cpp
   /**
    * @brief Brief description
    * @param param Parameter description
    * @return Return value description
    * @details Detailed explanation
    */
   ```

2. **Rebuild documentation** to see changes:
   ```bash
   make clean html
   ```

### Style Guidelines

- **Use RST syntax** for formatting
- **Include code examples** where helpful
- **Cross-reference** related sections using `:doc:` and `:ref:`
- **Add math equations** using LaTeX syntax
- **Use admonitions** for important notes and warnings

## Development Workflow

### Local Development

1. **Start development server:**
   ```bash
   make serve
   ```

2. **Edit files** and see changes automatically

3. **Build for production:**
   ```bash
   make all
   ```

### Continuous Integration

The documentation is automatically built and deployed when changes are pushed to the main branch.

## Troubleshooting

### Common Issues

1. **Doxygen not found:**
   - Install Doxygen using your package manager
   - Ensure it's in your PATH

2. **Breathe import error:**
   - Install breathe: `pip install breathe`
   - Check Python environment

3. **Build errors:**
   - Clean and rebuild: `make clean all`
   - Check RST syntax
   - Verify file paths

### Getting Help

- Check the Sphinx documentation: https://www.sphinx-doc.org/
- Review Doxygen manual: https://www.doxygen.nl/
- Browse Breathe documentation: https://breathe.readthedocs.io/

## License

This documentation is licensed under the same terms as PWDFT. 