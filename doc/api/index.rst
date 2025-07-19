API Reference
=============

.. toctree::
   :maxdepth: 2

   keywords
   functions
   modules

Input Keywords Reference
------------------------

This section provides comprehensive documentation of all PWDFT input keywords, automatically generated from the source code comments using Doxygen and Breathe.

**Keyword Categories:**

* **Control Parameters** - Calculation type, convergence, and output options
* **System Definition** - Atomic structure, cell parameters, and pseudopotentials
* **Electronic Structure** - SCF parameters, mixing, and smearing
* **Advanced Features** - DFT+U, spin-orbit coupling, and more

**Keyword Format:**

Each keyword entry includes:

* **Description**: What the parameter controls
* **Data Type**: Integer, float, string, logical
* **Default Value**: Standard setting if not specified
* **Valid Range**: Acceptable values and constraints
* **Physical Meaning**: Connection to underlying physics
* **Examples**: Sample input files and usage

**Auto-Generated Content:**

The keyword documentation is automatically extracted from source code comments using the following format:

.. code-block:: cpp

   /**
    * @brief Sets the kinetic energy cutoff for plane-wave basis
    * @param cutoff Energy cutoff in Rydberg
    * @details The cutoff determines the maximum kinetic energy of 
    *          plane waves used to expand the electronic wavefunctions.
    *          Higher values increase accuracy but computational cost.
    * @default 40.0 Ry
    * @range 20.0 - 200.0 Ry
    */
   void set_cutoff(double cutoff);

**Function Reference:**

Complete API documentation for all public functions, including:

* **Parameter descriptions**
* **Return values**
* **Error conditions**
* **Usage examples**

**Module Documentation:**

Detailed documentation of major code modules:

* **Core SCF routines**
* **Pseudopotential handling**
* **Parallelization utilities**
* **I/O and file formats**

**Search and Navigation:**

* **Full-text search** across all documentation
* **Cross-references** between related keywords
* **Index** of all functions and keywords
* **Call graphs** showing function relationships

**Contributing:**

To improve the API documentation:

1. **Add Doxygen comments** to source code functions
2. **Use consistent formatting** for parameter descriptions
3. **Include examples** in comments where helpful
4. **Update this documentation** when adding new features

**Building the Documentation:**

.. code-block:: bash

   # Generate Doxygen XML
   cd doc
   doxygen Doxyfile
   
   # Build Sphinx documentation
   make html
   make latexpdf

The generated documentation will be available in `doc/_build/html/` and `doc/_build/latex/`. 