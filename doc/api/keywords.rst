Input Keywords Reference
========================

This section provides a comprehensive analysis of PWDFT input keywords, their functions, data types, defaults, and underlying physics.

Keyword Analysis Table
---------------------

.. list-table:: PWDFT Input Keywords
   :widths: 20 25 15 15 25
   :header-rows: 1

   * - Keyword/Concept
     - Function/Purpose
     - Data Type & Example
     - Default Value & Rationale
     - Associated Physics

   * - calculation
     - Specifies calculation type
     - String: 'scf', 'relax', 'vc-relax'
     - 'scf'. Single-point SCF is fundamental
     - Kohn-Sham DFT, Structural Optimization

   * - ecutwfc
     - Kinetic energy cutoff for wavefunctions
     - Float: 40.0, 60.0 (Rydberg)
     - 40.0 Ry. Balance of cost and accuracy
     - Plane-Wave Basis Set, Fourier Expansion

   * - ecutrho
     - Kinetic energy cutoff for charge density
     - Float: 160.0, 240.0 (Rydberg)
     - 4 * ecutwfc. Higher frequency components
     - Charge Density Representation, Poisson's Equation

   * - xc
     - Exchange-correlation functional
     - String: 'pbe96', 'hse06', 'beef-vdw'
     - 'pbe96'. Robust GGA for solids
     - Exchange-Correlation Energy, Jacob's Ladder of DFT

   * - monkhorst-pack
     - K-point grid dimensions
     - Integer Array: [1,1,1], [8,8,8]
     - [1,1,1]. Gamma point for molecules
     - k-point Sampling, Brillouin Zone Integration

   * - smear
     - Occupation number smearing
     - String: 'methfessel-paxton', 'gaussian'
     - None. Only needed for metals
     - Fermi-Dirac Statistics, Electronic Temperature

   * - temperature
     - Smearing width (electronic temperature)
     - Float: 5000, 8000 (Kelvin)
     - 5000 K. Stabilize convergence in metals
     - Fermi-Dirac Distribution

   * - scf
     - SCF algorithm and mixing method
     - String: 'ks-grassmann-cg simple alpha 0.02'
     - 'ks-grassmann-cg pulay alpha 0.1'
     - Self-Consistent Field (SCF) Convergence

   * - loop
     - Maximum SCF iterations
     - String: '50 1', '20 20'
     - '50 1'. Allow sufficient iterations
     - SCF Convergence Criteria

   * - initial_wavefunction_guess
     - Initial wavefunction strategy
     - String: 'random', 'superposition', 'atomic'
     - 'superposition'. Physically reasonable
     - Wavefunction Initialization

   * - dipole_correction
     - Apply dipole correction for slabs
     - Logical: .true., .false.
     - .false. Only needed for asymmetric slabs
     - Surface Dipole, Work Function

   * - nspin
     - Number of spin components
     - Integer: 1, 2
     - 1. Non-spin-polarized default
     - Electron Spin, Magnetism

   * - pseudo_dir
     - Pseudopotential directory
     - String: './pseudos/'
     - '.'. Current directory default
     - Pseudopotential Approximation

Detailed Keyword Descriptions
----------------------------

Calculation Type Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^

.. _keyword-calculation:

calculation
~~~~~~~~~~~

**Purpose**: Controls the type of calculation to be performed.

**Values**:
- ``scf``: Single-point Self-Consistent Field calculation
- ``relax``: Structural relaxation with fixed cell
- ``vc-relax``: Variable-cell relaxation (cell + ions)

**Default**: ``scf``

**Physics**: The calculation type determines which degrees of freedom are optimized during the calculation.

Electronic Structure Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _keyword-ecutwfc:

ecutwfc
~~~~~~~

**Purpose**: Sets the kinetic energy cutoff for the plane-wave basis set.

**Range**: 20.0 - 200.0 Rydberg

**Default**: 40.0 Rydberg

**Physics**: The cutoff determines the maximum kinetic energy of plane waves used to expand electronic wavefunctions. Higher values increase accuracy but computational cost scales as :math:`E_{cut}^{3/2}`.

**Convergence**: Must be systematically converged for production calculations.

.. _keyword-xc:

xc
~~

**Purpose**: Selects the exchange-correlation functional.

**Available Functionals**:

* **LDA**: ``slater``, ``vosko``
* **GGA**: ``pbe96``, ``pbesol``, ``revpbe``
* **Hybrid**: ``hse06``, ``pbe0``
* **vdW**: ``optb88-vdw``, ``beef-vdw``

**Default**: ``pbe96``

**Physics**: The XC functional approximates the complex many-body electron interactions. Different rungs of "Jacob's Ladder" provide increasing accuracy at higher computational cost.

SCF Convergence Keywords
^^^^^^^^^^^^^^^^^^^^^^^

.. _keyword-scf:

scf
~~~

**Purpose**: Controls the SCF algorithm and mixing parameters.

**Format**: ``algorithm mixer alpha value``

**Algorithms**:
- ``ks-grassmann-cg``: Conjugate gradient on Grassmann manifold
- ``davidson``: Davidson diagonalization

**Mixers**:
- ``simple``: Simple mixing (most stable)
- ``pulay``: Pulay DIIS (most efficient)
- ``anderson``: Anderson mixing

**Default**: ``ks-grassmann-cg pulay alpha 0.1``

**Physics**: The SCF procedure iteratively solves the Kohn-Sham equations until self-consistency is achieved.

.. _keyword-initial_wavefunction_guess:

initial_wavefunction_guess
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Purpose**: Specifies the initial wavefunction strategy.

**Options**:
- ``random``: Random initialization
- ``superposition``: Atomic orbital superposition
- ``atomic``: Atomic wavefunctions
- ``gaussian``: Gaussian wavepackets
- ``mixed``: Combination of methods

**Default**: ``superposition``

**Physics**: The initial guess affects SCF convergence. Physically motivated guesses (superposition) often converge better than random initialization.

K-Point Sampling Keywords
^^^^^^^^^^^^^^^^^^^^^^^^^

.. _keyword-monkhorst-pack:

monkhorst-pack
~~~~~~~~~~~~~~

**Purpose**: Defines the Monkhorst-Pack k-point grid.

**Format**: ``nx ny nz`` (integers)

**Examples**:
- ``1 1 1``: Gamma point only (molecules)
- ``4 4 4``: 64 k-points (bulk crystals)
- ``8 8 1``: 64 k-points in xy-plane (slabs)

**Default**: ``1 1 1``

**Physics**: K-point sampling integrates over the Brillouin zone. Denser grids provide more accurate results but increase computational cost.

Smearing Keywords
^^^^^^^^^^^^^^^^

.. _keyword-smear:

smear
~~~~~

**Purpose**: Enables occupation number smearing for metallic systems.

**Options**:
- ``methfessel-paxton``: MP smearing (recommended)
- ``gaussian``: Gaussian smearing
- ``fermi-dirac``: Fermi-Dirac distribution

**Default**: None (no smearing)

**Physics**: Smearing helps SCF convergence in metals by handling fractional occupations at the Fermi level.

.. _keyword-temperature:

temperature
~~~~~~~~~~~

**Purpose**: Sets the smearing width (effective electronic temperature).

**Range**: 1000 - 10000 Kelvin

**Default**: 5000 Kelvin

**Physics**: The temperature controls the width of the smearing function. Higher values improve convergence but may affect accuracy.

Advanced Features
^^^^^^^^^^^^^^^^

.. _keyword-dipole_correction:

dipole_correction
~~~~~~~~~~~~~~~~~

**Purpose**: Applies dipole correction for asymmetric slabs.

**Values**: ``.true.``, ``.false.``

**Default**: ``.false.``

**Physics**: Asymmetric slabs create artificial electric fields. Dipole correction cancels this field to obtain correct surface energies and work functions.

.. _keyword-nspin:

nspin
~~~~~

**Purpose**: Controls spin polarization.

**Values**:
- ``1``: Non-spin-polarized (default)
- ``2``: Collinear spin-polarized

**Physics**: Spin polarization is needed for magnetic materials and some open-shell systems.

Best Practices
-------------

**Convergence Protocol**:

1. **Energy Cutoff**: Start with 40 Ry, increase until energy converges
2. **K-Points**: Use appropriate grid for system type
3. **SCF**: Start with conservative mixing (alpha=0.02)
4. **Smearing**: Enable for metals, use 5000-8000 K

**System-Specific Recommendations**:

* **Molecules**: Gamma point, 60+ Ry cutoff
* **Bulk Crystals**: Dense k-point grid, 50+ Ry cutoff
* **Slabs**: Dense xy-grid, dipole correction if asymmetric
* **Metals**: Smearing enabled, higher temperature

**Troubleshooting**:

* **SCF divergence**: Reduce mixing parameter, increase smearing
* **Memory issues**: Reduce cutoff or k-point grid
* **Slow convergence**: Try different initial guess or mixing method 