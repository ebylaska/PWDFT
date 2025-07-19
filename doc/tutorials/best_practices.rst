Best Practices Guide
===================

This guide provides system-specific recommendations for setting up reliable PWDFT calculations. These practices are based on extensive testing and experience with different types of materials and systems.

Gas-Phase Molecule Calculations
-------------------------------

**Supercell Construction:**

For isolated molecules, place them in a large periodic box to simulate gas-phase conditions:

* **Minimum vacuum spacing**: 10-15 Å between molecule and box edge in all directions
* **Box size**: Typically 20-30 Å cubic for small molecules
* **Rationale**: Prevents spurious interactions between periodic images

**K-Point Sampling:**

* **Use Gamma point only**: `monkhorst-pack = 1 1 1`
* **Rationale**: Electronic interactions in reciprocal space are negligible for isolated molecules

**Recommended Functionals:**

* **Standard GGA**: PBE for initial screening
* **Hybrid functionals**: PBE0, B3LYP for higher accuracy
* **Rationale**: Hybrids reduce self-interaction error important for molecular properties

**Convergence Protocol:**

* **Primary parameter**: Plane-wave energy cutoff (`ecutwfc`)
* **Typical range**: 40-80 Rydberg
* **Convergence criterion**: Energy change < 1 meV/atom

**Example Input:**

.. code-block:: text

   &nwpw
     xc = 'pbe96'
     cutoff = 60.0
     scf = 'ks-grassmann-cg simple alpha 0.02'
     loop = '50 1'
     monkhorst-pack = '1 1 1'
     initial_wavefunction_guess = 'superposition'
   /

Bulk Crystal Calculations
------------------------

**System Classification:**

**Insulators/Semiconductors:**
* **Smearing**: Disabled (`use_smearing = false`)
* **Rationale**: Distinct band gap allows integer occupations

**Metals:**
* **Smearing**: Enabled with Methfessel-Paxton scheme
* **Temperature**: 0.01-0.02 eV (typically 5000-8000 K)
* **Rationale**: Fractional occupations at Fermi level require smearing

**Convergence Protocol:**

**Two-Step Process:**

1. **K-Point Convergence** (fixed `ecutwfc`):
   * Start with 2×2×2 grid
   * Increase to 4×4×4, 6×6×6, etc.
   * Target: Energy convergence < 1-5 meV/atom

2. **Energy Cutoff Convergence** (fixed k-points):
   * Start with 40 Rydberg
   * Increase in 10 Rydberg steps
   * Target: Energy convergence < 1 meV/atom

**Recommended Functionals:**

* **PBE**: Excellent starting point for most bulk solids
* **Hybrids**: HSE06 for accurate band gaps
* **vdW**: optB88-vdW for layered materials

**Example Input (Metal):**

.. code-block:: text

   &nwpw
     xc = 'pbe96'
     cutoff = 50.0
     scf = 'ks-grassmann-cg simple alpha 0.02'
     smear = 'methfessel-paxton'
     temperature = 5000
     loop = '50 1'
     monkhorst-pack = '6 6 6'
     initial_wavefunction_guess = 'superposition'
   /

Slab and Adsorbate Calculations
------------------------------

**Slab Model Construction:**

**Slab Thickness:**
* **Minimum**: 4-6 atomic layers
* **Bottom layers**: Fix 1-2 layers in bulk positions
* **Rationale**: Ensures bulk-like properties in center

**Vacuum Spacing:**
* **Minimum**: 15-20 Å in non-periodic direction
* **Rationale**: Prevents slab-slab interactions

**Dipole Correction:**
* **When needed**: Asymmetric slabs or adsorbates on one side
* **Implementation**: Apply correction in non-periodic direction
* **Rationale**: Cancels artificial electric field

**Adsorption Energy Calculation:**

.. math::

   E_{ads} = E_{slab+adsorbate} - (E_{slab} + E_{adsorbate})

**Critical: van der Waals Forces**

* **Default choice**: vdW-inclusive functional
* **Recommended**: optB88-vdW or PBE+D3
* **Rationale**: Standard GGAs miss long-range dispersion forces

**K-Point Sampling:**

* **Periodic directions**: Dense grid (e.g., 8×8×1)
* **Non-periodic direction**: Single point
* **Rationale**: No periodicity in z-direction

**Example Input:**

.. code-block:: text

   &nwpw
     xc = 'optb88-vdw'
     cutoff = 50.0
     scf = 'ks-grassmann-cg simple alpha 0.02'
     smear = 'methfessel-paxton'
     temperature = 5000
     loop = '50 1'
     monkhorst-pack = '8 8 1'
     dipole_correction = .true.
     initial_wavefunction_guess = 'superposition'
   /

General Guidelines
-----------------

**SCF Convergence:**

* **Mixing parameter**: Start with 0.02 (conservative)
* **Maximum iterations**: 50-100 for complex systems
* **Convergence threshold**: 1e-6 Hartree/atom

**Wavefunction Initialization:**

* **Default**: `superposition` (atomic superposition)
* **Fallback**: `random` if convergence fails
* **Metals**: `superposition` often works better than `random`

**Performance Optimization:**

* **Parallelization**: Use MPI for large systems
* **Memory**: Monitor memory usage for large calculations
* **I/O**: Use scratch directories for temporary files

**Troubleshooting:**

* **SCF divergence**: Reduce mixing parameter, increase smearing
* **NaN errors**: Check pseudopotentials, reduce cutoff
* **Memory issues**: Reduce parallelization, use smaller k-point grids

**Validation:**

* **Test calculations**: Compare with known results
* **Convergence studies**: Always perform systematic convergence
* **Physical checks**: Verify forces, energies make sense 