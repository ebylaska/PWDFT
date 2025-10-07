# Symmetry-Aware K-Point Reduction Examples

This document provides examples of how to use the new symmetry-aware k-point reduction system for different crystal systems, including bulk crystals, surfaces, and adsorbates.

## Overview

The new system automatically detects crystal symmetry and reduces k-points to the irreducible Brillouin zone (IBZ) using the appropriate point group operations. This can significantly reduce computational cost while maintaining accuracy.

## Supported Symmetries

### Bulk Crystals
- **Oh**: Cubic symmetry (e.g., FCC, BCC metals)
- **D4h**: Tetragonal symmetry
- **D2h**: Orthorhombic symmetry

### Surfaces
- **C4v**: Square surface (e.g., (100) face of cubic crystal)
- **C3v**: Hexagonal surface (e.g., (111) face of FCC crystal)
- **C2v**: Rectangular surface
- **Cs**: Oblique surface with mirror plane

### Lower Symmetry
- **C2**: 2-fold rotation axis
- **Ci**: Inversion center
- **C1**: No symmetry (default)

## Example 1: FCC Copper Bulk (Oh symmetry)

```nwchem
start Cu_bulk
permanent_dir ./perm
scratch_dir ./scratch

geometry "Cu FCC Bulk" units angstroms
 system crystal
  lat_a 3.615
  lat_b 3.615
  lat_c 3.615
  alpha 90.0
  beta  90.0
  gamma 90.0
 end
 Cu  0.0  0.0  0.0
 Cu  0.5  0.5  0.0
 Cu  0.5  0.0  0.5
 Cu  0.0  0.5  0.5
end

nwpw
 simulation_cell
  fcc 3.615
 end
 cutoff 30.0
 smear fermi 0.01
 monkhorst-pack 8 8 8
end

task pspw energy
```

**Expected Results:**
- Full grid: 512 k-points
- After Oh symmetry reduction: ~24 k-points (reduction factor ~21x)

## Example 2: Cu(111) Surface (C3v symmetry)

```nwchem
start Cu_111_surface
permanent_dir ./perm
scratch_dir ./scratch

geometry "Cu(111) Surface" units angstroms
 system crystal
  lat_a 2.556
  lat_b 2.556
  lat_c 20.0
  alpha 90.0
  beta  90.0
  gamma 120.0
 end
 # Surface layer
 Cu  0.0  0.0  0.0
 Cu  0.5  0.5  0.0
 Cu  1.0  0.0  0.0
 Cu  1.5  0.5  0.0
 # Second layer
 Cu  0.333  0.333  2.089
 Cu  0.833  0.833  2.089
 Cu  1.333  0.333  2.089
 Cu  1.833  0.833  2.089
 # Third layer
 Cu  0.667  0.667  4.178
 Cu  1.167  1.167  4.178
 Cu  1.667  0.667  4.178
 Cu  2.167  1.167  4.178
end

nwpw
 simulation_cell
  lattice_vectors
   2.556  0.0    0.0
   1.278  2.213  0.0
   0.0    0.0    20.0
  end
 cutoff 30.0
 smear fermi 0.01
 monkhorst-pack 3 3 1
end

task pspw energy
```

**Expected Results:**
- Full grid: 9 k-points
- After C3v symmetry reduction: ~3 k-points (reduction factor ~3x)

## Example 3: Cu(100) Surface with Adsorbate (C4v symmetry)

```nwchem
start Cu_100_adsorbate
permanent_dir ./perm
scratch_dir ./scratch

geometry "Cu(100) with CO adsorbate" units angstroms
 system crystal
  lat_a 3.615
  lat_b 3.615
  lat_c 20.0
  alpha 90.0
  beta  90.0
  gamma 90.0
 end
 # Surface Cu atoms (2x2 supercell)
 Cu  0.0  0.0  0.0
 Cu  0.5  0.0  0.0
 Cu  0.0  0.5  0.0
 Cu  0.5  0.5  0.0
 # Second layer
 Cu  0.25  0.25  1.808
 Cu  0.75  0.25  1.808
 Cu  0.25  0.75  1.808
 Cu  0.75  0.75  1.808
 # CO adsorbate
 C   0.25  0.25  2.5
 O   0.25  0.25  3.2
end

nwpw
 simulation_cell
  lattice_vectors
   3.615  0.0    0.0
   0.0    3.615  0.0
   0.0    0.0    20.0
  end
 cutoff 30.0
 smear fermi 0.01
 monkhorst-pack 4 4 1
end

task pspw energy
```

**Expected Results:**
- Full grid: 16 k-points
- After C4v symmetry reduction: ~4 k-points (reduction factor ~4x)

## Example 4: Hexagonal Surface (C6v symmetry)

```nwchem
start hexagonal_surface
permanent_dir ./perm
scratch_dir ./scratch

geometry "Hexagonal Surface" units angstroms
 system crystal
  lat_a 2.5
  lat_b 2.5
  lat_c 20.0
  alpha 90.0
  beta  90.0
  gamma 120.0
 end
 # Hexagonal pattern
 Cu  0.0  0.0  0.0
 Cu  0.5  0.5  0.0
 Cu  1.0  0.0  0.0
 Cu  1.5  0.5  0.0
 Cu  0.25  0.75  0.0
 Cu  0.75  0.25  0.0
end

nwpw
 simulation_cell
  lattice_vectors
   2.5    0.0    0.0
   1.25   2.165  0.0
   0.0    0.0    20.0
  end
 cutoff 30.0
 smear fermi 0.01
 monkhorst-pack 6 6 1
end

task pspw energy
```

**Expected Results:**
- Full grid: 36 k-points
- After C6v symmetry reduction: ~6 k-points (reduction factor ~6x)

## Example 5: Low Symmetry System (C1)

```nwchem
start low_symmetry
permanent_dir ./perm
scratch_dir ./scratch

geometry "Low Symmetry System" units angstroms
 system crystal
  lat_a 3.0
  lat_b 4.0
  lat_c 5.0
  alpha 85.0
  beta  90.0
  gamma 95.0
 end
 Cu  0.1  0.2  0.3
 Cu  0.6  0.7  0.8
 Cu  0.3  0.4  0.5
end

nwpw
 simulation_cell
  lattice_vectors
   3.0  0.0  0.0
   0.0  4.0  0.0
   0.0  0.0  5.0
  end
 cutoff 30.0
 smear fermi 0.01
 monkhorst-pack 4 4 4
end

task pspw energy
```

**Expected Results:**
- Full grid: 64 k-points
- After C1 symmetry reduction: 64 k-points (no reduction)

## Advanced Usage

### Manual Symmetry Specification

You can also manually specify the symmetry group if automatic detection fails:

```nwchem
nwpw
 simulation_cell
  fcc 3.615
 end
 cutoff 30.0
 smear fermi 0.01
 monkhorst-pack 8 8 8
 # Manual symmetry specification (if needed)
 symmetry_group Oh
end
```

### Surface-Specific Considerations

For surfaces, consider these guidelines:

1. **Vacuum layer**: Ensure sufficient vacuum (typically 10-20 Å)
2. **Surface termination**: Choose appropriate surface termination
3. **Adsorbate placement**: Consider symmetry of adsorbate binding sites
4. **k-point sampling**: Use denser sampling in surface plane, sparse in vacuum direction

### Performance Tips

1. **Symmetry detection**: The system automatically detects symmetry, but you can verify by checking the output
2. **Memory usage**: Reduced k-points significantly decrease memory requirements
3. **Convergence**: Monitor energy convergence with respect to k-point density
4. **Accuracy**: Verify that symmetry reduction doesn't affect accuracy

## Troubleshooting

### Common Issues

1. **No symmetry detected**: Check lattice parameters and atomic positions
2. **Unexpected reduction**: Verify the detected symmetry group
3. **Convergence issues**: Increase k-point density if needed

### Debugging

To debug symmetry detection, add verbose output:

```nwchem
nwpw
 simulation_cell
  fcc 3.615
 end
 cutoff 30.0
 smear fermi 0.01
 monkhorst-pack 8 8 8
 symmetry_tolerance 1e-6
end
```

## Implementation Notes

The symmetry detection and k-point reduction system:

1. **Automatic detection**: Analyzes lattice vectors and atomic positions
2. **Point group operations**: Generates appropriate symmetry matrices
3. **IBZ reduction**: Maps k-points to irreducible Brillouin zone
4. **Weight preservation**: Maintains correct k-point weights
5. **Backward compatibility**: Works with existing input files

This system provides significant computational savings while maintaining accuracy for a wide range of crystal systems and surface structures. 