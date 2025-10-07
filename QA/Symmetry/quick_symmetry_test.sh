#!/bin/bash

# Quick test to verify symmetry detection is working
echo "Testing symmetry detection with FCC Cu..."

# Create test directory
mkdir -p test_symmetry
cd test_symmetry

# Generate FCC Cu input
cat > input.nw << EOF
start test_fcc
permanent_dir ./perm
scratch_dir ./scratch

geometry "FCC Cu" units angstroms
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
 monkhorst-pack 4 4 4
 tolerances 1e-6 1e-6 1e-6
end

task pspw energy
EOF

# Run the calculation
echo "Running calculation..."
~/PWDFT/build/pwdft < input.nw > output.log 2>&1

# Check for symmetry detection output
echo "Checking for symmetry detection output..."
if grep -q "Detected symmetry" output.log; then
    echo "✓ Symmetry detection is working!"
    grep "Detected symmetry" output.log
    grep "Monkhorst-Pack grid" output.log
    grep "Time-reversal reduction" output.log
    grep "Symmetry reduction" output.log
    grep "Total reduction" output.log
else
    echo "✗ Symmetry detection output not found"
    echo "Full output:"
    cat output.log
fi

cd ..
echo "Test completed. Check test_symmetry/output.log for details." 