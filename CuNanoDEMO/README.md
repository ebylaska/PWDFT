# CuNanoDemo Timings

|machine          | software   | ncpus | ngpus | cputime/step | non-local | lagrange | fft  |
|-----------------|------------|-------|-------|--------------|-----------|----------|------|
|jlse iris-gen9   |     SYCL   | 1     | 1     | 0.72         | 0.06      | 0.05     | 0.40 |
|                 |            |       |       |              |           |          |      |
|jlse iris-gen9   | MPI-SYCL   | 1     | 1     | 2.87         | 0.34      | 0.26     | 2.21 |
|jlse iris-gen9   | MPI-SYCL   | 2     | 1     | 2.51         | 0.30      | 0.24     | 1.84 |
|jlse iris-gen9   | MPI-SYCL   | 4     | 1     | 2.73         | 0.30      | 0.24     | 2.01 |
|                 |            |       |       |              |           |          |      |
|jlse iris-gen9   | MPI-ONLY   | 1     | 1     | 3.90         | 1.17      | 0.65     | 1.72 |
|jlse iris-gen9   | MPI-ONLY   | 2     | 1     | 2.39         | 0.70      | 0.41     | 1.01 |
|jlse iris-gen9   | MPI-ONLY   | 4     | 1     | 1.53         | 0.54      | 0.34     | 0.48 |
|                 |            |       |       |              |           |          |      |
|home macbook-pro | MPI-OPENCL | 1     | 1     | 2.87         | 0.34      | 0.26     | 2.21 |
