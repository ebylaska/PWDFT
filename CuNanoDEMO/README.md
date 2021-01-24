# CuNanoDemo Timings

## System: Cu20 + CH3OH
 - nfft = 80x80x80
 - input = cunano-10c.nw
 - cutoff = 60 Ry
 - FCC=38
 - nfft = 80x80x80

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
|home macbook-pro | MPI-OPENCL | 2     | 1     | 5.03         | 0.83      | 0.46     | 3.05 |
|home macbook-pro | MPI-OPENCL | 3     | 1     | 4.09         | 0.72      | 0.53     | 2.33 |
|home macbook-pro | MPI-OPENCL | 4     | 1     | 3.06         | 0.51      | 0.36     | 1.71 |
|home macbook-pro | MPI-OPENCL | 5     | 1     | 2.69         | 0.56      | 0.39     | 1.47 |
|home macbook-pro | MPI-OPENCL | 6     | 1     | 2.74         | 0.60      | 0.50     | 1.38 |
|                 |            |       |       |              |           |          |      |
|home macbook-pro | MPI-ONLY   | 2     | 1     | 19.70        | 10.65     | 4.52     | 3.15 |
|home macbook-pro | MPI-ONLY   | 3     | 1     | 14.35        | 7.67      | 3.35     | 2.30 |
|home macbook-pro | MPI-ONLY   | 4     | 1     | 10.88        | 5.77      | 2.53     | 1.73 |
|home macbook-pro | MPI-ONLY   | 5     | 1     | 9.17         | 4.83      | 2.13     | 1.51 |
|home macbook-pro | MPI-ONLY   | 6     | 1     | 7.67         | 4.02      | 1.77     | 1.25 |
