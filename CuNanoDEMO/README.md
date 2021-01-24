# CuNanoDemo Timings

## System: Cu20 + CH3OH
 - input = cunano-10c.nw
 - cutoff = 60 Ry
 - FCC=38
 - nfft = 80x80x80

|machine          | software   | ncpus | ngpus | cputime/step | non-local | lagrange | fft  | date      |
|-----------------|------------|-------|-------|--------------|-----------|----------|------|-----------|
|jlse iris-gen9   |     SYCL   | 1     | 1     | 0.72         | 0.06      | 0.05     | 0.40 | 1-23-2021 |
|                 |            |       |       |              |           |          |      |           |
|jlse iris-gen9   | MPI-SYCL   | 1     | 1     | 2.87         | 0.34      | 0.26     | 2.21 | 1-23-2021 |
|jlse iris-gen9   | MPI-SYCL   | 2     | 1     | 2.51         | 0.30      | 0.24     | 1.84 | 1-23-2021 |
|jlse iris-gen9   | MPI-SYCL   | 4     | 1     | 2.73         | 0.30      | 0.24     | 2.01 | 1-23-2021 |
|                 |            |       |       |              |           |          |      |           |
|jlse iris-gen9   | MPI-ONLY   | 1     | 0     | 3.90         | 1.17      | 0.65     | 1.72 | 1-23-2021 |
|jlse iris-gen9   | MPI-ONLY   | 2     | 0     | 2.39         | 0.70      | 0.41     | 1.01 | 1-23-2021 |
|jlse iris-gen9   | MPI-ONLY   | 4     | 0     | 1.53         | 0.54      | 0.34     | 0.48 | 1-23-2021 |
|                 |            |       |       |              |           |          |      |           |
|home macbook-pro | MPI-OPENCL | 2     | 1     | 5.03         | 0.83      | 0.46     | 3.05 | 1-23-2021 |
|home macbook-pro | MPI-OPENCL | 3     | 1     | 4.09         | 0.72      | 0.53     | 2.33 | 1-23-2021 |
|home macbook-pro | MPI-OPENCL | 4     | 1     | 3.06         | 0.51      | 0.36     | 1.71 | 1-23-2021 |
|home macbook-pro | MPI-OPENCL | 5     | 1     | 2.69         | 0.56      | 0.39     | 1.47 | 1-23-2021 |
|home macbook-pro | MPI-OPENCL | 6     | 1     | 2.74         | 0.60      | 0.50     | 1.38 | 1-23-2021 |
|                 |            |       |       |              |           |          |      |           |
|home macbook-pro | MPI-ONLY   | 2     | 0     | 19.70        | 10.65     | 4.52     | 3.15 | 1-23-2021 |
|home macbook-pro | MPI-ONLY   | 3     | 0     | 14.35        | 7.67      | 3.35     | 2.30 | 1-23-2021 |
|home macbook-pro | MPI-ONLY   | 4     | 0     | 10.88        | 5.77      | 2.53     | 1.73 | 1-23-2021 |
|home macbook-pro | MPI-ONLY   | 5     | 0     | 9.17         | 4.83      | 2.13     | 1.51 | 1-23-2021 |
|home macbook-pro | MPI-ONLY   | 6     | 0     | 7.67         | 4.02      | 1.77     | 1.25 | 1-23-2021 |


## System: Cu20 + CH3OH
 - input = cunano-10d.nw
 - cutoff = 100 Ry
 - FCC=38
 - nfft = 100x100x100
 
 |machine         | software   | ncpus | ngpus | cputime/step | non-local | lagrange | fft  | date      |
|-----------------|------------|-------|-------|--------------|-----------|----------|------|-----------|
|home macbook-pro | MPI-OPENCL | 2     | 1     | 11.79        | 1.62      | 0.95     | 7.45 | 1-23-2021 | 
|home macbook-pro | MPI-OPENCL | 3     | 1     | 8.61         | 1.39      | 0.91     | 5.24 | 1-23-2021 |
|home macbook-pro | MPI-OPENCL | 4     | 1     | 7.26         | 1.31      | 0.81     | 4.31 | 1-23-2021 |
|home macbook-pro | MPI-OPENCL | 5     | 1     | 6.29         | 1.19      | 0.73     | 3.60 | 1-23-2021 |
|home macbook-pro | MPI-OPENCL | 6     | 1     | 5.27         | 1.30      | 0.70     | 3.07 | 1-23-2021 |
|home macbook-pro | MPI-OPENCL | 7     | 1     | 5.09         | 1.17      | 0.72     | 2.79 | 1-23-2021 |
|home macbook-pro | MPI-OPENCL | 8     | 1     | 4.75         | 1.22      | 0.78     | 2.64 | 1-23-2021 |
|                 |            |       |       |              |           |          |      |           |
|home macbook-pro | MPI-ONLY   | 2     | 0     | 46.66        | 25.17     | 10.61    | 7.71 | 1-23-2021 |
|home macbook-pro | MPI-ONLY   | 3     | 0     | 31.47        | 16.84     | 7.10     | 5.29 | 1-23-2021 |
|home macbook-pro | MPI-ONLY   | 4     | 0     | 24.08        | 12.73     | 5.47     | 4.19 | 1-23-2021 |
|home macbook-pro | MPI-ONLY   | 5     | 0     | 20.41        | 10.68     | 4.58     | 3.61 | 1-23-2021 |
|home macbook-pro | MPI-ONLY   | 6     | 0     | 16.94        | 8.74      | 3.72     | 3.08 | 1-23-2021 |
|home macbook-pro | MPI-ONLY   | 7     | 0     | 15.59        | 8.00      | 3.44     | 2.87 | 1-23-2021 |
|home macbook-pro | MPI-ONLY   | 8     | 0     | 15.18        | 7.73      | 3.33     | 2.85 | 1-23-2021 |

 
 
