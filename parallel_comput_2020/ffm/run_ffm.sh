array=( 15000 12000 10000 8000 6000 4000 2000 1000 800 600 400 200 150 100 80 60 40 )

for i in "${array[@]}"
do
    for j in "${array[@]}"
    do
        # For GPU Executable
        SYCL_BE=PI_OPENCL ./GPU_EXE $i $j

        #(or)

        # For CPU Executable
        ./CPU_EXE $i $j
    done
done
