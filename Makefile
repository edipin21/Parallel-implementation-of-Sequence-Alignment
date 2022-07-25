build:
	mpicxx -fopenmp -c func.c -o func.o
	mpicxx -fopenmp -c main.c -o main.o
	nvcc  -I./inc -c computeOnGpu.cu -o computeOnGpu.o 
	mpicxx -fopenmp -o mpiCudaOpemMP  main.o func.o computeOnGpu.o   /usr/local/cuda/lib64/libcudart_static.a -ldl -lrt

clean:
	rm -f *.o ./mpiCudaOpemMP

run:
	mpiexec -np 4 ./mpiCudaOpemMP

runOn2:
	mpiexec -np 2 -machinefile  mf  -map-by  node  ./mpiCudaOpemMP
