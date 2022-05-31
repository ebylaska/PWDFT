#include <cstdio>
#include <cstdlib>
#include <string>
#include <cmath>
#include <time.h>

#include <chrono>
using namespace std::chrono;


#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/opencl.h>

#else
#include <sycl/cl.h>
#endif


#define MAX_SOURCE_SIZE (0x100000) 

typedef struct {
   cl_context context;
   cl_command_queue command_queue;
   cl_program program;
   cl_uint nbuf;
   cl_uint * bufsize;
   cl_mem  * buf_mem_obj;

   cl_kernel kernel;

} NWPW_GPU_Program_Type;

typedef struct {
   cl_uint num_platforms;
   cl_platform_id * platform_id;
   cl_uint * num_devices;
   cl_device_id ** device_id;
   cl_bool   ** avail;
   cl_bool   ** has_cl_khr_fp64;
   cl_uint   ** num_cores;
   cl_uint   ** freq;
   cl_uint   ** wdouble;
   cl_uint   ** wfloat;
   cl_ulong  ** mem;
   //Device_Type ** device;
} NWPW_GPU_Type;

/***************************************
 *                                     *
 *           nwpw_gpu_init             *
 *                                     *
 ***************************************/
void nwpw_gpu_init(NWPW_GPU_Type *gpu)
{
   size_t size;
   char str[1000];

   // Get platforms 
   cl_int ret = clGetPlatformIDs(0, NULL, &(gpu->num_platforms));
   gpu->platform_id = (cl_platform_id *) malloc(sizeof(cl_platform_id)*gpu->num_platforms);
   ret = clGetPlatformIDs(gpu->num_platforms,gpu->platform_id,NULL);

   gpu->num_devices = (cl_uint *) malloc(sizeof(cl_uint)*gpu->num_platforms);
   gpu->device_id = (cl_device_id **) malloc(sizeof(cl_device_id *)*gpu->num_platforms);
   gpu->avail     = (cl_bool **) malloc(sizeof(cl_bool *)*gpu->num_platforms);
   gpu->has_cl_khr_fp64 = (cl_bool **) malloc(sizeof(cl_bool *)*gpu->num_platforms);
   gpu->num_cores = (cl_uint **) malloc(sizeof(cl_uint *)*gpu->num_platforms);
   gpu->freq      = (cl_uint **) malloc(sizeof(cl_uint *)*gpu->num_platforms);
   gpu->wdouble   = (cl_uint **) malloc(sizeof(cl_uint *)*gpu->num_platforms);
   gpu->wfloat    = (cl_uint **) malloc(sizeof(cl_uint *)*gpu->num_platforms);
   gpu->mem       = (cl_ulong **) malloc(sizeof(cl_ulong *)*gpu->num_platforms);
   for (cl_uint i=0; i<gpu->num_platforms; ++i)
   {
      ret = clGetDeviceIDs(gpu->platform_id[i], CL_DEVICE_TYPE_ALL, 0,NULL,&(gpu->num_devices[i]));

      gpu->device_id[i] = (cl_device_id *) malloc(sizeof(cl_device_id)*gpu->num_devices[i]);
      gpu->avail[i]     = (cl_bool *) malloc(sizeof(cl_bool)*gpu->num_devices[i]);
      gpu->has_cl_khr_fp64[i] = (cl_bool *) malloc(sizeof(cl_bool)*gpu->num_devices[i]);
      gpu->num_cores[i] = (cl_uint *) malloc(sizeof(cl_uint)*gpu->num_devices[i]);
      gpu->freq[i]      = (cl_uint *) malloc(sizeof(cl_uint)*gpu->num_devices[i]);
      gpu->wdouble[i]   = (cl_uint *) malloc(sizeof(cl_uint)*gpu->num_devices[i]);
      gpu->wfloat[i]    = (cl_uint *) malloc(sizeof(cl_uint)*gpu->num_devices[i]);
      gpu->mem[i]       = (cl_ulong *) malloc(sizeof(cl_ulong)*gpu->num_devices[i]);

      ret = clGetDeviceIDs(gpu->platform_id[i],CL_DEVICE_TYPE_ALL,gpu->num_devices[i],gpu->device_id[i],NULL);

      for (cl_uint j=0; j<gpu->num_devices[i]; ++j)
      {
         ret = clGetDeviceInfo(gpu->device_id[i][j],CL_DEVICE_AVAILABLE,sizeof(cl_bool),&(gpu->avail[i][j]),&size);
         ret = clGetDeviceInfo(gpu->device_id[i][j],CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(cl_uint),&(gpu->num_cores[i][j]),&size);
         ret = clGetDeviceInfo(gpu->device_id[i][j],CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_uint),&(gpu->freq[i][j]),&size);
         ret = clGetDeviceInfo(gpu->device_id[i][j],CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,sizeof(cl_uint),&(gpu->wdouble[i][j]),&size);
         ret = clGetDeviceInfo(gpu->device_id[i][j],CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,sizeof(cl_uint),&(gpu->wfloat[i][j]),&size);
         ret = clGetDeviceInfo(gpu->device_id[i][j],CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&(gpu->mem[i][j]),&size);
         ret = clGetDeviceInfo(gpu->device_id[i][j],CL_DEVICE_EXTENSIONS,1000*sizeof(char),str,&size);
         gpu->has_cl_khr_fp64[i][j] = (strstr(str, "cl_khr_fp64") != NULL);
      }
   }
}

/***************************************
 *                                     *
 *           nwpw_gpu_end              *
 *                                     *
 ***************************************/
void nwpw_gpu_end(NWPW_GPU_Type *gpu)
{

   for (int i=0; i<gpu->num_platforms; ++i)
   {
      free(gpu->device_id[i]);
      free(gpu->avail[i]);
      free(gpu->has_cl_khr_fp64[i]);
      free(gpu->num_cores[i]);
      free(gpu->freq[i]);
      free(gpu->wdouble[i]);
      free(gpu->wfloat[i]);
      free(gpu->mem[i]);
   }
   free(gpu->platform_id); 
   free(gpu->num_devices);
   free(gpu->device_id);
   free(gpu->avail);
   free(gpu->has_cl_khr_fp64);
   free(gpu->num_cores);
   free(gpu->freq);
   free(gpu->wdouble);
   free(gpu->wfloat);
   free(gpu->mem);
}

/***************************************
 *                                     *
 *    nwpw_gpu_program_generate        *
 *                                     *
 ***************************************/
void nwpw_gpu_program_generate(NWPW_GPU_Program_Type *gpuprogram, 
                               cl_device_id device_id,
                               const char * source_str,
                               const size_t source_size,
                               const int nbuf,
                               const int bufsize[],
                               const int iswrite[])
{
   cl_int ret;

   // Create an OpenCL context
   gpuprogram->context = clCreateContext(NULL,1, &(device_id), NULL, NULL, &ret);

   // Create a command queue
   gpuprogram->command_queue = clCreateCommandQueue(gpuprogram->context, device_id, 0, &ret);

   // Create memory buffers on the device for each vector 
   gpuprogram->nbuf = (cl_uint) nbuf;
   gpuprogram->bufsize = (cl_uint *) bufsize;
   gpuprogram->buf_mem_obj = (cl_mem *) malloc(sizeof(cl_mem)*nbuf);
   for (int i=0; i<nbuf; ++i)
   {
      printf("id=%d bufsize=%d\n",i,gpuprogram->bufsize[i]);
      if (iswrite[i]) 
      {
         gpuprogram->buf_mem_obj[i] = clCreateBuffer(gpuprogram->context, 
                                                     CL_MEM_WRITE_ONLY, 
                                                     gpuprogram->bufsize[i]*sizeof(double),
                                                     NULL, &ret);
      }
      else 
      {
         gpuprogram->buf_mem_obj[i] = clCreateBuffer(gpuprogram->context,
                                                     CL_MEM_READ_ONLY,
                                                     gpuprogram->bufsize[i]*sizeof(double), 
                                                     NULL, &ret);
      }
   }

   // Create a program from the kernel source
   gpuprogram->program = clCreateProgramWithSource(gpuprogram->context, 1,
            (const char **)&source_str, (const size_t *)&source_size, &ret);

   // Build the program
   ret = clBuildProgram(gpuprogram->program, 1, &device_id, NULL, NULL, NULL);

    // Check for compilation errors
    size_t logSize;
    clGetProgramBuildInfo(gpuprogram->program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
    char* messages = (char*)malloc((1+logSize)*sizeof(char));
    clGetProgramBuildInfo(gpuprogram->program, device_id, CL_PROGRAM_BUILD_LOG, logSize, messages, NULL);
    messages[logSize] = '\0';
    if (logSize > 10) { printf(">>> Compiler message: %s\n", messages); }
    free(messages);

   // Create the OpenCL kernel
   gpuprogram->kernel = clCreateKernel(gpuprogram->program, "matmul", &ret);

}

/***************************************
 *                                     *
 *    nwpw_gpu_program_host2gpu        *
 *                                     *
 ***************************************/
void nwpw_gpu_program_host2gpu(NWPW_GPU_Program_Type *gpuprogram, 
                               const double *A, const int buf_id,
                               cl_event *event) 
{
   // Copy the lists A and B to their respective memory buffers
   cl_int ret = clEnqueueWriteBuffer(gpuprogram->command_queue, 
                                     gpuprogram->buf_mem_obj[buf_id], CL_FALSE, 0,
                                     gpuprogram->bufsize[buf_id]*sizeof(double), A, 0, NULL, event);
}

/***************************************
 *                                     *
 *    nwpw_gpu_program_gpu2host        *
 *                                     *
 ***************************************/
void nwpw_gpu_program_gpu2host(NWPW_GPU_Program_Type *gpuprogram, 
                               const int buf_id, double *C,
                               cl_event *event) 
{
   cl_uint ret = clEnqueueReadBuffer(gpuprogram->command_queue, 
                                     gpuprogram->buf_mem_obj[buf_id], CL_FALSE, 0,
                                     gpuprogram->bufsize[buf_id]*sizeof(double), 
                                     C, 0, NULL, event);
}

/***************************************
 *                                     *
 *    nwpw_gpu_program_run             *
 *                                     *
 ***************************************/
void nwpw_gpu_program_run(NWPW_GPU_Program_Type *gpuprogram, 
                          const int global_size,
                          const int local_size,
                          const cl_uint nevents,
                          cl_event  *events)
{
   cl_uint ret;

   for (int buf_id=0; buf_id<gpuprogram->nbuf; ++buf_id)
      ret = clSetKernelArg(gpuprogram->kernel, buf_id, sizeof(cl_mem), (void *)&(gpuprogram->buf_mem_obj[buf_id]));

   // Execute the OpenCL kernel on the list
   size_t global_item_size = (size_t) global_size;
   size_t local_item_size = (size_t) local_size; 
   ret = clEnqueueNDRangeKernel(gpuprogram->command_queue, 
                                gpuprogram->kernel, 1, NULL,
                                &global_item_size, 
                                &local_item_size, nevents,events, NULL);
}

/***************************************
 *                                     *             
 *    nwpw_gpu_program_run2            *             
 *                                     *             
 ***************************************/
void nwpw_gpu_program_run2(NWPW_GPU_Program_Type *gpuprogram,
                          const int niargs, const int *iargs,
                          const int *global_size,
                          const int *local_size,
                          const cl_uint nevents,
                          cl_event  *events)
{
   int ii = 0;
   cl_uint ret;
   
   // Set the int arguments of the kernel
   for (int i=0; i<niargs; ++i)
   {
      ret = clSetKernelArg(gpuprogram->kernel,ii,sizeof(int),(void*)&iargs[ii]);
      ++ii;
   }

   // Set the buffer arguments of the kernel
   for (int buf_id=0; buf_id<gpuprogram->nbuf; ++buf_id)
   {
      ret = clSetKernelArg(gpuprogram->kernel,ii,sizeof(cl_mem),(void *)&(gpuprogram->buf_mem_obj[buf_id]));
      ++ii;
   }

   // Execute the OpenCL kernel on the list
   const size_t global_item_size[2] = {global_size[0],global_size[1]};
   const size_t local_item_size[2]  = {local_size[0],local_size[1]};
   ret = clEnqueueNDRangeKernel(gpuprogram->command_queue,
                                gpuprogram->kernel, 2, NULL,
                                global_item_size,
                                local_item_size, nevents,events, NULL);
   //ret = clEnqueueNDRangeKernel(gpuprogram->command_queue,
   //                             gpuprogram->kernel, 2, NULL,
   //                             global_item_size,
   //                             local_item_size, 0,NULL,events);
}



void DisplayPlatformInfo(cl_platform_id id, cl_platform_info name, char  str[])
{
   cl_int err;
   size_t paramsize;

   err = clGetPlatformInfo(id,name,0,NULL,&paramsize);

   char * info = (char *) malloc(sizeof(char)*paramsize);

   err = clGetPlatformInfo(id,name,paramsize,info,NULL);
   printf(" - %s : %s\n",str,info);

   free(info);
}

void DisplayDevicesInfo(cl_platform_id id)
{
   size_t size,group;
   char str[1000];
   cl_ulong mem;
   cl_bool  avail;
   cl_uint num_devices,num_cores,freq,wdouble,wfloat;
   cl_int ret = clGetDeviceIDs(id, CL_DEVICE_TYPE_ALL, 0,NULL,&num_devices);
   printf("\n");
   printf(" - number of devices = %d\n", num_devices);

   cl_device_id * device_ids = (cl_device_id *) malloc(sizeof(cl_device_id)*num_devices);
   ret = clGetDeviceIDs(id,CL_DEVICE_TYPE_ALL,num_devices,device_ids,NULL);
   for (cl_uint i=0; i<num_devices; ++i)
   {
      printf("\n");
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(cl_uint),&num_cores,&size);
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_uint),&freq,&size);
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,sizeof(cl_uint),&wdouble,&size);
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,sizeof(cl_uint),&wfloat,&size);
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&mem,&size);
      printf("   -- device %d id = %ld, number of cores=%d @ %d MHz memory=%ld bytes wfloat=%d wdouble=%d\n",i,(size_t) device_ids[i],num_cores,freq,(long) mem,wfloat,wdouble);
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&group,&size);
      printf("   -- device group size %d\n",(int) group);
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_AVAILABLE,sizeof(cl_bool),&avail,&size);
      if (avail) {
         printf("   -- device available\n");
      } else {
         printf("   -- device not available\n");
      } 
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_NAME,1000*sizeof(char),str,&size);
      printf("   -- device name %s\n",str);
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_VERSION,1000*sizeof(char),str,&size);
      printf("   -- device version %s\n",str);
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_VENDOR,1000*sizeof(char),str,&size);
      printf("   -- device vendor %s\n",str);
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_PROFILE,1000*sizeof(char),str,&size);
      printf("   -- device profile %s\n",str);
      ret = clGetDeviceInfo(device_ids[i],CL_DEVICE_EXTENSIONS,1000*sizeof(char),str,&size);
      printf("   -- device extensions %s\n",str);
   }
}


/*
void nwpw_gpu_init()
{
   // Get platform and device information
    cl_platform_id  * platform_ids;
    cl_uint num_platforms;

    cl_int ret = clGetPlatformIDs(0, NULL, &num_platforms);
    printf("Number of gpu platforms = %d\n", num_platforms);

    platform_ids = (cl_platform_id *) malloc(sizeof(cl_platform_id)*num_platforms);

    // query the platform_ids
    ret = clGetPlatformIDs(num_platforms,platform_ids,NULL);
    for (cl_uint i=0; i<num_platforms; ++i)
    {
       printf(" - platform %d id = %d\n",i,(int) platform_ids[i]);
       DisplayPlatformInfo(platform_ids[i],CL_PLATFORM_PROFILE,"CL_PLATFORM_PROFILE");
       DisplayPlatformInfo(platform_ids[i],CL_PLATFORM_VERSION,"CL_PLATFORM_VERSION");
       DisplayPlatformInfo(platform_ids[i],CL_PLATFORM_VENDOR,"CL_PLATFORM_VENDOR");
       DisplayPlatformInfo(platform_ids[i],CL_PLATFORM_EXTENSIONS,"CL_PLATFORM_EXTENSIONS");
       DisplayDevicesInfo(platform_ids[i]);
    }

}
*/

                 

static NWPW_GPU_Type gpu;
static NWPW_GPU_Program_Type gpuprogram;
static int bufsize[6];
static int iswrite[6];


void matmul_start(int Msize, int Nsize, int Ksize)
{


   // Load the kernel source code into the array source_str
   FILE *fp;
   char *source_str;
   size_t source_size;

   fp = fopen("matmul_kernel.cl", "r");
   if (!fp) {
       fprintf(stderr, "Failed to load kernel.\n");
       exit(1);
   }
   source_str = (char*)malloc(MAX_SOURCE_SIZE);
   source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
   fclose(fp);
   printf("source_size=%d\n",source_size);

   // Allocate buffers on gpus
   bufsize[0] = Msize*Ksize; iswrite[0] = 0;
   bufsize[1] = Ksize*Nsize; iswrite[1] = 0;
   bufsize[2] = Msize*Nsize; iswrite[2] = 1;
   printf("buffsize[%d]=%d\n",0,bufsize[0]);
   printf("buffsize[%d]=%d\n",1,bufsize[1]);
   printf("buffsize[%d]=%d\n",2,bufsize[2]);
   nwpw_gpu_program_generate(&gpuprogram,
                             gpu.device_id[0][2],
                             source_str,
                             source_size,
                             3,bufsize,iswrite);
   printf("device_id=%d\n",gpu.device_id[0][2]);
}


int main() 
{
   //const int LIST_SIZE = 64*64*64*128*2;
   const int Msize = 1024;
   const int Nsize = 1024;
   const int Ksize = 1024;

   //NWPW_GPU_Type gpu;
   //NWPW_GPU_Program_Type gpuprogram;

   //int bufsize[6];
   //int iswrite[6];
   cl_uint nevents = 2;
   cl_event events[5];
   cl_event wevent;

   printf("+----------------------------------------------------+\n");
   printf("| OpenCL MatMul-like GPU timings on 2015 Macbook pro |\n");
   printf("+----------------------------------------------------+\n\n\n");

   nwpw_gpu_init(&gpu);

   printf("Number of platforms = %d\n",gpu.num_platforms);
   for (int i=0; i<gpu.num_platforms; ++i) 
   {
      printf(" - %d patform_id= %ld num_devices= %d\n",i,(size_t) gpu.platform_id[i],gpu.num_devices[i]);
      for (int j=0; j<gpu.num_devices[i]; ++j)
      {
         printf("   -- %d device_id= %ld num_cores=%3d mem=%12ld  %4d MHz wfloat=%d wdouble=%d avail=%d has_cl_khr_fp64=%d\n",j,(size_t) gpu.device_id[i][j],
                                                                gpu.num_cores[i][j],
                                                                (long) gpu.mem[i][j],
                                                                gpu.freq[i][j],
                                                                gpu.wfloat[i][j],
                                                                gpu.wdouble[i][j],
                                                                gpu.avail[i][j], 
                                                                gpu.has_cl_khr_fp64[i][j]);
      }
   }

   //clock_t time1,time2,time3,time4,time5,time6,time7,time8;
   //time1 = clock();
   auto ctime1 = high_resolution_clock::now();

   matmul_start(Msize,Nsize,Ksize);

   // Load the kernel source code into the array source_str
/*
   FILE *fp;
   char *source_str;
   size_t source_size;

   fp = fopen("matmul_kernel.cl", "r");
   if (!fp) {
       fprintf(stderr, "Failed to load kernel.\n");
       exit(1);
   }
   source_str = (char*)malloc(MAX_SOURCE_SIZE);
   source_size = fread( source_str, 1, MAX_SOURCE_SIZE, fp);
   fclose(fp);
   printf("source_size=%d\n",source_size);

   //bufsize[0] = 1; iswrite[0] = 2;
   //bufsize[1] = 1; iswrite[1] = 2;
   //bufsize[2] = 1; iswrite[2] = 2;
   // Allocate buffers on gpus
   bufsize[0] = Msize*Ksize; iswrite[0] = 0;
   bufsize[1] = Ksize*Nsize; iswrite[1] = 0;
   bufsize[2] = Msize*Nsize; iswrite[2] = 1;
   printf("buffsize[%d]=%d\n",0,bufsize[0]);
   printf("buffsize[%d]=%d\n",1,bufsize[1]);
   printf("buffsize[%d]=%d\n",2,bufsize[2]);
   nwpw_gpu_program_generate(&gpuprogram, 
                             gpu.device_id[0][2],
                             source_str,
                             source_size,
                             3,bufsize,iswrite);
   printf("device_id=%d\n",gpu.device_id[0][2]);
*/


   // Allocate buffers on host
   double *A = (double*)malloc(sizeof(double)*Msize*Ksize);
   double *B = (double*)malloc(sizeof(double)*Ksize*Nsize);
   double *C = (double*)malloc(sizeof(double)*Msize*Nsize);
   double *CC = (double*)malloc(sizeof(double)*Msize*Nsize);
   for (int i=0; i<Msize*Ksize; ++i) A[i] = ((double) i);
   for (int i=0; i<Ksize*Nsize; ++i) B[i] = sin(((double) i));
   for (int i=0; i<Msize*Nsize; ++i) C[i] = 0.0;

   //time2 = clock();
   auto ctime2 = high_resolution_clock::now();
   //nwpw_gpu_program_host2gpu(&gpuprogram,&Msize,0,&events[0]);
   //nwpw_gpu_program_host2gpu(&gpuprogram,&Nsize,1,&events[1]);
   //nwpw_gpu_program_host2gpu(&gpuprogram,&Ksize,2,&events[2]);
   nwpw_gpu_program_host2gpu(&gpuprogram,A,0,&events[0]);
   nwpw_gpu_program_host2gpu(&gpuprogram,B,1,&events[1]);
   //time3 = clock();
   auto ctime3 = high_resolution_clock::now();
   const int global[2] = {Msize,Nsize};
   const int local[2]  = {16,16};
   const int MNK[3] = {Msize,Nsize,Ksize};

   nwpw_gpu_program_run2(&gpuprogram,3,MNK,global,local,nevents,events);
   //nwpw_gpu_program_run(&gpuprogram,64*64,64,nevents,events);
   //time4 = clock();
   auto ctime4 = high_resolution_clock::now();

   nwpw_gpu_program_gpu2host(&gpuprogram,2,C,&wevent);
   //time5 = clock();
   auto ctime5 = high_resolution_clock::now();

   cl_uint ret = clWaitForEvents(1,&wevent);
   //time6 = clock();
   auto ctime6 = high_resolution_clock::now();


   // Display the result to the screen
/*
   for(int i=0; i<LIST_SIZE; ++i)
      printf("%lf + %lf = %lf\n",A[i],B[i],C[i]);
*/

   duration<double> deltatime21 = ctime2-ctime1;
   duration<double> deltatime32 = ctime3-ctime2;
   duration<double> deltatime43 = ctime4-ctime3;
   duration<double> deltatime54 = ctime5-ctime4;
   duration<double> deltatime65 = ctime6-ctime5;
   double dt21 = (double) deltatime21.count();
   double dt32 = (double) deltatime32.count();
   double dt43 = (double) deltatime43.count();
   double dt54 = (double) deltatime54.count();
   double dt65 = (double) deltatime65.count();

   printf("\n\nUsing gpu platform=%d device=%d\n",0,2);
   printf("\n");
   //printf("GPU compile and setup time           = time2-time1 = %lf ms %lf\n",1000.0*(double)(time2-time1)/CLOCKS_PER_SEC,dt21*1000.0);
   //printf("Host2gpu copy time                   = time3-time2 = %lf ms\n",1000.0*(double)(time3-time2)/CLOCKS_PER_SEC);
   //printf("Gpu time (includes host2gpu barrier) = time4-time3 = %lf ms\n",1000.0*(double)(time4-time3)/CLOCKS_PER_SEC);
   //printf("Gpu2host copy time                   = time5-time4 = %lf ms\n",1000.0*(double)(time5-time4)/CLOCKS_PER_SEC);
   //printf("Gpu2host barrier time                = time6-time5 = %lf ms\n",1000.0*(double)(time6-time5)/CLOCKS_PER_SEC);

   printf("GPU compile and setup time           = time2-time1 = %lf ms\n",dt21*1000.0);
   printf("Host2gpu copy time                   = time3-time2 = %lf ms\n",1000.0*dt32);
   printf("Gpu time (includes host2gpu barrier) = time4-time3 = %lf ms\n",1000.0*dt43);
   printf("Gpu2host copy time                   = time5-time4 = %lf ms\n",1000.0*dt54);
   printf("Gpu2host barrier time                = time6-time5 = %lf ms\n",1000.0*dt65);


   //time7 = clock();
   auto ctime7 = high_resolution_clock::now();
   for(int i=0; i<Msize; ++i)
   for(int j=0; j<Nsize; ++j)
   {
      double acc = 0.0;
      for(int k=0; k<Ksize; ++k)
         acc += A[i+k*Msize]*B[k+j*Nsize];
      CC[i+j*Msize] = acc;
   }
   //time8 = clock();
   auto ctime8 = high_resolution_clock::now();
   duration<double> deltatime87 = ctime8-ctime7;
   duration<double> deltatime52 = ctime5-ctime2;
   duration<double> deltatime62 = ctime6-ctime2;
   double dt87 = (double) deltatime87.count();
   double dt52 = (double) deltatime52.count();
   double dt62 = (double) deltatime62.count();

   printf("\nSerial time = %lf ms\n",1000.0*dt87);
   printf("\nspeedup (gpu only)                                = %lf\n",(dt87/dt43));
   printf("speedup1 (host2gpu+gpu+gpu2host)                  = %lf\n",(dt87/dt52));
   printf("speedup0 (host2gpu+gpu+gpu2host+gpu2host_barrier) = %lf\n",(dt87/dt62));
   for(int i=0; i<2; ++i)
      printf("i=%d C=%le %le\n",i,C[i],CC[i]);
  
   for(int i=10000; i<10002; ++i)
      printf("i=%d C=%le %le\n",i,C[i],CC[i]);


   nwpw_gpu_end(&gpu);
}
