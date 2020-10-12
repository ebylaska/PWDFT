#include <stdio.h>
#include <stdlib.h>

#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/cl.h>
#endif

#define MAX_SOURCE_SIZE (0x100000) 

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
      printf("   -- device %d id = %d, number of cores=%d @ %d MHz memory=%ld bytes wfloat=%d wdouble=%d\n",i,(int) device_ids[i],num_cores,freq,(long) mem,wfloat,wdouble);
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
/*
    ret = clGetDeviceIDs( platform_id, CL_DEVICE_TYPE_ALL, 1,
            &device_id, &ret_num_devices);
*/

}


int main() 
{
   nwpw_gpu_init();
}
