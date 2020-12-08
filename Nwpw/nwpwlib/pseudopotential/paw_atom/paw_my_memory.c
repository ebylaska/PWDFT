/*
   $Id$
*/

#include  <stdlib.h>
#include  <stdio.h>

/******************************


******************************/
float*  paw_alloc_1d_array(int n)
{
    float* array;
    int     i;

    array = (float*) malloc(n*sizeof(float));

    if (array==0)
    {
        printf("memory allocation problem in alloc_1d_array\n");
        abort();
    }

    for (i=0;i<n;i++)
        array[i]=0;

    return array;

}


/******************************


******************************/
void   paw_dealloc_1d_array(float* array)
{
    free(array);
}


/******************************


******************************/
float** paw_alloc_2d_array(int n1,int n2)
{
    float** array;
    int     i;
    int     j;

    array = (float**) malloc(n1*sizeof(float*));

    if (array==0)
    {
        printf("memory allocation problem in alloc_2d_array\n");
        abort();
    }

    for (i=0;i<n1;i++)
    {
        array[i] = (float*) malloc(n2*sizeof(float));

        if (array==0)
        {
            printf("memory allocation problem in alloc_1d_array\n");
            abort();
        }

    }

    for (i=0;i<n1;i++)
    {
        for (j=0;j<n2;j++)
            array[i][j] = 0.0;

    }

    return array;

}


/******************************


******************************/
void   paw_dealloc_2d_array(int n1,int n2,float** array)
{
    int     i;

    for (i=0;i<n1;i++)
        free(array[i]);

    free(array);

}


