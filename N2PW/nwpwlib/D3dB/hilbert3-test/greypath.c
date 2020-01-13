#include	<math.h>
#include	<stdio.h>

greypath2(int i)
{
   int j,n,m,tn;
   if (i<2)
      return i;
   else
   {
      n  = (int) floor(log(1.0*i)/log(2.0));
      tn = 1; for (j=0; j<n; ++j) tn *= 2;
      m  = 2*tn-i-1;
      return tn+greypath2(m);
   }
}

greypath3(int i, int d)
{
   int j,n,m,t2n,t3n,maxi;

   maxi = 1; for (j=0; j<d; ++j) maxi *= 2;

   if (i<2)
      return i;
   else if (i==maxi)
   {
       t3n = 1; for (j=1; j<d; ++j) t3n *= 3; 
       return(2*t3n);
   }
   else
   {
      n  = (int) floor(log(1.0*i)/log(2.0));
      t2n = 1; for (j=0; j<n; ++j) t2n *= 2;
      t3n = 1; for (j=0; j<n; ++j) t3n *= 3;
      m  = i-t2n;
      if (m==0)
         return(   t3n + greypath3(t2n-1,d) );
      else
         return( 2*t3n + greypath3(t2n-m,d) );
   }
}
      

      

main()
{
   int i;
   for (i=0; i<=32; ++i)
      printf("i = %d path2 = %d, path3 = %d \n",i,greypath2(i),greypath3(i,5));

}

