
#include "Parallel.hpp"

/***********************************
 *                                 *
 *         DMatrix_dgemm2          *
 *                                 *
 ***********************************/
void DMatrix_dgemm2(int m, int n, int k, int nblock,
                    double alpha,
                    double *A, int lda, int *ma, int *na,
                    double *B, int ldb, int *mb, int *nb,
                    double beta,
                    double *C, int ldc, int *mc, int *nc,
                    int  taskid_i,int taskid_j,
                    int  np_i, int np_j,
                    MPI_comm comm_i, MPI_comm comm_j,
                    double  *work1, double *work2)
{

   for (auto ij=0; ij<(nc[taskid_j]*mc[taskid_i]); ++ij)
      C[ij] = beta*C[ij];

   int iwrk;
   int ii=0;
   int jj=0;
   int kk=0;
   int icur=0;
   int jcur=0;

   /* loop over all row pannels of C */
   while (kk<m) 
   {
      iwrk = min(nblock, mc(icur+1)-ii)
      iwrk = min(iwrk,   na(jcur+1)-jj)

      /* iwrk*nc(taskid_j) submatrix !=0 */
      if (ma[taskid_i]>0)
      {
         /* pack current iwrk columns of A into work1 */
         if (taskid_j==jcur)
         {
            dlacpy("G",ma[taskid_i],iwrk,
                       &A[jj*lda],lda,
                       work1, ma[taskid_i]);
         }

         /* broadcast work1  within my row */
         ierr = MPI_Bcast(work1,iwrk*ma[taskid_i],MPI_DOUBLE_PRECISION,jcur,comm_j);

         if ((iwrk>0) && (nb[taskid_j]>0))
            DGEMM_OMP('T','N',iwrk,nb[taskid_j],ma[taskid_i],
                      alpha,
                      work1, ma[taskid_i],
                      B, ldb,
                      0.0,
                      work2, iwrk);
      }

      /* iwrk*nc(taskid_j+1) submatrix ==0 */
      else
         Parallel_shared_vector_zero(true,nc[taskid_j]*iwrk,work2);

      /* summ to node that holds current rows of C */
      ierr = MPI_Reduce(work2,work1,nc[taskid_j]*iwrk,
                          MPI_DOUBLE_PRECISION,MPI_SUM,icur,comm_i);

      /* add to current rows of C */
      if (taskid_i==icur) 
      {
         shift = 0;
         for (auto i=ii, i<(ii+iwrk-1)
         {
            DAXPY_OMP(nc[taskid_j],1.0,&work1[shift],iwrk,&C[i],mc[taskid_i]);
            ++shift;
         }
      }
      ii += iwrk;
      jj += iwrk;
      kk += iwrk;

      if (jj>=na[jcur]) 
      {
        ++jcur;
        jj=0;
      }
      if (ii>=mc[icur]) then
      {
        ++icur;
        ii=0;
      }

   }
}





c        **** summ to node that holds current rows of C ****
!$OMP MASTER
#ifdef MPI4
         stupid_msglen = nc(taskid_j+1)*iwrk
         stupid_taskid = icur
         call MPI_Reduce(work2,work1,stupid_msglen,
     >                   stupid_double,stupid_sum,
     >                   stupid_taskid,stupid_comm_i,stupid_ierr)
#else
         call MPI_Reduce(work2,work1,nc(taskid_j+1)*iwrk,
     >                   MPI_DOUBLE_PRECISION,MPI_SUM,icur,comm_i,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


c        **** add to current rows of C ****
         if (taskid_i.eq.icur) then
            shift = 1
            do i=ii,(ii+iwrk-1)
               call DAXPY_OMP(nc(taskid_j+1),1.0d0,work1(shift),iwrk,
     >                                    C(i+1,1),mc(taskid_i+1))
               shift = shift + 1
            end do
         end if

         ii = ii + iwrk
      }
   }

      ii = 0
      jj = 0
      kk = 0
      icur = 0
      jcur = 0
c     **** loop over all row pannels of C ***
      do while (kk.lt.m)
         iwrk = min(nblock, mc(icur+1)-ii)
         iwrk = min(iwrk,   na(jcur+1)-jj)


*        **** iwrk*nc(taskid_j+1) submatrix !=0 ****
         if (ma(taskid_i+1).gt.0) then

!$OMP MASTER
*           **** pack current iwrk columns of A into work1 ***
            if (taskid_j.eq.jcur) then
               call dlacpy("G", ma(taskid_i+1),iwrk,
     >                   A(1,jj+1), lda,
     >                   work1,     ma(taskid_i+1))
            end if

c           **** broadcast work1  within my row ***
#ifdef MPI4
            stupid_msglen = iwrk*ma(taskid_i+1)
            stupid_taskid = jcur
            call MPI_Bcast(work1,stupid_msglen,
     >                     stupid_double,stupid_taskid,
     >                     stupid_comm_j,stupid_ierr)
#else
            call MPI_Bcast(work1,iwrk*ma(taskid_i+1),
     >                     MPI_DOUBLE_PRECISION,jcur,comm_j,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


c            if ((iwrk.gt.0)          .and.
c     >          (nb(taskid_j+1).gt.0).and.
c     >          (ma(taskid_i+1).gt.0))
            if ((iwrk.gt.0)          .and.
     >          (nb(taskid_j+1).gt.0))
     >        call DGEMM_OMP('T','N',iwrk,nb(taskid_j+1),ma(taskid_i+1),
     >                   alpha,
     >                   work1, ma(taskid_i+1),
     >                   B, ldb,
     >                   0.0d0,
     >                   work2, iwrk)

*        **** iwrk*nc(taskid_j+1) submatrix ==0 ****
         else
            !call dcopy(nc(taskid_j+1)*iwrk,0.0d0,0,work2,1)
            call Parallel_shared_vector_zero(.true.,
     >                                       nc(taskid_j+1)*iwrk,work2)
         end if


c        **** summ to node that holds current rows of C ****
!$OMP MASTER
#ifdef MPI4
         stupid_msglen = nc(taskid_j+1)*iwrk
         stupid_taskid = icur
         call MPI_Reduce(work2,work1,stupid_msglen,
     >                   stupid_double,stupid_sum,
     >                   stupid_taskid,stupid_comm_i,stupid_ierr)
#else
         call MPI_Reduce(work2,work1,nc(taskid_j+1)*iwrk,
     >                   MPI_DOUBLE_PRECISION,MPI_SUM,icur,comm_i,ierr)
#endif
!$OMP END MASTER
!$OMP BARRIER


c        **** add to current rows of C ****
         if (taskid_i.eq.icur) then
            shift = 1
            do i=ii,(ii+iwrk-1)
               call DAXPY_OMP(nc(taskid_j+1),1.0d0,work1(shift),iwrk,
     >                                    C(i+1,1),mc(taskid_i+1))
               shift = shift + 1
            end do
         end if

         ii = ii + iwrk
         jj = jj + iwrk
         kk = kk + iwrk

         if (jj.ge.na(jcur+1)) then
           jcur = jcur + 1
           jj   = 0
         end if
         if (ii.ge.mc(icur+1)) then
           icur = icur + 1
           ii   = 0
         end if

      end do
   }
   }
