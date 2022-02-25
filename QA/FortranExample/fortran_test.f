      program stevedemo
      implicit none

      include 'mpif.h'

      integer MASTER,taskid,np,mpierr
      parameter (MASTER=0)
      real*8 Energy0,fion(3,2),rion(3,2)
      real*8 qion(2),uion(2)

*     **** MPI initiializer *****
      call MPI_INIT(mpierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,mpierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,np,mpierr)

*     *** needs to be written over tasks ***
      write(*,*) "taskid=",taskid
      rion(1,1) = 0.0d0
      rion(2,1) = 0.0d0
      rion(3,1) = -0.7d0
      rion(1,2) = 0.0d0
      rion(2,2) = 0.0d0
      rion(3,2) = 0.7d0
      uion(1) = 0.0d0
      uion(2) = 0.0d0

      call pspw_fortran_input("h2.nw",len("h2.nw"))
      call pspw_fortran_minimizer(MPI_COMM_WORLD,rion,uion,
     >                            Energy0,fion,qion)

      if (taskid.eq.MASTER) then
         write(*,*) Energy0
         write(*,*) rion(1,1),rion(2,1),rion(3,1)
         write(*,*) rion(1,2),rion(2,2),rion(3,2)
         write(*,*) fion(1,1),fion(2,1),fion(3,1)
         write(*,*) fion(1,2),fion(2,2),fion(3,2)
         write(*,*) qion(1)
         write(*,*) qion(2)
      end if

      call MPI_FINALIZE(mpierr)
      stop
      end
