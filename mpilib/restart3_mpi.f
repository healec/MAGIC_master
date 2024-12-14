c
c     =====================================================
      subroutine restart(meqn,mbc,mx,my,mz,
     &      xlower,ylower,zlower,dx,dy,dz,q)
c     =====================================================
c
c     --------------------------------------------------------
c     Modified in 2018-19 for 5.X/MAGIC compatibility by 
c             Jonathan B. Snively
c     --------------------------------------------------------
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
      include 'mpif.h'
c
c
      dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &      1-mbc:mz+mbc)
      character*10 fname1, fname2
      logical outt0
      common /restrt_block/ tinitial, iframe, outt0

c
c     # MPI: get number of processors and id of this process.
c
      call mpi_comm_rank(mpi_comm_world,id,ierr)

      if (id.eq.0) then
         write(*,*) '*** Restart is not supported from fort.qXXXX ',
     &        'output files from MPICLAW'
         write(*,*) '*** The HDF output version of MPICLAW supports ',
     &        'restarting'
      end if
      call mpi_finalize(ierr)
      STOP 

      return
      end
