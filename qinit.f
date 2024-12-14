c
c     =====================================================
      subroutine qinit(meqn,mbc,mx,my,mz,
     &      xlower,ylower,zlower,dx,dy,dz,q,maux,aux)
c     =====================================================
c
c
c  --------------------------------------------------------
c
c     Initialization / Array Copy for MAGIC Extensions
c        by J. B. Snively, 2003-2019
c        Subject to terms in the LICENSE file.
c
c  --------------------------------------------------------
c
c
c     # Set initial conditions for q.
c
c     # Eliminates a crash from the random number generator
c      use ifport
c
      implicit double precision (a-h,o-z)
      external ohsteady, restart
c
      dimension q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &       1-mbc:mz+mbc)
      dimension aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &       1-mbc:mz+mbc)
c
c     # MPI : Added so that we know which processor we are running on.
      common /mpi_proc_info/ np, id
c
      logical loadrestart
      logical usenoise
      data loadrestart /.false./
      data usenoise /.false./
      tweakfactor=1.d-5
c
      if (loadrestart) then
c        # Comment out because obsolete; use Inchin's HDF5 instead
         call restart(meqn,mbc,mx,my,mz,
     &                xlower,ylower,zlower,dx,dy,dz,q)
         dtoh=-1.d0
         goto 10
      else
         dtoh=0.d0
      endif
c
c     # Do you want 3D evolutions from 2D? Then add noise...
c
      if (usenoise) then
         do k = 1-mbc,mz+mbc
            call srand(k*(id+1))
            do j = 1-mbc,my+mbc
               do i = 1-mbc,mx+mbc
                  tweak = tweakfactor*2*(rand()-0.5d0)
                  q(1,i,j,k) = aux(1,i,j,k)*(1+tweak)
                  tweakmomentum = q(1,i,j,k)/aux(1,i,j,k)
                  q(2,i,j,k) = tweakmomentum*aux(2,i,j,k)
                  q(3,i,j,k) = tweakmomentum*aux(3,i,j,k)
                  q(4,i,j,k) = tweakmomentum*aux(4,i,j,k)
                  q(5,i,j,k) = aux(5,i,j,k)-
     &               0.5d0*(aux(2,i,j,k)**2+aux(3,i,j,k)**2+
     &               aux(4,i,j,k)**2)/aux(1,i,j,k)+
     &               0.5d0*(q(2,i,j,k)**2+q(3,i,j,k)**2+
     &               q(4,i,j,k)**2)/q(1,i,j,k)
                  do m = 6, meqn
                     q(m,i,j,k) = aux(m,i,j,k)
                  enddo
               enddo
            enddo
         enddo
      else
         do k = 1-mbc,mz+mbc
            do j = 1-mbc,my+mbc
               do i = 1-mbc,mx+mbc
                  do m = 1,meqn
                     q(m,i,j,k) = aux(m,i,j,k)
                  enddo
               enddo
            enddo
         enddo
      endif
c
c     # Initialize Photochemistry
 10    continue
! 10   call ohsteady(meqn,mbc,mx,my,mz,
!     &     xlower,ylower,zlower,dx,dy,dz,q,maux,aux,dtoh)
c
c
c
      return
      end
