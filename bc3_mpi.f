
c
c
c     =====================================================
      subroutine bc3(meqn,mbc,mx,my,mz,xlower,
     &               ylower,zlower,dx,dy,dz,q,maux,aux,t,dt,mthbc)
c     =====================================================
c
c
c  --------------------------------------------------------
c
c     Boundary Conditions for MAGIC Extensions
c        by J. B. Snively, 2003-2019
c        Subject to terms in the LICENSE file.
c
c     Based closely on the original MPI Clawpack code.
c
c  --------------------------------------------------------
c
c
c     # Standard boundary condition choices for claw3
c
c     # At each boundary  k = 1 (xlower),  2 (xupper),
c     #                       3 (ylower),  4 (yupper),
c     #                       5 (zlower),  6 (zupper):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  or 4'th (for k=5,6) component of q.
c     ------------------------------------------------
c
c     # Extend the data from the interior cells (1:mx, 1:my, 1:mz)
c     # to a layer of mbc ghost cells outside the region.
c
      implicit double precision (a-h,o-z)
      include 'mpif.h'

      dimension    q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc)
      dimension  aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc)

c     # MPI : Create arrays for communication with neighboring processors
c     # These will be automatically allocated here, and will be deallocated
c     # when we exit this routine.

      dimension ql(1:meqn, 1:mbc, 1:my, 1:mz)
      dimension qr(1:meqn, 1:mbc, 1:my, 1:mz)
      dimension qt(1:meqn, 1-mbc:mx+mbc, 1:mbc, 1:mz)
      dimension qb(1:meqn, 1-mbc:mx+mbc, 1:mbc, 1:mz)
      dimension qk(1:meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1:mbc)
      dimension qf(1:meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1:mbc)

      dimension ql_send(1:meqn, 1:mbc, 1:my, 1:mz)
      dimension qr_send(1:meqn, 1:mbc, 1:my, 1:mz)
      dimension qt_send(1:meqn, 1-mbc:mx+mbc, 1:mbc, 1:mz)
      dimension qb_send(1:meqn, 1-mbc:mx+mbc, 1:mbc, 1:mz)
      dimension qk_send(1:meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1:mbc)
      dimension qf_send(1:meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1:mbc)

      parameter (ileft = 1, iright = 2, ibottom = 3, itop = 4,
     &      iback = 5, ifront = 6)

      dimension mthbc(6), ireq(6), ireqs(6) 
      dimension MPIstatus(MPI_STATUS_SIZE,6)
      dimension id_coords(3), id_next(3), id_last(3)
      dimension mtotal(3), mstart(3)

c     # mpi : this communicator gives us information about the mapping 
c     # of nodes to grid subdomains.
      common /mpicomm/ mpi_comm_3d, lx, ly, lz, mtotal, mstart
      common /mpi_proc_info/ np, id
c
c
c     # mpi : initialize number of outstanding communication requests
      nreq = 0
      nreqs = 0

c     # Get the coordinates of the current processor in the processor array.
c     # id_coords(1) is the i coordinate of this processor in the array.
c     # id_coords(2) is the j coordinate of this processor in the array.
c     # id_coords(3) is the k coordinate of this processor in the array.
      call mpi_cart_coords(mpi_comm_3d,id,3,id_coords,ierr)
      imap = id_coords(1)
      jmap = id_coords(2)
      kmap = id_coords(3)
c
c     -------------------------------------------------------
c     # left boundary (xlower):
c     -------------------------------------------------------
c     # MPI : added option 4 : internal boundary
c
c     # Skip application of boundary if mthbc(i) <1
      if (mthbc(ileft).LT.0) go to 199
c
      go to (100,110,120,130,140) mthbc(ileft) + 1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc3'
      call mpi_finalize(ierr)
      stop
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do 115 k = 1-mbc, mz+mbc
         do 115 j = 1-mbc, my+mbc
            do 115 ibc=1,mbc
               do 115 m=1,meqn
                  q(m,1-ibc,j,k) = q(m,1,j,k)
  115 continue
      go to 199

  120 continue
c     # MPI : periodic
      if (lx == 1) then
         do 125 k=1,mz
            do 125 j=1,my
               do 125 i=1,mbc
                  do 125 m=1,meqn
                     ql(m,i,j,k) = q(m,mx-mbc+i,j,k)
  125    continue
      else
         print*,'Something happened...'
         STOP 'periodic boundaries are treated as internal for lx>1'
      endif
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 k = 1-mbc, mz+mbc
         do 135 j = 1-mbc, my+mbc
            do 135 ibc=1,mbc
               do 135 m=1,meqn
                  q(m,1-ibc,j,k) = q(m,ibc,j,k)
  135 continue
c     # negate the normal velocity:
      do 136 k = 1-mbc, mz+mbc
         do 136 j = 1-mbc, my+mbc
            do 136 ibc=1,mbc
               q(2,1-ibc,j,k) = -q(2,ibc,j,k)
  136 continue
      go to 199

  140 continue
c     # MPI : Internal boundary

c     # Send left side of domain
      nbcdata = mbc*my*mz*meqn
      id_last(1) = imap-1
      id_last(2) = jmap
      id_last(3) = kmap
      call mpi_cart_rank(mpi_comm_3d,id_last,idcomm,ierr)

c     # Receive left ghost cells
c     # Post non-blocking receive and continue with send
      nreq = nreq + 1
      call mpi_irecv(ql, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm, ileft, MPI_COMM_WORLD, ireq(nreq), ierr)
c
      do 190 k=1,mz
         do 190 j=1,my
            do 190 i=1,mbc
               do 190 m=1,meqn
                  ql_send(m,i,j,k) = q(m,i,j,k)
  190 continue
      nreqs = nreqs + 1
      call mpi_isend(ql_send, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm, iright, MPI_COMM_WORLD, ireqs(nreqs), ierr)

      go to 199

  199 continue
c
c     -------------------------------------------------------
c     # right boundary (xupper):
c     -------------------------------------------------------
c     # MPI : added option 4 : internal boundary
c
c     # Skip application of boundary if mthbc(i) <1
      if (mthbc(iright).LT.0) go to 299
c
      go to (200,210,220,230,240) mthbc(iright)+1
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc3'
      call mpi_finalize(ierr)
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do 215 k = 1-mbc, mz+mbc
         do 215 j = 1-mbc, my+mbc
            do 215 ibc=1,mbc
               do 215 m=1,meqn
                  q(m,mx+ibc,j,k) = q(m,mx,j,k)
  215 continue
      go to 299

  220 continue
c     # MPI : periodic:
      if (lx == 1) THEN
c        # All data is local. Apply periodic boundary condition.
         do 225 k=1,mz
            do 225 j=1,my
               do 225 i=1,mbc
                  do 225 m=1,meqn
                     qr(m,i,j,k) = q(m,i,j,k)
  225    continue
      else
         print*,'Something happened...'
         STOP 'periodic boundaries are treated as internal for lx>1'
      endif
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 k = 1-mbc, mz+mbc
         do 235 j = 1-mbc, my+mbc
            do 235 ibc=1,mbc
               do 235 m=1,meqn
                  q(m,mx+ibc,j,k) = q(m,mx+1-ibc,j,k)
  235 continue
c     # negate the normal velocity:
      do 236 k = 1-mbc, mz+mbc
         do 236 j= 1-mbc, my+mbc
            do 236 ibc=1,mbc
               q(2,mx+ibc,j,k) = -q(2,mx+1-ibc,j,k)
  236 continue
      go to 299


  240 continue
c     # MPI : Internal boundary condition

c     # Send right side of domain
      nbcdata = mbc*my*mz*meqn
      id_next(1) = imap+1
      id_next(2) = jmap
      id_next(3) = kmap
      call mpi_cart_rank(mpi_comm_3d,id_next,idcomm,ierr)

c     # Receive right ghost cells
c     # Post non-blocking receive and continue with send
      nreq = nreq + 1
      call mpi_irecv(qr, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm, iright, MPI_COMM_WORLD, ireq(nreq), ierr)

      do 290 k=1,mz
         do 290 j=1,my
            do 290 i=1,mbc
               do 290 m=1,meqn
                  qr_send(m,i,j,k) = q(m,mx-mbc+i,j,k)
  290 continue
      nreqs = nreqs + 1
      call mpi_isend(qr_send,nbcdata,MPI_DOUBLE_PRECISION,
     &      idcomm, ileft, MPI_COMM_WORLD, ireqs(nreqs), ierr)

      go to 299

  299 continue
c
c     -------------------------------------------------------
c     # MPI: If information has been passed from other processors, 
c     # record it in q so that corner ghost cells in y- and z-
c     # directions will have correct information.
c     -------------------------------------------------------
c     # MPI : Wait for all processes to complete
      do i = 1,nreq
         call mpi_wait(ireq(i),MPIstatus(1,i),ierr)
         call mpi_wait(ireqs(i),MPIstatus(1,i),ierr)
      enddo
c     # mpi : re-initialize number of outstanding communication requests
      nreq = 0
      nreqs = 0

c     # Copy values from communication buffers back into main array
c     # Left boundary
      select case (mthbc(ileft))
      case (2,4)
         do k=1,mz
            do j=1,my
               do i=1,mbc
                  do m=1,meqn
                     q(m,i-mbc,j,k) = ql(m,i,j,k)
                  end do
               end do
            end do
         end do
      case default
      end select

c     # Right boundary
      select case (mthbc(iright))
      case (2,4)
         do k=1,mz
            do j=1,my
               do i=1,mbc
                  do m=1,meqn
                     q(m,mx+i,j,k) = qr(m,i,j,k)
                  end do
               end do
            end do
         end do
      case default
      end select
c
c     -------------------------------------------------------
c     # back boundary (ylower):
c     -------------------------------------------------------
c     # MPI : added option 4 : internal boundary
c
c     # Skip application of boundary if mthbc(i) <1
      if (mthbc(ibottom).LT.0) go to 399
c
      go to (300,310,320,330,340) mthbc(ibottom)+1
c
  300 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc3'
      call mpi_finalize(ierr)
      stop
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 k = 1-mbc, mz+mbc
         do 315 jbc=1,mbc
            do 315 i = 1-mbc, mx+mbc
               do 315 m=1,meqn
                  q(m,i,1-jbc,k) = q(m,i,1,k)
  315 continue
      go to 399

  320 continue
c     # MPI : periodic:
      if (ly == 1) THEN
c        # All data is local. Apply periodic boundary condition.
         do 325 k=1,mz
            do 325 j=1,mbc
               do 325 i = 1-mbc, mx+mbc
                  do 325 m=1,meqn
                     qb(m,i,j,k) = q(m,i,my-mbc+j,k)
  325    continue
      else
         print*,'Something happened...'
         STOP 'periodic boundaries are treated as internal for ly>1'
      endif
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 k = 1-mbc, mz+mbc
         do 335 jbc=1,mbc
            do 335 i = 1-mbc, mx+mbc
               do 335 m=1,meqn
                  q(m,i,1-jbc,k) = q(m,i,jbc,k)
  335 continue
c     # negate the normal velocity:
      do 336 jbc=1,mbc
         do 336 i = 1-mbc, mx+mbc
            do 336 k = 1-mbc, mz+mbc
               q(i,1-jbc,k,3) = -q(i,jbc,k,3)
  336    continue
      go to 399

  340 continue
c     # MPI : Internal boundary condition

c     # Send bottom side of domain
      nbcdata = (mx+2*mbc)*mbc*mz*meqn
      id_last(1) = imap
      id_last(2) = jmap-1
      id_last(3) = kmap
      call mpi_cart_rank(mpi_comm_3d,id_last,idcomm,ierr)

c     # Receive bottom ghost cells
c     # Post non-blocking receive and continue with send
      nreq = nreq + 1
      CALL mpi_irecv(qb, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm, ibottom, MPI_COMM_WORLD, ireq(nreq), ierr)

      do 390 k = 1,mz
         do 390 j = 1,mbc
            do 390 i = 1-mbc,mx+mbc
               do 390 m=1,meqn
                  qb_send(m,i,j,k) = q(m,i,j,k)
  390 continue
      nreqs = nreqs + 1
      call mpi_isend(qb_send, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm,itop, MPI_COMM_WORLD, ireqs(nreqs), ierr)

      go to 399

  399 continue
c
c     -------------------------------------------------------
c     # front boundary (yupper):
c     -------------------------------------------------------
c     # MPI : added option 4 : internal boundary
c
c     # Skip application of boundary if mthbc(i) <1
      if (mthbc(itop).LT.0) go to 499
c
      go to (400,410,420,430,440) mthbc(itop)+1
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc3'
      call mpi_finalize(ierr)
      stop
      go to 499

  410 continue
c     # zero-order extrapolation:
      do 415 k = 1-mbc, mz+mbc
         do 415 jbc=1,mbc
            do 415 i = 1-mbc, mx+mbc
               do 415 m=1,meqn
                  q(m,i,my+jbc,k) = q(m,i,my,k)
  415 continue
      go to 499

  420 continue
c     # MPI : periodic:
      if(ly == 1) THEN
c        # All data is local. Apply periodic boundary condition.
         do 425 k=1,mz
            do 425 j=1,mbc
               do 425 i = 1-mbc,mx+mbc
                  do 425 m=1,meqn
                     qt(m,i,j,k) = q(m,i,j,k)
  425    continue 
      else
         print*,'Something happened...'
         STOP 'periodic boundaries are treated as internal for ly>1'
      endif
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 k = 1-mbc, mz+mbc
         do 435 jbc=1,mbc
            do 435 i = 1-mbc, mx+mbc
               do 435 m=1,meqn
                  q(m,i,my+jbc,k) = q(m,i,my+1-jbc,k)
  435 continue
c     # negate the normal velocity:
      do 436 k = 1-mbc, mz+mbc
         do 436 jbc=1,mbc
            do 436 i = 1-mbc, mx+mbc
               q(3,i,my+jbc,k) = -q(3,i,my+1-jbc,k)
  436    continue
      go to 499

  440 continue
c     # MPI : Internal boundary condition

c     # Send top of domain
      nbcdata = (mx+2*mbc)*mbc*mz*meqn
      id_next(1) = imap
      id_next(2) = jmap+1
      id_next(3) = kmap
      call mpi_cart_rank(mpi_comm_3d,id_next,idcomm,ierr)

c     # Receive top ghost cells
c     # Post non-blocking receive and continue with send
      nreq = nreq + 1
      call mpi_irecv(qt, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm,itop, MPI_COMM_WORLD, ireq(nreq), ierr)

      do 490  k=1,mz
         do 490 j=1,mbc
            do 490 i = 1-mbc,mx+mbc
               do 490 m=1,meqn
                  qt_send(m,i,j,k) = q(m,i,my-mbc+j,k)
  490 continue
      nreqs = nreqs + 1
      call mpi_isend(qt_send,nbcdata,MPI_DOUBLE_PRECISION,
     &      idcomm,ibottom, MPI_COMM_WORLD, ireqs(nreqs), ierr)

      go to 499

  499 continue

c
c     -------------------------------------------------------
c     # MPI: If information has been passed from other processors, 
c     # record it in q so that corner ghost cells in z-direction
c     # will have correct information.
c     -------------------------------------------------------
c     # MPI : Wait for all processes to complete
      do i = 1,nreq
         call mpi_wait(ireq(i),MPIstatus(1,i),ierr)
         call mpi_wait(ireqs(i),MPIstatus(1,i),ierr)
      enddo
c     # mpi : re-initialize number of outstanding communication requests
      nreq = 0
      nreqs = 0

c     # Copy values from communication buffers back into main array
c     #   Bottom boundary
      select case (mthbc(ibottom))
      case (2,4)
         do k=1,mz
            do j=1,mbc
               do i = 1-mbc, mx+mbc
                  do m=1,meqn
                     q(m,i,j-mbc,k) = qb(m,i,j,k)
                  end do
               end do
            end do
         end do
      case default
      end select

c     # Top boundary
      select case (mthbc(itop))
      case (2,4)
         do k=1,mz
            do j=1,mbc
               do i = 1-mbc, mx+mbc
                  do m=1,meqn
                     q(m,i,my+j,k) = qt(m,i,j,k)
                  end do
               end do
            end do
         end do
      case default
      end select
c
c     -------------------------------------------------------
c     # ground boundary (zlower):
c     -------------------------------------------------------
c     # MPI : added option 4 : internal boundary
c
c     # Skip application of boundary if mthbc(i) <1
      if (mthbc(iback).LT.0) go to 599
c
      go to (500,510,520,530,540) mthbc(iback)+1
c
  500 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(5)=0 and no BCs specified in bc3'
      call mpi_finalize(ierr)
      stop
      go to 599
c
  510 continue
c     # zero-order extrapolation:
      do 515 kbc=1,mbc
         do 515 j = 1-mbc, my+mbc
            do 515 i = 1-mbc, mx+mbc
            q(1,i,j,1-kbc)=aux(1,i,j,1-kbc)+
     &           (q(1,i,j,1)-aux(1,i,j,1))*
     &           (aux(1,i,j,1-kbc)/aux(1,i,j,1))**(.5)
            q(2,i,j,1-kbc)=aux(2,i,j,1-kbc)+
     &           (q(2,i,j,1)-aux(2,i,j,1))*
     &           (aux(1,i,j,1-kbc)/aux(1,i,j,1))**(.5)
            q(3,i,j,1-kbc)=aux(3,i,j,1-kbc)+
     &           (q(3,i,j,1)-aux(3,i,j,1))*
     &           (aux(1,i,j,1-kbc)/aux(1,i,j,1))**(.5)
            q(4,i,j,1-kbc)=q(4,i,j,1)*
     &           (aux(1,i,j,1-kbc)/aux(1,i,j,1))**(.5)
            q(5,i,j,1-kbc)=aux(5,i,j,1-kbc)+
     &           (q(5,i,j,1)-aux(5,i,j,1))*
     &           (aux(1,i,j,1-kbc)/aux(1,i,j,1))**(.5)
            do 515 m=6,meqn
               q(m,i,j,1-kbc)=q(m,i,j,1)
  515       continue
      go to 599

  520 continue
c     # MPI : periodic:
      if (lz == 1) then
c        # All data is local. Apply periodic boundary condition.
         do 525 k = 1,mbc
            do 525 j = 1-mbc,my+mbc
               do 525 i = 1-mbc,mx+mbc
                  do 525 m = 1,meqn
                     qb(m,i,j,k) = q(m,i,j,mz-mbc+k)
  525    continue
      else
         print*,'Something happened...'
         STOP 'periodic boundaries are treated as internal for lz>1'
      endif
      go to 599

  530 continue
c     # solid wall (assumes 4'rd component is velocity or momentum in y):
      do 535 kbc=1,mbc
         do 535 j = 1-mbc, my+mbc
            do 535 i = 1-mbc, mx+mbc
               q(1,i,j,1-kbc)=aux(1,i,j,1-kbc)+
     &           (q(1,i,j,1)-aux(1,i,j,1))*
     &           (aux(1,i,j,1-kbc)/aux(1,i,j,1))**(.5)
               q(2,i,j,1-kbc)=aux(2,i,j,1-kbc)+
     &           (q(2,i,j,1)-aux(2,i,j,1))*
     &           (aux(1,i,j,1-kbc)/aux(1,i,j,1))**(.5)
               q(3,i,j,1-kbc)=aux(3,i,j,1-kbc)+
     &           (q(3,i,j,1)-aux(3,i,j,1))*
     &           (aux(1,i,j,1-kbc)/aux(1,i,j,1))**(.5)
               q(4,i,j,1-kbc)=q(4,i,j,1)*
     &           (aux(1,i,j,1-kbc)/aux(1,i,j,1))**(.5)
               q(5,i,j,1-kbc)=aux(5,i,j,1-kbc)+
     &           (q(5,i,j,1)-aux(5,i,j,1))*
     &           (aux(1,i,j,1-kbc)/aux(1,i,j,1))**(.5)
            do 535 m = 6, meqn
               q(m,i,j,1-kbc)=q(m,i,j,1)
  535 continue
c     # negate the normal velocity:
      do 536 kbc=1,mbc
         do 536 j = 1-mbc, my+mbc
            do 536 i = 1-mbc, mx+mbc
               q(4,i,j,1-kbc) = -q(4,i,j,kbc)
  536 continue
      go to 599

  540 continue
c     # MPI : Internal boundary

c     # Send back side of domain
      nbcdata = (mx+2*mbc)*(my+2*mbc)*mbc*meqn
      id_last(1) = imap
      id_last(2) = jmap
      id_last(3) = kmap-1
      call mpi_cart_rank(mpi_comm_3d,id_last,idcomm,ierr)

c     # Receive back ghost cells
      nreq = nreq + 1
      call mpi_irecv(qk, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm, iback, MPI_COMM_WORLD, ireq(nreq), ierr)

      do 590 k = 1,mbc
         do 590 j = 1-mbc,my+mbc
            do 590 i = 1-mbc,mx+mbc
                do 590 m = 1,meqn
                  qk_send(m,i,j,k) = q(m,i,j,k)
  590 continue
      nreqs = nreqs + 1
      call mpi_isend(qk_send, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm, ifront, MPI_COMM_WORLD, ireqs(nreqs), ierr)

      go to 599

  599 continue
c
c     -------------------------------------------------------
c     # top boundary (zupper):
c     -------------------------------------------------------
c     # MPI : added option 4 : internal boundary
c
c     # Skip application of boundary if mthbc(i) <1
      if (mthbc(ifront).LT.0) go to 699
c
      go to (600,610,620,630,640) mthbc(ifront)+1
c
  600 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(6)=0 and no BCs specified in bc3'
      call mpi_finalize(ierr)
      stop
      go to 699

  610 continue
c     # zero-order extrapolation:
      do 615 kbc=1,mbc
         do 615 j = 1-mbc, my+mbc
            do 615 i = 1-mbc, mx+mbc
               q(1,i,j,mz+kbc)=aux(1,i,j,mz+kbc)+
     &           (q(1,i,j,mz)-aux(1,i,j,mz))*
     &           (aux(1,i,j,mz+kbc)/aux(1,i,j,mz))**(.5)
               q(2,i,j,mz+kbc)=aux(2,i,j,mz+kbc)+
     &           (q(2,i,j,mz)-aux(2,i,j,mz))*
     &           (aux(1,i,j,mz+kbc)/aux(1,i,j,mz))**(.5)
               q(3,i,j,mz+kbc)=aux(3,i,j,mz+kbc)+
     &           (q(3,i,j,mz)-aux(3,i,j,mz))*
     &           (aux(1,i,j,mz+kbc)/aux(1,i,j,mz))**(.5)
               q(4,i,j,mz+kbc)=q(4,i,j,mz)*
     &           (aux(1,i,j,mz+kbc)/aux(1,i,j,mz))**(.5)
               q(5,i,j,mz+kbc)=aux(5,i,j,mz+kbc)+
     &           (q(5,i,j,mz)-aux(5,i,j,mz))*
     &           (aux(1,i,j,mz+kbc)/aux(1,i,j,mz))**(.5)
            do 615 m=6,meqn
               q(m,i,j,mz+kbc)=q(m,i,j,mz)
  615 continue
      go to 699

  620 continue
c     # MPI : periodic:
      if (lz == 1) then
c        # All data is local. Apply periodic boundary condition.
         do 625 k = 1,mbc
            do 625 j = 1-mbc,my+mbc
               do 625 i = 1-mbc,mx+mbc
                  do 625 m = 1,meqn
                     qf(m,i,j,k) = q(m,i,j,k)
  625    continue
      else
         print*,'Something happened...'
         STOP 'periodic boundaries are treated as internal for lz>1'
      endif
      go to 699

  630 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 635 kbc=1,mbc
         do 635 j = 1-mbc, my+mbc
            do 635 i = 1-mbc, mx+mbc
               q(1,i,j,mz+kbc)=aux(1,i,j,mz+kbc)+
     &           (q(1,i,j,mz)-aux(1,i,j,mz))*
     &           (aux(1,i,j,mz+kbc)/aux(1,i,j,mz))**(.5)
               q(2,i,j,mz+kbc)=aux(2,i,j,mz+kbc)+
     &           (q(2,i,j,mz)-aux(2,i,j,mz))*
     &           (aux(1,i,j,mz+kbc)/aux(1,i,j,mz))**(.5)
               q(3,i,j,mz+kbc)=aux(3,i,j,mz+kbc)+
     &           (q(3,i,j,mz)-aux(3,i,j,mz))*
     &           (aux(1,i,j,mz+kbc)/aux(1,i,j,mz))**(.5)
               q(4,i,j,mz+kbc)=q(4,i,j,mz)*
     &           (aux(1,i,j,mz+kbc)/aux(1,i,j,mz))**(.5)
               q(5,i,j,mz+kbc)=aux(5,i,j,mz+kbc)+
     &           (q(5,i,j,mz)-aux(5,i,j,mz))*
     &           (aux(1,i,j,mz+kbc)/aux(1,i,j,mz))**(.5)
            do 635 m=6,meqn
               q(m,i,j,mz+kbc)=q(m,i,j,mz)
  635 continue
c     # negate the normal velocity:
      do 636 kbc=1,mbc
         do 636 j = 1-mbc, my+mbc
            do 636 i = 1-mbc, mx+mbc
               q(4,i,j,mz+kbc) = -q(4,i,j,mz+1-kbc)
  636 continue
      go to 699

  640 continue
c     # MPI : Internal boundary

c     # Send front side of domain
      nbcdata = (mx+2*mbc)*(my+2*mbc)*mbc*meqn
      id_next(1) = imap
      id_next(2) = jmap
      id_next(3) = kmap+1
      call mpi_cart_rank(mpi_comm_3d,id_next,idcomm,ierr)

c     # Receive front ghost cells
c     # Post non-blocking receive and continue with send
      nreq = nreq + 1
      call mpi_irecv(qf, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm,ifront, MPI_COMM_WORLD, ireq(nreq), ierr)

      do 690 k=1,mbc
         do 690 j = 1-mbc,my+mbc
            do 690 i = 1-mbc,mx+mbc
               do 690 m=1,meqn
                  qf_send(m,i,j,k) = q(m,i,j,mz-mbc+k)
  690 continue
      nreqs = nreqs + 1
      call mpi_isend(qf_send, nbcdata, MPI_DOUBLE_PRECISION,
     &      idcomm,iback, MPI_COMM_WORLD, ireqs(nreqs),  ierr)

      go to 699

  699 continue

c     # MPI : Wait for all processes to complete
      do i = 1,nreq
         call mpi_wait(ireq(i),MPIstatus(1,i),ierr)
         call mpi_wait(ireqs(i),MPIstatus(1,i),ierr)
      enddo
c     # mpi : re-initialize number of outstanding communication requests
      nreq = 0
      nreqs = 0

c     # Copy values from communication buffers back into main array
c     # Back boundary
      select case (mthbc(iback))
      case (2,4)
         do k=1,mbc
            do j = 1-mbc,my+mbc
               do i = 1-mbc,mx+mbc
                  do m=1,meqn
                     q(m,i,j,k-mbc) = qk(m,i,j,k)
                  end do
               end do
            end do
         end do
      case default
      end select

c     # Front boundary
      select case (mthbc(ifront))
      case (2,4)
         do k=1,mbc
            do j = 1-mbc,my+mbc
               do i = 1-mbc,mx+mbc
                  do m=1,meqn
                     q(m,i,j,mz+k) = qf(m,i,j,k)
                  end do
               end do
            end do
         end do
      case default
      end select

      return
      end
