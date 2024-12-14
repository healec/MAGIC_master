c     =======================================================
      subroutine conduct(mbc,mx,my,mz,
     &                   xlower,ylower,zlower,dx,dy,dz,meqn,
     &                   q,mdiffq,maux,aux,mcoefaux,dtr,ndtrc)
c     =======================================================
c
      implicit double precision (a-h,o-z)
c
      dimension    q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc)
      dimension  aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc)
c     # Local Storage
      dimension    qt(1-mbc:mx+mbc, 1-mbc:my+mbc,
     &                1-mbc:mz+mbc)
      dimension     u(1-mbc:mx+mbc, 1-mbc:my+mbc,
     &                1-mbc:mz+mbc)
      dimension coefs(1-mbc:mx+mbc, 1-mbc:my+mbc,
     &                1-mbc:mz+mbc)
c
c     --------------------------------------------------------
c
c     Conduction Solution for 3D MAGIC Extensions
c        by J. B. Snively, 2003-2019
c        Subject to terms in the LICENSE file.
c
c     --------------------------------------------------------
c
c     Explicit First Order Conduction Solver
c
      dtdx2 = (dtr) / (dx*dx)
      dtdy2 = (dtr) / (dy*dy)
      dtdz2 = (dtr) / (dz*dz) 
c
c
c     # Set Coefs, copy q to improve caching
      do 5 k=1-mbc,mz+mbc
         do 5 j=1-mbc,my+mbc
            do 5 i=1-mbc,mx+mbc
               coefs(i,j,k)=aux(mcoefaux,i,j,k)
               u(i,j,k)=q(mdiffq,i,j,k)
  5   continue
c
      do 50 n=1,ndtrc
c 
c     # Apply conduction or diffusion (x)
c
      do 10 k=1,mz
         do 10 j=1,my
            do 10 i=1,mx
               am=u(i-1,j,k)
               aq=u(i,j,k)
               ap=u(i+1,j,k)
               qt(i,j,k) = (ap - aq - aq + am)
  10  continue
c     # Make assignments
      do 15 k=1,mz
         do 15 j=1,my
            do 15 i=1,mx
               u(i,j,k) = u(i,j,k) +
     &         dtdx2 * coefs(i,j,k) * qt(i,j,k)
  15  continue
c
c     # Apply conduction or diffusion (y)
c
      do 20 k=1,mz
         do 20 j=1,my
            do 20 i=1,mx
               ab=u(i,j-1,k)
               aq=u(i,j,k)
               af=u(i,j+1,k)
               qt(i,j,k) = (af - aq - aq + ab)
  20  continue
c     # Make assignments
      do 25 k=1,mz
         do 25 j=1,my
            do 25 i=1,mx
               u(i,j,k) = u(i,j,k) +
     &         dtdy2 * coefs(i,j,k) * qt(i,j,k)
  25  continue
c
c     # Apply simplified z boundary conditions to u
c     # Note that these were not set by bc3, and 
c     # never cross CPU domain boundaries.
c     
      do 27 j = 1-mbc, my+mbc
         do 27 i = 1-mbc, mz+mbc
            do 27 kbc=1,mbc
               u(i,j,1-kbc)=u(i,j,1)
               u(i,j,mz+kbc)=u(i,j,mz)
  27  continue
c
c     # Apply conduction or diffusion (z)
c
      do 30 k=1,mz
         do 30 j=1,my
            do 30 i=1,mx
               au=u(i,j,k+1)
               aq=u(i,j,k)
               ad=u(i,j,k-1)
               qt(i,j,k) = (ad - aq - aq + au)
  30  continue
c     # Make assignments to u
      do 35 k=1,mz
         do 35 j=1,my
            do 35 i=1,mx
               u(i,j,k) = u(i,j,k) +
     &                    dtdz2 * coefs(i,j,k) * qt(i,j,k)
  35  continue
c
c
  50  continue
c
      do 60 k=1,mz
         do 60 j=1,my
            do 60 i=1,mx
               q(mdiffq,i,j,k) = u(i,j,k)
  60  continue
c
      return
      end
