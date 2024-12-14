c     =======================================================
      subroutine diffuse(mbc,mx,my,mz,
     &                   xlower,ylower,zlower,dx,dy,dz,
     &                   meqn,q,maux,aux,mcoefaux,dtr)
c     =======================================================
c
      implicit double precision (a-h,o-z)
c
      dimension     q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &                1-mbc:mz+mbc)
      dimension   aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &                1-mbc:mz+mbc)
c     # Local Storage
      dimension    qt(1-mbc:mx+mbc, 1-mbc:my+mbc,
     &                1-mbc:mz+mbc)
      dimension     u(1-mbc:mx+mbc, 1-mbc:my+mbc,
     &                1-mbc:mz+mbc, 3)
      dimension coefs(1-mbc:mx+mbc, 1-mbc:my+mbc,
     &                1-mbc:mz+mbc)
c
c     --------------------------------------------------------
c
c     Momentum Diffusion Solution for 3D MAGIC Extensions
c        by J. B. Snively, 2003-2019
c        Subject to terms in the LICENSE file.
c
c     --------------------------------------------------------
c
c     Explicit Euler Viscous Diffusion Solver
c     (Assumes Stokes' Hypothesis and Small Variations of mu)
c
c     # Three dimensions, one algorithm... split solution,
c     # optimized for fast performance, not readability nor
c     # higher-order accuracy in time. (Sorry.)
c     # Pending is a 3D ODE-based solution; it will cost CPU time.
c
      fourthirds=4.d0/3.d0
c
c     # Second derivative factors
c
      dtdx2 = (dtr) / (dx*dx)
      dtdy2 = (dtr) / (dy*dy)
      dtdz2 = (dtr) / (dz*dz)
c
c     # Apply lower/upper boundary conditions to q;
c     # rest were done by bc3 in prior call.
c
      do 1 kbc=1,mbc
         do 1 j = 1-mbc, my+mbc
            do 1 i = 1-mbc, mx+mbc
               do 1 m=2,4
                  q(m,i,j,1-kbc)=q(m,i,j,1)
                  q(m,i,j,mz+kbc)=q(m,i,j,mz)
  1   continue
c
c     # Set Coefs, copy q to improve memory access
c
      do 2 k = 1-mbc, mz+mbc
         do 2 j = 1-mbc, my+mbc
            do 2 i = 1-mbc, mx+mbc
               coefs(i,j,k)=aux(mcoefaux,i,j,k)
               u(i,j,k,1)=q(2,i,j,k)
               u(i,j,k,2)=q(3,i,j,k)
               u(i,j,k,3)=q(4,i,j,k)
  2   continue
c
c     # Apply 3D Cross-Term Solution (not split, using new u)
c       (Again, note that this shares the same stability constraint.)
c     # Cross derivative initial factors
      dtr32 = dtr/6
      dt2dxdy = dtr32 / (dx*dy)
      dt2dxdz = dtr32 / (dx*dz)
      dt2dydz = dtr32 / (dy*dz)
c
      do 5 k=1,mz
        do 5 j=1,my
          do 5 i=1,mx
            ad1a=u(i+1,j+1,k,2)
            ad2a=u(i-1,j+1,k,2)
            ad3a=u(i-1,j-1,k,2)
            ad4a=u(i+1,j-1,k,2)
            ctxxy=(ad1a-ad2a+ad3a-ad4a)
            ad1b=u(i+1,j,k+1,3)
            ad2b=u(i-1,j,k+1,3)
            ad3b=u(i-1,j,k-1,3)
            ad4b=u(i+1,j,k-1,3)
            ctxxz=(ad1b-ad2b+ad3b-ad4b)
            ad1a=u(i+1,j+1,k,1)
            ad2a=u(i-1,j+1,k,1)
            ad3a=u(i-1,j-1,k,1)
            ad4a=u(i+1,j-1,k,1)
            ctyxy=(ad1a-ad2a+ad3a-ad4a)
            ad1b=u(i,j+1,k+1,3)
            ad2b=u(i,j-1,k+1,3)
            ad3b=u(i,j-1,k-1,3)
            ad4b=u(i,j+1,k-1,3)
            ctyyz=(ad1b-ad2b+ad3b-ad4b)
            ad1a=u(i+1,j,k+1,1)
            ad2a=u(i-1,j,k+1,1)
            ad3a=u(i-1,j,k-1,1)
            ad4a=u(i+1,j,k-1,1)
            ctzxz=(ad1a-ad2a+ad3a-ad4a)
            ad1b=u(i,j+1,k+1,2)
            ad2b=u(i,j-1,k+1,2)
            ad3b=u(i,j-1,k-1,2)
            ad4b=u(i,j+1,k-1,2)
            ctzyz=(ad1b-ad2b+ad3b-ad4b)
c
c         # Cross terms, copying solution to q array
c
            q(2,i,j,k) = u(i,j,k,1) + coefs(i,j,k) *
     &                        (dt2dxdy * ctxxy + dt2dxdz * ctxxz)
            q(3,i,j,k) = u(i,j,k,2) + coefs(i,j,k) *
     &                        (dt2dxdy * ctyxy + dt2dydz * ctyyz)
            q(4,i,j,k) = u(i,j,k,3) + coefs(i,j,k) *
     &                        (dt2dxdz * ctzxz + dt2dydz * ctzyz)
c
  5   continue
c
c     # Copy q to u
c
      do 6 k = 1, mz
         do 6 j = 1, my
            do 6 i = 1, mx
               do 6 mdqm1 = 1,3
                  u(i,j,k,mdqm1)=q(mdqm1+1,i,j,k)
  6   continue
c
c     --------------------------------------------------------
c
c     # Solve Heat Equation for Each 3D Direction & Component
c
c     --------------------------------------------------------
c
      do 40 mdqm1=1,3   !mdqm1=mdiffq-1
c
c     # Select Appropriate Factors
c
      if (mdqm1 .eq. 1) then
         dtdx2 = fourthirds * dtdx2
      else if (mdqm1 .eq. 2) then
         dtdy2 = fourthirds * dtdy2
         dtdx2 = dtdx2 / fourthirds
      else
         dtdz2 = fourthirds * dtdz2
         dtdy2 = dtdy2 / fourthirds
      endif
c
c     --------------------------------------------------------
c
c     # Apply diffusion to u,v,w (x)
c
      do 10 k=1,mz
        do 10 j=1,my
          do 10 i=1,mx
            am=u(i-1,j,k,mdqm1)
            aq=u(i,j,k,mdqm1)
            ap=u(i+1,j,k,mdqm1)
            qt(i,j,k) = (ap - aq - aq + am)     
  10  continue
c     # Make assignments
      do 15 k=1,mz
        do 15 j=1,my
          do 15 i=1,mx
            u(i,j,k,mdqm1) = u(i,j,k,mdqm1) + 
     &      dtdx2 * coefs(i,j,k) * qt(i,j,k)
  15  continue
c
c     --------------------------------------------------------
c
c     # Apply diffusion to u,v,w (y)
c
      do 20 k=1,mz
        do 20 j=1,my
          do 20 i=1,mx
            ab=u(i,j-1,k,mdqm1)
            aq=u(i,j,k,mdqm1)
            af=u(i,j+1,k,mdqm1)
            qt(i,j,k) = (af - aq - aq + ab)
  20  continue
c     # Make assignments
      do 25 k=1,mz
        do 25 j=1,my
          do 25 i=1,mx
            u(i,j,k,mdqm1) = u(i,j,k,mdqm1) +
     &      dtdy2 * coefs(i,j,k) * qt(i,j,k)
  25  continue
c
c     --------------------------------------------------------
c
c     # Apply diffusion to u,v,w (z)
c
      do 30 k=1,mz
        do 30 j=1,my
          do 30 i=1,mx
            au=u(i,j,k+1,mdqm1)
            aq=u(i,j,k,mdqm1)
            ad=u(i,j,k-1,mdqm1)
            qt(i,j,k) = (ad - aq - aq + au)
  30  continue
c     # Make assignments
      do 35 k=1,mz
        do 35 j=1,my
          do 35 i=1,mx
c           # Only update q(mdqm1+1)
            q(mdqm1+1,i,j,k) = u(i,j,k,mdqm1) +
     &      dtdz2 * coefs(i,j,k) * qt(i,j,k)
  35  continue 
c
  40  continue
c
c
      return
      end
