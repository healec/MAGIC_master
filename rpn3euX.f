c
c    ==================================================================
      subroutine rpn3(ixyz,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,
     &		      auxl,auxr,fwave,s,amdq,apdq)
c     ==================================================================
c
c     --------------------------------------------------------
c
c     Normal Riemann Solution for MAGIC Extensions
c        by J. B. Snively, 2003-2019,
c        based on work by LeVeque, Langseth, et al.
c        using limiters of Kemm [2011, 2012] applied to
c        wave strengths (use without Clawpack limiters).
c
c        Subject to terms in the LICENSE file.
c
c     --------------------------------------------------------
c
c     # Roe solver for the Euler equations, f-wave, with gravity and tracers.
c     # Flux limiting is applied within this routine, so that only 3 waves 
c     # are needed in total.
c
c     # Tracers are transported based on edge (Roe) velocities.
c     # The number of tracers = meqn-5 (Each lumped into wave 2)
c     # The first tracers nspecadv are solved via simple advection,
c     # the next tracers meqn-5-nspecadv are solved via continuity eqns.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixyz=1
c     #                               the y-direction if ixyz=2.
c     #                               the z-direction if ixyz=3.
c     # On output, fwave contains the f-waves, s the speeds,
c     # and amdq, apdq the left-going and right-going flux differences,
c     # respectively.  
c
c     # A limiter is applied for each tracer species individually;
c     # this is because limiting them together does not guarantee stability.
c     # Method parameters are specified via logical parameters; see below.
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      dimension fwave(meqn, mwaves, 1-mbc:maxm+mbc)
      dimension     s(mwaves, 1-mbc:maxm+mbc)
      dimension    ql(meqn, 1-mbc:maxm+mbc)
      dimension    qr(meqn, 1-mbc:maxm+mbc)
      dimension  amdq(meqn, 1-mbc:maxm+mbc)
      dimension  apdq(meqn, 1-mbc:maxm+mbc)
      dimension  auxl(maux, 1-mbc:maxm+mbc)
      dimension  auxr(maux, 1-mbc:maxm+mbc)
c
c
c     local arrays and common block (comroe is passed to rpt3eu)
c     ------------
      parameter (maxmrp = 3000)
      dimension delta(5,-1:maxmrp)
c     # Now, to enable CFL-dependent limiting... (modified Clawpack)
      dimension  dtdx(  -1:maxmrp)
      dimension    ID(3,-1:maxmrp)
      dimension   CFL(3,-1:maxmrp)
c
      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom
      common /comroe/ u2v2w2(-1:maxmrp),
     &       u(-1:maxmrp),v(-1:maxmrp),w(-1:maxmrp),enth(-1:maxmrp),
     &       a(-1:maxmrp),g1a2(-1:maxmrp),euv(-1:maxmrp)
c
      common /options/ loadprof, nspecadv
      logical efix, laxliu, thirdorder, nonlocal
c
      data efix /.true./       !# use entropy fix for transonic rarefactions
      data laxliu /.true./     !# use more-expensive Lax-Liu limiting
      data thirdorder /.true./ !# use third order CFL-dependent method
      data nonlocal /.false./   !# use nonlocal CFL [Kemm, 2012]
c
c     # CFL Superbee beta-theta limiter parameters
c     # (Theta = 1, Beta = 2/3 Recommended)
c     
      theta = 1.0d0
      beta  = 0.666666666666666d0
c
c
c     # check hard-coded limits and settings
c
      if (-1.GT.1-mbc .OR. maxmrp .LT. maxm+mbc) then
	 write(6,*) 'need to increase maxmrp in rpn3eu, '
         write(6,*) 'maxmrp = ',maxmrp
         write(6,*) '  maxm = ',maxm
         write(6,*) '   mbc = ',mbc
	 stop
	 endif
c
c     # check that mwaves=3, because this solver combines waves at u(i)
c
      if (mwaves .NE. 3) then
         write(6,*) 'need mwaves=3 for this Riemann solver'
         stop
         endif
c
c     # set mu to point to  the component of the system that corresponds
c     # to momentum in the direction of this slice, mv and mw to the
c     # orthogonal momentums:
c
      if (ixyz .EQ. 1) then
	  mu = 2
	  mv = 3
          mw = 4
      else if (ixyz .EQ. 2) then
	  mu = 3
	  mv = 4
          mw = 2
      else
          mu = 4
          mv = 2
          mw = 3
      endif
c
c     # note that notation for u,v, and w reflects assumption that the
c     # Riemann problems are in the x-direction with u in the normal
c     # direction and v and w in the orthogonal directions, but with the
c     # above definitions of mu, mv, and mw the routine also works with
c     # ixyz=2 and ixyz = 3
c     # and returns, for example, f0 as the Godunov flux g0 for the
c     # Riemann problems u_t + g(u)_y = 0 in the y-direction.
c
c     # Compute the Roe-averaged variables needed in the Roe solver.
c     # These are stored in the common block comroe since they are
c     # later used in routine rpt3eu to do the transverse wave
c     # splitting.
c
c     # Loop over rest ...
c
      do 20 i = 2-mbc, mx+mbc
c        # Check for negative densities, however unlikely...
         if (qr(1,i-1) .le. 0.d0 .or. ql(1,i) .le. 0.d0) then
            write(*,*) i, mu, mv, mw
            write(*,990) (qr(j,i-1),j=1,5)
            write(*,990) (ql(j,i),j=1,5)
 990        format(5e12.4)
            if (ixyz .eq. 1) 
     &         write(6,*) '*** rho .le. 0 in x-sweep at ',i,jcom,kcom
            if (ixyz .eq. 2) 
     &         write(6,*) '*** rho .le. 0 in y-sweep at ',icom,i,kcom
            if (ixyz .eq. 3) 
     &         write(6,*) '*** rho .le. 0 in z-sweep at ',icom,jcom,i
            write(6,*) 'stopped with rho < 0...'
            stop
         endif
c        # Also set dtdx for flux limiters
         dx = dxcom*merge(1.d0,0.d0,(ixyz.EQ.1)) +
     &        dycom*merge(1.d0,0.d0,(ixyz.EQ.2)) +
     &        dzcom*merge(1.d0,0.d0,(ixyz.EQ.3))
         dtdx(i) = (dtcom/dx) ! Assume no capacity function, uniform
	 rhsqrtl = dsqrt(qr(1,i-1))
	 rhsqrtr = dsqrt(ql(1,i))
	 rhsq2   = rhsqrtl + rhsqrtr
         gamma   = (rhsqrtr*auxl(19,i)+rhsqrtl*auxr(19,i-1))/rhsq2
         gamma1  = gamma-1.d0
	 pl = (auxr(19,i-1)-1.d0)*(qr(5,i-1) - 0.5d0*(qr(mu,i-1)**2 +
     &	       qr(mv,i-1)**2 + qr(mw,i-1)**2)/qr(1,i-1))
	 pr = (auxl(19,i)-1.d0)*(ql(5,i) - 0.5d0*(ql(mu,i)**2 +
     &	       ql(mv,i)**2 + ql(mw,i)**2)/ql(1,i))
	 u(i) = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2 
	 v(i) = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
	 w(i) = (qr(mw,i-1)/rhsqrtl + ql(mw,i)/rhsqrtr) / rhsq2
	 enth(i) = (((qr(5,i-1)+pl)/rhsqrtl
     &	         + (ql(5,i)+pr)/rhsqrtr)) / rhsq2
	 u2v2w2(i) = u(i)**2 + v(i)**2 + w(i)**2
         a2 = gamma1*(enth(i) - .5d0*u2v2w2(i))
c        # Further error checking, also unlikely...
         if (a2 .le. 0.d0) then
            if (ixyz .eq. 1) 
     &         write(6,*) '*** a2 .le. 0 in x-sweep at ',i,jcom,kcom
            if (ixyz .eq. 2) 
     &         write(6,*) '*** a2 .le. 0 in y-sweep at ',icom,i,kcom
            if (ixyz .eq. 3) 
     &         write(6,*) '*** a2 .le. 0 in z-sweep at ',icom,jcom,i
            write(6,*) 'stopped with a2 < 0...'
            stop
         endif
         a(i) = dsqrt(a2)
	 g1a2(i) = gamma1 / a2
	 euv(i) = enth(i) - u2v2w2(i)
c
c
c     # Now split the jump in q1d at each interface into waves
c
c     # Find b1 thru b5, the coefficients of the 5 eigenvectors:
         ur = ql(mu,i)/ql(1,i)
         ul = qr(mu,i-1)/qr(1,i-1)
         delta(1,i) = ql(mu,i) - qr(mu,i-1)
         delta(2,i) = (ql(mu,i)*ur + pr) - (qr(mu,i-1)*ul + pl)
         delta(3,i) = ur*ql(mv,i) - ul*qr(mv,i-1)
         delta(4,i) = ur*ql(mw,i) - ul*qr(mw,i-1)
         delta(5,i) = (ql(5,i)+pr)*ur - (qr(5,i-1)+pl)*ul
c     # Apply gravity for vertical direction ...
         if (ixyz.eq.3) then
c        # Based on cell centers
            delta(2,i) = delta(2,i) - 
     &                0.5d0*(ql(1,i)+qr(1,i-1))*auxl(17,i)
            delta(5,i) = delta(5,i) -
     &                0.5d0*(ql(mu,i)+qr(mu,i-1))*auxl(17,i)
         endif
         b4 = g1a2(i) * (euv(i)*delta(1,i)
     &      + u(i)*delta(2,i) + v(i)*delta(3,i) + w(i)*delta(4,i)
     &      - delta(5,i))
         b2 = delta(3,i) - v(i)*delta(1,i)
         b3 = delta(4,i) - w(i)*delta(1,i)
         b5 = (delta(2,i) + (a(i)-u(i))*delta(1,i)  
     &      - a(i)*b4) / (2.d0*a(i))
         b1 = delta(1,i) - b4 - b5
c
c     # Compute the waves.
c     # Note that the 2-, 3- and 4-wave travel at the same speed,
c     # and are combined. They are limited separately, however -- 
c     # This is a critical correction for stability.
c     # The 5-wave is stored in fwave(.,.,3). 
c     # Tracers are in fwave(.,.,2).
c
c     # Acoustic:
         fwave(1,1,i)  = b1
         fwave(mu,1,i) = b1*(u(i)-a(i))
         fwave(mv,1,i) = b1*v(i)
         fwave(mw,1,i) = b1*w(i)
         fwave(5,1,i)  = b1*(enth(i) - u(i)*a(i))
         s(1,i) = u(i)-a(i)
c
c     # Shear + Entropy (+Tracers):
         fwave(1,2,i)  = b4
         fwave(mu,2,i) = b4*u(i)
         fwave(mv,2,i) = b2+b4*v(i)
         fwave(mw,2,i) = b3+b4*w(i)
         fwave(5,2,i)  = b2*v(i)+b3*w(i)+b4*0.5d0*u2v2w2(i)
         s(2,i) = u(i)
c
c     # Acoustic:
         fwave(1,3,i)  = b5
         fwave(mu,3,i) = b5*(u(i)+a(i))
         fwave(mv,3,i) = b5*v(i)
         fwave(mw,3,i) = b5*w(i)
         fwave(5,3,i)  = b5*(enth(i)+u(i)*a(i))
         s(3,i) = u(i)+a(i)
c
c     # Additional nspec waves added for tracers or mixing ratios;
c     # assume conserved/advected species coupled at Roe velocity.
c     # but, if not...
         if (meqn.EQ.5) goto 20
c     # Specify fluxes/fwaves for meqn-5 species:
         do 15 m = 6, meqn
c        # Total tracers = meqn-5
c        # Total advected tracers = nspecadv
c        # Total conserved tracers = meqn-5-nspecadv
c        # Set neighboring fwaves to zero
            fwave(m,1,i) = 0.d0
            fwave(m,3,i) = 0.d0
c        # Select conservative vs. pure advection
            if ((m-5) .LE. nspecadv) then
c           # Set fwave without velocity differencing
               fwave(m,2,i) = s(2,i)*(ql(m,i) - qr(m,i-1))
            else
c           # Set flux difference at cell-center velocity:
               fwave(m,2,i) = ur*ql(m,i)-ul*qr(m,i-1)
            endif
   15    continue
c
   20 continue
c
c
c     # compute flux differences amdq and apdq.
c     ---------------------------------------
c
      if (efix) goto 110
c
c     # no entropy fix
c     ----------------
c
c     # amdq = SUM fwave   over left-going waves
c     # apdq = SUM fwave   over right-going waves
c
      do 100 m = 1, meqn
         do 100 i = 2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do 90 mws = 1, 3
               if (s(mws,i) .LT. 0.d0) then
                  amdq(m,i) = amdq(m,i) + fwave(m,mws,i)
               else
                  apdq(m,i) = apdq(m,i) + fwave(m,mws,i)
               endif
   90       continue
  100 continue
      go to 300
c
c     -----------------------------------------------------
c
  110 continue
c
c     # With entropy fix
c     ------------------
c
c    # compute flux differences amdq and apdq.
c    # First compute amdq as sum of fwave for left going waves.
c    # Incorporate entropy fix by adding a modified fraction of wave
c    # if s should change sign.
c
      do 200 i = 2-mbc, mx+mbc
c
c        # check 1-wave:
c        ---------------
c
         gamma  = auxr(19,i-1)
         gamma1 = gamma-1.d0
	 rhoim1 = qr(1,i-1)
	 pim1   = gamma1*(qr(5,i-1) - 0.5d0*(qr(mu,i-1)**2
     &          + qr(mv,i-1)**2 + qr(mw,i-1)**2) / rhoim1)
	 cim1   = dsqrt(gamma*pim1/rhoim1)
	 s0     = qr(mu,i-1)/rhoim1 - cim1     !# u-c in left state (cell i-1)
c
c
c     # check for fully supersonic case:
         if (s0.GE.0.d0 .AND. s(1,i).GT.0.d0)then
c        # everything is right-going
            do 115 m = 1, meqn
               amdq(m,i) = 0.d0
  115       continue
            go to 200
         endif
c
         rho1   = qr(1,i-1) + fwave(1,1,i)/s(1,i)
         rhou1  = qr(mu,i-1) + fwave(mu,1,i)/s(1,i)
         rhov1  = qr(mv,i-1) + fwave(mv,1,i)/s(1,i)
         rhow1  = qr(mw,i-1) + fwave(mw,1,i)/s(1,i)
         en1    = qr(5,i-1) + fwave(5,1,i)/s(1,i)
         p1     = gamma1*(en1 - 0.5d0*(rhou1**2 + rhov1**2 +
     &                rhow1**2)/rho1)
         c1     = dsqrt(gamma*p1/rho1)
         s1     = rhou1/rho1 - c1  !# u-c to right of 1-wave
         if (s0.LT.0.d0 .AND. s1.GT.0.d0) then
c        # transonic rarefaction in the 1-wave
	    sfract = s0 * (s1-s(1,i)) / (s1-s0)
	 else if (s(1,i) .LT. 0.d0) then
c	 # 1-wave is leftgoing
	    sfract = s(1,i)
         else
c	 # 1-wave is rightgoing
            sfract = 0.d0   !# this should not happen since s0 < 0
         endif
         do 120 m = 1, meqn
            amdq(m,i) = sfract*fwave(m,1,i)/s(1,i)
  120    continue
c
c        # check 2-wave:
c        ---------------
c
         if (s(2,i) .GE. 0.d0) go to 200  !# 2-,3- and 4- waves are rightgoing
         do 140 m = 1, meqn
            amdq(m,i) = amdq(m,i) + fwave(m,2,i)
  140    continue
c
c        # check 5-wave:
c        ---------------
c
         gamma  = auxl(19,i)
         gamma1 = gamma-1.d0
         rhoi   = ql(1,i)
         pi     = gamma1*(ql(5,i) - 0.5d0*(ql(mu,i)**2
     &          + ql(mv,i)**2 + ql(mw,i)**2) / rhoi)
         ci     = dsqrt(gamma*pi/rhoi)
         s3     = ql(mu,i)/rhoi + ci     !# u+c in right state  (cell i)
c
         rho2   = ql(1,i) - fwave(1,3,i)/s(3,i)
         rhou2  = ql(mu,i) - fwave(mu,3,i)/s(3,i)
         rhov2  = ql(mv,i) - fwave(mv,3,i)/s(3,i)
         rhow2  = ql(mw,i) - fwave(mw,3,i)/s(3,i)
         en2    = ql(5,i) - fwave(5,3,i)/s(3,i)
         p2     = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2 +
     &                    rhow2**2)/rho2)
         c2     = dsqrt(gamma*p2/rho2)
         s2     = rhou2/rho2 + c2   !# u+c to left of 5-wave
         if (s2 .LT. 0.d0 .AND. s3.GT.0.d0 ) then
c        # transonic rarefaction in the 5-wave
            sfract = s2 * (s3-s(3,i)) / (s3-s2) 
         else if (s(3,i) .LT. 0.d0) then
c        # 5-wave is leftgoing
            sfract = s(3,i)
         else
c        # 5-wave is rightgoing
            go to 200
         endif
c
         do 160 m = 1, meqn
            amdq(m,i) = amdq(m,i) + sfract*fwave(m,3,i)/s(3,i)
  160    continue
  200 continue
c
c     # compute the rightgoing flux differences:
c     # df = SUM s*wave   is the total flux difference and apdq = df - amdq
c
      do 220 m = 1, meqn
         do 220 i = 2-mbc, mx+mbc
            df = 0.d0
            do 210 mws = 1, 3
               df = df + fwave(m,mws,i)
  210       continue
         apdq(m,i) = df - amdq(m,i)
  220 continue
c
c
c     # Calculate local and nonlocal (upstream) CFL and upwind Id
c     # (... may as well do both, even if CFL not needed)
c
      do 240 i = 0, mx+1
         do 240 mws = 1, 3
            Id(mws,i) = int(i-dsign(1.d0,s(mws,i)))
            CFL(mws,i) = dabs(s(mws,i)*dtdx(max(i,Id(mws,i))))
  240 continue
c
c
  300 continue
c
      if (meqn.EQ.5) goto 325
c
c     # Apply flux limiters to species individually:
c     # Reduce compute time by sharing wave 3; ensure
c     # stability by separating limiter -- must also set
c     # mthlim = 0 or limiters will be applied twice.
c     ---------------------------------------------------
c
      wlimitr = 0.d0
      do 320 m = 6, meqn
         fwavec = 0.d0
         do 310 i = 0, mx+1
            fwavel = fwavec
            fwavec = fwave(m,2,i)
            fwaver = fwave(m,2,i+1)
            if (i .EQ. 0) goto 310
            if (fwavec .EQ. 0.d0) goto 310
            if (s(2,i) .GT. 0.d0) then
               r = fwavel/fwavec
            else
               r = fwaver/fwavec
            endif
c           # Apply limiter to fwave
            if ((m-5).LE.nspecadv) then
c           # Advection Equation gets same limiter as Entropy Wave
               if (thirdorder) then
c              # CFL-Superbee Theta = 1/2, Beta = 2/3
                  wlimitr = cflsupbee(0.5d0,beta,r,
     &                                CFL(2,i),CFL(2,Id(2,i)))
               else
c              # MC Limiter (CFL independent, Good Choice)
                  c = (1.d0 + 1.d0*r)/2.d0
                  wlimitr = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r))
               endif
            else
c           # Continuity Equation gets minmod, as a safer choice
               wlimitr = dmax1(0.d0,dmin1(1.d0, r))
            endif
c           # Apply limiter to fwave
            fwave(m,2,i) = wlimitr * fwave(m,2,i)
  310    continue
  320 continue
c
  325 continue
c
c     # Apply Lax-Liu style limiters to other fwaves:
c     ---------------------------------------------------
c
c     Taking full advantage of the generosity of F77 compilers,
c     we use if-then statements to define the limiting scheme.
c     As these depend only on local parameter data, the compiler
c     will eliminate them automatically as "pre-processing".
c     So, there is NO PENALTY to this approach in performance,
c     and it obviates the need for a real pre-processor!
c     (Just be sure to enable optimization; -O3 works well.)
c
      wlimiter1 = 0.d0
      wlimiter2 = 0.d0
      wlimiter3 = 0.d0
      wlimiter4 = 0.d0
      wlimiter5 = 0.d0
      b1 = 0.d0
      b2 = 0.d0
      b3 = 0.d0
      b4 = 0.d0
      b5 = 0.d0
      do 330 i = 0, mx+1
         if (.NOT.(laxliu)) then
c        # Store left coefficients from previous
            b1l = b1
            b2l = b2
            b3l = b3
            b4l = b4
            b5l = b5
c        # Store right coefficients from fwaves
            b1r = fwave(1,1,i+1)
            b4r = fwave(1,2,i+1)
            b2r = fwave(mv,2,i+1)-b4r*v(i+1)
            b3r = fwave(mw,2,i+1)-b4r*w(i+1)
            b5r = fwave(1,3,i+1)
         endif
c     # Obtain coefficients from fwaves
         b1 = fwave(1,1,i)
         b4 = fwave(1,2,i)
         b2 = fwave(mv,2,i)-b4*v(i)
         b3 = fwave(mw,2,i)-b4*w(i)
         b5 = fwave(1,3,i)
         if (i .EQ. 0) goto 330
c
c     # Acoustic Wave 1
         if (b1 .NE. 0.d0) then
c        # Solve coefficient, find ratio r1
            if (laxliu) then
               b4Id = g1a2(i) * (euv(i)*delta(1,Id(1,i))
     &              + u(i)*delta(2,Id(1,i)) + v(i)*delta(3,Id(1,i)) 
     &              + w(i)*delta(4,Id(1,i)) - delta(5,Id(1,i)))
               b5Id = (delta(2,Id(1,i)) + (a(i)-u(i))*delta(1,Id(1,i)) 
     &              - a(i)*b4Id) / (2.d0*a(i))
               b1Id = delta(1,Id(1,i)) - b4Id - b5Id
            else  !# Classic Clawpack Limiter
               if (s(1,i) .GT. 0.d0) then
                  b1Id = b1l
               else
                  b1Id = b1r
               endif
            endif
            r1 = b1Id/b1
            if (thirdorder) then
c           # SuperPower Limiter for Wave 1:
               if (nonlocal) then
                  CFLup=CFL(1,Id(1,i))
               else
                  CFLup=CFL(1,i)
               endif
               wlimiter1 = suppow(r1,CFL(1,i),CFLup)
            else
c           # MC Limiter (CFL independent, Good Choice)
               c = (1.d0 + 1.d0*r1)/2.d0
               wlimiter1 = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r1))
            endif
         endif
c
c     # Shear Wave 2
         if (b2 .NE. 0.d0) then
c        # Solve coefficient, find ratio r2
            if (laxliu) then
               b2Id = delta(3,Id(2,i)) - v(i)*delta(1,Id(2,i))
            else  !# Classic Clawpack Limiter
               if (s(2,i) .GT. 0.d0) then
                  b2Id = b2l
               else
                  b2Id = b2r
               endif
            endif
            r2 = b2Id/b2
            if (thirdorder) then
c           # CFL-Superbee beta=2/3 Limiter for Wave 2:
               if (nonlocal) then
                  CFLup=CFL(2,Id(2,i))
               else
                  CFLup=CFL(2,i)
               endif
               wlimiter2 = cflsupbee(theta,beta,r2,
     &                               CFL(2,i),CFLup)
            else
c           # MC Limiter (CFL independent, Good Choice)
               c = (1.d0 + 1.d0*r2)/2.d0
               wlimiter2 = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r2))
            endif
         endif
c
c     # Shear Wave 3
         if (b3 .NE. 0.d0) then
c        # Solve coefficient, find ratio r3
            if (laxliu) then
               b3Id = delta(4,Id(2,i)) - w(i)*delta(1,Id(2,i))
            else  !# Classic Clawpack Limiter
               if (s(2,i) .GT. 0.d0) then
                  b3Id = b3l
               else
                  b3Id = b3r
               endif
            endif
            r3 = b3Id/b3
            if (thirdorder) then
c           # CFL-Superbee beta=2/3 Limiter for Wave 3:
               wlimiter3 = cflsupbee(theta,beta,r3,
     &                               CFL(2,i),CFLup) !CFLup already set
            else
c           # MC Limiter (CFL independent, Good Choice)
               c = (1.d0 + 1.d0*r3)/2.d0
               wlimiter3 = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r3))
            endif
         endif
c
c     # Entropy Wave 4
         if (b4 .NE. 0.d0) then
c        # Solve coefficient, find ratio r4
            if (laxliu) then
               b4Id = g1a2(i) * (euv(i)*delta(1,Id(2,i))
     &              + u(i)*delta(2,Id(2,i)) + v(i)*delta(3,Id(2,i))
     &              + w(i)*delta(4,Id(2,i)) - delta(5,Id(2,i)))
            else  !# Classic Clawpack Limiter
               if (s(2,i) .GT. 0.d0) then
                  b4Id = b4l
               else
                  b4Id = b4r
               endif
            endif
            r4 = b4Id/b4
            if (thirdorder) then
c           # CFL-Superbee beta=2/3 Limiter for Wave 4:
               wlimiter4 = cflsupbee(theta,beta,r4,
     &                               CFL(2,i),CFLup) !CFLup already set
            else
c           # MC Limiter (CFL independent, Good Choice)
               c = (1.d0 + 1.d0*r4)/2.d0
               wlimiter4 = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r4))
            endif
         endif
c
c     # Acoustic Wave 5
         if (b5 .NE. 0.d0) then
c        # Solve coefficient, find ratio r5
            if (laxliu) then
               b4Id = g1a2(i) * (euv(i)*delta(1,Id(3,i))
     &              + u(i)*delta(2,Id(3,i)) + v(i)*delta(3,Id(3,i))
     &              + w(i)*delta(4,Id(3,i)) - delta(5,Id(3,i)))
               b5Id = (delta(2,Id(3,i)) + (a(i)-u(i))*delta(1,Id(3,i)) 
     &              - a(i)*b4Id) / (2.d0*a(i))
            else  !# Classic Clawpack Limiter
               if (s(3,i) .GT. 0.d0) then
                  b5Id = b5l
               else
                  b5Id = b5r
               endif
            endif      
            r5 = b5Id/b5
            if (thirdorder) then
c           # CFL-Superbee beta=2/3 Limiter for Wave 5:
               if (nonlocal) then
                  CFLup=CFL(3,Id(3,i))
               else
                  CFLup=CFL(3,i)
               endif
               wlimiter5 = suppow(r5,CFL(3,i),CFLup)
            else
c           # MC Limiter (CFL independent, Good Choice)
               c = (1.d0 + 1.d0*r5)/2.d0
               wlimiter5 = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r5))
            endif
         endif
c        
c     # Recompute the limited fwaves, collapsing into three ...
c
c     # Acoustic Waves 1 and 5:
         do m = 1, 5
            fwave(m,1,i) = fwave(m,1,i)*wlimiter1
            fwave(m,3,i) = fwave(m,3,i)*wlimiter5
         enddo
c     # Shear + Entropy Waves 2, 3, and 4:
         b2ltd = wlimiter2*b2
         b3ltd = wlimiter3*b3
         b4ltd = wlimiter4*b4
         fwave(1,2,i)  = b4ltd
         fwave(mu,2,i) = b4ltd*u(i)
         fwave(mv,2,i) = b2ltd+b4ltd*v(i)
         fwave(mw,2,i) = b3ltd+b4ltd*w(i)
         fwave(5,2,i)  = b2ltd*v(i)+b3ltd*w(i) + b4ltd*0.5d0*u2v2w2(i)
c
  330 continue

  500 continue

      return
      end
c
c     # Limiter functions are pasted below for convenience
c
c     =====================================================
      double precision function cflsupbee(theta,beta,r,
     &                                    CFL,CFLup)
c     =====================================================

      implicit double precision (a-h,o-z)

c     # CFL-Superbee Limiter:
      CFLe = dmax1(0.001d0,dmin1(0.999d0,CFL))
      CFLupe = dmax1(0.001d0,dmin1(0.999d0,CFLup))
      weight = (1.d0 - CFLupe) / (1.d0 - CFLe)
      s1 = theta*weight*2.d0/CFLe
      s2 = (1.d0 + CFLe) / 3.d0
      phimax = theta * 2.d0 / (1.d0 - CFLe)
      ultra = dmax1(0.d0,dmin1(s1*r,phimax))
      c1 = 1.d0 + (s2 - beta/2.d0) * (r-1.d0)
      c2 = 1.d0 + (s2 + beta/2.d0) * (r-1.d0)
      cflsupbee = dmax1(0.d0, dmin1(ultra,dmax1(c1,c2)))

      return
      end
c
c
c     =====================================================
      double precision function suppow(r,CFL,CFLup)
c     =====================================================

      implicit double precision (a-h,o-z)

c     # Superpower Limiter:
      CFLe = dmax1(0.001d0,dmin1(0.999d0,CFL))
      CFLupe = dmax1(0.001d0,dmin1(0.999d0,CFLup))
      s2 = (1.d0 + CFLe) / 3.d0
      s3 = 1.d0 - s2
      weight = (1.d0 - CFLupe) / (1.d0 - CFLe)
      if (r.LE.1.d0) then
         pp = weight*(2.d0/CFLe)*2.d0*s3
      else
         pp = dabs(2.d0/(1.d0-CFLe))*2.d0*s2
      endif
      rabs = dabs(r)
      rfrac = dabs((1.d0-rabs)/(1.d0+rabs))
      signum = dmax1(0.d0,dsign(1.d0,r))
      suppow = signum * (s3+s2*r) * (1.d0-rfrac**pp)

      return
      end
