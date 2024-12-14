c     ============================================
      subroutine b4step3(mbc,mx,my,mz,meqn,q,
     &            xlower,ylower,zlower,dx,dy,dz,t,dt,maux,aux)
c     ============================================
c
c
c  --------------------------------------------------------
c
c     MAGIC state updates, J. B. Snively;
c     based on work by Pineyro [MS Thesis, 2018]
c        2018-2019 
c
c        Subject to terms in the LICENSE file.
c
c  --------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
      dimension    q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc)
      dimension  aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 
     &               1-mbc:mz+mbc)
c
      logical varyR
      logical varyGamma
      logical varyMu
      data varyR     /.false./
      data varyGamma /.false./
      data varyMu    /.false./
c
c     # Escape this routine if there are no tracers ...
c
      if (meqn.EQ.5) goto 999
c     
      do 15 k = 1-mbc, mz+mbc
        do 15 j = 1-mbc, my+mbc
          do 15 i = 1-mbc, mx+mbc
c
c        # Constrain Mass Fractions
c        # Yes*, we force these to be bounded and nonzero;
c        # sometimes 0.d0 is okay, though. 
c        # (* for theta=0.5 or MC limiter, the error is negligible)
            do 5 m=6,meqn
               q(m,i,j,k) = dmax1(0.d0,dmin1(q(m,i,j,k),1.d0))
    5       continue
c
c         # Calculate Number Densities
            if (varyR.OR.varyGamma.OR.varyMu) then
              dox   = 6.022d23*(q(6,i,j,k)/16.d0)*q(1,i,j,k)*1.d-3
              dox2  = 6.022d23*(q(7,i,j,k)/32.d0)*q(1,i,j,k)*1.d-3
c             # Hope the following is never needed, but... just in case
              dnit2 = dmax1(0.d0,(1.d0-q(6,i,j,k)-q(7,i,j,k)))
              dnit2 = 6.022d23*(dnit2/28.d0)
     &              * q(1,i,j,k)*1.d-3
c             # Alternate approach: Rewrite to combine O2+N2=Major Gas
c             # This alternative is robust and straightforward to code.
            end if
c
c         # Specific Gas Contant R
            if (varyR) then
              aux(18,i,j,k) = 8.31d3 /
     &          ((dox*16.d0+dnit2*28.d0+dox2*32.d0) /
     &           (dox+dnit2+dox2))
            end if
c
c         # Ratio of Specific Heats Gamma
c         # Calculate similar to Hickey and Walterscheid ...
            if (varyGamma) then
              aux(19,i,j,k) = (7*(dox2+dnit2)+5*dox)/
     &                        (5*(dox2+dnit2)+3*dox)
            end if
c     
c         # Molecular Diffusivities / Conductivities
c         # Calculate from MSIS Profile
c         # [Rees, 1989]
            if (varyMu) then
              temp = (aux(19,i,j,k)-1.d0)*(q(5,i,j,k)-
     &          0.5d0*(q(2,i,j,k)**2+q(3,i,j,k)**2+q(4,i,j,k)**2)
     &          /q(1,i,j,k)) / (q(1,i,j,k)*aux(18,i,j,k))
              Cp = 1.d0/(1.d0-1.d0/aux(19,i,j,k))*aux(18,i,j,k)
              dmolecu = (3.43*dnit2+4.03*dox2+3.90*dox)*1.d-7
              dmolecu = dmolecu*(temp**(0.69))/(dox+dnit2+dox2)
              dtherm = (56.0*(dnit2+dox2)+75.9*dox)*1.d-5
              dtherm = dtherm*(temp**(0.69))/(dox+dnit2+dox2)
              aux(11,i,j,k) = dmolecu*(1.d0/q(1,i,j,k))
              aux(12,i,j,k) = (dtherm*(1.d0/(q(1,i,j,k)*Cp)))
     &                 * aux(19,i,j,k)
            end if
c
  15  continue
c
  999 return
      end
