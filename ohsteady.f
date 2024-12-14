c
c
c     =======================================================
      subroutine ohsteady(meqn,mbc,mx,my,mz,
     &     xlower,ylower,zlower,dx,dy,dz,q,maux,aux,dtoh)
c     =======================================================
c
c  --------------------------------------------------------
c
c     Hydroxyl Solution for Atmospheric Model Extensions
c        by J. B. Snively, 2003-2019
c        Subject to terms in the LICENSE file.
c
c  --------------------------------------------------------
c
c
      implicit double precision (a-h,o-z)
c
      dimension    q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc)
      dimension    aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc)
      dimension b(0:9),a8(0:9),a9(1:9,0:8),suma9(1:9),
     &          a10(0:9),aa(1:9,1:6)
c
c     # Coefficients below are from Adler-Golden [1997],
c     # and Turnbull and Lowe [1989].
c
      data b/0,0,0,0,0,0,0.08,0.17,0.27,0.48/
      data a8/3.9,10.5,25,25,25,25,25,25,25,25/
      data ((a9(i,j),j=0,8),i=1,9)/
     *2,0,0,0,0,0,0,0,0,
     *0,4,0,0,0,0,0,0,0,
     *0,1,7,0,0,0,0,0,0,
     *0,1,2,10,0,0,0,0,0,
     *0,1,2,6,16,0,0,0,0,
     *1,1,3,6,11,22,0,0,0,
     *4,6,9,12,16,23,32,0,0,
     *4,6,8,10,14,19,25,33,0,
     *28,29,31,32,34,36,38,40,42/
      data suma9/2,4,8,13,25,44,102,119,310/
      data a10/0,0.58,1.0,1.7,3.0,5.2,9.1,16,70,48/
      data ((aa(i,j),j=1,6),i=1,9)/
     *22.74,0,0,0,0,0,
     *30.43,15.42,0,0,0,0,
     *28.12,40.33,2.032,0,0,0,
     *20.30,69.77,7.191,0.299,0,0,
     *11.05,99.42,15.88,1.315,0.051,0,
     *4,125.6,27.94,3.479,0.274,0.010,
     *2.34,145.1,42.91,7.165,0.847,0.063,
     *8.6,154.3,59.98,12.68,2.007,0.230,
     *23.72,148.9,78.64,19.94,4.053,0.620/
c
c     Find partial steady-state solution for
c     OH Airglow vibrational states v=0-9
c     [e.g., Adler-Golden, 1997]
c
c     Balance is imposed for H and O3, although
c     they are evolved time-dependently. This is 
c     appropriate for typical time-scales ~minutes.
c     [Snively et al., 2010]
c
      if (dtoh.LE.0) then ! Negative for restart, otherwise set [O3]
        hydmax = 0.d0
        kmin = int(80*(1000/dz))
        kmax = int(90*(1000/dz))
        do k=kmin,kmax
          dox   = 6.022d23*(aux(6,1,1,k)/16.d0)*aux(1,1,1,k)*1.d-3
          dnit2 = 6.022d23*((1.d0-aux(6,1,1,k)-aux(7,1,1,k))
     &            /28.d0)*aux(1,1,1,k)*1.d-3
          dox2  = 6.022d23*(aux(7,1,1,k)/32.d0)*aux(1,1,1,k)*1.d-3
          dhyd  = aux(9,1,1,k)*(dox+dox2+dnit2) 
          if (dhyd.GT.hydmax) then 
            hydmax  = dhyd
            hydzmax = zlower + (k-0.5d0)*dz
          endif
        enddo
      endif
      kmin = int(75*(1000/dz)) 
      kmax = int(125*(1000/dz))
      do k=kmin,kmax
        zcell = zlower + (k-0.5d0)*dz
        do j=1,my
          do i=1,mx
            if (dtoh.GT.0) then ! Proceed with a time step...
              do m=6,10         ! Bound the tracer variables
                 q(m,i,j,k) = dmax1(0.d0,dmin1(q(m,i,j,k),1.d0))
              enddo
              dox   = 6.022d23*(q(6,i,j,k)/16.d0)*q(1,i,j,k)*1.d-3
              dnit2 = 6.022d23*((1.d0-q(6,i,j,k)-q(7,i,j,k))
     &                /28.d0)*q(1,i,j,k)*1.d-3
              dox2  = 6.022d23*(q(7,i,j,k)/32.d0)*q(1,i,j,k)*1.d-3   
              dhyd  = q(9,i,j,k)*(dox+dox2+dnit2)
              dox3  = q(10,i,j,k)*(dox+dox2+dnit2)
              gamma = aux(19,i,j,k)
              gamma1= gamma-1.d0
c           # Convert energy to temperature
              tempk = gamma1*(q(5,i,j,k)-0.5d0*(q(2,i,j,k)**2+
     &           q(3,i,j,k)**2+q(4,i,j,k)**2)/q(1,i,j,k))/
     &           (q(1,i,j,k)*aux(18,i,j,k))
            else ! Use aux values to update data, in case restart...
              dox   = 6.022d23*(aux(6,i,j,k)/16.d0)*aux(1,i,j,k)*1.d-3
              dnit2 = 6.022d23*((1.d0-aux(6,i,j,k)-aux(7,i,j,k))
     &                /28.d0)*aux(1,i,j,k)*1.d-3
              dox2  = 6.022d23*(aux(7,i,j,k)/32.d0)*aux(1,i,j,k)*1.d-3
              dhyd  = aux(9,i,j,k)*(dox+dox2+dnit2)
              dox3  = aux(10,i,j,k)*(dox+dox2+dnit2)
              gamma = aux(19,i,j,k)
              gamma1= gamma-1.d0
c           # Convert energy to temperature
              tempk = gamma1*(aux(5,i,j,k)-0.5d0*(aux(2,i,j,k)**2+
     &           aux(3,i,j,k)**2+aux(4,i,j,k)**2)/aux(1,i,j,k))/
     &           (aux(1,i,j,k)*aux(18,i,j,k))
            endif
            rk1   = (1.4d-10)*dexp(-470/tempk)
            rk3N2 = (5.7d-34)*((300/tempk)**(2.62))
            rk3O2 = (5.96d-34)*((300/tempk)**(2.37))
            rk8 = 1d-11
            rk9 = 1d-13
            rk10= 1d-14
c           # First Step only...
            if (dtoh.LE.0) then ! Initializing or restarting...
c           # Approx. Ozone Steady-State Profile
              dhyd1 = hydmax*dexp(-(zcell-hydzmax)/8620)
              rL = dhyd1*rk1
              rP = dox*dox2*(dox2*rk3O2+dnit2*rk3N2)
              dox3 = rP/rL
              if (dtoh.EQ.0) then ! Set q because not a restart...
                aux(10,i,j,k) = dox3/(dox+dox2+dnit2)
                q(10,i,j,k)   = dox3/(dox+dox2+dnit2)
              else                ! Don't set q, because restart...
                aux(10,i,j,k) = dox3/(dox+dox2+dnit2) 
              endif
            endif
c        # Calculate [OH(v=9)]
            rP=rk1*b(9)*dhyd*dox3
            rL=rk8*a8(9)*dox+rk9*suma9(9)*dox2+rk10*a10(9)*dnit2
     &         +aa(9,1)+aa(9,2)+aa(9,3)+aa(9,4)+aa(9,5)+aa(9,6)
            aux(29,i,j,k)=rP/rL
c        # Calculate [OH(v=8)]
            rP=rk1*b(8)*dhyd*dox3+aux(29,i,j,k)*aa(9,1)
     &         +rk9*(a9(9,8)*aux(29,i,j,k))*dox2
     &         +rk10*a10(9)*dnit2*aux(29,i,j,k)
            rL=rk8*a8(8)*dox+rk9*suma9(8)*dox2+rk10*a10(8)*dnit2
     &         +aa(8,1)+aa(8,2)+aa(8,3)+aa(8,4)+aa(8,5)+aa(8,6)
            aux(28,i,j,k)=rP/rL
c        # Calculate [OH(v=7)]
            rP=rk1*b(7)*dhyd*dox3+aux(29,i,j,k)*aa(9,2)
     &         +aux(28,i,j,k)*aa(8,1)+rk9*(a9(9,7)*aux(29,i,j,k)
     &         +a9(8,7)*aux(28,i,j,k))*dox2
     &         +rk10*a10(8)*dnit2*aux(28,i,j,k)
            rL=rk8*a8(7)*dox+rk9*suma9(7)*dox2+rk10*a10(7)*dnit2
     &         +aa(7,1)+aa(7,2)+aa(7,3)+aa(7,4)+aa(7,5)+aa(7,6)
            aux(27,i,j,k)=rP/rL
c        # Calculate [OH(v=6)]
            rP=rk1*b(6)*dhyd*dox3+aux(29,i,j,k)*aa(9,3)
     &         +aux(28,i,j,k)*aa(8,2)+aux(27,i,j,k)*aa(7,1)
     &         +rk9*(a9(9,6)*aux(29,i,j,k)+a9(8,6)*aux(28,i,j,k)
     &         +a9(7,6)*aux(27,i,j,k))*dox2
     &         +rk10*a10(7)*dnit2*aux(27,i,j,k)
            rL=rk8*a8(6)*dox+rk9*suma9(6)*dox2+rk10*a10(6)*dnit2
     &         +aa(6,1)+aa(6,2)+aa(6,3)+aa(6,4)+aa(6,5)+aa(6,6)
            aux(26,i,j,k)=rP/rL
c        # Calculate [OH(v=5)]
            rP=aux(29,i,j,k)*aa(9,4)+aux(28,i,j,k)*aa(8,3)
     &         +aux(27,i,j,k)*aa(7,2)+aux(26,i,j,k)*aa(6,1)
     &         +rk9*(a9(9,5)*aux(29,i,j,k)+a9(8,5)*aux(28,i,j,k)
     &         +a9(7,5)*aux(27,i,j,k)+a9(6,5)*aux(26,i,j,k))*dox2
     &         +rk10*a10(6)*dnit2*aux(26,i,j,k)  
            rL=rk8*a8(5)*dox+rk9*suma9(5)*dox2+rk10*a10(5)*dnit2
     &         +aa(5,1)+aa(5,2)+aa(5,3)+aa(5,4)+aa(5,5)
            aux(25,i,j,k)=rP/rL
c        # Calculate [OH(v=4)]
            rP=aux(29,i,j,k)*aa(9,5)+aux(28,i,j,k)*aa(8,4)
     &         +aux(27,i,j,k)*aa(7,3)+aux(26,i,j,k)*aa(6,2)
     &         +aux(25,i,j,k)*aa(5,1)+rk9*(a9(9,4)*aux(29,i,j,k)
     &         +a9(8,4)*aux(28,i,j,k)+a9(7,4)*aux(27,i,j,k)
     &         +a9(6,4)*aux(26,i,j,k)+a9(5,4+1)*aux(25,i,j,k))*dox2
     &         +rk10*a10(5)*dnit2*aux(25,i,j,k)
            rL=rk8*a8(4)*dox+rk9*suma9(4)*dox2+rk10*a10(4)*dnit2
     &         +aa(4,1)+aa(4,2)+aa(4,3)+aa(4,4)
            aux(24,i,j,k)=rP/rL
c        # Calculate [OH(v=3)]
            rP=aux(29,i,j,k)*aa(9,6)+aux(28,i,j,k)*aa(8,5)
     &         +aux(27,i,j,k)*aa(7,4)+aux(26,i,j,k)*aa(6,3)
     &         +aux(25,i,j,k)*aa(5,2)+aux(24,i,j,k)*aa(4,1)
     &         +rk9*(a9(9,3)*aux(29,i,j,k)+a9(8,3)*aux(28,i,j,k)
     &         +a9(7,3)*aux(27,i,j,k)+a9(6,3)*aux(26,i,j,k)
     &         +a9(5,3)*aux(25,i,j,k)+a9(4,3)*aux(24,i,j,k))*dox2
     &         +rk10*a10(4)*dnit2*aux(24,i,j,k)
            rL=rk8*a8(3)*dox+rk9*suma9(3)*dox2+rk10*a10(3)*dnit2
     &         +aa(3,1)+aa(3,2)+aa(3,3)
            aux(23,i,j,k)=rP/rL
c        # Calculate [OH(v=2)]
            rP=aux(28,i,j,k)*aa(8,6)+aux(27,i,j,k)*aa(7,5)
     &         +aux(26,i,j,k)*aa(6,4)+aux(25,i,j,k)*aa(5,3)
     &         +aux(24,i,j,k)*aa(4,2)+aux(23,i,j,k)*aa(3,1)
     &         +rk9*(a9(9,2)*aux(29,i,j,k)+a9(8,2)*aux(28,i,j,k)
     &         +a9(7,2)*aux(27,i,j,k)+a9(6,2)*aux(26,i,j,k)
     &         +a9(5,2)*aux(25,i,j,k)+a9(4,2)*aux(24,i,j,k)
     &         +a9(3,2)*aux(23,i,j,k))*dox2
     &         +rk10*a10(3)*dnit2*aux(23,i,j,k)
            rL=rk8*a8(2)*dox+rk9*suma9(2)*dox2
     &         +rk10*a10(2)*dnit2+(aa(2,1)+aa(2,2))
            aux(22,i,j,k)=rP/rL
c        # Calculate [OH(v=1)]
            rP=aux(27,i,j,k)*aa(7,6)+aux(26,i,j,k)*aa(6,5)
     &         +aux(25,i,j,k)*aa(5,4)+aux(24,i,j,k)*aa(4,3)
     &         +aux(23,i,j,k)*aa(3,2)+aux(22,i,j,k)*aa(2,1)
     &         +rk9*(a9(9,1)*aux(29,i,j,k)+a9(8,1)*aux(28,i,j,k)
     &         +a9(7,1)*aux(27,i,j,k)+a9(6,1)*aux(26,i,j,k)
     &         +a9(5,1)*aux(25,i,j,k)+a9(4,1)*aux(24,i,j,k)
     &         +a9(3,1)*aux(23,i,j,k)+a9(2,1)*aux(22,i,j,k))*dox2
     &         +rk10*a10(2)*dnit2*aux(22,i,j,k)
            rL=rk8*a8(1)*dox+rk9*suma9(1)*dox2
     &         +rk10*a10(1)*dnit2+aa(1,1)
            aux(21,i,j,k)=rP/rL
c        # Calculate [OH(v=0)]
            rP=aux(26,i,j,k)*aa(6,6)+aux(25,i,j,k)*aa(5,5)
     &         +aux(24,i,j,k)*aa(4,4)+aux(23,i,j,k)*aa(3,3)
     &         +aux(22,i,j,k)*aa(2,2)+aux(21,i,j,k)*aa(1,1)
     &         +rk9*(a9(9,0)*aux(29,i,j,k)+a9(8,0)*aux(28,i,j,k)
     &         +a9(7,0)*aux(27,i,j,k)+a9(6,0)*aux(26,i,j,k)
     &         +a9(5,0)*aux(25,i,j,k)+a9(4,0)*aux(24,i,j,k)
     &         +a9(3,0)*aux(23,i,j,k)+a9(2,0)*aux(22,i,j,k)
     &         +a9(1,0)*aux(21,i,j,k))*dox2
     &         +rk10*a10(1)*dnit2*aux(21,i,j,k)
            rL=rk8*a8(0)*dox
            aux(20,i,j,k)=rP/rL
c        # Set Hydrogen Steady-State
            rP = 0.d0
c        # Production and loss rates
            do mvib=0,9
               rP = rP+aux(20+mvib,i,j,k)*dox*rk8*a8(mvib)
            enddo
            rL = dhyd*dox3*rk1
            if (dtoh.LE.0) then
c        #  ... define initial tendency ...
c        # Hydrogen
              aux(30,i,j,k) = rP-rL
c        # Ozone
              rP = dox*dox2*(dox2*rk3O2+dnit2*rk3N2)
c        #  ... define initial tendency ...
              aux(31,i,j,k) = rP-rL
            else
c        #  ... time step P and L for H and O3 ...
c        # Hydrogen
              dhyd = dhyd+dtoh*((rP-rL)-aux(30,i,j,k))
              q(9,i,j,k) = dhyd/(dox+dox2+dnit2)
c        # Ozone
              rP   = dox*dox2*(dox2*rk3O2+dnit2*rk3N2)
              dox3 = dox3+dtoh*((rP-rL)-aux(31,i,j,k))
              q(10,i,j,k) = dox3/(dox+dox2+dnit2)
            endif 
          enddo
        enddo
      enddo
c
c
c
      return
      end

