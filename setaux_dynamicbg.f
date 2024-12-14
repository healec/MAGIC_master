c
c
c     ==================================================================
      subroutine setaux(mbc,mx,my,mz,xlower,ylower,
     &                  zlower,dx,dy,dz,maux,iframe,aux)
c     ==================================================================
c
c     --------------------------------------------------------
c
c     Aux. Initialization for MAGIC Extensions
c        by J. B. Snively, 2003-2019
c        Subject to terms in the LICENSE file.
c
c     --------------------------------------------------------
c
      implicit double precision (a-h,o-z)
      dimension  aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc)
c
c     # Provide option to use CGS for compatible profile.data
      logical cgs
      data cgs /.true./
c
      common /param/ gamma, gamma1, grav, dens, pres,
     &                  speed, dtherm0, dmolec0
      common /wind/ omega1, amp0, amp1, prop1z, z1width, z1pos,
     &                  z1Gpos, phi, windscale1
      common /gauss/ xpos, ypos, zpos, xwidth, ywidth,
     &               zwidth, tcenter, twidth
      common /options/ loadprof, nspecadv

      character(16) currentprof
      character(4) lpp
      integer popl
c
c     # MPI : Added so that we know which processor we are running on.
      common /mpi_proc_info/ np, id
c
c
c     # Open profile if desired
c
       popl=9000+iframe
       write(lpp,'(i4)') popl
       currentprof = 'profile'//lpp//'.data'

      if (loadprof.EQ.1) then
         open(unit=7,file=currentprof,
     &            status='old',form='formatted')
      elseif (loadprof.LE.0) then
         goto 100
      endif
c
c     # Initialize auxiliary variables
      do 15 i = 1-mbc, mx+mbc
         xcell = xlower + (i-0.5d0)*dx
         do 15 j = 1-mbc, my+mbc
         ycell = ylower + (j-0.5d0)*dy
         rewind(7)
            do 15 k = 1-mbc, mz+mbc
            zcell = zlower + (k-0.5d0)*dz
c        # Isothermal Atmosphere
c
            if (loadprof.EQ.0) then
               aux(1,i,j,k) = density(zcell)
c           # Analytical wind function (if desired)
               if ((amp0 .EQ. 0) .AND. (amp1 .EQ. 0)) then
                  aux(2,i,j,k) = 0
               else
                  scale=(aux(1,i,j,1)/aux(1,i,j,k))**0.5d0
                  aux(2,i,j,k) = aux(1,i,j,k)*
     &               windfnx(xcell,ycell,zcell,0,scale)
               endif
               aux(3,i,j,k) = 0
               aux(4,i,j,k) = 0
               pres1= pressure(zcell)
               velx = aux(2,i,j,k)/aux(1,i,j,k)
               vely = 0
               velz = 0
               aux(5,i,j,k) = pres1/gamma1+
     &            .5*(velx**2+vely**2+velz**2)*aux(1,i,j,k)
c           # Reserve the following variables
               aux(6,i,j,k) = 0
               aux(7,i,j,k) = 0
               aux(8,i,j,k) = 0
               aux(9,i,j,k) = 0
               aux(10,i,j,k)= 0
c           # Default thermodynamic parameters
               aux(19,i,j,k) = gamma
               aux(18,i,j,k) = 287
c           # Set viscosity / conductivity
               dmolecular=dmolec0*(aux(1,i,j,1)/aux(1,i,j,k))
               dthermal=dtherm0*(aux(1,i,j,1)/aux(1,i,j,k))
               aux(11,i,j,k)=dmolecular
               aux(12,i,j,k)=dthermal*aux(19,i,j,k)
c
            else
c
c        # Otherwise, load vertical profile from profile.data
c
c           # Require exact height alignment
 10            read (7,*) height,dox,dnit2,dox2,aux(1,i,j,k),
     &             aux(5,i,j,k),dhyd,dne,windx,windy
               if (height.lt.(zcell/1.d3)) then
c                 # Commented out to suppress
c                  print*,height,'.LT.',(zcell/1000)
                  goto 10
               elseif (height.gt.(zcell/1.d3)) then
                  print*,height,'.GT.',(zcell/1.d3)
                  stop
               endif
c           # R / M    (approximately diatomic, gamma=constant)
               aux(18,i,j,k) = 8.31d3/
     &              ((dox*16.d0+dnit2*28.d0+dox2*32.d0)/
     &               (dox+dnit2+dox2))
c           # Gamma
               aux(19,i,j,k) = (7*(dox2+dnit2)+5*dox)/
     &                         (5*(dox2+dnit2)+3*dox)
               temp  =aux(5,i,j,k)
               RoverM=aux(18,i,j,k)
               gamma =aux(19,i,j,k)
               gamma1=gamma-1.d0
c
c           # Convert to MKS from CGS
               if (cgs) then
                  cnv=1.d3
               else
                  cnv=1.d-3
               endif
c
               aux(1,i,j,k) = (dox*16.d0+dox2*32.d0+dnit2*28.d0)*cnv
               aux(1,i,j,k) = aux(1,i,j,k)/(6.022d23)
c           # Analytical wind function in x (if desired)
               if ((amp0 .EQ. 0) .AND. (amp1 .EQ. 0)) then
                   aux(2,i,j,k) = aux(1,i,j,k)*windx
                   aux(3,i,j,k) = aux(1,i,j,k)*windy
               else
                  windx=amp0
                  windy=0.d0
                  scale=(aux(1,i,j,1)/aux(1,i,j,k))**0.5d0
                  aux(2,i,j,k) = aux(1,i,j,k)
     &             * windfnx(xcell,ycell,zcell,0,scale)
                  aux(2,i,j,k) = aux(2,i,j,k)+windx*aux(1,i,j,k)
                  aux(3,i,j,k) = windy
               endif
               aux(4,i,j,k) = 0.d0
               aux(5,i,j,k) = aux(1,i,j,k)*aux(5,i,j,k)*
     &                aux(18,i,j,k)/gamma1+
     &                .5*(aux(2,i,j,k)**2)/aux(1,i,j,k)+
     &                .5*(aux(3,i,j,k)**2)/aux(1,i,j,k)
c           # Set tracer number density
c               tracer = dne   ! Just use the Ne/Na data from profile
c           # Apply ampules at multiple heights
               tracer = ampule(xcell,ycell,zcell,1d4,5d3,9d4) +
     &                  ampule(xcell,ycell,zcell,1d4,5d3,1d5) +
     &                  ampule(xcell,ycell,zcell,1d4,5d3,1.1d5) +
     &                  ampule(xcell,ycell,zcell,1d4,5d3,1.2d5)
c           # Reserve the following variables
               aux(6,i,j,k) = dox*16.d0*cnv/(6.022d23*aux(1,i,j,k))
               aux(7,i,j,k) = dox2*32.d0*cnv/(6.022d23*aux(1,i,j,k))
c               aux(i,j,k,8) = 1.d8*tracer/(dox+dox2+dnit2)
               aux(8,i,j,k) = dne/(dox+dox2+dnit2)
               aux(9,i,j,k) = dhyd/(dox+dox2+dnit2)
c           # Reserve for Ozone
               aux(10,i,j,k)= 0.d0
c           # Calculate diffusivities from MSIS Profile
               Cp=1/(1-1/gamma)*RoverM
               dmolec=(3.43*dnit2+4.03*dox2+3.90*dox)*1d-7
               dmolec=dmolec*(temp**(0.69))/(dox+dnit2+dox2)
               dtherm=(56.0*(dnit2+dox2)+75.9*dox)*1d-5
               dtherm=dtherm*(temp**(0.69))/(dox+dnit2+dox2)
c           # Assign values of diffusivity
               dmolecular=dmolec*(1/aux(1,i,j,k))
               dthermal=dtherm*(1/(aux(1,i,j,k)*Cp))
               aux(11,i,j,k)=dmolecular
               aux(12,i,j,k)=dthermal*aux(19,i,j,k)
               if (k .gt. 375) THEN
               aux(11,i,j,k)=aux(11,i,j,375)
               aux(12,i,j,k)=aux(12,i,j,375)
               end if

            endif
c
c
c           # Set "virtual" gravity force
c
            if (k.GT.(1-mbc)) then
               auxpl = (aux(19,i,j,k-1)-1.d0)
     &          * (aux(5,i,j,k-1)-0.5d0*(aux(2,i,j,k-1)**2
     &          + aux(3,i,j,k-1)**2 + aux(4,i,j,k-1)**2)/aux(1,i,j,k-1))
               auxpr = (aux(19,i,j,k)-1.d0)
     &          * (aux(5,i,j,k)-0.5d0*(aux(2,i,j,k)**2
     &          + aux(3,i,j,k)**2 + aux(4,i,j,k)**2)/aux(1,i,j,k))
               aux(17,i,j,k) = (auxpr-auxpl)/
     &           (.5d0*(aux(1,i,j,k)+aux(1,i,j,k-1)))
            endif
c
c
  15       continue
c
      if (id.EQ.0) then
         print*,'Max Diffusivities (molec,therm): ',
     &           aux(11,1,1,mz), aux(12,1,1,mz)
c        # Calculate fractional time step for typical diff. CFL
         gamma=aux(19,1,1,mz)
         gamma1=gamma-1
         pres = gamma1*(aux(5,1,1,mz)-0.5d0*(aux(2,1,1,mz)**2
     &       + aux(3,1,1,mz)**2 + aux(4,1,1,mz)**2)/aux(1,1,1,mz))
         dt=min(dx,dy,dz)/sqrt(gamma*pres/aux(1,1,1,mz))
         dxyzmin=min(dx,dy,dz)
         difmax=aux(11,1,1,mz)*(4.d0/3.d0)
         difcfl=0.45d0
         ndtr=ceiling((dt*difmax/(dxyzmin*dxyzmin))/difcfl)
         dtr=dt/ndtr
         conmax=aux(12,1,1,mz)
         concfl=0.45d0
         ndtrc=ceiling((dtr*conmax/(dxyzmin*dxyzmin))/concfl)
         dtrc=dtr/ndtrc
         print*,'Diffusive Steps Per Advective Step: ', ndtr
         print*,'Conductive Steps Per Diffusive Step: ', ndtrc
      endif
c
c     # Close profile.data
c
      if (loadprof.EQ.1) then
         close(7)
      endif
c
c
      goto 999
c
c     # Special Case for qinit-defined tests
  100 continue
      do 150 i = 1-mbc, mx+mbc
         do 150 j = 1-mbc,my+mbc
            do 150 k = 1-mbc, mz+mbc
               aux(1,i,j,k)  = 1.d0
               aux(17,i,j,k) = 0.d0   !Mustn't be updated
  150 continue
c
  999 return
      end
c
c
c     ======================================================
      double precision function ampule(x,y,z,xpos,ypos,zpos)
c     ======================================================
c
      implicit double precision (a-h,o-z)
      swidth = 1d2
      ampule = dexp(-.5*((x-xpos)**2/swidth**2))*
     &         dexp(-.5*((y-ypos)**2/swidth**2))*
     &         dexp(-.5*((z-zpos)**2/swidth**2))
c
      return
      end
