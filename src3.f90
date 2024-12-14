
      module thunderstorm

      implicit double precision (a-h,o-z)
      real, allocatable :: preciparray(:,:)
      integer :: minx,miny,maxx,maxy
      double precision :: dxp, dyp, dtp

      contains
!     =======================================================
      subroutine readprecip2array(t,dtp,dxp,dyp,preciparray,minx,miny,maxx,maxy)
!     =======================================================
!     Takes the inputs t (the current time) and dtp (the precip file time step) in seconds and loads the correct precipitation file into an
!     array (preciparray) of length (nlen,4): returns the preciparray, and the minimum and maximum x and y values for which there is precipdata.
!      precip columns are x (m), y (m), precip rate at (t/dtp) and precip rate at (t/dtp)+dtp
!     -----------------------------------------------------------------------------------
      ! Declarations
      implicit double precision (a-h,o-z)
      integer Reason,nlen
      double precision, intent(in) :: t
      character(len=100) :: ioerrmsg
      character*20 precipname
      character*5 aString
      logical :: file_exists, opn
      real, allocatable, intent(out) :: preciparray(:,:)
      integer, intent(out) :: minx,miny,maxx,maxy
      double precision, intent(out) :: dxp, dyp, dtp

      ! Define tlog to build the character string------
      tlog = int(t/dtp)*dtp

      ! Build the correct character string to open the correct file(s)------------------
      write(aString, 10) int(tlog)
  10  format (I5.5)
      precipname='precip'//aString//'.txt'

      ! Check if this filename exists and exit if it does not-----------------------------
      INQUIRE(FILE=precipname, EXIST=file_exists)

      IF (.not. file_exists) THEN
       print *, ' Error: file doesnt exist'
       stop

      ENDIF

      !Open file(s)-----------------------------------------------------------------------

      ! If already allocated, then deallocate
       IF (allocated(preciparray)) then
          deallocate(preciparray)
       END IF

      OPEN(UNIT=8, FILE=precipname,IOSTAT=Reason)

      ! Read header line to skip it
      READ(8,*,IOSTAT=Reason)
      ! Read precipitation array resolution
      READ(8,*,IOSTAT=Reason) dxp, dyp, dtp
      ! Read header line to skip it
      READ(8,*,IOSTAT=Reason)

      ! Perform a Do loop to cycle through the precip file to count lines so we can allocate
      !the array size
       nlen=0
      DO
       READ(8,*,IOSTAT=Reason) arbvar
       if (Reason /= 0) exit
       nlen=nlen+1
      END DO

      ! Allocate the array size based upon the lines counted
       allocate(preciparray(nlen,4))


      ! rewind the files and read the header lines to skip them again
      rewind(8)
      READ(8,*,IOSTAT=Reason)
      READ(8,*,IOSTAT=Reason)
      READ(8,*,IOSTAT=Reason)

      ! Perform the read loop into an array
      DO i=1,nlen
       READ(8,*,IOSTAT=Reason)(preciparray(i,j), j=1,4)

      END DO

      minx=MINVAL(preciparray(1:nlen,1))
      miny=MINVAL(preciparray(1:nlen,2))
      maxx=MAXVAL(preciparray(1:nlen,1))
      maxy=MAXVAL(preciparray(1:nlen,2))


       !----------------------------------------------------------------------------------------
      CLOSE(8)


      RETURN
      END subroutine



!     ===================================================================================
      subroutine preciparray2heating(x,y,z,t,dxp,dyp,dtp,dx,dy,heat,preciparray)
!     ===================================================================================

!
!     Takes the inputs of x, y, and z in m, t in seconds and the precipitation array (preciparray)
!     precip file resolutions dxp,dyp,dtp, length of the file  and return the latent heating rate (K/s) for this location and time.
!     -----------------------------------------------------------------------------------
      ! Declarations
      implicit double precision (a-h,o-z)
      real, allocatable :: preciparray(:,:)
      integer n

      ! Exit if altitude if greater than 15 km (no forcing here)
      IF (z .gt. 15000) then

      RETURN
      END IF

       ! convert the precip to latent heating-----------------------------------------
    ! NEXT SEARCH THE ARRAY FOR A PRECIP VALUE AT x,y,z

     ! Define logical comparison for if the input x, y, and t fall in xp+dtp, yp+dtp, tp+dtp------
      tlog = int(t/dtp)*dtp
      xlog = int(x/dxp)*dxp
      ylog = int(y/dyp)*dyp
      ! define the length of the preciparray
      n=size(preciparray,1)

      ! Perform a Do loop to cycle through the precip file and look for matches between
      !input co-ordinates and precipitation values in the file
      DO l=1,n

            xp=preciparray(l,1)
            yp=preciparray(l,2)
            precipt=preciparray(l,3)
            precipdt=preciparray(l,4)

       IF ((((xp .eq. xlog)) .AND. &
            (yp .eq. ylog))) THEN

            ! Take weighted average
             w1=((t/dtp) - int(t/dtp))
             precipw=((1.d0-w1)*precipt + w1*precipdt)*((dx*dy)/(dxp*dyp))
!             precipw=precipt

      !Define heating parameters from precipitation rate using Stephan and Alexander [2015]
      th=10.79d0+0.05d0*precipw !top of heating profile
      ph=6.93d0-0.03d0*precipw !peak of heating profile
      bh=3.43d0-0.14d0*precipw !bottom of heating profile
      hamp=-0.0001d0+0.0025d0*precipw !heating amplitude

      tc=1.12d0+0.10d0*precipw ! top of cooling profile
      pc=0.63d0+0.02d0*precipw ! peak of cooling profile
      bc=0.37d0+0.00d0*precipw ! bottom of cooling profile
      camp=-0.0003d0-0.0005d0*precipw !cooling amplitude


         !Implement four quarter sines for each section
         !pi=3.141592653589793238
         pi=4.D0*DATAN(1.D0)
         z2=z/1000.d0
          !peak of cooling to top of cooling
         IF ((bc .LE. z2) .AND. (z2 .LE. pc)) THEN
           heat=camp*sin((2.d0*pi/(4.d0*(pc-bc))) &
           *(z2-bc))

          ! peak to top of cooling
         ELSEIF ((pc .LE. z2) .AND. (z2 .LE. tc)) THEN
           heat=-camp*sin((2.d0*pi/(4.d0*(tc-pc))) &
           *(z2-tc))

         ! bottom to peak heating
         ELSEIF ((pc .LE. z2) .AND. (z2 .LE. ph)) THEN
           heat=hamp*sin((2.d0*pi/(4.d0*(ph-tc))) &
           *(z2-tc))

         ! peak to top of heating
         ELSEIF ((ph .LE. z2) .AND. (z2 .LE. th)) THEN
           heat=-hamp*sin((2.d0*pi/(4.d0*(th-ph))) &
           *(z2-th))

         ELSE
          heat=0.d0

         ENDIF

         !If we have a match then no point to continue looping through file so exit
         RETURN

        END IF

      END DO


      RETURN
      END subroutine

      end module thunderstorm

!    =======================================================
      subroutine src3(meqn,mbc,mx,my,mz, &
                     xlower,ylower,zlower,dx,dy,dz,q, &
                     maux,aux,t,dt,mthbc)
!     =======================================================
      use thunderstorm
      implicit double precision (a-h,o-z)
      real :: start, finish
      external gaussian, forcefnz, viscous,&
              conduct, bc3, ohsteady

      common /param/ gamma, gamma1, grav, dens, pres, &
                       speed, dtherm0, dmolec0
      common /forcing/ omega, amplitude, propx, propz, &
                       w_n, w_o, vsrcx, forcemth
      common /wind/ omega1, amp0, amp1, prop1z, z1width, z1pos, &
                       z1Gpos, phi, windscale1
      common /gauss/ xpos, ypos, zpos, xwidth, ywidth, &
                    zwidth, tcenter, twidth
      common /options/ loadprof, nspecadv
!     # MPI : Added so that we know which processor we are running on.
      common /mpi_proc_info/ np, id

      dimension    q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc, &
                    1-mbc:mz+mbc)
      dimension    aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, &
                    1-mbc:mz+mbc)

      dimension mthbc(6)
      dimension  az(1-mbc:mz+mbc)


!  --------------------------------------------------------

!     Source Terms for Atmospheric Model Extensions
!        by J. B. Snively, 2003-2018

!  --------------------------------------------------------
!     # Apply viscous diffusion and thermal conduction

      do k = 1-mbc, mz+mbc
         do  j = 1-mbc, my+mbc
            do  i = 1-mbc, mx+mbc
            gamma=aux(19,i,j,k)
            gamma1=gamma-1
!        # Convert energy to temperature
            qtemp = gamma1* &
              (q(5,i,j,k)-0.5d0*(1/q(1,i,j,k))*(q(2,i,j,k)**2 &
             + q(3,i,j,k)**2 + q(4,i,j,k)**2))&
             /(q(1,i,j,k)*aux(18,i,j,k))
            auxtemp = gamma1* &
              (aux(5,i,j,k)-0.5d0*(1/aux(1,i,j,k))&
             *(aux(2,i,j,k)**2 + aux(3,i,j,k)**2&
             + aux(4,i,j,k)**2))&
             /(aux(1,i,j,k)*aux(18,i,j,k))
            q(5,i,j,k) = qtemp-auxtemp
!        # Subtract Wind from u, convert momentum to velocity
            q(2,i,j,k) = q(2,i,j,k)/q(1,i,j,k) &
             -aux(2,i,j,k)/aux(1,i,j,k)
!        # Convert momentum to velocity
            q(3,i,j,k) = q(3,i,j,k)/q(1,i,j,k) &
             -aux(3,i,j,k)/aux(1,i,j,k)
!        # Convert momentum to velocity
            q(4,i,j,k) = q(4,i,j,k)/q(1,i,j,k) &
             -aux(4,i,j,k)/aux(1,i,j,k)
!        # Define a sponge layer
            zcell = zlower + (k-0.5d0)*dz
            ycell = ylower + (j-0.5d0)*dy
            xcell = xlower + (i-0.5d0)*dx
          ! az(k)=0.09d0*dexp(-log(2.d0)*((zcell-300d3)/20d3)**2)
          ! q(2,i,j,k)=q(2,i,j,k)-dt*az(k)*q(2,i,j,k)
          ! q(3,i,j,k)=q(3,i,j,k)-dt*az(k)*q(3,i,j,k)
          ! q(4,i,j,k)=q(4,i,j,k)-dt*az(k)*q(4,i,j,k)
          ! q(5,i,j,k)=q(5,i,j,k)-dt*az(k)*q(5,i,j,k)

          ! funfu=0.01d0*dexp(-log(2.d0)*((xcell-2400d3)/5d3)**2)
          ! q(2,i,j,k)=q(2,i,j,k)-dt*funfu*q(2,i,j,k)
          ! q(3,i,j,k)=q(3,i,j,k)-dt*funfu*q(3,i,j,k)
          ! q(4,i,j,k)=q(4,i,j,k)-dt*funfu*q(4,i,j,k)
          ! q(5,i,j,k)=q(5,i,j,k)-dt*funfu*q(5,i,j,k)

          ! funfu2=0.01d0*dexp(-log(2.d0)*((xcell)/5d3)**2)
          ! q(2,i,j,k)=q(2,i,j,k)-dt*funfu2*q(2,i,j,k)
          ! q(3,i,j,k)=q(3,i,j,k)-dt*funfu2*q(3,i,j,k)
          ! q(4,i,j,k)=q(4,i,j,k)-dt*funfu2*q(4,i,j,k)
          ! q(5,i,j,k)=q(5,i,j,k)-dt*funfu2*q(5,i,j,k)

          ! funfv=0.01d0*dexp(-log(2.d0)*((ycell-1600d3)/5d3)**2)
          ! q(2,i,j,k)=q(2,i,j,k)-dt*funfv*q(2,i,j,k)
          ! q(3,i,j,k)=q(3,i,j,k)-dt*funfv*q(3,i,j,k)
          ! q(4,i,j,k)=q(4,i,j,k)-dt*funfv*q(4,i,j,k)
          ! q(5,i,j,k)=q(5,i,j,k)-dt*funfv*q(5,i,j,k)

          ! funfv2=0.01d0*dexp(-log(2.d0)*((ycell)/5d3)**2)
          ! q(2,i,j,k)=q(2,i,j,k)-dt*funfv2*q(2,i,j,k)
          ! q(3,i,j,k)=q(3,i,j,k)-dt*funfv2*q(3,i,j,k)
          ! q(4,i,j,k)=q(4,i,j,k)-dt*funfv2*q(4,i,j,k)
          ! q(5,i,j,k)=q(5,i,j,k)-dt*funfv2*q(5,i,j,k)


           end do
         end do
       end do

!     # Apply viscous diffusion of momentum

!     Hint: conduct3 and viscous3 are unsplit; use with cfl 0.1
!           (will be faster in nearly inviscid problems)
!           conduct and viscous are split; use with cfl 0.4
!           (may or may not be faster, TBD; probably more accurate)

!     # Calculate fractional time step
      dxyzmin2=min(dx,dy,dz)**2
      difmax=(4.d0/3.d0)*aux(11,1,1,mz)
      difcfl=0.45d0
      ndtr=ceiling((dt*difmax/(dxyzmin2))/difcfl)
      dtr=dt/ndtr
      conmax=aux(12,1,1,mz)
      concfl=0.45d0
      ndtrc=ceiling((dtr*conmax/(dxyzmin2))/concfl)
      dtrc=dtr/ndtrc

!     # Compute actual
      difcfl=difmax*dtr/(dxyzmin2)
      concfl=conmax*dtrc/(dxyzmin2)

!     # Store Boundary conditions
!     # Note that u and T are subject to same side wall boundary
!     # conditions, however tops and bottoms will differ from
!     # momentum and internal energy. Just disable these...

      mthbc_save_5=mthbc(5)
      mthbc_save_6=mthbc(6)
      mthbc(5)=-1
      mthbc(6)=-1


!     # Substep...

      do n=1,ndtr

!     # Set time
       	 tr = t

!        # Apply boundary conditions
         call bc3(meqn,mbc,mx,my,mz, &
              xlower,ylower,zlower,dx,dy,dz, &
              q,maux,aux,tr,dtr,mthbc)

!        # Molecular Viscosity (Diffusion of Momentum)
         call diffuse(mbc,mx,my,mz,xlower, &
           ylower,zlower,dx,dy,dz,meqn,q,maux,aux,11,dtr)

!        # Thermal Conduction (Diffusion of Temperature)
         call conduct(mbc,mx,my,mz,xlower, &
             ylower,zlower,dx,dy,dz,meqn,q,5,maux,aux,12,dtrc,ndtrc)

!        # Calculate intermediate time
         tr=t+n*dtr


       end do

!    # Restore boundary conditions

      mthbc(5)=mthbc_save_5
      mthbc(6)=mthbc_save_6

      if (id.eq.0) then
         write(6,601) ndtr,ndtrc,difcfl,concfl
  601    format('SRC3...  Vis. Steps',i5,'  Con. Sub-steps',i5, &
                        ' VisCFL =',f6.3,' ConCFL =',f6.3)
      endif

!    # Add sponge layer here

!     # Reassemble q variables, but not on BCs

      do k = 1, mz
         zcell = zlower + (k-0.5d0)*dz
         do j = 1, my
            ycell = ylower + (j-0.5d0)*dy
            do i = 1, mx
            xcell = xlower + (i-0.5d0)*dx
            gamma = aux(19,i,j,k)
            gamma1= gamma-1
!        # Update time-dependent wind field
            if (omega1 .NE. 0) then
               aux(5,i,j,k) = aux(5,i,j,k)- &
                 0.5d0*((aux(2,i,j,k))**2)/aux(1,i,j,k)
               scale = (aux(1,i,j,1)/aux(1,i,j,k))**0.5d0
               aux(2,i,j,k) = aux(2,i,j,k)-aux(1,i,j,k) &
                 *windfnx(xcell,ycell,zcell,t-dt,scale)
               aux(2,i,j,k) = aux(2,i,j,k)+aux(1,i,j,k) &
                 *windfnx(xcell,ycell,zcell,t,scale)
               aux(5,i,j,k) = aux(5,i,j,k)+ &
                 0.5d0*((aux(2,i,j,k))**2)/aux(1,i,j,k)
            endif
!        # Convert velocity to momentum
            q(2,i,j,k) = q(2,i,j,k)*q(1,i,j,k) + &
                   (aux(2,i,j,k)/aux(1,i,j,k))*q(1,i,j,k)
            q(3,i,j,k) = q(3,i,j,k)*q(1,i,j,k) + &
                   (aux(3,i,j,k)/aux(1,i,j,k))*q(1,i,j,k)
            q(4,i,j,k) = q(4,i,j,k)*q(1,i,j,k) + &
                   (aux(4,i,j,k)/aux(1,i,j,k))*q(1,i,j,k)
!        # Convert temperature perturbation to energy density
            auxtemp = gamma1* &
              (aux(5,i,j,k)-0.5d0*(1/aux(1,i,j,k))* &
              (aux(2,i,j,k)**2 + aux(3,i,j,k)**2 + &
               aux(4,i,j,k)**2))/(aux(1,i,j,k)*aux(18,i,j,k))
            q(5,i,j,k) = (auxtemp + q(5,i,j,k))* &
               q(1,i,j,k)*aux(18,i,j,k)/gamma1 + &
               0.5d0*(1/q(1,i,j,k))*(q(2,i,j,k)**2 + &
               q(3,i,j,k)**2 + q(4,i,j,k)**2)
            end do
          end do
        end do

!     # OH Airglow Chemistry
!      call ohsteady(meqn,mbc,mx,my,mz, &
!          xlower,ylower,zlower,dx,dy,dz,q,maux,aux,dt)






!----------# Oscillatory Forcing Function-------------------------------------------------
if (forcemth.EQ.1) then

do i=1,mx
  xcell = xlower + (i-0.5d0)*dx - t*vsrcx
   do j=1,my
    ycell = ylower + (j-0.5d0)*dy
     do k=1,mz
       zcell = zlower + (k-0.5d0)*dz

              bf=dt*amplitude* &
                   gaussian(xcell,ycell,zcell,t)* &
                   forcefnz(xcell,ycell,zcell,t)
               q(4,i,j,k) = q(4,i,j,k)+bf*q(1,i,j,k)
      enddo
    enddo
  enddo
  !-----------Thunderstorm forcing with precip files---------------------------------------------------------
  elseif (forcemth.EQ.2) then

  ! Set the resolution of the input precipitation dataset
              dtp=60
              dxp=2000
              dyp=2000

   ! Call the readprecip array if (t/dtf) is a whole number.
              IF ((t .eq. 0 ) .OR. (MOD(t,dtp) .eq. 0 )) THEN
              print*, 'readprecip2array has been called'
              CALL readprecip2array(t,dtp,dxp,dyp,preciparray,minx,miny,maxx,maxy)
              ENDIF

do i=1,mx
  xcell = xlower + (i-0.5d0)*dx - t*vsrcx
   do j=1,my
    ycell = ylower + (j-0.5d0)*dy
     do k=1,mz
       zcell = zlower + (k-0.5d0)*dz
         if (zcell .gt. 15000) exit

                 if (xcell .ge. minx) then
                 if (xcell .le. maxx) then
                 if (ycell .ge. miny) then
                 if (ycell .le. maxy) then
  
           CALL preciparray2heating(xcell,ycell,zcell,t,dxp,dyp,dtp,dx,dy,heat,preciparray)

                        gamma=aux(19,i,j,k)
                        gamma1=gamma-1
              !        # Convert energy to temperature

                        qtemp = gamma1* &
                          (q(5,i,j,k)-0.5d0*(1/q(1,i,j,k))*(q(2,i,j,k)**2 &
                         + q(3,i,j,k)**2 + q(4,i,j,k)**2)) &
                         /(q(1,i,j,k)*aux(18,i,j,k))
                        auxtemp = gamma1* &
                          (aux(5,i,j,k)-0.5d0*(1/aux(1,i,j,k)) &
                         *(aux(2,i,j,k)**2 + aux(3,i,j,k)**2 &
                         + aux(4,i,j,k)**2)) &
                         /(aux(1,i,j,k)*aux(18,i,j,k))
                        q(5,i,j,k) = qtemp-auxtemp

              !        # Add latent heating to the temperature
                     q(5,i,j,k)=q(5,i,j,k)+dt*heat

              !        # Convert temperature perturbation to energy density
                        auxtemp = gamma1* &
                          (aux(5,i,j,k)-0.5d0*(1/aux(1,i,j,k))* &
                          (aux(2,i,j,k)**2 + aux(3,i,j,k)**2 + &
                           aux(4,i,j,k)**2))/(aux(1,i,j,k)*aux(18,i,j,k))
                        q(5,i,j,k) = (auxtemp + q(5,i,j,k))* &
                           q(1,i,j,k)*aux(18,i,j,k)/gamma1 + &
                           0.5d0*(1/q(1,i,j,k))*(q(2,i,j,k)**2 + &
                           q(3,i,j,k)**2 + q(4,i,j,k)**2)


                         end if
                         end if
                         end if
                         end if

                           end do
                       end do
                     end do

!-----------Thunderstorm forcing with Analytical precip---------------------------------------------------------
elseif (forcemth.EQ.3) then

do i=1,mx
  xcell = xlower + (i-0.5d0)*dx - t*vsrcx
   do j=1,my
    ycell = ylower + (j-0.5d0)*dy
     do k=1,mz
       zcell = zlower + (k-0.5d0)*dz
         if (zcell .gt. 15000) exit

         gaussian2=dexp(-.5*((xcell-xpos)**2/xwidth**2))* &
          dexp(-.5*((ycell-ypos)**2/ywidth**2))* &
          dexp(-.5*((t-tcenter)**2/twidth**2))


         precipw=amplitude*gaussian2

         if (precipw .gt. 0.1d0) then

         !Define heating parameters from precipitation rate using Stephan and Alexander [2015]
         th=10.79d0+0.05d0*precipw !top of heating profile
         ph=6.93d0-0.03d0*precipw !peak of heating profile
         bh=3.43d0-0.14d0*precipw !bottom of heating profile
         hamp=-0.0001d0+0.0025d0*precipw !heating amplitude

         tc=1.12d0+0.10d0*precipw ! top of cooling profile
         pc=0.63d0+0.02d0*precipw ! peak of cooling profile
         bc=0.37d0+0.00d0*precipw ! bottom of cooling profile
         camp=-0.0003d0-0.0005d0*precipw !cooling amplitude


           !Implement four quarter sines for each section
           !pi=3.141592653589793238
           pi=4.D0*DATAN(1.D0)
           z2=zcell/1000.d0
            !peak of cooling to top of cooling
           IF ((bc .LE. z2) .AND. (z2 .LE. pc)) THEN
             heat=camp*sin((2.d0*pi/(4.d0*(pc-bc))) &
             *(z2-bc))

            ! peak to top of cooling
           ELSEIF ((pc .LE. z2) .AND. (z2 .LE. tc)) THEN
             heat=-camp*sin((2.d0*pi/(4.d0*(tc-pc))) &
             *(z2-tc))

           ! bottom to peak heating
           ELSEIF ((pc .LE. z2) .AND. (z2 .LE. ph)) THEN
             heat=hamp*sin((2.d0*pi/(4.d0*(ph-tc))) &
             *(z2-tc))

           ! peak to top of heating
           ELSEIF ((ph .LE. z2) .AND. (z2 .LE. th)) THEN
             heat=-hamp*sin((2.d0*pi/(4.d0*(th-ph))) &
             *(z2-th))

           ELSE
            heat=0.d0

           ENDIF


              gamma=aux(19,i,j,k)
              gamma1=gamma-1
  !        # Convert energy to temperature

              qtemp = gamma1* &
                (q(5,i,j,k)-0.5d0*(1/q(1,i,j,k))*(q(2,i,j,k)**2 &
               + q(3,i,j,k)**2 + q(4,i,j,k)**2)) &
               /(q(1,i,j,k)*aux(18,i,j,k))
              auxtemp = gamma1* &
                (aux(5,i,j,k)-0.5d0*(1/aux(1,i,j,k)) &
               *(aux(2,i,j,k)**2 + aux(3,i,j,k)**2 &
               + aux(4,i,j,k)**2)) &
               /(aux(1,i,j,k)*aux(18,i,j,k))
              q(5,i,j,k) = qtemp-auxtemp

  !        # Add latent heating to the temperature
           q(5,i,j,k)=q(5,i,j,k)+dt*heat

  !        # Convert temperature perturbation to energy density
              auxtemp = gamma1* &
                (aux(5,i,j,k)-0.5d0*(1/aux(1,i,j,k))* &
                (aux(2,i,j,k)**2 + aux(3,i,j,k)**2 + &
                 aux(4,i,j,k)**2))/(aux(1,i,j,k)*aux(18,i,j,k))
              q(5,i,j,k) = (auxtemp + q(5,i,j,k))* &
                 q(1,i,j,k)*aux(18,i,j,k)/gamma1 + &
                 0.5d0*(1/q(1,i,j,k))*(q(2,i,j,k)**2 + &
                 q(3,i,j,k)**2 + q(4,i,j,k)**2)

           endif

                 end do
             end do
           end do

  endif

!--------------------END OF Forcing section-----------------------------------

!     # Correct Negative Densities
!     # (This shouldn't happen, but with layers it might.)

      if (meqn.EQ.5) goto 500

      do k = 1, mz
         do j = 1, my
            do i = 1, mx
               do m=6,meqn
                  q(m,i,j,k) = dmax1(0.d0,dmin1(q(m,i,j,k),1.d0))
               end do
            end do
          end do
      end do

!     NOTE THAT SRC3 HAS RUINED ALL BOUNDARY CELLS!!!
!        # Apply boundary conditions once more,
!        # Only necessary with Strang Splitting!
!        # Hence, we've modified claw3 to do this.

  500 return
      end
