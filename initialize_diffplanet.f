      subroutine diffplanetgen
c     Creates a gaseous planet from the data file yrec output
      include 'sphu.h'
      include 'mpif.h'
      integer numlines,i,maxTempI
      integer idumb,ip,ix,iy,iz
      real*8 anumden,rhotry,rhoex,rtry,rhomax,hc,xcm,ycm,zcm,amtot,
     $     ammin,ammax,xtry,ytry,ztry,ri,rhoi,ran1
      integer irtry
      real*8 amass,masscgs,radius
      real*8 tem(kdm),pres(kdm),rarray(kdm),
     &     rhoarray(kdm),rhoarray1(kdm),rhoarray2(kdm),rhoarray3(kdm),
     $     uarray(kdm),uarray1(kdm),uarray2(kdm),uarray3(kdm),
     $     muarray(kdm),muarray1(kdm),muarray2(kdm),muarray3(kdm)
c     Added by jason
      real*8 u_g,u_c,u_m,x_m,x_c,xm,xc
      real*8 rhoarrayi(kdm),uarrayi(kdm),muarrayi(kdm)
      real*8 rhoarraym(kdm),uarraym(kdm),muarraym(kdm)
      real*8 rhoarraye(kdm),uarraye(kdm),muarraye(kdm)
      real*8 farray1(kdm),farray2(kdm),farray(kdm),farrayc(kdm)
      real*8 egsol
c     Astronomical constants:
      parameter(egsol=1.9891d+33)
c     Derived constants:
      real*8 integratednum
      common/splinestuff/rarray,uarray,muarray,rhoarray,
     $     uarray1,uarray2,uarray3,
     $     muarray1,muarray2,muarray3,
     $     rhoarray1,rhoarray2,rhoarray3,
     $     uarraye,uarraym,uarrayi,
     $     muarraye,muarraym,muarrayi,
     $     rhoarraye,rhoarraym,rhoarrayi,
     $     amass,radius,integratednum,maxmu,minmu,numlines,maxTempI,
     $     farray1,farray2,farray,farrayc
      real*8 utot2,wtota2
      integer ixmax,iymax,izmax,corepts
      double precision cellvolume,a1
      real*8 integral
      real*8 maxmu,minmu,drhodhi
      double precision utottest
      real*8 epot
      real*8 redge1,redge2
      real*8 hpguess,xacc,dxmax
      real*8 mci
      integer nmin
      real*8 amass1,amass2
      common/forcompbest/ amass1,amass2
      common/presarray/ pres,i
      common/hack/tem,redge1,masscgs,utot2,wtota2
      integer mygravlength, ierr
      integer comm_worker, irank
      common/gravworkers/comm_worker
      integer status(mpi_status_size)

      firstflag = 1

      call splinesetupdiffplanet

      idumb=-2391
      rhomax=rhoarray(1)
      if(myrank.eq.0)write(69,*)'maximum central density=',rhomax

      if(myrank.eq.0)write(69,*)'85 percent radius redge1=',redge1

      if (n.lt.0) then
         n=abs(n)*masscgs/egsol
      endif

      if (gflag.eq.0) then
         hc=0.5d0*radius*(1.8d0*NNOPT/n)**(1.d0/3.d0) ! neighbor number is about 1.8*NNOPT, roughly
      else
         hc=0.5d0*radius*(1.3d0*NNOPT/n)**(1.d0/3.d0) !Trying to better estimate the nearest neighbor number... about 1.3*NNOPT
      endif

      redge2=radius-3.d0*hc
      if(myrank.eq.0)write(69,*)'3h away from surface is redge2=',redge2

      redge=max(redge1,redge2)
      if (myrank.eq.0) write(69,*)'redge=max(redge1,redge2)=',redge

      if (myrank.eq.0) write(69,*)'hc=',hc
      
      nmin=max(156,int(12.d0*(0.5d0*redge/hc)**3.d0))
      if (myrank.eq.0) write(69,*) 'N should be at least',nmin

      n = max(n,nmin)

      if (myrank.eq.0) write(69,*) 'will try for n=',n

      mci=rhomax*4.d0/3.d0*pi*redge**3/n
      if (myrank.eq.0) write(69,*)'a central particle would have mass',mci

      if (nnopt.le.0) then
         nnopt=max(nint(n*(hc/redge)**3.d0*8.d0),13)
      endif
      if (myrank.eq.0) write(69,*) 'Using nnopt=',nnopt

      if (4*mci.ge.amass) then
c     If corepts=1 then the core point will be particle 1... first sph particle will be particle 2
         corepts=1
         if(myrank.eq.0)write(69,*)'We will use a core point'
      else
c     If corepts=0 then there is no core particle and the first sph particle will be particle 1
         corepts=0
         if (myrank.eq.0) write(69,*)'We will not use a core point'
      endif
      ip=corepts

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     make an hcp lattice
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(myrank.eq.0)write(69,*)'keeping particles up to a distance',radius-redge,
     $     'less than the full radius',radius

c     The fraction of particles at a region of density rhoex that
c     will be kept is (rhoex/rhomax)**equalmass.  so if the number
c     density of lattice points that will be tried is n, then we expect
c     the number of particles to be N=n*integrate[4*pi*r**2*
c     (rhoex(r)/rhomax)**equalmass,{r,0,redge}].  we can solve this for
c     n, and then use that the cell volume is 2/n (as there are two
c     particles per cell)
      
      if (equalmass.ne.0.d0) then
         i=2
         integral=0
         do while (rarray(i).lt.redge)
            integral=integral+pi*(rarray(i)+rarray(i-1))**2*
     $         (rarray(i)-rarray(i-1))*
     $         (0.5d0*(rhoarray(i)+rhoarray(i-1))/rhomax)**equalmass
            i=i+1
         enddo
         if(myrank.eq.0)write(69,*)'integral=',integral,4.d0/3.d0*pi*redge**3.d0
      else
         integral=4.d0/3.d0*pi*redge**3.d0
         if(myrank.eq.0)write(69,*)'volume integral=',integral,4.d0/3.d0*pi*redge**3.d0
      endif

      cellvolume=2.d0*integral/n
      a1=(cellvolume/2.d0**0.5d0)**(1.d0/3.d0)
      if(myrank.eq.0)write(69,*)'a1=',a1,cellvolume,integral
      
c     looking at figure 9(b) and page 18 of kittel
c     (a1 vector)=-0.5*a1*(x hat)-3^0.5/2*a1*(y hat)
c     (a2 vector)=a1*(x hat)
c     (a3 vector)=(8/3)^0.5*a1*(z hat)
      ixmax=int(redge/a1)+2
      iymax=int(redge/(3.d0**0.5d0/2.d0*a1))+2
      izmax=int(redge/(0.5d0*(8.d0/3.d0)**0.5d0*a1))+2
      if(myrank.eq.0)write(69,*) 'coreRad1=',coreRad1
      if(myrank.eq.0)write(69,*) 'coreRad2=',coreRad2
      if(myrank.eq.0)write(69,*) 'redge=',redge

      do ix=-ixmax,ixmax
         do iy=-iymax,iymax
            do iz=-izmax,izmax
               xtry=(ix-0.5d0)*a1+mod(abs(iy),2)*0.5d0*a1
               ytry=iy*3.d0**0.5d0/2.d0*a1
     $            -(mod(abs(iz),2)-0.5d0)*1.d0/3.d0**0.5d0*a1
               ztry=(iz-0.5d0)*0.5d0*(8.d0/3.d0)**0.5d0*a1
               rtry=sqrt(xtry**2.d0+ytry**2.d0+ztry**2.d0)
               if(rtry.lt.redge) then
c     using two arrays to interpolate density
                  if(rtry.ge.coreRad2) then
                    call sph_splint(rarray,rhoarraye,rhoarray1,numlines,
     $                    rtry,rhoex)
                  else
                    if(rtry.ge.coreRad1) then
                      call sph_splint(rarray,rhoarraym,rhoarray2,numlines,
     $                      rtry,rhoex)
                    else
                      call sph_splint(rarray,rhoarrayi,rhoarray3,numlines,
     $                      rtry,rhoex)
                    endif
                  endif
                  rhotry=rhomax**equalmass*ran1(idumb)
c     Determine if particle is accepted
                  if (rhotry.le.rhoex**equalmass) then
                     ip=ip+1                 
                     if(rtry.le.a1 .and. myrank.eq.0) then
                         write(69,'(4i5,4e11.4)')ip,ix,iy,iz,
     $                        xtry,ytry,ztry,rtry
                     endif
!                     if (rtry.gt.coreRad2) then
!                        cc(ip) = 1
!                     else
!                        if (rtry.gt.coreRad1) then
!                          cc(ip) = 2
!                        else
!                          cc(ip) = 2
!                        endif
!                     endif
                     cc(ip) = 2
                     if (rtry.gt.(coreRad2+coreRad1)/2.) cc(ip) = 1
                     x(ip)=xtry        
                     y(ip)=ytry            
                     z(ip)=ztry
                  endif
               endif
            enddo
         enddo
      enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     (note that n will be slightly changed!)
      n=ip-corepts
      if(myrank.eq.0) write (69,*) 'egrgparent: n=',n
      amass1=n
      amass2=0
      if (n+corepts.gt.nmax) stop 'parent: n>nmax ???'
c     Assign particle masses (to represent density):
      amtot=0.d0
      xcm=0.d0
      ycm=0.d0
      zcm=0.d0

c     Start following loop at the first sph particle (with index 2, not 1, if there is a core particle)
      do i=1+corepts,n+corepts
        ri=sqrt(x(i)**2+y(i)**2+z(i)**2)
c using two arrays to interpolate
        if(ri.ge.coreRad2) then
          call sph_splint(rarray,rhoarraye,rhoarray1,numlines,ri,rhoi)
          call sph_splint(rarray,uarraye,  uarray1,  numlines,ri,u(i))
          call sph_splint(rarray,muarraye, muarray1, numlines,ri,
     $         meanmolecular(i))
          call sph_splint(rarray,farray, farray1,numlines,ri,ufactor(i))
        else
          call sph_splint(rarray,farrayc,farray2,numlines,ri,ufactor(i))
          if(ri.ge.coreRad1) then
            call sph_splint(rarray,rhoarraym,rhoarray2,numlines,ri,rhoi)
            call sph_splint(rarray,uarraym,  uarray2,  numlines,ri,u(i))
            call sph_splint(rarray,muarraym, muarray2, numlines,ri,
     $           meanmolecular(i))
          else
            call sph_splint(rarray,rhoarrayi,rhoarray3,numlines,ri,rhoi)
            call sph_splint(rarray,uarrayi,  uarray3,  numlines,ri,u(i))
            call sph_splint(rarray,muarrayi, muarray3, numlines,ri,
     $           meanmolecular(i))
          endif
        endif

        if(rhoi.le.0.d0) then
          if(myrank.eq.0)
     $      write(69,*)'warning: rho(',i,')<=0 at r=',ri,'???'
        endif
        rho0(i) = rhoi+0.0d0
        T0(i) = tem(maxTempI)
        if (ri.lt.coreRad1) then
           if (u_cc.gt.u(i)) then
              u_cc  = u(i)
              rhocc = rhoi
           endif
        elseif (ri.lt.coreRad2) then
           if (u_cm.lt.u(i)) then
              u_cm  = u(i)
              rhocm = rhoi
           endif
           if (u_ci.gt.u(i)) then
              u_ci  = u(i)
              rhoci = rhoi
           endif
        else
           if (u_ei.lt.u(i)) then
              u_ei  = u(i)
              rhoei = rhoi
           endif
!          if (ri.ge.coreRad2) then
           T0(i) = (rho0(i)*munit/runit**3.d0)**(envGam-1.)/(envGam-1.)
           T0(i) = u(i)*gravconst*munit/runit-envConst*T0(i)
           T0(i) = T0(i)/(boltz*ufactor(i))*meanmolecular(i)
        endif
!        T0(i) = min(tem(maxTempI),T0(i))
        am(i)=amass/n*(integral*rhoi/amass)**(1.d0-equalmass)
        xcm=xcm+am(i)*x(i)
        ycm=ycm+am(i)*y(i)
        zcm=zcm+am(i)*z(i)
        amtot=amtot+am(i)
        anumden=rhoi/am(i)
c     The actual number of neighbors is closer to 1.9*nnopt
        hp(i)=(3.d0/32.d0/3.1415926535897932384626d0*
     $        1.9d0*nnopt/anumden)**(1.d0/3.d0)

      enddo

      if (myrank.eq.0) then
         open(13,file='sph.rhoparams',status='unknown')
         write(13,'(10g15.7)')rhoei,rhoci,rhocm,rhocc,envConst,envGam,
     $      u_ei,u_ci,u_cm,u_cc
         close(13)
      endif
      
      xcm=xcm/amtot
      ycm=ycm/amtot
      zcm=zcm/amtot
      if(myrank.eq.0) then
        write(69,'(a,3g13.4)') 'center of mass at',xcm,ycm,zcm
        write(69,'(a,g11.4,a,g11.4,a)')
     $       'PARENT: total mass SPH fluid was M=',amtot
     $       ,
     $       ', total mass of star will equal', amass
      endif

      hmin=1.d30
      hmax=0.d0
      do i=1+corepts,n+corepts
        hmin=min(hmin,hp(i))
        hmax=max(hmax,hp(i))
      enddo
      if(myrank.eq.0)then
        write(69,*)'hmin=',hmin
        write(69,*)'hmax=',hmax
      endif

c     Above here, n=number of SPH particles
      n=n+corepts
c     Below here, n=total number of particles

      ntot=n
      if(myrank.eq.0) write(69,*) 'PARENT: N=',n,' NTOT=', ntot
      if (n.gt.nmax) stop 'PARENT: N>NMAX ???'
      if(ntot.gt.nmax) then
        if(myrank.eq.0)write(69,*) 'N is too large'
        stop
      endif

      if(corepts.gt.0) then
        x(1)=0.d0
        y(1)=0.d0
        z(1)=0.d0
        vx(1)=0.d0
        vy(1)=0.d0
        vz(1)=0.d0
        am(1)=amass-amtot
        if(myrank.eq.0)write(69,*) 'Mass of core = ',am(1),'Msun'
        hp(1)=hmin/1.5d0
        if(myrank.eq.0)write(69,*)'hp(core mass)',hp(1)
        cc(1)=int(2.d0*log(1.35d0*a1*
     $    (integral/(4.d0/3.d0*pi*redge**3))**(1.d0/3.d0))
     $    /log(2.d0)+15.d0)
        cc(1)=max(cc(1),0)
        cc(1)=min(cc(1),ntypes-1)
        u(1)=0.d0
        call sph_splint(rarray,muarray,muarray2,numlines,0.d0,
     $       meanmolecular(1))
        if (am(1).le.0) then
          if(myrank.eq.0)write(69,*)'CORE MASS < 0'
          stop
        endif
      else
        ammin=1.d30
        ammax=0.d0
        do i=1,n
          am(i)=am(i)/amtot*amass
          if(am(i).gt.0.01d0*amass) then
            if(myrank.eq.0)write(69,*)'warning: particle',i,
     $        'has mass',am(i)
            if(myrank.eq.0)write(69,*)'x,y,z=',x(i),y(i),z(i)
            if(myrank.eq.0)write(69,*)'r,rho=',ri,rhoi
          endif
          ammin=min(ammin,am(i))
          ammax=max(ammax,am(i))
          rtry=sqrt(x(i)**2+y(i)**2+z(i)**2)
c using two arrays to interpolate
          if(rtry.ge.coreRad2) then
            call sph_splint(rarray,rhoarraye,rhoarray1,numlines,
     $           rtry,rhoex)
          else
            if(rtry.ge.coreRad1) then
              call sph_splint(rarray,rhoarraym,rhoarray2,numlines,
     $             rtry,rhoex)
            else
              call sph_splint(rarray,rhoarrayi,rhoarray3,numlines,
     $             rtry,rhoex)
            endif
          endif
          anumden=rhoex/am(i)
          hp(i)=(3.d0/32.d0/3.1415926535897932384626d0*
     $          1.9d0*nnopt/anumden)**(1.d0/3.d0)
        enddo

     
        if(myrank.eq.0)then
          write (69,*) 'parent: min particle mass=',ammin,ammin*munit/egsol
          write (69,*) 'parent: max particle mass=',ammax,ammax*munit/egsol
        endif
      endif

c     Because the particles were distributed uniformly in the sphere, the
c     the number density is 3*n/(4*pi*radius**3) and we need to choose
c     the smoothing length hc such that 
c     4*pi*(2*hc)**3/3 * (number density)=nnopt.
c     This gives 8*hc**3 * n/radius**3= nnopt, or:
c     hc=(nnopt/(8.d0*n))**(1.d0/3.d0)*redge
c     (should give a number of nearest neighbors close to nnopt)
         
      call stride_setup
      if(ntot.gt.nmax) then
        if(myrank.eq.0)write(69,*) 'n is too large'
        stop
      endif

      if(myrank.eq.0)write(69,*)'corepts=',corepts,ntot

      if(myrank.eq.0)write(69,*) masscgs/egsol,'solar masses has cc=',
     $     int(10000*masscgs/egsol)

      if(corepts.eq.0 .and. treloff.le.0.d0) then
         if(myrank.eq.0) write(69,*)'WILL TRY TO GET CORRECT U AND W...'

         utottest=0.d0
         do i=1,n
            utottest=utottest+u(i)*am(i)
         enddo
         if(myrank.eq.0) write (69,'(a,g11.4,a,g11.4,a)')
     $        'parent: total internal energy of star was',utottest
     $        ,
     $        ', which should equal', utot2/(gravconst*munit**2/runit),
     $        '... renormalizing...'
         do i=1,n
            u(i)=u(i)*utot2/utottest/(gravconst*munit**2/runit)
         enddo

c     Do loop to get total gravitational potential energy:
         if(ngr.ne.0)then
            call gravquant      ! must call for inital set up
            call rho_and_h
            call gravforce
            if(myrank.lt.ngravprocs)then
               if(nusegpus.eq.1)then
                  call lasthalf_grav_forces(ntot, gx, gy, gz, grpot)
               else
                  call get_gravity_using_cpus
               endif
               if(ngravprocs.gt.1) then
                  mygravlength=ngrav_upper-ngrav_lower+1
                  if(myrank.ne.0)then
                     call MPI_GATHERV(grpot(ngrav_lower), mygravlength, MPI_DOUBLE_PRECISION,
     $                    grpot, gravRECVCOUNTS, gravDISPLS, MPI_DOUBLE_PRECISION, 0,
     $                    comm_worker, ierr)
                  else
                     call MPI_GATHERV(MPI_IN_PLACE, mygravlength, MPI_DOUBLE_PRECISION,
     $                    grpot, gravRECVCOUNTS, gravDISPLS, MPI_DOUBLE_PRECISION, 0,
     $                    comm_worker, ierr)
                  endif
               endif
            endif

            if(myrank.eq.0) then
               epot=0.d0
               if(ngr.ne.0)then
                  do i=1,ntot
                     epot=epot+am(i)*grpot(i)
                  enddo
                  epot=0.5d0*epot
               endif
               do irank=1,nprocs-1
                  call MPI_Send(epot, 1, MPI_DOUBLE_PRECISION,
     $                 irank,irank,MPI_COMM_WORLD,ierr)
               enddo
            else
               call MPI_Recv(epot, 1, MPI_DOUBLE_PRECISION,
     $              0,myrank,MPI_COMM_WORLD,status,ierr)
            endif
            
            if(myrank.eq.0)write (69,'(a,g11.4,a,g11.4,a)')
     $           'parent: gravitational potential energy of star was',epot
     $           ,
     $           ', which should equal', wtota2/(gravconst*munit**2/runit),
     $           '... rescaling positions...'
            do i=1,n
               x(i)=x(i)*epot/(wtota2/(gravconst*munit**2/runit))
               y(i)=y(i)*epot/(wtota2/(gravconst*munit**2/runit))
               z(i)=z(i)*epot/(wtota2/(gravconst*munit**2/runit))
               hp(i)=hp(i)*epot/(wtota2/(gravconst*munit**2/runit))
               ri=sqrt(x(i)**2+y(i)**2+z(i)**2)
               call sph_splint(rarray,muarray,muarray2,numlines,ri,
     $              meanmolecular(i))
               call sph_splint(rarray,farray, farray2,numlines,ri,
     $              ufactor(i))
               meanmolecular(i)=min(maxmu,meanmolecular(i))
               meanmolecular(i)=max(minmu,meanmolecular(i))
               if(meanmolecular(i).lt.1d-25 .or.
     $              meanmolecular(i).gt.1.d-23) then
                  if(myrank.eq.0)write(69,*)'mean molecular problem...',i,meanmolecular(i)
                  stop
               endif
            enddo
         endif
      endif

      if(myrank.eq.0 .and. corepts.gt.0)
     $     write(69,*)'properties of core:',x(1),y(1),z(1),vx(1),
     $     vy(1),vz(1),am(1),hp(1),cc(1),u(1)

      if(myrank.eq.0) then
        write(69,*) 'radius   density   u'
        do i=1,numlines,numlines/10
           write(69,*) rarray(i),rhoarray(i),uarray(i)
        enddo
      endif
      
      if(myrank.eq.0)then
         write(69,*)'rarray(1) and rarray(2)=',rarray(1),rarray(2)
         write(69,*)'rhoarray(1) and rhoarray(2)=',
     $     rhoarray(1),rhoarray(2)
         write(69,*)'rarray(numlines) =',rarray(numlines)
         write(69,*)'rhoarray(numlines) =',rhoarray(numlines)
         write(69,*)'near origin:'
         write(69,*)' r    rho'
      endif
      do irtry=0,nint(1000*rarray(2)),nint(100*rarray(2)+0.5)
         rtry=irtry/1000.d0
c using two arrays to interpolate
        if(rtry.ge.coreRad2) then
          call sph_splint(rarray,rhoarraye,rhoarray1,numlines,
     $         rtry,rhoi)
        else
          if(rtry.ge.coreRad1) then
            call sph_splint(rarray,rhoarraym,rhoarray2,numlines,
     $           rtry,rhoi)
          else
            call sph_splint(rarray,rhoarrayi,rhoarray3,numlines,
     $           rtry,rhoi)
          endif
        endif
c       if(myrank.eq.0)write(69,*) rtry,rhoi
      enddo
      if(myrank.eq.0)then
        write(69,*)'near surface:'
        write(69,*)' r    rho'
      endif
      do irtry=nint(50*radius),nint(100*radius),nint(5*radius+0.5)
        rtry=min(irtry/100.d0,radius)
c     using two arrays to interpolate
        if(rtry.ge.coreRad1) then
          call sph_splint(rarray,rhoarraye,rhoarray1,numlines,
     $         rtry,rhoi)
        else
          if(rtry.ge.coreRad2) then
            call sph_splint(rarray,rhoarraym,rhoarray2,numlines,
     $           rtry,rhoi)
          else
            call sph_splint(rarray,rhoarrayi,rhoarray3,numlines,
     $           rtry,rhoi)
          endif
        endif

        if(myrank.eq.0)write(69,*) rtry,rhoi
      enddo

c     Assign velocities (all zero)
      do i=1,ntot
         vx(i)=0.d0
         vy(i)=0.d0
         vz(i)=0.d0
         vxdot(i)=0.d0
         vydot(i)=0.d0
         vzdot(i)=0.d0
         udot(i)=0.d0
      enddo


c     Prepare leap-frog scheme for first iteration:
      call lfstart

      open(13,file='sph.particle_comp',form='unformatted')
      do i=1,ntot
         write(13) rho0(i),xmix(i),ufactor(i),T0(i),umix(i),dmix(i)
      enddo
      close(13)

      do i=1,ntot
         vy(i)=0.d0
         vz(i)=0.d0
         vxdot(i)=0.d0
         vydot(i)=0.d0
         vzdot(i)=0.d0
         udot(i)=0.d0
      enddo

      hmin=1.d30
      hmax=0.d0
      do i=1+corepts,ntot
         hmin=min(hmin,hp(i))
         hmax=max(hmax,hp(i))
      enddo

      if(myrank.eq.0)then
         write(69,*)'hmin=',hmin
         write(69,*)'hmax=',hmax
         write(69,*)'EGRGPARENT: TF=',tf,dtout
         write(69,*)'PARENT: EXITING EGRGPARENT'
      endif
      firstflag = 0

      return

      end

      subroutine splinesetupdiffplanet
c     Creates a gaseous planet from the data file yrec output
cFIXME-remove psi, maxTemp, and coreRad from common
      include 'sphu.h'
      integer numlines,i,j,coreRad1I,coreRad2I,maxTempI
      real*8 amass,radiuscgs,masscgs,radius
      real*8 integratedmass2
      integer*4 io
      parameter(io=13)
      real*8 tem(kdm),pres(kdm),rarray(kdm),xm(kdm)
      real*8 rhoarray(kdm),rhoarray1(kdm),rhoarray2(kdm),rhoarray3(kdm)
      real*8 uarray(kdm),uarray1(kdm),uarray2(kdm),uarray3(kdm)
      real*8 muarray(kdm),muarray1(kdm),muarray2(kdm),muarray3(kdm)
      real*8 psi(kdm),pdeg(kdm),pgas(kdm),prad(kdm)
      real*8 egsol,solrad,mue
      real*8 rhoarrayi(kdm),uarrayi(kdm),muarrayi(kdm)
      real*8 rhoarraym(kdm),uarraym(kdm),muarraym(kdm)
      real*8 rhoarraye(kdm),uarraye(kdm),muarraye(kdm)
      real*8 farray1(kdm),farray2(kdm),farray(kdm),farrayc(kdm)
c     Astronomical constants:
      parameter(egsol=1.9891d+33,solrad=6.9599d10,mue=2.0d0)
c     Derived constants:
      integer ndim
      parameter(ndim=9)
      real*8 xx(ndim,kdm)
      real*8 integratednum,integratednum2
      common/splinestuff/rarray,uarray,muarray,rhoarray,
     $     uarray1,uarray2,uarray3,
     $     muarray1,muarray2,muarray3,
     $     rhoarray1,rhoarray2,rhoarray3,
     $     uarraye,uarraym,uarrayi,
     $     muarraye,muarraym,muarrayi,
     $     rhoarraye,rhoarraym,rhoarrayi,
     $     amass,radius,integratednum,maxmu,minmu,numlines,maxTempI,
     $     farray1,farray2,farray,farrayc
      real*8 utot2,wtota2,sum
      real*8 maxmu,minmu
      real*8 zeroin
      real*8 redge1
      real*8 amass1,amass2
      common/forcompbest/ amass1,amass2
      real*8 temperaturefunction
      real*8 useeostable,temupperlimit,uupperlimit
      common/presarray/ pres,i
      common/hack/tem,redge1,masscgs,utot2,wtota2
      integer k
      real*8 mH
      parameter(mH = 1.67262158d-24)

c     Get profiles:
      open(io,file='sph.profile',status='old')
      i=0
      do k=1,kdm
         i=i+1
         read(io,*, end=21) xm(i),rarray(i),pres(i),rhoarray(i),
     $       tem(i),muarray(i),uarray(i)
         if(i.gt.1 .and. rarray(i).le.rarray(i-1)) then
            if(myrank.eq.0) write(69,*)'Ignoring shell',k
            i=i-1
         endif
      enddo
      if(myrank.eq.0)write(69,*)'need to increase kdm'
      stop
 21   close(io)

c     Assign outermost gridpoint
      numlines=i-1
      xm(numlines+1)=xm(numlines)
      rarray(numlines+1)=rarray(numlines)
      pres(numlines+1)=0.d0
      rhoarray(numlines+1)=0.d0
      if(myrank.eq.0)write(69,*)'number of shells=',numlines

c     Assign mean molecular weight abundances
      maxmu=0.d0
      minmu=1.d30
      coreRad1I=1
      coreRad2I=1
      do i=1,numlines-1
c       Find location of core-envelope interface
        muarray(i)=muarray(i)*mH
        if (abs(rhoarray(i)-rhoarray(i+1))/rhoarray(i+1).gt.0.2) then
           if (coreRad1I.eq.1) then
              coreRad1I = i+1
              coreRad1 = rarray(i)
           else
              if (coreRad2I.eq.1) then
                 coreRad2I = i+1
                 coreRad2  = rarray(i)
              else
                if(myrank.eq.0)write(69,*)'error with detecting core'
              endif
           endif
        endif
        maxmu=max(maxmu,muarray(i))
        minmu=min(minmu,muarray(i))
      enddo
      muarray(numlines)=muarray(numlines)*mH

c      if(myrank.eq.0)
c     $     write(69,*)'done reading eg.last1.muse_s2mm',muarray(1)
      
c     This routine matches the pressure and density profiles from the stellar evolution code.
c     This way, hydrostatic equilbrium dP/dr=-g*rho is maintained even though the 
c     stellar evolution code's EoS is different than SPH's EoS.
c     The code below solves for the temperature profile that is necessary 
c     to give the desired pressure and density profiles.
c     Then it uses this temperature to get the internal energy.
c     The equation of state that we choose to use in this routine is
c     generated by fitting the pressure, density, temperature, and
c     internal energy profiles from MESA for a range of planet masses

c     Calculate the polytrope constants using the boundary conditions
      i = coreRad2I
      envConst= pres(i)-boltz*rhoarray(i)*tem(i)/muarray(i)
      envConst= envConst/(rhoarray(i)**envGam)
c      envConst= 2.75d12
c      envConst= 2.740893d12

c     Calculate f, temperature, pressure, and internal energy
c     for each gas particle
      do i=coreRad2I,numlines-1
c        Calculate Temperature for each particle
         tem(i)=pres(i)-envConst*rhoarray(i)**envGam
         tem(i)=tem(i)/(boltz*rhoarray(i)/muarray(i))
c        Calculate f (degrees of freedom) for each particle
         farray(i)=envConst*rhoarray(i)**(envGam-1.0)/(envGam-1.0)
         farray(i)=(uarray(i)-farray(i))*muarray(i)/(boltz*tem(i))
         farrayc(i)=farray(coreRad2I)
c        Calculate internal energy
         uarray(i)=envConst*rhoarray(i)**(envGam-1.0)/(envGam-1.0)
         uarray(i)=uarray(i)+farray(i)*boltz*tem(i)/muarray(i)
      enddo


c     Calculate temperature, pressure, and internal energy
c     for each core particle
      do i=1,coreRad2I-1
         farray(i)=farray(coreRad2I)
         farrayc(i)=farray(coreRad2I)
         tem(i) = tem(coreRad2I)
      enddo
      maxTempI = coreRad2I

      if(myrank.eq.0)then
        write(69,*)pres(i-1),rhoarray(i-1),coreGam,tem(i-1)
        write(69,*)pres(i),rhoarray(i),envGam,tem(i)
        write(69,*)'Inner Core Index=',coreRad1I
        write(69,*)'Outer Core Index=',coreRad2I
        write(69,*)'Envelope Polytropic Index=',envGam
        write(69,*)'Envelope Polytropic Constant=',envConst
      endif


c     Integrate mass, number, W, and internal energy from grid points
      integratednum=0.d0
      integratedmass2=0.d0
      integratednum2=0.d0
      utot2=0.d0
      wtota2=0.d0

      do i=1,numlines
         integratedmass2=integratedmass2+rhoarray(i)*4.d0*pi*
     $      rarray(i)**2*0.5d0*(rarray(i+1)-rarray(i-1))
         integratednum2=integratednum2+rhoarray(i)*4.d0*pi*
     $      rarray(i)**2*0.5d0*(rarray(i+1)-rarray(i-1))/muarray(i)
         wtota2=wtota2-3.d0*pres(i)*4.d0*pi*
     $      rarray(i)**2*0.5d0*(rarray(i+1)-rarray(i-1))
         utot2=utot2+uarray(i)*rhoarray(i)*4.d0*pi*
     $      rarray(i)**2*0.5d0*(rarray(i+1)-rarray(i-1))
      enddo
      if(myrank.eq.0)then
         write(69,*)'central u=',uarray(1)
         write(69,*)'central temperature=',tem(1)
         write(69,*)'Core Radius = ',rarray(i)
      endif

      if (myrank.eq.0) then
         write(69,*)'mass from integrating rho profile=',
     $     integratedmass2/egsol,'msun'
         write(69,*)'number from integrating=',integratednum2
         write(69,*)
     $     'utot (in sph units)=',utot2/(gravconst*munit**2/runit)
         write(69,*)
     $     'wtot (in sph units)=',wtota2/(gravconst*munit**2/runit)
      endif

      masscgs=xm(numlines)
      if (myrank.eq.0)
     $   write(69,*)'total mass=',xm(numlines),xm(numlines)/egsol

      do i=1,numlines
         if (xm(i).gt.0.85d0*masscgs) then
            redge1=rarray(i)/runit
            goto 143
         endif
      enddo
 143  continue

      radiuscgs=rarray(numlines)
      radius=radiuscgs/runit
      amass=masscgs/munit
      if (myrank.eq.0) then
         write(69,'(3(a,g11.4),a)')
     $      '  mass =',masscgs,  'grams =',masscgs/egsol,'msun =',
     $      amass,'code units'
         write(69,'(3(a,g11.4),a)')
     $      'radius =',radiuscgs,'cm    =',radiuscgs/solrad,'rsun =',
     $      radius,'code units'
         write(69,*)'   i   rarray(i)  rhoarray(i) pres(i) ',
     $      '     muarray(i) tem(i)     uarray(i)'
         do i=1,numlines,numlines/10
            write(69,'(i5,9g12.5)') i,
     $      rarray(i), rhoarray(i), pres(i),
     $      muarray(i),tem(i), uarray(i)
        enddo
      endif

      open(21, file='pre-parent.sph',status='unknown')
      do i=1,numlines
         write(21,'(9g15.7)') rarray(i),pres(i),
     $        pres(i)*(radiuscgs**2/masscgs)**2/gravconst/
     $        rhoarray(i),tem(i),muarray(i),uarray(i),
     $        pdeg(i),pgas(i),prad(i)
      enddo
      close(21)

c     Convert arrays and values to code units
      coreRad1 = coreRad1/runit
      coreRad2 = coreRad2/runit
      do i=1,numlines
         rarray(i)=rarray(i)/runit
         rhoarray(i)=rhoarray(i)*(runit**3/munit)
         uarray(i)=uarray(i)/(gravconst*munit)*runit
         if (i.gt.1) then
            if(rarray(i).lt.rarray(i-1) .and. myrank.eq.0) then
               write(69,*) 'radius warning'
               write(69,*) i-1,rarray(i-1),rhoarray(i-1),uarray(i-1),
     $            pres(i-1)
               write(69,*) i,rarray(i),rhoarray(i),uarray(i),pres(i)
            endif
            if (rhoarray(i).gt.rhoarray(i-1) .and. myrank.eq.0)
     $         write(69,'(a,3g11.4)')
     $            'rho increases outward warning near shell',
     $            i,rhoarray(i-1),rhoarray(i)
            if (pres(i).gt.pres(i-1) .and. myrank.eq.0)
     $         write(69,*) 'pressure warning'
         endif
      enddo

      open(21, file='parent.sph',status='unknown')
      do i=1,numlines
         write(21,'(6g15.7)') rarray(i),pres(i)/punit,
     $        rhoarray(i),tem(i),muarray(i),uarray(i)
      enddo
      close(21)

      CALL FLUSH(69)

c     Set envelope and core densities at the interface (code units)
      rhoei = rhoarray(coreRad2I)
      u_ei  = uarray(coreRad2I)
      rhoci = rhoarray(coreRad2I-1)
      u_ci  = uarray(coreRad2I-1)
c      rhoci = rhoei+(rhoarray(coreRad2I-1)-rhoei)*0.75d0
      rhocm = rhoarray(coreRad1I)
      u_cm  = uarray(coreRad1I)
      rhocc = rhoarray(coreRad1I-1)
      u_cc  = uarray(coreRad1I-1)

c     output the interface densities and polytropic constants and indices
!      if (myrank.eq.0) then
!         open(13,file='sph.rhoparams',status='unknown')
!         write(13,'(10g15.7)')rhoei,rhoci,rhocm,rhocc,envConst,envGam,
!     $      u_ei,u_ci,u_cm,u_cc
!         close(13)
!      endif

c     Generate arrays used for splines
      do i=1,numlines
         if (i.lt.coreRad1I) then
c           Iron core profile
            rhoarrayi(i) = rhoarray(i)
            uarrayi(i)   = uarray(i)
            muarrayi(i)  = muarray(coreRad2I)

c           Silicate mantle profile
            rhoarraym(i) = rhoarray(coreRad1I)
            uarraym(i)   = uarray(coreRad1I)
            muarraym(i)  = muarray(coreRad2I)

c           Gas envelope profile
            rhoarraye(i) = rhoarray(coreRad2I)
            uarraye(i)   = uarray(coreRad2I)
            muarraye(i)  = muarray(coreRad2I)
         else
            rhoarrayi(i) = rhoarray(coreRad1I-1)
            uarrayi(i)   = uarray(coreRad1I-1)
            muarrayi(i)  = muarray(coreRad1I-1)

            if (i.lt.coreRad2I) then
              rhoarraym(i) = rhoarray(i)
              uarraym(i)   = uarray(i)
              muarraym(i)  = muarray(coreRad2I)

              rhoarraye(i) = rhoarray(coreRad2I)
              uarraye(i)   = uarray(coreRad2I)
              muarraye(i)  = muarray(coreRad2I)
            else
              rhoarraym(i) = rhoarray(coreRad2I-1)
              uarraym(i)   = uarray(coreRad2I-1)
              muarraym(i)  = muarray(coreRad2I-1)

              rhoarraye(i) = rhoarray(i)
              uarraye(i)   = uarray(i)
              muarraye(i)  = muarray(i)
              
            endif

c           Gas profile
         endif
      enddo
      call sph_spline(rarray,rhoarraye,numlines,1.d30,1.d30,rhoarray1)
      call sph_spline(rarray,rhoarraym,numlines,1.d30,1.d30,rhoarray2)
      call sph_spline(rarray,rhoarrayi,numlines,1.d30,1.d30,rhoarray3)
      call sph_spline(rarray,uarraye,  numlines,1.d30,1.d30,uarray1)
      call sph_spline(rarray,uarraym,  numlines,1.d30,1.d30,uarray2)
      call sph_spline(rarray,uarrayi,  numlines,1.d30,1.d30,uarray3)
      call sph_spline(rarray,muarraye, numlines,1.d30,1.d30,muarray1)
      call sph_spline(rarray,muarraym, numlines,1.d30,1.d30,muarray2)
      call sph_spline(rarray,muarrayi, numlines,1.d30,1.d30,muarray3)
      call sph_spline(rarray,farray,   numlines,1.d30,1.d30,farray1)
      call sph_spline(rarray,farrayc,  numlines,1.d30,1.d30,farray2)


      CALL FLUSH(69)
      end
