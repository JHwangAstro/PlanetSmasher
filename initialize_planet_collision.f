      SUBROUTINE jupiters                                         
********************************************************
C     Sets up a system of two planets using previously
c     calculated relaxed planet models. The planets are
c     made up of a core approximated as a polytrope and
c     an envelope with user-generated profile
c     Called by INIT
c     Calls GRAVQUANT,LFSTART
**********************************************************
      INCLUDE 'sphu.h'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      integer n2,i,nchk,corepts
      integer nnoptold,noutold,nitold,navold,ngrold,nrelaxold
      real*8 hminold,hmaxold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,tskipahead,trelaxold,dtold
      integer istart,corepts1,corepts2
      common/core/corepts1,corepts2
      logical resetsep0,twofiles,planets
      real*8 am1chk,am2chk,realdummy1,realdummy2,realdummy3
      real*8 deltax1,deltay1,deltaz1,deltax2,deltay2,deltaz2
      real*8 egsol,solrad
      parameter(egsol=1.9891d+33,solrad=6.9599d10)
      real*8 deltaxbin,deltaybin,deltazbin
      character*7 dummy
      integer*4 io,iostatus
      real*8 vxorb1,vyorb1,vzorb1
      real*8 vxorb2,vyorb2,vzorb2
      real*8 xorb1, yorb1, zorb1
      real*8 xorb2, yorb2, zorb2
      real*8 mass1,mass2
      parameter(io=13)
      real*8 rhoei0,rhoci0,rhocm0,rhocc0,rhocgs,rhocgs0,ucgs
      real*8 gam1,gam2,rho1,rho2,const1,const2
      real*8 rho1max,rho2max
      parameter(rho1max=20.0)
      parameter(rho2max=6.5)
      parameter(gam1 = 1./0.528)
      parameter(gam2 = 1./0.549)
      parameter(rho1 = 8.300)
      parameter(rho2 = 4.260)

      const1 = (1.0/(3.49*10.**(-(1./gam1+6.))))**gam1
      const2 = (1.0/(1.27*10.**(-(1./gam2+6.))))**gam2

      if(myrank.eq.0) write(69,*) 'STARTING ROUTINE'

c     Get core-envelope interface parameters
      if (coreGam.ne.0) then
         open(io,file='sph.rhoparams',status='old',action='read')
         read(io,'(10g15.7)')rhoei,rhoci,rhocm,rhocc,envConst,envGam,
     $      u_ei,u_ci,u_cm,u_cc
!         read(io,'(6g15.7)')
!     $       rhoei,rhoci,rhocm,rhocc,envConst,envGam
         close(io)
         open(io,file='sph.rhoparams2',status='old',action='read')
         read(io,'(10g15.7)')rhoei2,rhoci2,rhocm2,rhocc2,envConst2,
     $      envGam2,u_ei2,u_ci2,u_cm2,u_cc2
!         read(io,'(6g15.7)')
!     $       rhoei2,rhoci2,rhocm2,rhocc2,envConst2,envGam2
         close(io)
      endif

c     Get orbital elements
      open(io,file='sph.orbits',status='old',action='read')
      read(io,*,end=21) mass1,xorb1,yorb1,zorb1,vxorb1,vyorb1,vzorb1
      read(io,*,end=21) mass2,xorb2,yorb2,zorb2,vxorb2,vyorb2,vzorb2
 21   close(io)

      if(myrank.eq.0) write(69,*) 'p1 coordinates: ',xorb1,yorb1,zorb1
      if(myrank.eq.0) write(69,*) 'p1 velocities: ',vxorb1,vyorb1,vzorb1
      if(myrank.eq.0) write(69,*) 'p2 coordinates: ',xorb2,yorb2,zorb2
      if(myrank.eq.0) write(69,*) 'p2 velocities: ',vxorb2,vyorb2,vzorb2

      corepts1=0
      corepts2=0
      if(myrank.eq.0) write (69,*) 'JUPITERS: READING START FILES ...'
      
      open(12,file='sph.start1u',form='unformatted')
c     (The following READ sequence must match exactly the WRITE sequence
c     used in subroutine DUMP)
      read(12) n1,nnoptold,hminold,hmaxold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,
     $     navold,alphaold,betaold,tskipahead,
     $     ngrold,
     $     nrelaxold,trelaxold,dtold,omega2
      am1=0.d0
      do i=1,n1
         read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $        vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        grpot(i),meanmolecular(i),
     $        cc(i), realdummy1

c        Assign particle type         
!         cc(i) = cc(i) + 2
!         if (rho0(i).gt.(rhoci+rhocm)/2.) cc(i) = 2

         am1=am1+am(i)
         if(hp(i).le.0.d0 .or. u(i).eq.0.d0) then
            if(myrank.eq.0) write(69,*)'particle',i,
     $           'is a corepoint of mass',am(i)
         endif
      enddo
      read(12) nchk
      close(12)

      rhoei0 = rhoei*munit/runit**3.d0
      rhoci0 = rhoci*munit/runit**3.d0
      rhocm0 = rhocm*munit/runit**3.d0
      rhocc0 = rhocc*munit/runit**3.d0
      open(12,file='sph.particle_comp',form='unformatted')
      do i=1,n1
         read(12) rho0(i),xmix(i),ufactor(i),T0(i),umix(i),dmix(i)
         rhocgs  = rho(i)*munit/runit**3.d0
         rhocgs0 = rho0(i)*munit/runit**3.d0
         ucgs = u(i)*gravconst*munit/runit
         ! for mantle-gas particlees
         if (rhocgs0.le.(rhoci0+rhocm0)/2.) then
            call getDensityComponents(xmix(i),rhocgs,ucgs,const2,
     $           envConst,gam2,envGam,rho2,ufactor(i),
     $           meanmolecular(i),T0(i),rho2max,rhocgs0,rhoei0,
     $           rhoci0,0,rhoh(i),rhol(i))
         ! for core-mantle particles
         else
            call getDensityComponents(xmix(i),rhocgs,ucgs,const1,
     $           const2,gam1,gam2,rho1,rho2,meanmolecular(i),T0(i),
     $           rho1max,rhocgs0,rhocm0,rhocc0,1,rhoh(i),rhol(i))
         endif
      enddo
      close(12)

      if(myrank.eq.0) write(69,*) 'mass1= ', am1

      if (nchk.ne.n1) stop 'JUPITERS: PROBLEM WITH FILE sph.start1u'

      if(myrank.eq.0) write(69,*)'n1=',n1

      inquire(file='sph.start2u',exist=twofiles)
      if(twofiles) then
         open(12,file='sph.start2u',form='unformatted')
c     (The following READ sequence must match exactly the WRITE sequence
c     used in subroutine DUMP)
         read(12) n2,nnoptold,hminold,hmaxold,sep0old,
     $        tfold,dtoutold,noutold,nitold,told,navold,
     $        alphaold,betaold,tskipahead,ngrold,
     $        nrelaxold,trelaxold,dtold
         if(myrank.eq.0) write(69,*)'n2=',n2
         am2=0.d0
         istart=n1+1-corepts1
         do i=istart,istart+n2-1
            read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $           vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $           grpot(i),meanmolecular(i),cc(i),realdummy1
            cc(i) = cc(i)+2!3
!            if (rho0(i).gt.(rhoci2+rhocm2)/2.) cc(i) = 4

            am2=am2+am(i)
            if(hp(i).le.0.d0 .or. u(i).eq.0.d0) then
               if(myrank.eq.0) write(69,*)'particle',i,
     $              'is a corepoint of mass',am(i)
            endif
         enddo
         read(12) nchk
         close(12)

         rhoei0 = rhoei2*munit/runit**3.d0
         rhoci0 = rhoci2*munit/runit**3.d0
         rhocm0 = rhocm2*munit/runit**3.d0
         rhocc0 = rhocc2*munit/runit**3.d0
         open(12,file='sph.particle_comp2',form='unformatted')
         do i=istart,istart+n2-1
            read(12) rho0(i),xmix(i),ufactor(i),T0(i),umix(i),dmix(i)
            rhocgs  = rho(i)*munit/runit**3.d0
            rhocgs0 = rho0(i)*munit/runit**3.d0
            ucgs = u(i)*gravconst*munit/runit
            ! for mantle-gas particlees
            if (rhocgs0.le.(rhoci0+rhocm0)/2.) then
               call getDensityComponents(xmix(i),rhocgs,ucgs,const2,
     $              envConst2,gam2,envGam2,rho2,ufactor(i),
     $              meanmolecular(i),T0(i),rho2max,rhocgs0,rhoei0,
     $              rhoci0,0,rhoh(i),rhol(i))
            ! for core-mantle particles
            else
               call getDensityComponents(xmix(i),rhocgs,ucgs,const1,
     $              const2,gam1,gam2,rho1,rho2,meanmolecular(i),T0(i),
     $              rho1max,rhocgs0,rhocm0,rhocc0,1,rhoh(i),rhol(i))
            endif
         enddo
         close(12)

         if (nchk.ne.n2) then
            write(69,*) 'JUPITERS: PROBLEM WITH FILE sph.start2u'
            STOP
         endif
      endif
     
      if(myrank.eq.0) write(69,*) 'mass2= ', am2

c     check for existence of planets file:
      inquire(file='sph.planets',exist=planets)
      
!      i = n1+n2
c     if dump file exists, read it and restart from it:
      if (planets) then
        open(12,file='sph.planets',status='old',action='read')
        do i = n1+n2+1, n1+n2+nplanets
          gx(i)=0.d0
          gy(i)=0.d0
          gz(i)=0.d0
          vxdot(i)=0.d0
          vydot(i)=0.d0
          vzdot(i)=0.d0

          u(i)=0.d0
          udot(i)=0.d0
          if(hmin.gt.0.d0) then
             hp(i)=hmin
          else
             hp(i)=hp(n1)
          endif
!          read(12,'(7g15.7)') am(i),x(i),y(i),z(i),vx(i),vy(i),vz(i)
          read(12,*,end=22) am(i),x(i),y(i),z(i),vx(i),vy(i),vz(i)
c         Create planets
        enddo
 22     close(12)

!        ntot=n1+n2+1
!        ntot = i
        if (ntot.gt.nmax) then
           write(69,*) 'MUST INCREASE NMAX...'
           stop
        endif
      endif

c     Create the central star
!      i = i + n_planets+1
      i = n1+n2+nplanets+1
      if(myrank.eq.0) write(69,*)'solar_mass',egsol
      am(i)=starmass!1.071
      x(i)=0.d0
      y(i)=0.d0
      z(i)=0.d0
      vx(i)=0.d0
      vy(i)=0.d0
      vz(i)=0.d0
      gx(i)=0.d0
      gy(i)=0.d0
      gz(i)=0.d0
      vxdot(i)=0.d0
      vydot(i)=0.d0
      vzdot(i)=0.d0
      u(i)=0.d0
      udot(i)=0.d0
      if(hmin.gt.0.d0) then
         hp(i)=hmin
      else
         hp(i)=hp(n1)
      endif
      corepts=corepts1+corepts2
      ntot=n1+n2+1+nplanets
!      ntot=ntot-1!+nplanets
      if (ntot.gt.nmax) then
         write(69,*) 'MUST INCREASE NMAX...'
         stop
      endif

      n=ntot-corepts

      if(resetsep0) then
         if(dabs(am1*munit/egsol/am1chk-1.d0).gt.1.d-8 .or.
     $        dabs(am2*munit/egsol/am2chk-1.d0).gt.1.d-8)then
            write(69,*)'mass(es) in sph.start?u does not match with'
            write(69,*)'mass(es) in input.bs'
            stop
         endif
      endif

      if(ntot.gt.nmax) then
         write(69,*)'ERROR: N1+N2>NMAX'
         stop
      endif

      if(myrank.eq.0) write (69,*)'JUPITERS: N=',n,'NTOT=',ntot
C     Place planets on their orbits
      if(myrank.eq.0) write(69,*)'SEP0=',sep0
      DO I=1,n1-corepts1
         X(I)=X(I)+xorb1
         Y(I)=Y(I)+yorb1
         Z(I)=Z(I)+zorb1
         vx(I)=vx(I)+vxorb1
         vy(I)=vy(I)+vyorb1
         vz(I)=vz(I)+vzorb1
      ENDDO
      DO I=n1+1-corepts1,ntot-corepts1-1-nplanets
         X(I)=X(I)+xorb2
         Y(I)=Y(I)+yorb2
         Z(I)=Z(I)+zorb2
         vx(I)=vx(I)+vxorb2
         vy(I)=vy(I)+vyorb2
         vz(I)=vz(I)+vzorb2
      ENDDO

      open(12,file='sph.particle_comp0',form='unformatted')
!      do i=1,ntot-corepts1-1
      do i=1,n1+n2
         write(12) rho0(i),xmix(i),ufactor(i),T0(i),umix(i),dmix(i)
      enddo
      close(12)

      if(myrank.eq.0) then
         write(69,*)'Jupiter 1 has mass',am1,n1-corepts1,ntot-corepts1+1,ntot
         write(69,*)'Jupiter 1 has nominal mass',mass1,(mass1-am1)/mass1
         write(69,*)'Jupiter 2 has mass',am2,n1+1-corepts1,ntot-corepts1
         write(69,*)'Jupiter 2 has nominal mass',mass2,(mass2-am2)/mass2
      endif


      call stride_setup

c     Prepare leap-frog scheme for first iteration:
      call lfstart

      return
      end
