      subroutine init
c     initialization of run (new or restart)
      include 'sphu.h'
      include 'mpif.h'
      real*8 xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2
      common/centersofmass/xcm1,ycm1,zcm1,xcm2,ycm2,zcm2,am1,am2

      integer i,j,nchk,nrelaxold,nnoptold,navold,ngrold
      integer corepts, istart, n2,nitold,noutold
      real*8 hminold,hmaxold,tfold,dtoutold,sep0old,told,
     $     alphaold,betaold,trelaxold,dtold,realdummy1
      logical restart,twofiles
      character*3 iname
      namelist /initt/ iname
      common /inittcom/ iname
      integer corepts1
      common/core/corepts1
      real*8 divv(nmax)
      common/commdivv/divv
      real*8 amass1,amass2
      common/forcompbest/ amass1,amass2
      real*8 tskipahead
      common /skipcomm/ tskipahead
      real*8 xcm,ycm,zcm,vxcm,vycm,vzcm,amtot
      integer*4 io
      real*8 utab2(10000000)
      real*8 utab3(10000000)
      parameter(io=13)
      real*8 rhocgs,rhocgs0,ucgs,xm,xc,temp_rho
      real*8 rhoei0,rhoci0,rhocm0,rhocc0
      real*8 poly2,poly1,ucore2,ucore1
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

      firstflag = 0

      call energy_setup

c     read in hosts, find gravhost:
c      call rank_setup
c     set up hydro communicator:
c      call hydrocomm_setup
c     initialize "my" (parallel) arrays:
c      call myinit

c     get parameters from input file:
      call get_input
c     compute look-up tables:
      call tabulinit

      if(myrank.eq.0) write(69,*)'nav, nintvar, neos, ncooling=',
     $     nav,nintvar,neos,ncooling
c
      punit=gravconst*(munit/runit**2)**2

      if(myrank.eq.0) then
        write(69,*) 'mass unit=',munit
        write(69,*) 'radius unit=',runit
        write(69,*) 'pressure unit=',punit
      endif  

c     check for existence of restart dump file:
      inquire(file='restartrad.sph',exist=restart)

c     if dump file exists, read it and restart from it:
      if (restart) then
c        Get core-envelope interface parameters
         if (coreGam.ne.0) then
            open(io,file='sph.rhoparams',status='old',action='read')
            read(io,'(10g15.7)')rhoei,rhoci,rhocm,rhocc,envConst,
     $             envGam,u_ei,u_ci,u_cm,u_cc
            close(io)
            if (NRELAX.eq.0) then
               open(io,file='sph.rhoparams2',status='old',action='read')
               read(io,'(10g15.7)')rhoei2,rhoci2,rhocm2,rhocc2,
     $             envConst2,envGam2,u_ei2,u_ci2,u_cm2,u_cc2
               close(io)
            endif
         endif

         if(myrank.eq.0) then
            write(69,*)'init: continuing run'
            write(69,*)'init: reading restart dump file ...'
         endif
         corepts=0

         open(12,file='restartrad.sph',form='unformatted')
c     (the following read sequence must match exactly the write sequence
c     used in subroutine dump)
         read(12) ntot,nnoptold,hminold,hmaxold,sep0,tfold,
     $        dtoutold,nout,nit,t,navold,alphaold,betaold,tskipahead,
     $        ngrold,nrelaxold,trelaxold,dt,omega2

         if(dtoutold.ne.dtout) then
            if(myrank.eq.0)write(69,*)'NOUT from restartrad.sph is',NOUT
            NOUT = T/DTOUT-1
            if(myrank.eq.0)write(69,*)'NOUT changed to',NOUT
         endif

         if(tskipahead.gt.0) tskipahead=-tskipahead
         tf=sign(max(abs(tfold),abs(tf)),tf)

         if(myrank.eq.0) then
            write(69,*) myrank,'tskipahead=',tskipahead,'tf=',tf
         endif

         omeg=sqrt(omega2)
         if(myrank.eq.0) write(69,*)'restartrad: omeg=',omeg
         do i=1,ntot
            read(12) x(i),y(i),z(i),am(i),hp(i),rho(i),
     $           vx(i),vy(i),vz(i),vxdot(i),vydot(i),vzdot(i),
     $           u(i),udot(i),!gx(i),gy(i),gz(i),
     $           grpot(i),meanmolecular(i),cc(i),divv(i)
         enddo
         read(12) nchk
         close(12)

         if (NRELAX.eq.1) then
            open(12,file='sph.particle_comp',form='unformatted')
            do i=1,ntot
               read(12) rho0(i),xmix(i),ufactor(i),T0(i),umix(i),dmix(i)
               u(i)   = u(i)*gravconst*munit/runit
               rhocgs = rho(i)*munit/runit**3.d0
               rhocgs0= rho0(i)*munit/runit**3.d0
               rhoei0 = rhoei*munit/runit**3.d0
               rhoci0 = rhoci*munit/runit**3.d0
               rhocm0 = rhocm*munit/runit**3.d0
               rhocc0 = rhocc*munit/runit**3.d0
               if (cc(i).eq.1) then
!                  if (xmix(i).ne.0.and.xmix(i).ne.1) then
!                  call getCompositions2(rhocgs,const2,envConst,gam2,
!     $                 envGam,rho2,ufactor(i),meanmolecular(i),
!     $                 T0(i),rhocgs0,rhoei0,rhoci0,0,xmix(i),u(i))
!                  endif
                  call getDensityComponents(xmix(i),rhocgs,u(i),
     $                 const2,envConst,gam2,envGam,rho2,ufactor(i),
     $                 meanmolecular(i),T0(i),rho2max,rhocgs0,
     $                 rhoei0,rhoci0,0,rhoh(i),rhol(i))
               else
!                  if (xmix(i).ne.0.and.xmix(i).ne.1) then
!                  call getCompositions2(rhocgs,const1,const2,gam1,
!     $                 gam2,rho1,rho2,meanmolecular(i),
!     $                 T0(i),rhocgs0,rhocm0,rhocc0,1,xmix(i),u(i))
!                  endif
                  call getDensityComponents(xmix(i),rhocgs,u(i),
     $                 const1,const2,gam1,gam2,rho1,rho2,
     $                 meanmolecular(i),T0(i),rho1max,rhocgs0,
     $                 rhocm0,rhocc0,1,rhoh(i),rhol(i))
               endif
               u(i)  = u(i)/gravconst/munit*runit
               dmix(i) = rho(i)
            enddo
            close(12)
            open(13,file='sph.particle_comp_mod',form='unformatted')
            do i=1,ntot-1-nplanets
               write(13)rho0(i),xmix(i),ufactor(i),T0(i),u(i),dmix(i)
            enddo
            close(13)
         endif

         if (NRELAX.eq.0)then
            open(12,file='sph.particle_comp0',form='unformatted')
            do i=1,ntot-1-nplanets
               ucgs = u(i)*gravconst*munit/runit
               read(12) rho0(i),xmix(i),ufactor(i),T0(i),umix(i),dmix(i)
               ! read in planet 1
               if (cc(i).le.2) then
                  rhoei0 = rhoei*munit/runit**3.d0
                  rhoci0 = rhoci*munit/runit**3.d0
                  rhocm0 = rhocm*munit/runit**3.d0
                  rhocc0 = rhocc*munit/runit**3.d0
                  if (cc(i).eq.1) then
                      call getDensityComponents(xmix(i),rhocgs,ucgs,
     $                     const2,envConst,gam2,envGam,rho2,ufactor(i),
     $                     meanmolecular(i),T0(i),rho2max,rhocgs0,
     $                     rhoei0,rhoci0,0,rhoh(i),rhol(i))
                  else
                      call getDensityComponents(xmix(i),rhocgs,ucgs,
     $                     const1,const2,gam1,gam2,rho1,rho2,
     $                     meanmolecular(i),T0(i),rho1max,rhocgs0,
     $                     rhocm0,rhocc0,1,rhoh(i),rhol(i))
                  endif
               ! read in planet 2
               else
                  rhoei0 = rhoei2*munit/runit**3.d0
                  rhoci0 = rhoci2*munit/runit**3.d0
                  rhocm0 = rhocm2*munit/runit**3.d0
                  rhocc0 = rhocc2*munit/runit**3.d0
                  if (cc(i).eq.3) then
                      call getDensityComponents(xmix(i),rhocgs,ucgs,
     $                    const2,envConst2,gam2,envGam2,rho2,ufactor(i),
     $                     meanmolecular(i),T0(i),rho2max,rhocgs0,
     $                     rhoei0,rhoci0,0,rhoh(i),rhol(i))
                  else
                      call getDensityComponents(xmix(i),rhocgs,ucgs,
     $                     const1,const2,gam1,gam2,rho1,rho2,
     $                     meanmolecular(i),T0(i),rho1max,rhocgs0,
     $                     rhocm0,rhocc0,1,rhoh(i),rhol(i))
                  endif
               endif
            enddo
            close(12)
         endif

         corepts1=(corepts+1)/2
         n=ntot-corepts
         if(nrelax.ge.2) then

c     Use these lines of code if doing binaries for triple collisions:
            n1=1
            do while(cc(n1).eq.cc(1))
               n1=n1+1
               if(u(n1).eq.0.d0) goto 168
            enddo
 168        continue
            n1=n1-1
            if(myrank.eq.0) write(69,*)
     $         'number of particles in star 1 is n1=',n1
         else
            n1=ntot
         endif
         
c         if(nrelax.eq.1 .and. tresplintmuoff.gt.t) then
c            if(myrank.eq.0) write(69,*)'calling splinesetup...'
c            call splinesetup
c         endif

c         if(trelax.ne.trelaxold) then
c            if(myrank.eq.0) write(69,*)'***using trelax=',trelaxold
c            trelax=trelaxold
c         endif
c         if(nrelax.ne.nrelaxold) then
c            if(myrank.eq.0) write(69,*)'***using nrelax=',nrelaxold
c            nrelax=nrelaxold                  
c         endif
         if(nrelax.ne.2 .and. nrelaxold.eq.2) then
            nout=0
            nit=0
            t=0.d0
            if(treloff.le.0) then
               if(myrank.eq.0) then
                  write(69,*)
     $                 'Assuming a dynamical calculation is starting.'
                  write(69,*)'Resetting nout=nit=t=0'
                  write(69,*)'Setting trelax=1.d30'
               endif
               trelax=1.d30
               if(alpha.le.0.d0 .and. beta.le.0.d0) then
                  alpha=1.d0
                  beta=2.d0
                  if(myrank.eq.0) then
                     write(69,*)'Resetting alpha,beta=',alpha,beta
                  endif
               endif
            endif
         endif

         if(omega2.eq.0.d0 .or. (nrelax.eq.3 .and. treloff.le.0))
     $        gonedynamic=.true.
         if(myrank.eq.0) write(69,*)'gonedynamic=',gonedynamic
         if(myrank.eq.0) write(69,*)'n=',n,'ntot=',ntot
         
         if(nchk.ne.ntot) stop 'init: problem with dump file ???'
         
         if(myrank.eq.0) write(69,*) 'about to call stride_setup'
         
         call stride_setup
         
         dth=0.5d0*dt

         call gravquant         ! must call for inital set up
         
c         time1=rtc()
c         call linkedlists
c         time2=rtc()
c         if(myrank.eq.0) write(69,*) 'linkedlists lasted',time2-time1
c         do i=n_lower,n_upper
c            call newnene(i)
c         enddo
c         time1=rtc()
c         if(myrank.eq.0) write(69,*) 'newnene calls lasted',time1-time2

         call pressure
      
         if(nrelaxold.ge.2) then
            if(nrelax.lt.2) then
               if(myrank.eq.0) write(69,*) 
     $               'did you mean to have nrelax>=2?'
c               stop
            endif
            call getcoms
            if(.not.gonedynamic) then
               if(myrank.eq.0) write(69,*) 'sep0=',sep0
               sepfinal=sep0
               if(myrank.eq.0) write(69,*) 'sepfinal=',sepfinal
c     use the following lines if want to restart a scanning run into a fixed
c     separation run that relaxes until t=treloff before going
c     dynamic:
ccccccccccccccc               treloff=dble(nint(t))
               sep0=dabs(xcm2-xcm1)
               sepfinal=sep0
               if(myrank.eq.0) write(69,*) 'sep0 & sepfinal set to',sep0
            endif

            hmax=hmaxold
            if(myrank.eq.0) write(69,*) 'set hmax=hmaxold=',hmax            

         endif
   
         if(myrank.eq.0) write(69,*)'init:         ... done(ish)'
         
c     otherwise create new data set:
      else
         nout=0
         nit=0
c         write(69,*) "setting t=0", myrank
         t=0.d0
c     initialize output parameters:
c     get 3-letter code for type of initial condition from init file
         open(12,file='sph.init',err=100)
         read(12,initt)
         close(12)
         if(myrank.eq.0) write(69,*)'init: new run, iname=',iname
         if(iname.eq.'1es') then
            call polyes
         elseif(iname.eq.'1mc')then
            call polymces
         elseif(iname.eq.'2cr')then
            call corotating
         elseif(iname.eq.'chj')then
            call corothotjup
         elseif(iname.eq.'ppp')then
            call jupiters
         elseif(iname.eq.'bps')then
            call bps            !binary plus single
         elseif(iname.eq.'bph')then
            call bpbh           !binary plus black hole
         elseif(iname.eq.'hyp')then
            call hyperbolic
         elseif(iname.eq.'hbs')then
            call hyperbolic_binary_single
         elseif(iname.eq.'erg')then
            call egrgparent
            if(myrank.eq.0) write(69,*)'init: tf=',tf,dtout
         elseif(iname.eq.'plt')then
            call planetgen
         elseif(iname.eq.'dpl')then
            call diffplanetgen
         elseif(iname.eq.'tri')then
            call triple
         elseif(iname.eq.'bhe')then
            call smbh
         elseif(iname.eq.'rin')then
            call readin
         elseif(iname.eq.'grs')then
            call grsph
         elseif(iname.eq.'txt')then
            call asciiimage
         else
            if(myrank.eq.0) write(69,*)'init: unknown name',iname
            stop
         endif
      endif

      if(n.gt.nmax)then
         if(myrank.eq.0) write(69,*)'init: increase nmax to',n
         stop
      endif

      if(myrank.eq.0)then
         xcm=0d0
         ycm=0d0
         zcm=0d0
         vxcm=0d0
         vycm=0d0
         vzcm=0d0
         amtot=0d0
         do i=1,ntot
            xcm=xcm+am(i)*x(i)
            ycm=ycm+am(i)*y(i)
            zcm=zcm+am(i)*z(i)
            vxcm=vxcm+am(i)*vx(i)
            vycm=vycm+am(i)*vy(i)
            vzcm=vzcm+am(i)*vz(i)
            amtot=amtot+am(i)
         enddo
         xcm=xcm/amtot
         ycm=ycm/amtot
         zcm=zcm/amtot
         vxcm=vxcm/amtot
         vycm=vycm/amtot
         vzcm=vzcm/amtot
         write(69,*)'Total mass=',amtot
         write(69,'(a,9g15.7)')'Center of mass position=',xcm,ycm,zcm
         write(69,'(a,9g15.7)')'Center of mass velocity=',vxcm,vycm,vzcm
      endif
         
      do i=n+1,ntot
         nn(i)=0
      enddo

c     write run parameters:
c HERE !!!!!
      if(myrank.eq.0) then
         write(69,*)'init: t=',t,' nit=',nit
         write(69,101) n,nnopt,hmin,hmax,dtout,nout,tf,sep0,nav,
     $     alpha,beta,ngr,nrelax,trelax,dt,ntot-n
 101     format (' init: parameters for this ideal + radiation run:',/,
     $        ' n=',i7,' nnopt=',i4,' hmin=',g12.4,' hmax=',g12.4,/,
     $        ' dtout=',g10.3,' nout=',i4,' tf=',g10.3,/,
     $        ' sep0=',g12.4,/,
     $        ' nav=',i2,' alpha=',f6.2,' beta=',f6.2,/,
     $        ' ngr=',i3,/,
     $        ' nrelax=',i2,' trelax=',g12.4,
     $        ' dt=',g12.4,
     $        ' corepts=',i2,/)
         if(nintvar.eq.1) then
            write(69,*)'integrating entropic variable a'
         elseif(nintvar.eq.2)then
            write(69,*)'integrating energy density u'
         else
            write(69,*)'must integrate either a or u'
            stop
         endif
         if(neos.eq.0) then
            write(69,*)'polytropic equation of state'
         elseif(neos.eq.1)then
            write(69,*)'ideal gas + radiation pressure EOS'
         elseif(neos.eq.2)then
            write(69,*)'tabulated equation of state'
            if(nintvar.ne.2) then
               write(69,*)'must integrate u'
               stop
            endif
         else
            write(69,*)'invalid neos=',neos
            stop
         endif
         if(nusegpus.eq.0)then
            write(69,*)'CPUs will be used for any gravity'
         else
            write(69,*)'GPUs will be used for gravity'
         endif
         if(nselfgravity.eq.1)then
            write(69,*)'The gas is self-gravitating'
         else
            write(69,*)'The gas is *not* self-gravitating'
         endif
         write(69,*)'Courant numbers=',cn1,cn2,cn3,cn4,cn5,cn6,cn7
         if(reat.ge.0.d0) write(69,*)'Eating radius reat=',reat
         write(69,*)'init:            ... done'
         if(ncooling.eq.0)then
            write(69,*) 'No radiative cooling'
         else
            write(69,*) 'Radiative cooling will be implemented'
         endif
      endif

      return

c     error condition:
 100  stop 'init:  error reading input file ???'
      end
************************************************************************
c      subroutine rank_setup
c      include 'sphu.h'
c      integer i,j
c      logical strideexists,notfound,notinstride
c      character*15 allhosts(nhmax),hosts(nhmax)
c      real*8 allstride(nhmax)
c
c      gravhost=-1
c      inquire(file='stride.in',exist=strideexists)
c      if(strideexists)then
c         open(23,file='stride.in',status='old')
c         do i=1,nhmax
c            read(23,*,end=33) allhosts(i),allstride(i)
c         enddo
c 33      continue
c         close(23)
c         if(myrank.eq.0) write(69,*) 'read stride.in file',nprocs
c         scalesum=0.d0
c         open(24,file='hosts',status='old')
c         do i=1,nprocs
c            read(24,*) hosts(i)
c            notfound=.true.
c            notinstride=.true.
c            j=1
c            do while(notfound.and.j.le.nhmax)
c               if(hosts(i).eq.allhosts(j))then
c                  stride(i)=allstride(j)
c                  scalesum=scalesum+stride(i)
c                  if(stride(i).gt.0) lasthydroproc=i-1
c                  notfound=.false.
c                  notinstride=.false.
c               else
c                  j=j+1
c               endif
c            enddo
c            if(notinstride)then
c               write(69,*)'need to put ',hosts(i),' in stride.in'
c               stop
c            endif
c            if(hosts(i).eq.'toe') gravhost=i-1
c         enddo
c         close(24)
c      else
c         scalesum=0.d0
c         open(24,file='hosts',status='old')
c         do i=1,nprocs
c            read(24,*) hosts(i)
c            if(hosts(i).eq.'toe') then
c               stride(i)=0
c               gravhost=i-1
c            else
c               stride(i)=1
c               scalesum=scalesum+stride(i)
c               lasthydroproc=i-1
c            endif
c         enddo
c         close(24)
c      endif
c      
c      return
c      end
************************************************************************
c      subroutine hydrocomm_setup
cc     make new mpi_comm for machines doing hydro calculations
cc     (exclude toe if it's stride value is 0)
c      include 'sphu.h'
c      integer color
c      if(gravhost.gt.0 .and. stride(gravhost+1).eq.0)then
c         color=1
c      endif
c      
c      end
************************************************************************
      subroutine get_input
      include 'sphu.h'
      logical autotf
      common/autotfblock/autotf
      common/orbitalelements/e0,semimajoraxis

      semimajoraxis=0.d0
      bimpact=-1.d30
      e0=-1.d30
      vinf2=1.d30

C     Set some default values, so that they don't necessarily have to be set in the sph.input file:
      TF=50000                 ! desired final time to stop simulation
      DTOUT=100                ! how often an out*.sph files should be dumped
      N=100000                 ! desired number of particles.  If N<0 then |N|=number of particles *per solar mass*.  Used only if making a new star.
      Gflag=1                   ! set to 0 for G function from appendix of Gaburov et al. (2010); set to 1 for a G function that works better when there are black holes
      NNOPT=22+Gflag                 ! controls neighbor number.  Leave it at 22 to get almost 40 neighbors.
      NAV=3                    ! artificial viscosity (AV) flag.  Leave it at 3 to get a hybrid Balsara-Monaghan AV.
      ALPHA=1                  ! AV coefficient for term linear in mu
      BETA=2                   ! AV coefficient for mu^2 term
      NGR=3                    ! gravity flag.  Leave it at 3.  If your want no gravity, NGR=0 might still work.
      HMIN=-.5                 ! minimum smoothing length allowed.  Leave it at 0 or a negative number to be safe.  No longer used.
      HMAX=1.d30               ! maximum smoothing length allowed.  Leave it very large to be safe.  No longer used.
      NRELAX=1                 ! Relaxation flag.  0=dynamical calculation, 1=relaxation of single star, 2=relaxation of binary in corotating frame with centrifugal force, 3=calculation rotating frame with centrifugal and Coriolis forces
      TRELAX=1.d30             ! timescale for artificial drag force.  Keep it very large to turn off the drag force, which seems best even in relaxation runs (as the AV can do the relaxation).
      SEP0=200                 ! initial separation of two stars in a binary or collision calculation
      EQUALMASS=0              ! particle mass is proportional to rho^(1-EQUALMASS), so EQUALMASS=1 has equal mass particles and EQUALMASS=0 is for constant number density.
      TRELOFF=0                ! time to end a relaxation and switch to a dynamical calculation
      TRESPLINTMUOFF=0.        ! time to stop resplinting the mean molecular weight.  Leave this at 0.
      NITPOT=1                 ! Number of iterations between evaluation of the gravitational potential energy.
      TSCANON=0                ! time that the scan of a binary starts
      SEPFINAL=1.d30           ! final separation for the scan of a binary
      NINTVAR=2                ! 1=integrate entropic variable A, 2=integrate internal energy u
      ngravprocs=-2             ! The number of gravity processors (must be <= min(nprocs,ngravprocsmax))
      qthreads=0                ! Number of GPU threads per particle. Typically set to 1, 2, 4, or 8.  Set to a negative value to optimize the number of threads by timing.  Set to 0 to guess the best number of threads without timing.
      mbh=4.0                   !Mass of black hole
      runit=6.9599d10          ! number of cm in the unit of length.  Use 6.9599d10 if want solar radius.
      munit=1.9891d33          ! number of g in unit of mass.  Use 1.9891d33 if want solar mass.
!     The Courant numbers cn1, cn2, cn3, and cn4 are for SPH particles:
!     dt_sph=1/(1/dt1 + 1/dt2 + 1/dt3 + 1/dt4)
      cn1=.6d0                 ! dt1=cn1*h/v_signal
      cn2=0.06d0               ! dt2=cn2*(h/|a-a_smoothed|)^0.5
      cn3=0.06d0               ! dt3=cn3*u/|du/dt|
      cn4=1.d30                ! dt4=cn4*v_signal/|a-a_smoothed|
!     The Courant numbers cn5, cn6, and cn7 are for a particle i that is a compact object (co):
!     dt_co=1/(1/dt5 + 1/dt6)
      cn5=0.02d0               ! dt5=cn5*r_ij/v_ij (minimized over all other particles j)
      cn6=0.02d0               ! dt6=cn6*(r_ij/a_ij)^.5 (minimized over all other particles j)
      cn7=4.d0                 ! r_ij=(x_ij^2+y_ij^2+z_ij^2+cn7*h_i^2)^.5
!     The final timestep dt is the minimum of dt_sph and dt_co for all particles i
      COMPUTEEXCLUSIVEMODE=0   ! Set this to 1 if on machine like grapefree with GPUs in compute exclusive mode; set this to 0 on supercomputers like lincoln
      omega_spin=0.d0 ! Angular rotation rate of star, used in nrelax=1 relaxations to give a rigidly rotating model
      ppn=8
      neos=1 ! 0 for polytropic equation of state (EOS), 1 for ideal gas + radiation pressure, 2 for tabulated EOS
      nselfgravity=0 ! 0 if just do gravity to point particles, 1 if self-gravitating
      gam=5.d0/3.d0 ! leave this set at a reasonable value even if using NEOS=1 or 2 (because the value of gam is used in estimating the local sound speed in balAV3.f)
      reat=-1.d0
      starmass=1d0
      starradius=1d0
      ncooling=0 ! 0 if no cooling, otherwise radiative cooling
      hclamp=1.0d30
      rhobreak=0.5d0
      coregam=0.0d0
      sigmoid1=0.0d0
      sigmoid2=0.0d0
      envGam  =1.333333333d0
      nplanets=0
 
      open(12,file='sph.input',err=100)
      read(12,input)
      close(12)

      call set_nusegpus         ! If using GPUs, this sets nusegpus=1 *and* nselfgravity=1

      if(cn1.lt.0.d0 .or. cn2.lt.0.d0 .or. cn3.lt.0.d0 .or.
     $     cn4.lt.0.d0 .or. cn5.lt.0.d0 .or. cn6.lt.0.d0 .or.
     $     cn7.lt.0.d0) then
         write(69,*) 'Strange cn:',cn1,cn2,cn3,cn4,cn5,cn6,cn7
         stop
      endif

      if(ngr.ne.0 .and. nusegpus.eq.0) then
         ngravprocs=nprocs
         if(myrank.eq.0) then
            write(69,*)'Using CPUs to calculate gravity w/ ngravprocs=',
     $           ngravprocs
         endif
      endif

c      if(myrank.eq.0 .and. ngr.gt.0 .and. ngravprocs.eq.0) then
c         write(69,*) 'Need to have at least one gravity process'
c         write(69,*) 'Make sure ngravprocs is set in sph.input'
c         stop
c      endif

      if(ngr.eq.0 .and. ngravprocs.ne.0) then
         ngravprocs=0
         if(myrank.eq.0)
     $        write(69,*)'Because ngr=0, we are setting ngravprocs=0'
      endif

      autotf=.false.
      if(tf.lt.0.d0) then
         autotf=.true.
         tf=abs(tf)
         open(23,file='ecc.sph')
         open(34,file='skipahead.sph')
      endif

      return

 100  stop 'init: error reading input file'

      end
************************************************************************
      subroutine stride_setup
c     assign n_lower and n_upper to each processor
      include 'sphu.h'
c      integer i,stridesum,intstride(nprocs)
      integer i

c      write(69,*)'nprocs=',nprocs
c      write(69,*)'lasthydroproc=',lasthydroproc

c      if(nprocs.eq.1)then
c         n_lower=1
c         n_upper=n
c      else
c         stridesum=0
c         do i=1,myrank+1
c            n_lower=stridesum+1
c            if(i.eq.lasthydroproc+1.or.nprocs.eq.1)then
c               n_upper=n
c               intstride(i)=n_upper-n_lower+1
c            else
c               intstride(i)=int(stride(i)/scalesum*n+0.5d0)
c               n_upper=min(n_lower+intstride(i)-1,n)
c            endif
c            stridesum=stridesum+intstride(i)
c         enddo
c      endif

      if(nrelax.ne.1)then
         n=ntot
      else
         n=n
      endif

c      if(myrank.eq.0) write(69,*)'n_lower,n_upper,n',n_lower,n_upper,n

      do i=1,n
         vxdotsm(i)=0.d0
         vydotsm(i)=0.d0
         vzdotsm(i)=0.d0
      enddo

      return
      end
************************************************************************
      subroutine energy_setup
      include 'sphu.h'
      logical energyfilealreadyexists!,strideexists
      integer filenum
      character*11 energyfile
      character*8 logfile

      if(myrank.eq.0) then 
         energyfilealreadyexists=.true.
         filenum=-1
         do while(energyfilealreadyexists)
            filenum=filenum+1
            if(myrank.eq.0) write(energyfile,101)filenum
 101        format('energy',i1.1,'.sph')
            inquire(file=energyfile,exist=energyfilealreadyexists)
         enddo
      
         write(logfile,103)filenum
 103     format('log',i1.1,'.sph')
         open(22,file=energyfile,status='unknown')
         open(69,file=logfile,status='unknown')
         write(69,*)'writing energy data to ',energyfile
         write(69,*)'writing log data to ',logfile
      endif
      end
************************************************************************
      subroutine lfstart
c     prepare for first leap-frog iteration
      include 'sphu.h'                                          
      integer i!,status(mpi_status_size)
      real*8 rhocgs,rhocgs0,const1,const2,rho1,rho2,gam1,gam2
      real*8 rhoei0,rhoci0,rhocm0,rhocc0,ucgs
      real*8 rho1max,rho2max,ep
      parameter(rho1 = 8.300)
      parameter(rho2 = 4.260)
      parameter(gam1 = 1./0.528)
      parameter(gam2 = 1./0.549)
      parameter(rho1max=20.0)
      parameter(rho2max=6.5)
      integer j

      const1 = (1.0/(3.49*10.**(-(1./gam1+6.))))**gam1
      const2 = (1.0/(1.27*10.**(-(1./gam2+6.))))**gam2

c     omega_spin>0 gives counterclockwise spin:
      if(omega_spin.ne.0.d0 .and. nrelax.eq.1) then
         if(myrank.eq.0) write(69,*) 'omega_spin=',omega_spin
         do i=1,ntot
            vx(i)=vx(i)-omega_spin*y(i)
            vy(i)=vy(i)+omega_spin*x(i)
         enddo
      endif

      call gravquant   ! must call for inital set up
      call rho_and_h

      if(ngr.ne.0) call gravforce

      if (NRELAX.eq.1) then
         rhoei0 = rhoei*munit/runit**3.d0
         rhoci0 = rhoci*munit/runit**3.d0
         rhocm0 = rhocm*munit/runit**3.d0
         rhocc0 = rhocc*munit/runit**3.d0
         do i=n_lower,n_upper
            rhocgs  = rho(i)*munit/runit**3.d0
            rhocgs0 = rho0(i)*munit/runit**3.d0
            ! for mantle-gas particlees
            if (rhocgs0.le.(rhoci0+rhocm0)/2.) then
               cc(i) = 1
               call getCompositions(rhocgs,const2,envConst,gam2,envGam,
     $              rho2,ufactor(i),meanmolecular(i),T0(i),rhocgs0,
     $              rhoei0,rhoci0,0,xmix(i),u(i))
               call getDensityComponents(xmix(i),rhocgs,u(i),const2,
     $              envConst,gam2,envGam,rho2,ufactor(i),
     $              meanmolecular(i),T0(i),rho2max,rhocgs0,rhoei0,
     $              rhoci0,0,rhoh(i),rhol(i))
            ! for core-mantle particles
            else
               cc(i) = 2
               call getCompositions(rhocgs,const1,const2,gam1,gam2,
     $              rho1,rho2,meanmolecular(i),T0(i),rhocgs0,rhocm0,
     $              rhocc0,1,xmix(i),u(i))
               call getDensityComponents(xmix(i),rhocgs,u(i),const1,
     $              const2,gam1,gam2,rho1,rho2,meanmolecular(i),T0(i),
     $              rho1max,rhocgs0,rhocm0,rhocc0,1,rhoh(i),rhol(i))
            endif
            u(i)  = u(i)/gravconst/munit*runit
            umix(i) = u(i)
            dmix(i) = rho(i)
         enddo
      endif

      call uvdots                ! evaluate udot and accelerations at half-timestep

      call enout(.true.)

c     advance velocities to half-timestep:
      call tstep

      dth=0.5d0*dt
c      do i=n_lower,n
c      do i=n_lower,n_upper
      do i=1,n
         vx(i)=vx(i)+vxdot(i)*dth
         vy(i)=vy(i)+vydot(i)*dth
         vz(i)=vz(i)+vzdot(i)*dth
      enddo
c      do i=n_lower,n_upper
      do i=1,n
         u(i)=u(i)+dth*udot(i)  ! update u to half-timestep
      enddo
      t=t+dth

      return
      end

************************************************************************
      subroutine readin
c     6/22/94 - all initial parameters except for n come from sph.input
      include 'sphu.h' 
      integer i,nchk,nrelaxold,nnoptold,noutold,nitold,navold,
     $     ngrold,corepts
      real*8 hminold,hmaxold,sep0old,tfold,dtoutold,told,
     $     alphaold,betaold,tskipahead,trelaxold

c     read in data from previous run:      
      corepts=0
      if(myrank.eq.0) write(69,*)'READIN: reading start file ...' 
      open(12,file='startu.sph',form='unformatted')
c     (The following READ sequence must match exactly the WRITE sequence
c     used in subroutine DUMP)
      read(12) ntot,nnoptold,hminold,hmaxold,sep0old,
     $     tfold,dtoutold,noutold,nitold,told,navold,
     $     alphaold,betaold,tskipahead,ngrold,nrelaxold,
     $     trelaxold,dt
      do i=1,ntot
         read (12) x(i),y(i),z(i),am(i),hp(i),rho(i),vx(i),vy(i),
     $        vz(i),vxdot(i),vydot(i),vzdot(i),u(i),udot(i),
     $        gx(i),gy(i),gz(i),grpot(i),meanmolecular(i),
     $        cc(i)
         if(u(i).eq.0.d0) then
            corepts=corepts+1
            if(myrank.eq.0) write(69,*)
     $         'Should corepts really be set properly?'
            sToP
         endif
      enddo
      read(12) nchk
      close(12)
      n=ntot-corepts
      
      if (nchk.ne.ntot) stop 'READIN: PROBLEM WITH START FILE ???'
      call stride_setup
      if(myrank.eq.0) write(69,*)'READIN:          ... DONE'
      return                                                           
      end
