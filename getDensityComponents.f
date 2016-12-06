
      subroutine getDensityComponents2(x,d,u,K1,K2,g1,g2,b1,b2,m,T,dmax,
     $  dl,dh,u0,d0,core,df1,df2)
      ! Subroutine to solve for densities of mixed-composition particles
      ! Inputs:
      ! x - mass-fraction of heavier composition
      ! d - density of particle
      ! u - internal energy of particle
      ! K1 - polytropic constant of heavier material
      ! K2 - polytropic constant of lighter material
      ! g1 - polytropic index of heavier material
      ! g2 - polytropic index of lighter material
      ! b1 - degrees of freedom of gas-molecule
      ! m - mean molecular weight
      ! T - temperature
      ! dmax - maximum value of density to consider
      ! dl - density of lighter material
      ! dh - density of heavier material
      ! u0 - initial internal energy of particle
      ! d0 - initial density of particle
      ! core - if 1: core-mantle interface, if 0: mantle-envelope interface
      ! Written by Jason Hwang 11/01/16

      implicit none
      real*8, intent(in) :: x,d,u,K1,K2,g1,g2,b1,b2,m,T,dmax,dl,dh,u0,d0
      integer*4, intent(in) :: core
      real*8, intent(out) :: df1,df2

      integer*4 nsection,i
      real*8 d_array(6)
      real*8 ep,e6,uh0,ul0,duh,dul,du,dum,u1,u2
      real*8 uh,ul,ph,pl,dudrho,golden_section,yold,ynew
      external uh,ul,ph,pl,dudrho,golden_section
      parameter(ep=0.00001)

      write(193,*) 'diagnosis: ', x,d,u

      du = (u-u0)/(d-d0)
      uh0 = 0.
      ul0 = 0.
      if (x.gt.0.) uh0 = x*uh(dh,K1,g1,b1,m,T,core)
      if (x.lt.1.) ul0 = (1.-x)*ul(dl,K2,g2,b2,m,T,core)

      duh = (u-u0)*uh0/u0
      dul = (u-u0)*ul0/u0

      ! calculate rho_min of heavier element
      d_array(1) = max(d,b1)

      ! find left and right boundaries
      ! value of dh when dl = dl_initial
      dum = x*dl*d/(dl+(x-1.)*d)
      if (dum.gt.dh) then
         d_array(2) = dh-ep
         d_array(3) = dh+ep
         d_array(4) = dum-ep
         d_array(5) = dum+ep
      else
         d_array(2) = dum-ep
         d_array(3) = dum+ep
         d_array(4) = dh-ep
         d_array(5) = dh+ep
      endif

      ! calculate u(rhom_max)
      e6 = (1.-x)*d+ep
      if (core==1) e6 = max(b2+ep,e6)
      d_array(6) = min(dmax,x*e6*d/(e6+(x-1.)*d))

      ! find number of sections
      if (d_array(2)+ep.gt.d_array(6)) then
         nsection = 2
         d_array(2) = d_array(6)
      elseif (d_array(4)+ep.gt.d_array(6)) then
         nsection = 4
         d_array(4) = d_array(6)
      else
         nsection = 6
      endif

      ! solve for density of heavier component
      df1 = (d_array(1)+d_array(6))/2.
      yold = dudrho(x,df1,d,dh,dl,duh,dul)-du

      do i=1,nsection-1,2
         if (d_array(i).gt.d_array(i+1)) then 
            if (i.lt.5) d_array(i+2) = d_array(i)
         else
            u1 = dudrho(x,d_array(i),d,dh,dl,duh,dul)
            u2 = dudrho(x,d_array(i+1),d,dh,dl,duh,dul)
            dum = golden_section(x,d,du,d_array(i),d_array(i+1),
     $            u1,u2,dh,dl,duh,dul)
            ynew = dudrho(x,dum,d,dh,dl,duh,dul)-du
            if (abs(ynew).lt.abs(yold)) then
               df1 = dum
               yold = ynew
            endif
         endif
      enddo
      df2 = (1.-x)*df1*d/(df1-x*d)

      end



      subroutine getDensityComponents3(x,d,u,K1,K2,g1,g2,b1,b2,m,T,dmax,
     $  dl,dh,u0,d0,lb,ub,core,df1,df2)
      ! Subroutine to solve for densities of mixed-composition particles
      ! Inputs:
      ! x - mass-fraction of heavier composition
      ! d - density of particle
      ! u - internal energy of particle
      ! K1 - polytropic constant of heavier material
      ! K2 - polytropic constant of lighter material
      ! g1 - polytropic index of heavier material
      ! g2 - polytropic index of lighter material
      ! b1 - degrees of freedom of gas-molecule
      ! m - mean molecular weight
      ! T - temperature
      ! dmax - maximum value of density to consider
      ! dl - density of lighter material
      ! dh - density of heavier material
      ! u0 - initial internal energy of particle
      ! d0 - initial density of particle
      ! lb - lower-bound of density at interface
      ! ub - upper-bound of density at interface
      ! core - if 1: core-mantle interface, if 0: mantle-envelope interface
      ! Written by Jason Hwang 11/01/16

      implicit none
      real*8, intent(in) :: x,d,u,K1,K2,g1,g2,b1,b2,m,T,dmax,dl,dh,u0
      real*8, intent(in), target :: lb,ub,d0
      integer*4, intent(in) :: core
      real*8, intent(out) :: df1,df2

      integer*4 ncount,maxcount
      real*8 ep,uerror
      real*8, pointer :: dc,uc,dn,un,dt
      real*8, target :: u1,u2,u3,u4,u5
      real*8, target :: d1,d2,d3,d4,d5
      real*8 e1,e2,e3,e4,e5
      real*8 uh0,ul0,du
      real*8 uh,ul,ph,pl
      external uh,ul,ph,pl
      parameter(ep=0.001)
      parameter(uerror = 0.000000001)
      parameter(maxcount = 20)

      du = u-u0
      uh0 = 0.
      ul0 = 0.
      if (x.gt.0.) uh0 = x*uh(dh,K1,g1,b1,m,T,core)
      if (x.lt.1.) ul0 = (1.-x)*ul(dl,K2,g2,b2,m,T,core)

      ! calculate rho_min of heavier element
      d1 = max(d,b1)
      e1 = (1.-x)*d1*d/(d1-x*d)
      u1 = x*uh(d1,K1,g1,b1,m,T,core)+
     $     (1.-x)*ul(e1,K2,g2,b2,m,T,core)-uh0-ul0

      ! calculate u(rhom_max)
      e2 = (1.-x)*d+ep
      d2 = min(dmax,x*e2*d/(e2+(x-1.)*d))
      e2 = (1.-x)*d2*d/(d2-x*d)
      u2 = x*uh(d2,K1,g1,b1,m,T,core)+
     $     (1.-x)*ul(e2,K2,g2,b2,m,T,core)-uh0-ul0

      ! conduct golden section search until solution is bracketed by higher-density branch
      d3 = (d1+d2)/2.
      e3 = (1.-x)*d3*d/(d3-x*d)
      u3 = x*uh(d3,K1,g1,b1,m,T,core)+
     $     (1.-x)*ul(e3,K2,g2,b2,m,T,core)-uh0-ul0

      ncount = 0
      if ((u1-du)*(u2-du)>=0.) then
         ! set preference based on initial particle location
         dc = dh
         if (d3.gt.dh) then
            uc => u1
            dt => d1
            un => u2
            dn => d2
         else
            uc => u2
            dt => d2
            un => u1
            dn => d1
         endif
         write(195,*) 'target density: ', dc,u1,u2,du
         do while(ncount < maxcount.and.(u1-du)*(u3-du) >= 0.
     $            .and.(u2-du)*(u3-du)>=0.)
            d4 = (d1+d3)/2.
            e4 = (1.-x)*d4*d/(d4-x*d)
            u4 = x*uh(d4,K1,g1,b1,m,T,core)+
     $           (1.-x)*ul(e4,K2,g2,b2,m,T,core)-uh0-ul0

            d5 = (d2+d3)/2.
            e5 = (1.-x)*d5*d/(d5-x*d)
            u5 = x*uh(d5,K1,g1,b1,m,T,core)+
     $           (1.-x)*ul(e5,K2,g2,b2,m,T,core)-uh0-ul0

            if (u4.lt.u3) then
               d2 = d3
               u2 = u3
               d3 = d4
               u3 = u4
            elseif (u5.lt.u3) then
               d1 = d3
               u1 = u3
               d3 = d5
               u3 = u5
            else
               d1 = d4
               u1 = u4
               d2 = d5
               u2 = u5
            endif

            if (d3.gt.dc) then
               uc => u1
               dt => d1
               un => u2
               dn => d2
            else
               uc => u2
               dt => d2
               un => u1
               dn => d1
            endif
            ncount = ncount + 1
         end do
         ! assign values for bisection search                
         if ((uc-du)*(u3-du) <= 0.) then
            d2 = dt
            u2 = uc
            write(195,*) ncount, 'preferable region',u3,u2,(u3+u2)/2.,du
         elseif ((un-du)*(u3-du) <= 0.) then
            d2 = dn
            u2 = un
            write(195,*) ncount, 'unpreferable region',u3,u2,(u3+u2)/2.,du
         else
            d2 = d3
            u2 = u3
            write(195,*) ncount, 'no region',u3,u2,(u3+u2)/2.,du
         endif
         d1 = d3
         u1 = u3
      else
         write(195,*) 'skipping search'
      endif
      ncount = 0

      ! perform bisection search
      if ((u1-du)*(u2-du)<=0.) then
         d3 = (d1+d2)/2.
         e3 = (1.-x)*d3*d/(d3-x*d)
         u3 = x*uh(d3,K1,g1,b1,m,T,core)+
     $        (1.-x)*ul(e3,K2,g2,b2,m,T,core)-uh0-ul0
         do while(ncount < maxcount)! .and. abs(u3-u)/u >= uerror)
            if (u3.gt.du) then
               d2 = d3
               u2 = u3
            else
               d1 = d3
               u1 = u3
            endif
            d3 = (d1+d2)/2.
            e3 = (1.-x)*d3*d/(d3-x*d)
            u3 = x*uh(d3,K1,g1,b1,m,T,core)+
     $           (1.-x)*ul(e3,K2,g2,b2,m,T,core)-uh0-ul0
            ncount = ncount + 1
         end do
      else
         write(195,*) 'cannot perform bisection'
      endif
      write(195,*) 'after bisection: ', u1,u2,u3,du

      df1 = d3
      if (core==0) then
        if (d0.gt.(ub+lb)/2.) write(195,*)'rhom: ', ncount,x,d,df1,d0
        if (d0.le.(ub+lb)/2.) write(195,*)'rhom: ', ncount,x,d,df1,ub
      else
        if (d0.gt.(ub+lb)/2.) write(195,*)'rhoc: ', ncount,x,d,df1,d0
        if (d0.le.(ub+lb)/2.) write(195,*)'rhoc: ', ncount,x,d,df1,ub
      endif

      df2 = (1.-x)*df1*d/(df1-x*d)
      if (core==0) then
        if (d0.gt.(ub+lb)/2.) write(195,*)'rhoe: ', ncount,x,d,df2,lb
        if (d0.le.(ub+lb)/2.) write(195,*)'rhoe: ', ncount,x,d,df2,d0
      else
        if (d0.gt.(ub+lb)/2.) write(195,*)'rhom: ', ncount,x,d,df2,lb
        if (d0.le.(ub+lb)/2.) write(195,*)'rhom: ', ncount,x,d,df2,d0
      endif

      if (u == u0) then
         df1 = dh
         df2 = dl
      endif

      end



      subroutine getDensityComponents4(x,d,u,K1,K2,g1,g2,b1,b2,m,T,dmax,
     $  dl,dh,u0,d0,lb,ub,core,df1,df2)
      ! Subroutine to solve for densities of mixed-composition particles
      ! Inputs:
      ! x - mass-fraction of heavier composition
      ! d - density of particle
      ! u - internal energy of particle
      ! K1 - polytropic constant of heavier material
      ! K2 - polytropic constant of lighter material
      ! g1 - polytropic index of heavier material
      ! g2 - polytropic index of lighter material
      ! b1 - degrees of freedom of gas-molecule
      ! m - mean molecular weight
      ! T - temperature
      ! dmax - maximum value of density to consider
      ! dl - density of lighter material
      ! dh - density of heavier material
      ! u0 - initial internal energy of particle
      ! d0 - initial density of particle
      ! lb - lower-bound of density at interface
      ! ub - upper-bound of density at interface
      ! core - if 1: core-mantle interface, if 0: mantle-envelope interface
      ! Written by Jason Hwang 11/01/16

      implicit none
      real*8, intent(in) :: x,d,u,K1,K2,g1,g2,b1,b2,m,T,dmax,dl,dh,u0,d0
      real*8, intent(in), target :: lb,ub
      integer*4, intent(in) :: core
      real*8, intent(out) :: df1,df2

      integer*4 ncount,maxcount
      real*8 ep,uerror
      real*8, pointer :: uc,dn,un,dt
      real*8, target :: u1,u2,u3,u4,u5
      real*8, target :: d1,d2,d3,d4,d5
      real*8 e1,e2,e3,e4,e5
      real*8 uh0,ul0,du
      real*8 uh,ul,ph,pl
      external uh,ul,ph,pl
      parameter(ep=0.001)
      parameter(uerror = 0.000000001)
      parameter(maxcount = 20)

      du = (u-u0)/(d-d0)
      uh0 = 0.
      ul0 = 0.
      if (x.gt.0.) uh0 = x*(uh(dh,K1,g1,b1,m,T,core)+(u-u0))
      if (x.lt.1.) ul0 = (1.-x)*(ul(dl,K2,g2,b2,m,T,core)+(u-u0))

!      write(189,*)core,u/u0,d,d0,du

      ! calculate rho_min of heavier element
      d1 = max(d,b1)
      e1 = (1.-x)*d1*d/(d1-x*d)
      u1 = (x*ph(uh0,d1,K1,g1,b1,m,core)+
     $     (1.-x)*pl(ul0,e1,K2,g2,b2,m,core))/d/d

      ! calculate u(rhom_max)
      e2 = (1.-x)*d+ep
      d2 = min(dmax,x*e2*d/(e2+(x-1.)*d))
      e2 = (1.-x)*d2*d/(d2-x*d)
      u2 = (x*ph(uh0,d2,K1,g1,b1,m,core)+
     $     (1.-x)*pl(ul0,e2,K2,g2,b2,m,core))/d/d

      ! conduct golden section search until solution is bracketed by higher-density branch
      d3 = (d1+d2)/2.
      e3 = (1.-x)*d3*d/(d3-x*d)
      u3 = (x*ph(uh0,d3,K1,g1,b1,m,core)+
     $     (1.-x)*pl(ul0,e3,K2,g2,b2,m,core))/d/d

      ncount = 0
      if ((u1-du)*(u2-du)>=0.) then
         ! set preference based on initial particle location
         if (d3.gt.dh) then
            uc => u1
            dt => d1
            un => u2
            dn => d2
         else
            uc => u2
            dt => d2
            un => u1
            dn => d1
         endif
         write(195,*) 'target density: ', dh,u1,u2,du
         do while(ncount < maxcount.and.(u1-du)*(u3-du) >= 0.
     $            .and.(u2-du)*(u3-du)>=0.)
            d4 = (d1+d3)/2.
            e4 = (1.-x)*d4*d/(d4-x*d)
            u4 = (x*ph(uh0,d4,K1,g1,b1,m,core)+
     $           (1.-x)*pl(ul0,e4,K2,g2,b2,m,core))/d/d

            d5 = (d2+d3)/2.
            e5 = (1.-x)*d5*d/(d5-x*d)
            u5 = (x*ph(uh0,d5,K1,g1,b1,m,core)+
     $           (1.-x)*pl(ul0,e5,K2,g2,b2,m,core))/d/d

            if (u4.lt.u3) then
               d2 = d3
               u2 = u3
               d3 = d4
               u3 = u4
            elseif (u5.lt.u3) then
               d1 = d3
               u1 = u3
               d3 = d5
               u3 = u5
            else
               d1 = d4
               u1 = u4
               d2 = d5
               u2 = u5
            endif

            if (d3.gt.dh) then
               uc => u1
               dt => d1
               un => u2
               dn => d2
            else
               uc => u2
               dt => d2
               un => u1
               dn => d1
            endif
            ncount = ncount + 1
         end do
         ! assign values for bisection search                
         if ((uc-du)*(u3-du) <= 0.) then
            d2 = dt
            u2 = uc
            write(195,*) ncount, 'preferable region'
         elseif ((un-du)*(u3-du) <= 0.) then
            d2 = dn
            u2 = un
            write(195,*) ncount, 'unpreferable region'
         else
            d2 = d3
            u2 = u3
            write(195,*) ncount, 'no region'
         endif
         d1 = d3
         u1 = u3
      else
         write(195,*) 'skipping search',(u1-du)*(u2-du)
      endif
      ncount = 0

      ! perform bisection search
      if ((u1-du)*(u2-du)<=0.) then
         write(195,*) 'beginning bisection: ', d1, d2
         d3 = (d1+d2)/2.
         e3 = (1.-x)*d3*d/(d3-x*d)
         u3 = (x*ph(uh0,d3,K1,g1,b1,m,core)+
     $        (1.-x)*pl(ul0,e3,K2,g2,b2,m,core))/d/d
         do while(ncount < maxcount)! .and. abs(u3-u)/u >= uerror)
            if (u3.gt.du) then
               d2 = d3
               u2 = u3
            else
               d1 = d3
               u1 = u3
            endif
            d3 = (d1+d2)/2.
            e3 = (1.-x)*d3*d/(d3-x*d)
            u3 = (x*ph(uh0,d3,K1,g1,b1,m,core)+
     $           (1.-x)*pl(ul0,e3,K2,g2,b2,m,core))/d/d
            ncount = ncount + 1
            write(195,*) (u3-du)/du
         end do
      else
         write(195,*) 'cannot perform bisection'
      endif

      df1 = d3
      df2 = (1.-x)*df1*d/(df1-x*d)

      write(195,*) ncount,x,d,df1,dh,df2,dl

      if (d == d0) then
         df1 = dh
         df2 = dl
      endif
      if (u == u0) then
         df1 = dh
         df2 = dl
      endif

      end



      subroutine getDensityComponents(x,d,u,K1,K2,g1,g2,b1,b2,m,T,dmax,
     $  d0,lb,ub,core,df1,df2)
      ! Subroutine to solve for densities of mixed-composition particles
      ! Inputs:
      ! x - mass-fraction of heavier composition
      ! d - density of particle
      ! u - internal energy of particle
      ! K1 - polytropic constant of heavier material
      ! K2 - polytropic constant of lighter material
      ! g1 - polytropic index of heavier material
      ! g2 - polytropic index of lighter material
      ! b1 - degrees of freedom of gas-molecule
      ! m - mean molecular weight
      ! T - temperature
      ! dmax - maximum value of density to consider
      ! d0 - initial density of particle
      ! lb - lower-bound of density at interface
      ! ub - upper-bound of density at interface
      ! core - if 1: core-mantle interface, if 0: mantle-envelope interface
      ! Written by Jason Hwang 11/01/16

      implicit none
      real*8, intent(in) :: x,d,u,K1,K2,g1,g2,b1,b2,m,T,dmax
      real*8, intent(in), target :: lb,ub,d0
      integer*4, intent(in) :: core
      real*8, intent(out) :: df1,df2

      integer*4 ncount,maxcount
      real*8 ep,uerror
      real*8 uh,ul
      real*8 dc,uc,dn,un,dt
      real*8 u1,u2,u3,u4,u5
      real*8 d1,d2,d3,d4,d5
      real*8 e1,e2,e3,e4,e5
      external uh,ul
      parameter(ep=0.001)
      parameter(uerror = 0.000000001)
      parameter(maxcount = 20)

      ! calculate rho_min of heavier element
      d1 = max(d,b1)
      e1 = (1.-x)*d1*d/(d1-x*d)
      u1 = x*uh(d1,K1,g1,b1,m,T,core)+(1.-x)*ul(e1,K2,g2,b2,m,T,core)

      ! calculate u(rhom_max)
      e2 = (1.-x)*d+ep
      if (core.eq.1) e2 = max(e2,b2+ep)
      d2 = min(dmax,x*e2*d/(e2+(x-1.)*d))
      e2 = (1.-x)*d2*d/(d2-x*d)
      u2 = x*uh(d2,K1,g1,b1,m,T,core)+(1.-x)*ul(e2,K2,g2,b2,m,T,core)

      ! conduct golden section search until solution is bracketed by higher-density branch
      d3 = (d1+d2)/2.
      e3 = (1.-x)*d3*d/(d3-x*d)
      u3 = x*uh(d3,K1,g1,b1,m,T,core)+(1.-x)*ul(e3,K2,g2,b2,m,T,core)

      ncount = 0
      if ((u1-u)*(u2-u)>=0.) then
         ! set preference based on initial particle location
         dc = d0
         if (d0.lt.(ub+lb)/2.) dc = ub

         do while(ncount < maxcount.and.(u1-u)*(u3-u) >= 0.
     $            .and.(u2-u)*(u3-u)>=0.)
            d4 = (d1+d3)/2.
            e4 = (1.-x)*d4*d/(d4-x*d)
            u4 = x*uh(d4,K1,g1,b1,m,T,core)+
     $           (1.-x)*ul(e4,K2,g2,b2,m,T,core)

            d5 = (d2+d3)/2.
            e5 = (1.-x)*d5*d/(d5-x*d)
            u5 = x*uh(d5,K1,g1,b1,m,T,core)+
     $           (1.-x)*ul(e5,K2,g2,b2,m,T,core)

            if (u4.lt.u3) then
               d2 = d3
               u2 = u3
               d3 = d4
               u3 = u4
            elseif (u5.lt.u3) then
               d1 = d3
               u1 = u3
               d3 = d5
               u3 = u5
            else
               d1 = d4
               u1 = u4
               d2 = d5
               u2 = u5
            endif
            ncount = ncount + 1
         end do

         if (d3.gt.dc) then
            uc = u1
            dt = d1
            un = u2
            dn = d2
         else
            uc = u2
            dt = d2
            un = u1
            dn = d1
         endif

         ! assign values for bisection search                
         if ((uc-u)*(u3-u) <= 0.) then
            d2 = dt
            u2 = uc
         elseif ((un-u)*(u3-u) <= 0.) then
            d2 = dn
            u2 = un
         else
            d2 = d3
            u2 = u3
         endif
         d1 = d3
         u1 = u3
      else
      ncount = 0

      ! perform bisection search
      if ((u1-u)*(u2-u)<=0.) then
         d3 = (d1+d2)/2.
         e3 = (1.-x)*d3*d/(d3-x*d)
         u3 = x*uh(d3,K1,g1,b1,m,T,core)+
     $        (1.-x)*ul(e3,K2,g2,b2,m,T,core)
         do while(ncount < maxcount)
            if (u3.gt.u) then
               d2 = d3
               u2 = u3
            else
               d1 = d3
               u1 = u3
            endif
            d3 = (d1+d2)/2.
            e3 = (1.-x)*d3*d/(d3-x*d)
            u3 = x*uh(d3,K1,g1,b1,m,T,core)+
     $       (1.-x)*ul(e3,K2,g2,b2,m,T,core)

            ncount = ncount + 1
         end do
      endif

      if (x.ne.0.) then
         df1 = d3
         df2 = (1.-x)*df1*d/(df1-x*d)
      else
         df1 = 0.
         df2 = 0.
         if (d0.gt.(ub+lb)/2.) df1 = d0
         if (d0.le.(ub+lb)/2.) df2 = d0
      endif

      end














