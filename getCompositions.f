

      subroutine getCompositions(d,K1,K2,g1,g2,b1,b2,m,T,
     $  d0,lb,ub,core,x,u)
      ! Subroutine to solve for densities of mixed-composition particles
      ! Inputs:
      ! d - density of particle
      ! K1 - polytropic constant of heavier material
      ! K2 - polytropic constant of lighter material
      ! g1 - polytropic index of heavier material
      ! g2 - polytropic index of lighter material
      ! b1 - degrees of freedom of gas-molecule
      ! m - mean molecular weight
      ! T - temperature
      ! d0 - initial density of particle
      ! lb - lower-bound of density at interface
      ! ub - upper-bound of density at interface
      ! core - if 1: core-mantle interface, if 0: mantle-envelope interface

      ! x - mass-fraction of heavier composition
      ! u - internal energy of particle
      ! Written by Jason Hwang 11/04/16

      implicit none
      real*8, intent(in) :: d,K1,K2,g1,g2,b1,b2,m,T,d0,lb,ub
      integer*4, intent(in) :: core
      real*8, intent(out) :: x,u

      real*8 uh,ul
      external uh,ul

      x = 0.0d0
      if (core == 1) then
         if (d0.le.(lb+ub)/2.) then
            !calculate composition, forcing 0<x<1
            if (d.gt..95*lb) x = ub*(d-d0)/(ub-d0)/d
            x=max(0.0d0,min(1.0d0,x))
            !calculate internal energy of mixed particle
            u=x*uh(ub,K1,g1,b1,m,T,core)+(1.-x)*ul(d0,K2,g2,b2,m,T,core)
         else
            !calculate composition, forcing 0<x<1
            x=d0*(d-lb)/(d0-lb)/d
            x=max(0.0d0,min(1.0d0,x))
            !calculate internal energy of mixed particle
            u=x*uh(d0,K1,g1,b1,m,T,core)+(1.-x)*ul(lb,K2,g2,b2,m,T,core)
         endif
         if (x.le.0.005) x = 0.0d0
         if (x.gt.0.995) x = 1.0d0
      else
         if (d0.le.(lb+ub)/2.) then
            !calculate composition, forcing 0<x<1
            if (d.gt.lb) x = ub*(d-d0)/(ub-d0)/d
            x=max(0.0d0,min(1.0d0,x))
            !calculate internal energy of mixed particle
            u=x*uh(ub,K1,g1,b1,m,T,core)+(1.-x)*ul(d0,K2,g2,b2,m,T,core)
         else
            !calculate composition, forcing 0<x<1
            x=d0*(d-lb)/(d0-lb)/d
            x=max(0.0d0,min(1.0d0,x))
            !calculate internal energy of mixed particle
            u=x*uh(d0,K1,g1,b1,m,T,core)+(1.-x)*ul(lb,K2,g2,b2,m,T,core)
         endif
         if (x.lt.0.005)  x = 0.0d0
         if (x.gt.0.9995) x = 1.0d0
      endif

      end





      subroutine getCompositions2(d,K1,K2,g1,g2,b1,b2,m,T,
     $  d0,lb,ub,core,x,u)
      ! Subroutine to solve for densities of mixed-composition particles
      ! in a restarted run
      ! Inputs:
      ! d - density of particle
      ! K1 - polytropic constant of heavier material
      ! K2 - polytropic constant of lighter material
      ! g1 - polytropic index of heavier material
      ! g2 - polytropic index of lighter material
      ! b1 - degrees of freedom of gas-molecule
      ! m - mean molecular weight
      ! T - temperature
      ! d0 - initial density of particle
      ! lb - lower-bound of density at interface
      ! ub - upper-bound of density at interface
      ! core - if 1: core-mantle interface, if 0: mantle-envelope interface

      ! x - mass-fraction of heavier composition
      ! u - internal energy of particle
      ! Written by Jason Hwang 11/04/16

      implicit none
      real*8, intent(in) :: d,K1,K2,g1,g2,b1,b2,m,T,d0,lb,ub!,uin
      integer*4, intent(in) :: core
      real*8, intent(inout) :: x,u

      real*8 uh,ul,K0
      external uh,ul
      !K0 is the recalculated entropic constant

      if (core == 1) then
         if (d0.le.(lb+ub)/2.) then
            !change in u from initial assignment
            u=u-x*uh(ub,K1,g1,b1,m,T,core)
            u=u-(1.-x)*ul(d0,K2,g2,b2,m,T,core)
            if (d.gt..95*lb) x = ub*(d-d0)/(ub-d0)/d
            x=max(0.0d0,min(1.0d0,x))
            !adjust for new x
            u=u+x*uh(ub,K1,g1,b1,m,T,core)
            u=u+(1.-x)*ul(d0,K2,g2,b2,m,T,core)
         else
            !change in u from initial assignment
            u=u-x*uh(d0,K1,g1,b1,m,T,core)
            u=u-(1.-x)*ul(lb,K2,g2,b2,m,T,core)
            x=d0*(d-lb)/(d0-lb)/d
            x=max(0.0d0,min(1.0d0,x))
            !adjust for new x
            u=u+x*uh(d0,K1,g1,b1,m,T,core)
            u=u+(1.-x)*ul(lb,K2,g2,b2,m,T,core)
         endif
         if (x.le.0.005) x = 0.0d0
         if (x.gt.0.995) x = 1.0d0
      else
         if (d0.le.(lb+ub)/2.) then
            !change in u from initial assignment
            u=u-x*uh(ub,K1,g1,b1,m,T,core)
            u=u-(1.-x)*ul(d0,K2,g2,b2,m,T,core)
            if (d.gt.lb) x = ub*(d-d0)/(ub-d0)/d
            x=max(0.0d0,min(1.0d0,x))
            !adjust for new x
            u=u+x*uh(ub,K1,g1,b1,m,T,core)
            u=u+(1.-x)*ul(d0,K2,g2,b2,m,T,core)
         else
            !change in u from initial assignment
            u=u-x*uh(d0,K1,g1,b1,m,T,core)
            u=u-(1.-x)*ul(lb,K2,g2,b2,m,T,core)
            x=d0*(d-lb)/(d0-lb)/d
            x=max(0.0d0,min(1.0d0,x))
            !adjust for new x
            u=u+x*uh(d0,K1,g1,b1,m,T,core)
            u=u+(1.-x)*ul(lb,K2,g2,b2,m,T,core)
         endif
         if (x.lt.0.005)  x = 0.0d0
         if (x.gt.0.9995) x = 1.0d0
      endif

      end










