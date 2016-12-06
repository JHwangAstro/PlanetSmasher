
      real*8 function dudrho(x,d1,d,dh,dl,duh,dul)
         real*8, intent (in) :: x,d1,d,dh,dl,duh,dul
         real*8 d2

         d2 = (1.-x)*d1*d/(d1-x*d)
         dudrho = x*duh/(d1-dh)+(1.-x)*dul/(d2-dl)
         return
      end


      real*8 function uh(d,K,gam,b,m,T,core)
         real*8, intent (in) :: d,K,gam,b,m,T
         integer, intent (in) :: core
         real*8 p

         p = 0.99535151+ 8.40976345*(b/d)-0.56542391*(b/d)**2.
         if(core==1) p=0.99663563+15.98296227*(b/d)-0.71640366*(b/d)**2.
         uh = K*d**(gam-1.)/(gam-1.)*p
         return
      end

      
      real*8 function ul(d,K,gam,b,m,T,core)
         real*8, intent (in) :: d,K,gam,b,m,T
         integer, intent (in) :: core
         real*8 p
         real*8 boltzcgs
         parameter(boltzcgs=1.380658d-16)

         if (core==1) then
            p = 0.99535151+ 8.40976345*(b/d)-0.56542391*(b/d)**2.
            ul = K*d**(gam-1.)/(gam-1.)*p
         else
            ul = K*d**(gam-1.)/(gam-1.)+b*boltzcgs*T/m
         endif

         return
      end


      real*8 function ph(u,d,K,gam,b,m,core)
         real*8, intent (in) :: u,d,K,gam,b,m
         integer, intent (in) :: core
         real*8 p

         p = 0.99535151+ 8.40976345*(b/d)-0.56542391*(b/d)**2.
         if(core==1) p=0.99663563+15.98296227*(b/d)-0.71640366*(b/d)**2.
         ph = u*(gam-1.)*d*(1.-b/d)**gam/p

         return
      end


      real*8 function pl(u,d,K,gam,b,m,core)
         real*8, intent (in) :: u,d,K,gam,b,m
         integer, intent (in) :: core
         real*8 p

         if (core==1) then
            p = 0.99535151+ 8.40976345*(b/d)-0.56542391*(b/d)**2.
            pl = u*(gam-1.)*d*(1.-b/d)**gam/p
         else
            pl = d*u/b+K*d**gam*(1.-1./b/(gam-1.))
         endif

         return
      end






      real*8 function pIron(b,d)
         real*8, intent (in) :: b,d
         pIron = 0.99663563+15.98296227*(b/d)-0.71640366*(b/d)**2.
         return
      end


      real*8 function pMantle(b,d)
         real*8, intent (in) :: b,d
         pMantle = 0.99535151+ 8.40976345*(b/d)-0.56542391*(b/d)**2.
         return
      end


      real*8 function uSeager(d,K,gam,b,m,T,p)
         real*8, intent (in) :: d,K,gam,b,m,T
         real*8 p
         uSeager = K*d**(gam-1.)/(gam-1.)*p(b,d)
         return
      end
      

      real*8 function uMESA(d,K,gam,b,m,T,p)
         real*8, intent (in) :: d,K,gam,b,m,T
         real*8 p
         uMESA = K*d**(gam-1.)/(gam-1.)+b*boltz*T/m
         return
      end
