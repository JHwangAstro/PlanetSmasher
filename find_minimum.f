
      ! x1 - min value of x
      ! x2 - max value of x
      real*8 function golden_section(xm,x,y,x10,x20,y10,y20,xh,xl,yh,yl)
        real*8, intent(in) :: xm,x,y,xh,xl,yh,yl
        real*8, intent(in) :: x10,x20,y10,y20

        integer*4 ncount,maxcount
        real*8 x1,y1,x2,y2
        real*8 xt,xn,yt,yn
        real*8 x3,y3,x4,y4,x5,y5
        real*8 dudrho
        external dudrho
        parameter(maxcount = 20)

        ! conduct golden section search until solution is bracketed by higher-density branch
        x1 = x10
        y1 = y10
        x2 = x20
        y2 = y20
        x3 = (x1+x2)/2.
        y3 = dudrho(xm,x3,x,xh,xl,yh,yl)

        ncount = 0
        if ((y1-y)*(y2-y)>=0.) then
           do while(ncount < maxcount.and.(y1-y)*(y3-y) >= 0.
     $              .and.(y2-y)*(y3-y)>=0.)
              x4 = (x1+x3)/2.
              y4 = dudrho(xm,x4,x,xh,xl,yh,yl)

              x5 = (x2+x3)/2.
              y5 = dudrho(xm,x5,x,xh,xl,yh,yl)

              if (y4.lt.y3) then
                 x2 = x3
                 y2 = y3
                 x3 = x4
                 y3 = y4
              elseif (y5.lt.y3) then
                 x1 = x3
                 y1 = y3
                 x3 = x5
                 y3 = y5
              else
                 x1 = x4
                 y1 = y4
                 x2 = x5
                 y2 = y5
              endif

              ncount = ncount + 1
           end do

           if (x3.gt.xh) then
              xt = x1
              yt = y1
              xn = x2
              yn = y2
           else
              xt = x2
              yt = y2
              xn = x1
              yn = y1
           endif

           ! assign values for bisection search                
           if ((yt-y)*(y3-y) <= 0.) then
              x2 = xt
              y2 = yt
           elseif ((yn-y)*(y3-y) <= 0.) then
              x2 = xn
              y2 = yn
           else
              x2 = x3
              y2 = y3
           endif
           x1 = x3
           y1 = y3
        endif

        ! perform bisection search
        ncount = 0
        if ((y1-y)*(y2-y)<=0.) then
           x3 = (x1+x2)/2.
           y3 = dudrho(xm,x3,x,xh,xl,yh,yl)
           do while(ncount < maxcount)
              if (y3.gt.y) then
                 if (y1.gt.y2) then
                    x1 = x3
                    y1 = y3
                 else
                    x2 = x3
                    y2 = y3
                 endif
              else
                 if (y2.gt.y1) then
                    x1 = x3
                    y1 = y3
                 else
                    x2 = x3
                    y2 = y3
                 endif
              endif
              x3 = (x1+x2)/2.
              y3 = dudrho(xm,x3,x,xh,xl,yh,yl)
              ncount = ncount + 1
           end do
        endif

        golden_section = x3
        return
      end
