      program example
      implicit none
      integer IX,IQ
      integer i
      real*8 xx(1:47),QS(1:30)
      real*8 DUV,DDV,DUBAR,DDBAR,DSTR,DGLU,x,Q2,Q
      real*8 DUP,DDN

      data xx /
     1           1.0D-5, 4.D-5, 6.7D-5, 1.0D-4, 1.4D-4, 2.0D-4,
     2           3.0D-4, 4.5D-4, 6.7D-4, 1.0D-3, 1.4D-3, 2.0D-3,
     3           3.0D-3, 4.5D-3, 6.7D-3, 1.0D-2, 1.4D-2, 2.0D-2,
     4           3.0D-2, 4.5D-2, 0.06, 0.08, 0.1, 0.125, 0.15,
     5           0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325,
     6           0.35, 0.375, 0.4,  0.45, 0.5, 0.55, 0.6,
     7           0.65,  0.7,  0.75,  0.8,  0.85, 0.9, 0.95, 1.0 / 
 
      data QS  / 0.8D0, 1.0D0, 1.25d0, 1.5D0, 2.d0, 2.5D0, 
     1     4.0D0, 6.4D0, 1.0D1, 1.5D1, 2.5D1, 4.0D1, 6.4D1,
     2     1.0D2, 1.8D2, 3.2D2, 5.8D2, 1.0D3, 1.8D3,
     3     3.2D3, 5.8D3, 1.0D4, 1.8D4, 3.2D4, 5.8D4, 
     4     1.0D5, 1.8D5, 3.2D5, 5.8D5, 1.0D6  /
    
       call DSSVINI

       do i = 1,100
          x = 0.01*i
          Q = dsqrt(2.4d0)
!          do Q=1d0
             Q2=Q*Q
             call DSSVFIT(X,Q2,DUV,DDV,DUBAR,DDBAR,DSTR,DGLU)

             DUP = DUV + DUBAR
             DDN = DDV + DDBAR


             write(21,80) dlog(x),dlog(Q2),DUP/x
             write(22,80) dlog(x),dlog(Q2),DUBAR/x
             write(23,80) dlog(x),dlog(Q2),DDN/x
             write(24,80) dlog(x),dlog(Q2),DDBAR/x
             write(25,80) dlog(x),dlog(Q2),DSTR/x
             write(31,80) dlog(x),dlog(Q2),DGLU/x

!          enddo
       enddo

 80    format(3(1pE12.4))



       stop
       end program example 

