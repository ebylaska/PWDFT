      subroutine sinti (n,wsave)
      dimension       wsave(*)
      data pi /3.14159265358979/
      if (n .le. 1) return
      ns2 = n/2
      np1 = n+1
      dt = pi/float(np1)
      do 101 k=1,ns2
         wsave(k) = 2.*sin(k*dt)
  101 continue
      call rffti (np1,wsave(ns2+1))
      return
      end
