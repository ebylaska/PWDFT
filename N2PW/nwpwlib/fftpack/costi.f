      subroutine costi (n,wsave)
      dimension       wsave(1)
      data pi /3.14159265358979/
      if (n .le. 3) return
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      dt = pi/float(nm1)
      fk = 0.
      do 101 k=2,ns2
         kc = np1-k
         fk = fk+1.
         wsave(k) = 2.*sin(fk*dt)
         wsave(kc) = 2.*cos(fk*dt)
  101 continue
      call rffti (nm1,wsave(n+1))
      return
      end
