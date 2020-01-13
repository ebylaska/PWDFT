      subroutine cosqi (n,wsave)
      dimension       wsave(*)
      data pih /1.57079632679491/
      dt = pih/float(n)
      fk = 0.
      do 101 k=1,n
         fk = fk+1.
         wsave(k) = cos(fk*dt)
  101 continue
      call rffti (n,wsave(n+1))
      return
      end
