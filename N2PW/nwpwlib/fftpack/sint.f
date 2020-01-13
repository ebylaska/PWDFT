      subroutine sint (n,x,wsave)
      dimension       x(*)       ,wsave(*)
      np1 = n+1
      iw1 = n/2+1
      iw2 = iw1+np1
      iw3 = iw2+np1
      call sint1(n,x,wsave,wsave(iw1),wsave(iw2),wsave(iw3))
      return
      end
