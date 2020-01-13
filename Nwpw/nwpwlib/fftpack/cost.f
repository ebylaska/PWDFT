      subroutine cost (n,x,wsave)
      dimension   x(*),wsave(*)
      nm1 = n-1
      np1 = n+1
      ns2 = n/2
      if (n-2) 106,101,102
  101 x1h = x(1)+x(2)
      x(2) = x(1)-x(2)
      x(1) = x1h
      return
  102 if (n .gt. 3) go to 103
      x1p3 = x(1)+x(3)
      tx2 = x(2)+x(2)
      x(2) = x(1)-x(3)
      x(1) = x1p3+tx2
      x(3) = x1p3-tx2
      return
  103 c1 = x(1)-x(n)
      x(1) = x(1)+x(n)
      do 104 k=2,ns2
         kc = np1-k
         t1 = x(k)+x(kc)
         t2 = x(k)-x(kc)
         c1 = c1+wsave(kc)*t2
         t2 = wsave(k)*t2
         x(k) = t1-t2
         x(kc) = t1+t2
  104 continue
      modn = mod(n,2)
      if (modn .ne. 0) x(ns2+1) = x(ns2+1)+x(ns2+1)
      call rfftf (nm1,x,wsave(n+1))
      xim2 = x(2)
      x(2) = c1
      do 105 i=4,n,2
         xi = x(i)
         x(i) = x(i-2)-x(i-1)
         x(i-1) = xim2
         xim2 = xi
  105 continue
      if (modn .ne. 0) x(n) = xim2
  106 return
      end
