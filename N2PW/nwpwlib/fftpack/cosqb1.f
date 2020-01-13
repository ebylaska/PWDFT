      subroutine cosqb1 (n,x,w,xh)
      dimension  x(*),w(*),xh(*)
      ns2 = (n+1)/2
      np2 = n+2
      do 101 i=3,n,2
         xim1 = x(i-1)+x(i)
         x(i) = x(i)-x(i-1)
         x(i-1) = xim1
  101 continue
      x(1) = x(1)+x(1)
      modn = mod(n,2)
      if (modn .eq. 0) x(n) = x(n)+x(n)
      call rfftb (n,x,xh)
      do 102 k=2,ns2
         kc = np2-k
         xh(k) = w(k-1)*x(kc)+w(kc-1)*x(k)
         xh(kc) = w(k-1)*x(k)-w(kc-1)*x(kc)
  102 continue
      if (modn .eq. 0) x(ns2+1) = w(ns2)*(x(ns2+1)+x(ns2+1))
      do 103 k=2,ns2
         kc = np2-k
         x(k) = xh(k)+xh(kc)
         x(kc) = xh(k)-xh(kc)
  103 continue
      x(1) = x(1)+x(1)
      return
      end
