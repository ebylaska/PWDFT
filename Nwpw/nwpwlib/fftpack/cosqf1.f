      subroutine cosqf1 (n,x,w,xh)
      dimension  x(*),w(*),xh(*)
      ns2 = (n+1)/2
      np2 = n+2
      do 101 k=2,ns2
         kc = np2-k
         xh(k) = x(k)+x(kc)
         xh(kc) = x(k)-x(kc)
  101 continue
      modn = mod(n,2)
      if (modn .eq. 0) xh(ns2+1) = x(ns2+1)+x(ns2+1)
      do 102 k=2,ns2
         kc = np2-k
         x(k) = w(k-1)*xh(kc)+w(kc-1)*xh(k)
         x(kc) = w(k-1)*xh(k)-w(kc-1)*xh(kc)
  102 continue
      if (modn .eq. 0) x(ns2+1) = w(ns2)*xh(ns2+1)
      call rfftf (n,x,xh)
      do 103 i=3,n,2
         xim1 = x(i-1)-x(i)
         x(i) = x(i-1)+x(i)
         x(i-1) = xim1
  103 continue
      return
      end
