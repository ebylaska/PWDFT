      subroutine sinqf (n,x,wsave)
      dimension       x(*)       ,wsave(*)
      if (n .eq. 1) return
      ns2 = n/2
      do 101 k=1,ns2
         kc = n-k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
  101 continue
      call cosqf (n,x,wsave)
      do 102 k=2,n,2
         x(k) = -x(k)
  102 continue
      return
      end
