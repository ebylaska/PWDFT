      subroutine sinqb (n,x,wsave)
      dimension       x(*)       ,wsave(*)
      if (n .gt. 1) go to 101
      x(1) = 4.*x(1)
      return
  101 ns2 = n/2
      do 102 k=2,n,2
         x(k) = -x(k)
  102 continue
      call cosqb (n,x,wsave)
      do 103 k=1,ns2
         kc = n-k
         xhold = x(k)
         x(k) = x(kc+1)
         x(kc+1) = xhold
  103 continue
      return
      end
