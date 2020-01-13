      subroutine cosqb (n,x,wsave)
      dimension  x(*),wsave(*)
      data tsqrt2 /2.82842712474619/
      if (n-2) 101,102,103
  101 x(1) = 4.*x(1)
      return
  102 x1 = 4.*(x(1)+x(2))
      x(2) = tsqrt2*(x(1)-x(2))
      x(1) = x1
      return
  103 call cosqb1 (n,x,wsave,wsave(n+1))
      return
      end
