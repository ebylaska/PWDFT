      subroutine cosqf (n,x,wsave)
      dimension   x(*),wsave(*)
      data sqrt2 /1.4142135623731/
      if (n-2) 102,101,103
  101 tsqx = sqrt2*x(2)
      x(2) = x(1)-tsqx
      x(1) = x(1)+tsqx
  102 return
  103 call cosqf1 (n,x,wsave,wsave(n+1))
      return
      end
