      subroutine rffti (n,wsave)
      real*8          wsave(*)
      if (n .eq. 1) return
      call rffti1 (n,wsave(n+1),wsave(2*n+1))
      return
      end
