      SUBROUTINE XERBLA(SRNAME,INFO)
      INTEGER INFO
      CHARACTER*6 SRNAME
c      WRITE (*,9999) SRNAME,INFO
c      STOP
       return
 9999 FORMAT (' ** On entry to ',A6,' parameter number ',I2,' had ',
     +       'an illegal value')
      END
