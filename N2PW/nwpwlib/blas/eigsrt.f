
      SUBROUTINE EIGSRT(D,V,N,NP)
*     -----------------------------------------------------------------
*     SORT EIGENVALUES INTO DESCENDING ORDER AND REARRANGE EIGENVECTORS
*     -----------------------------------------------------------------
      implicit double precision(a-h, o-z)
      double precision D(NP),V(NP,NP)
      DO 13 I=1,N-1
         K=I
         P=D(I)
         DO 11 J=I+1,N
            IF(D(J).GE.P) THEN
               K=J
               P=D(J)
            ENDIF
   11    CONTINUE
         IF(K.NE.I) THEN
            D(K)=D(I)
            D(I)=P
            DO 12 J=1,N
               P=V(J,I)
               V(J,I)=V(J,K)
               V(J,K)=P
   12       CONTINUE
         ENDIF
   13 CONTINUE
      RETURN
      END
