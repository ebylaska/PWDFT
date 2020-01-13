

      SUBROUTINE TQLI(D,E,N,NP,Z,IERR)
*     ==========================================
*     EIGENVALUES AND EIGENVECTORS
*     FOR A REAL SYMMETRIC TRIDIAGONAL MATRIX 
*     ==========================================
      implicit double precision(a-h, o-z)
      double precision D(NP),E(NP),Z(NP,NP)
      IF(N.GT.1) THEN
        DO 11 I=2,N
          E(I-1)=E(I)
   11   CONTINUE
        E(N)=0
        DO 15 L=1,N
          ITER=0
    1     DO 12 M=L,N-1
            DD=DABS(D(M))+DABS(D(M+1))
            IF(DABS(E(M))+DD.EQ.DD) GO TO 2
   12     CONTINUE
          M=N
    2     IF(M.NE.L) THEN
            IF(ITER.EQ.30) THEN
              IERR=1
              RETURN
            ENDIF
            ITER=ITER+1
            G=(D(L+1)-D(L))/(2.0d0*E(L))
            R=DSQRT(G**2+1.0d0)
            G=D(M)-D(L)+E(L)/(G+DSIGN(R,G))
            S=1.0d0
            C=1.0d0
            P=0.0d0
            DO 14 I=M-1,L,-1
              F=S*E(I)
              B=C*E(I)
              IF(DABS(F).GE.DABS(G)) THEN
                C=G/F
                R=DSQRT(C**2+1.0d0)
                E(I+1)=F*R
                S=1.0d0/R
                C=C*S
              ELSE
                S=F/G
                R=DSQRT(S**2+1.0d0)
                E(I+1)=G*R
                C=1.0d0/R
                S=S*C
              ENDIF
              G=D(I+1)-P
              R=(D(I)-G)*S+2.0d0*C*B
              P=S*R
              D(I+1)=G+P
              G=C*R-B
              DO 13 K=1,N
                F=Z(K,I+1)
                Z(K,I+1)=S*Z(K,I)+C*F
                Z(K,I)=C*Z(K,I)-S*F
   13         CONTINUE
   14       CONTINUE
            D(L)=D(L)-P
            E(L)=G
            E(M)=0.0d0
            GO TO 1
          ENDIF
   15   CONTINUE
      ENDIF
      RETURN
      END
