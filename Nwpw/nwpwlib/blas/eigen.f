

      SUBROUTINE EIGEN(NEMAX,NE,HML,EIG,WORK,IERR)
*     =======================================
*     DIAGONALIZATION OF HAMILTONIAN MATRIX 
*     =======================================
      implicit double precision(a-h, o-z)
      DIMENSION EIG(*),HML(NEMAX,*),WORK(*)

      CALL DCOPY(NEMAX,0.0d0,0,EIG,1)
      IF(NE.EQ.0) RETURN
      IF(NE.EQ.1) THEN
        EIG(1)=HML(1,1)
        HML(1,1)=1.0d0
        RETURN
      ENDIF
      CALL TRED2(HML,NE,NEMAX,EIG,WORK)
      CALL TQLI(EIG,WORK,NE,NEMAX,HML,IERR)

      CALL EIGSRT(EIG,HML,NE,NEMAX)

 1000 FORMAT(I8,E20.10)
      RETURN
      END
