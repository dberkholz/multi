      PROGRAM MULTICORE
C
C Reads multiple  PROTEIN DATA BANK files with orthogonal A coordinates 
C The best rotation is determined and carried out based on the
C Kabsch algorithm for every pairwise comparison using a user-given
C distance cutoff for atoms to be considered structurally equivalent.
C
      CHARACTER*80 FILE1,NAME
      CHARACTER*4 AT1,RES1,ATYP,HEAD
      CHARACTER*1 HFLAG
      DIMENSION XIN(3,5001,30),RMSUSE(31,31),RMSALL(31,31),NUSE(31,31)
      DIMENSION X(3,5001),Y(3,5001),WT(5000),ISITIN(5000)
      DIMENSION ROT(3,3),DIST(5000)
      DIMENSION RES1(5000),AT1(5000),B(5000),P(5000)
C
      OPEN(UNIT=80,FILE='multicore.log',STATUS='UNKNOWN')
C
c
      WRITE(*,*) 'Program Multicore version 1.0 - 26 July 2005)'
      WRITE(*,*) '**********************************'
      WRITE(*,*) 'All coordinate sets must have the exact same' 
      WRITE(*,*) 'number of atoms. Each model must start with a '
      WRITE(*,*) 'MODEL line and end with an ENDMDL line'
      WRITE(*,*) ' '
      WRITE(*,*) 'ALL PAIRS OF COORDINATE FILES WILL BE COMPARED'
      WRITE(*,*) ' '
      WRITE(*,*) 'TER, ANISOU and all other lines are ignored'
c
      WRITE(*,*) ' Name of the coordinate file?'
      READ(5,1000) FILE1
      OPEN(1,FILE=FILE1,STATUS='OLD',READONLY)
      I = 0
 5    READ(1,1000,END=30) NAME
      WRITE(*,*) name
      IF (NAME(1:5).EQ.'MODEL') GOTO 7
      GOTO 5
 7    I = I +1
      IPT = 1
 10   READ (1,1010,END=20,ERR=20)  HEAD,HFLAG,AT1(IPT),RES1(IPT),
     $   XIN(1,IPT,I),XIN(2,IPT,I),XIN(3,IPT,I)
      IF (HEAD(1:4).EQ.'ENDM') GOTO 20
      IF (HFLAG.EQ.'H'.OR.AT1(IPT)(1:1).EQ.'H') GOTO 10  !SKIP HYDROGENS
      IF (HEAD(1:4).NE.'ATOM'.AND.HEAD(1:4).NE.'HETA') GOTO 10
      IPT = IPT + 1
      GOTO 10
 20   NATOMS = IPT - 1
      WRITE(6,*) NATOMS,' non-H atoms read for set',I
      IF (I.EQ.1) IPT1 = NATOMS
      WRITE(80,1040) NATOMS,I
c
      IF (NATOMS.NE.IPT1) THEN
        WRITE(*,*) 'wrong number of atoms in file',I
        PAUSE 'To continue anyway enter "C"'
        IPT1 = MIN(IPT1,NATOMS)
        WRITE(*,*) 'Fitting based on the first',IPT1,' coordinates.'
        WRITE(80,1060) IPT1
      ENDIF
      GOTO 5
c
 30   CONTINUE
      NFILE = I
      DO 35 I = 1 , IPT1
        ISITIN(I) = 1
 35   CONTINUE
C
      WRITE(*,*) ' Name of the output flag file?'
      READ(5,1000) NAME
      OPEN(2,FILE=NAME,STATUS='NEW')
C
      WRITE(*,*) 'Max deviation (in A) for points used in superposition'
      READ(*,*) DMAX
      WRITE(80,*) ' Cutoff for acceptable atoms is',DMAX
      NCYC = 50
C
C Cycle over all pairwise combinations of the coordinates and
C accumulate in ISITIN a flag of whether a given atom was
C always within the allowed error or not. In this way generate 
C a consitent set of atoms for use in the final overlays
C
      DO 200 IFILE = 1 , NFILE
      DO 190 JFILE = IFILE+1 , NFILE
        DO 40 K = 1 , IPT1      ! fill the coordinate arrays with the 
        DO 40 L = 1 , 3         ! current pair
          X(L,K) = XIN(L,K,IFILE)
          Y(L,K) = XIN(L,K,JFILE)
          WT(K) = 1.0
 40     CONTINUE
        
        CALL ITERATE(NCYC,DMAX,IPT1,X,Y,WT,ROT,RMS)
C
C Write out final coordinates and information.
C
      IUSE = 0
      SUMD2 = 0.0
      DO 50 I = 1 , IPT1
        D2 = (X(1,I) - Y(1,I))**2 + 
     $       (X(2,I) - Y(2,I))**2 +
     $       (X(3,I) - Y(3,I))**2 
        SUMD2 = SUMD2 + D2
        DIST(I) = SQRT(D2)
        IF (WT(I).EQ.1.0) IUSE = IUSE + 1
        ISITIN(I) = ISITIN(I)*WT(I)
 50   CONTINUE
      NUSE(JFILE,IFILE) = IUSE
      RMSUSE(JFILE,IFILE) = RMS
      RMSALL(JFILE,IFILE) = SQRT(SUMD2/IPT1)
      WRITE(*,1150)  JFILE,IFILE,IUSE,RMS,RMSALL(JFILE,IFILE)
      WRITE(*,1110)  ROT

 190  CONTINUE
 200  CONTINUE
C
C Output overlay statistics
C
      IUSE = 0
      DO 225 I = 1 , IPT1
        WRITE(2,1015) RES1(I),AT1(I),ISITIN(I)
        IF (ISITIN(I).EQ.1) IUSE = IUSE + 1
 225  CONTINUE
C
      II = 0
      DO 250 IFILE = 1 , NFILE
      DO 250 JFILE = IFILE+1 , NFILE
        II = II + 1
        AVEUSE = AVEUSE + NUSE(JFILE,IFILE)
        AVERMS = AVERMS + RMSUSE(JFILE,IFILE)
        AVERMSALL = AVERMSALL + RMSALL(JFILE,IFILE)
 250  CONTINUE
      AVEUSE = AVEUSE/II
      AVERMS = AVERMS/II
      AVERMSALL = AVERMSALL/II
C
      WRITE(80,*) 'With cutoff=',DMAX,':',IUSE,' out of',IPT1,
     $ ' are atoms common to all superpositions.'
      WRITE(80,*) II,' pairs were compared'
      WRITE(80,*) 'Average common=',AVEUSE,';Average rms=',AVERMS,
     $            ';average all atom rms=',AVERMSALL
      WRITE(81,*) 'With cutoff=',DMAX,': atoms common to all=',IUSE
      WRITE(81,*) II,' pairs were compared'
      WRITE(81,*) 'Average common=',AVEUSE,';Average rms=',AVERMS,
     $            ';average all atom rms=',AVERMSALL
C
      WRITE(80,*) ' summary of # atoms, rms subset, rms all atoms'
      DO 300 j = 1, nfile
        WRITE(80,1155) (NUSE(i,j),i=j,nfile)
        WRITE(80,1160) (RMSUSE(i,j),i=j,nfile)
        WRITE(80,1160) (RMSALL(i,j),i=j,nfile)
 300  continue
 1155 FORMAT(<5*(j-1)>x,<nfile-j+1>I5)
 1160 FORMAT(<5*(j-1)>x,<nfile-j+1>F5.2)
      STOP
 1000 FORMAT(A)
 1010 FORMAT(A4,8X,A1,A4,5X,A4,4X,3F8.3)
 1015 FORMAT('ATOM',1X,A4,1X,A4,I5)
 1020 FORMAT(10X,' atom-1    atom-2     dist')
 1030 FORMAT(10X,2(2A4,2X),F7.2)
 1040 FORMAT(/,I6,' atomic coordinates read from set',I5) 
 1050 FORMAT(/,I6,' atomic coordinates read from ',A80) 
 1060 FORMAT(' !!!!!Fit based only on first',I5,' common coordinates',/)
 1070 FORMAT(/,I6,' atomic coordinates written to ',A80)
 1080 FORMAT(/' Best rotation of',I5,' equivalent atoms based on',
     $ I5,' atoms with',/,' final deviation less than',F6.2,' A.')
 1090 FORMAT(/,' The center of mass of set 1 before superposition is:',
     $ /,/,20X,3F10.4)
 1100 FORMAT(/,' The center of mass of set 2 is:',/,/,20X,3F10.4)
 1110 FORMAT(/,' The rotation matrix is:',/,3(/,20X,3F10.6))
 1120 FORMAT(/,' The final rms error for the',I5,' points is',F7.3)
 1130 FORMAT(/,' The final rms error for all',I5,' points is',F7.3)
 1140 FORMAT(/,' Here follows a listing of the deviations for the',
     $            ' individual atoms.',/)
 1150 FORMAT(/'**Set',I3,' onto',I3,
     $   ':nats=',i4,';rmsuse=',f6.2,';rmsall=',f6.2)
      END
      SUBROUTINE ITERATE(NCYC,DMAX,NPT,X,Y,WT,ROT,RMS)
      DIMENSION X(3,NPT+1),Y(3,NPT+1),WT(NPT),ROT(3,3)
      DIMENSION XINIT(3,5000)
C
C Program for steering superposition of two sets of points.
C NCYC is the number of iterations desired and DMAX is
C the maximum distance between equivalent points for which
C the point will be included in the next superposition attempt.
C The WT buffer must be filled ahead of time signaling which
C atoms should be used during the first superposition.
C X coordinates are rotated to match the Y coordinate set and
C the rotation matrix is returned. The center of mass for the
C two unrotated coordinate sets is also returned in X(NPT+1),
C Y(NPT+1).
C
C
C Put the original X-coordinates in a temporary array
C and find out how many are to be used in first pass.
C
      NUSE = 0
      DO 10 I = 1,NPT
        IF (WT(I).EQ.1.0) NUSE = NUSE + 1
        DO 5 J = 1,3
 5      XINIT(J,I) = X(J,I)
 10   CONTINUE
C
C Repeat process for desired number of cycles or convergence
C
      DO 100 ICYC = 1,NCYC
C
C      Refill X array so that rotation matrix always is the desired one
C
        DO 50 I = 1,NPT
          DO 40 J = 1,3
 40       X(J,I) = XINIT(J,I)
 50     CONTINUE
C
C      Do the job
C
        CALL SUPERIMPOSE(NPT,X,Y,WT,ROT,RMS)
C        WRITE(6,*)'*****CYCLE***************',ICYC
C        WRITE(6,*)' FOR',NUSE,' ATOMS USED THE RMS DEVIATION IS',RMS
        NSAME = 0
        NUSE = 0
        DO 30 I=1,NPT
          DIST2 = 0.
          DO 20 J=1,3
 20       DIST2 = DIST2 + (X(J,I) - Y(J,I))**2
          DIST = SQRT(DIST2)
          WTOLD = WT(I)
          WT(I) = 1.
          IF (DIST.GT.DMAX) WT(I) = 0.
          IF (WT(I).EQ.1.) NUSE = NUSE + 1
          IF (WT(I).EQ.WTOLD) NSAME = NSAME + 1
 30     CONTINUE
        TYPE *,NUSE,' ATOMS TO BE USED FOR THE NEXT TRY.'
        IF (NSAME.EQ.NPT) THEN
C          WRITE(6,*) 'PROCESS CONVERGED IN',ICYC,'CYCLES.'
          RETURN
        ENDIF
C
 100  CONTINUE
C
C All done
C
      WRITE(6,*) 'COMPLETED',NCYC,' ITERATIONS WITHOUT CONVERGENCE'
      RETURN
      END
C
C
      SUBROUTINE SUPERIMPOSE(NPT,X,Y,WT,ROT,RMS)
C***********************************************************************
C*****VERSION OF 7.FEB.1986 BY P.A. KARPLUS
C***********************************************************************
C*****Given two sets of coordinates and it will calculate the best
C*****rotation for 1 onto 2 and return the rotated coordinates
C*****as well as the rms error, and the transformation matrix.
C*****The program is based on the method of Kabsch: Acta Cryst. A32 
C*****(1976),p922 & A34(1978),p826.
C*****If RMS error is returned as -1. that means the points were on
C*****a line and this method does not work.
C***********************************************************************
      DIMENSION X(3,NPT+1),Y(3,NPT+1),WT(NPT),CENX(3),CENY(3)
      DIMENSION XREL(3,5000),YREL(3,5000)
      REAL  ROT(3,3),RTR(3,3),MU(3),ROTA(3,3)
      REAL  MU1,MU2,MU3
      DIMENSION R(3,3),A(3,3),B(3,3),B3OLD(3)
      COMMON A1(3),A2(3),A3(3),B1(3),B2(3),B3(3),R1(3),R2(3),R3(3)
      EQUIVALENCE (A1(1),A(1,1)),(B1(1),B(1,1)),(R1(1),ROTA(1,1))
C
C Zero arrays
C
      RMS = -1.
      EZERO = 0.
      SUMWT = 0.
      DO 5 I=1,3
        CENX(I) = 0.
        CENY(I) = 0.
      DO 5 J=1,3
        R(I,J) = 0.
        B(I,J) = 0.
        RTR(I,J) = 0.
 5      ROTA(I,J) = 0.
C
C translate coordinates such that center of mass is zero 
C
      DO 10 N = 1,NPT
        SUMWT = SUMWT + WT(N)
      DO 10 I = 1,3
        CENX(I) = CENX(I) + X(I,N)*WT(N)
 10     CENY(I) = CENY(I) + Y(I,N)*WT(N)
C
      DO 15 I = 1,3
        CENX(I) = CENX(I)/SUMWT
        X(I,NPT+1) = CENX(I)
        CENY(I) = CENY(I)/SUMWT
 15     Y(I,NPT+1) = CENY(I)
      DO 20 N = 1,NPT
      DO 20 I = 1,3
        XREL(I,N) = X(I,N) - CENX(I)
        YREL(I,N) = Y(I,N) - CENY(I)
 20     EZERO = EZERO + WT(N)*(XREL(I,N)**2 + YREL(I,N)**2)
      EZERO = EZERO/2.
C
C Fill the matrix R
C
      DO 30 N = 1,NPT
      DO 30 J = 1,3
      DO 30 I = 1,3
 30     R(I,J) = R(I,J) + WT(N)*YREL(I,N)*XREL(J,N)
C
C Generate R-tilda R
C
      DO 40 J = 1,3
      DO 40 I = 1,J
        DO 35 K = 1,3
 35       RTR(I,J) = RTR(I,J) + R(K,I)*R(K,J)
 40       RTR(J,I) = RTR(I,J)
C
C Get eigen-vectors and eigen-values of RTR.
C Reverse the order of the vectors and the values so 
C that they are in descending order.
C
      CALL EIGEN(RTR,A,MU)
      TMP=MU(1)
      MU(1) = MU(3)
      MU(3) = TMP
      IF (MU(2).EQ.0.) RETURN  ! degenerate
      A(1,1)=A(1,3)
      A(2,1)=A(2,3)
      A(3,1)=A(3,3)
C
C Set A3 = A1 cross A2 
C
      CALL CROSS(A1,A2,A3)
      CALL NORM(A3,A3LEN)
C
C Determine B(I) = R * A(I) 
C
       DO 41 I=1,3
       DO 41 J=1,3
       DO 41 K=1,3
 41      B(J,I) = B(J,I) + R(J,K)*A(K,I)
C
C Set B3 = B1 cross B2. Put original B3 in B3OLD for determining sigma3
C
      CALL NORM(B1,B1LEN)
      CALL NORM(B2,B2LEN)
      CALL NORM(B3,B3LEN)
      DO 42 I=1,3
 42   B3OLD(I) = B3(I)
      CALL CROSS(B1,B2,B3)
      CALL NORM(B3,B3LEN)
C
C Make rotation matrix (ROT)
C
      DO 50 J=1,3
      DO 50 I=1,3
      DO 50 K=1,3
 50     ROTA(J,I) = ROTA(J,I) + B(I,K)*A(J,K)
      CALL NORM(R1,R1LEN)
      CALL NORM(R2,R2LEN)
      CALL NORM(R3,R3LEN)
      DO 55 I=1,3
      DO 55 J=1,3
 55     ROT(I,J) = ROTA(I,J) 
C
C Apply rotation to X and translate X to YCEN
C
      DO 60 N = 1,NPT
      DO 60 I = 1,3
        X(I,N) = CENY(I)
      DO 60 J = 1,3
 60     X(I,N) = X(I,N) + ROT(J,I)*XREL(J,N)
C
C Calculate RMS deviation
C
      DOT = B3(1)*B3OLD(1) + B3(2)*B3OLD(2) + B3(3)*B3OLD(3)
      DO 65 I=1,3
 65   IF (MU(I).LT.0.) MU(I)=0.
      MU1 = SQRT(MU(1))
      MU2 = SQRT(MU(2))
      MU3 = SQRT(MU(3))
      IF (DOT.LT.0) MU3 = -MU3
      ERROR = EZERO - (MU1 + MU2 + MU3)
      AVGSQERROR = 2. * ERROR / SUMWT
      IF (AVGSQERROR.LT.0) AVGSQERROR = 0.
      RMS = SQRT(AVGSQERROR)
C
C All done. Write some stuff to unit 8 in case it is wanted.
C
      WRITE (8,*) 'NUMBER OF POINTS',NPT
      WRITE (8,*) 'CENTER OF MASS SET 1'
      WRITE (8,1000) CENX
      WRITE (8,*) 'CENTER OF MASS SET 2'
      WRITE (8,1000) CENY
      WRITE (8,*) 'RMS DEVIATION',RMS
      WRITE (8,*) 'ROTATION MATRIX'
      WRITE (8,1000) ROT
      RETURN
 1000 FORMAT(3F10.5)
      END

      SUBROUTINE NORM(A,ALEN)
C*****To normalize a vector a(3)
      REAL A(3),ALEN
      ALEN = SQRT(A(1)**2 + A(2)**2 + A(3)**2)
      DO 10 I=1,3
 10   A(I)=A(I)/ALEN
      RETURN
      END
C
C
      SUBROUTINE CROSS(A,B,C)
C*****To calculate C = A cross B
      REAL A(3),B(3),C(3)
      C(1) = A(2)*B(3) - A(3)*B(2)
      C(2) = A(3)*B(1) - A(1)*B(3)
      C(3) = A(1)*B(2) - A(2)*B(1)
      RETURN
      END
c***********************************************************************
      subroutine eigen ( a, b, c)
      dimension a(3,3), b(3,3), c(3), d(3)
c
c      three-by-three eigensystem routine called by
c      programs solving for principal planes.
c
c      a is the incoming symetric matrix.
c      b is the three eigen vectors returned for matrix a.
c      c is the eigen values of matrix a.
c
      logical conv
c
      DO 10 i = 1, 3
      DO 10 j = 1, 3
 10     b(i,j) = a(i,j)
C
      call tred2 (3, 3, b, c, d)
      call tql2 (3, 3, c, d, b, conv)
C
      IF (.not. conv) STOP 'Eigenvalue did not converge'
C
      return
      END
      subroutine tred2 (n, m, a, d, e)
      dimension  a(m,n), d(n), e(n)
c
c     Real symmetric matrix a(1:n,1:n) is reduced to tridiagonal
c     form by Householder transformations.
c
c     On return a contains the orthogonal matrix q such that qt=a,
c     where t is tridiagonal.  The diagonal of t is returned in
c     d(1:n) and the subdiagonal in e(2:n).
c
c     Tol is machine dependent and should be eta/epsilon, where
c     eta is the smallest representable real number and epsilon
c     is the smallest positive real number such that
c     (1.0 + epsilon .ne. 1.0) is .true.
c
c     Algorithm taken from Numerische Mathamatik 11, 181-195(1968)
c
      double precision ff, gg
c
c     Tol for the VAX 11/780
c
      data tol / '00000c00'x /
c
      IF (n .gt. 2) go to 10
      d(1) = a(1,1)
      e(1) = 0.0
      a(1,1) = 1.0
      IF (n .lt. 2) return
c
      d(2) = a(2,2)
      e(2) = a(2,1)
      a(2,1) = 0.0
      a(1,2) = 0.0
      a(2,2) = 1.0
      return
c
   10 continue
      DO 120 ii = 2, n
      l = n - ii
      i = l + 2
      f = a(i,i-1)
      gg = 0.0d0
      IF (l .lt. 1) go to 30
c
      DO 20 k = 1, l
   20 gg = gg + a(i,k)**2
c
   30 g = gg
      h = gg + f*f
c
c     IF g is too small, skip this transformation
c
      IF (g .le. tol) go to 100
c
      l = l + 1
      g = sqrt(h)
      IF (f .ge. 0.0) g = -g
      e(i) = g
      h = h - f*g
      a(i,i-1) = f - g
      ff = 0.0d0
      DO 70 j = 1, l
      a(j,i) = a(i,j)/h
      gg = 0.0d0
      DO 40 k = 1, j
   40 gg = gg + a(j,k)*a(i,k)
      IF (j .ge. l) go to 60
      jj = j + 1
      DO 50 k = jj, l
   50 gg = gg + a(k,j)*a(i,k)
   60 continue
      e(j) = gg/h
      ff = ff + gg*a(j,i)
   70 continue
      hh = ff/(h + h)
      DO 90 j = 1, l
      f = a(i,j)
      g = e(j) - hh*f
      e(j) = g
      DO 80 k = 1, j
   80 a(j,k) = a(j,k) - f*e(k) - g*a(i,k)
   90 continue
      go to 110
c
  100 continue
      e(i) = f
      h = 0.0
c
  110 continue
      d(i) = h
  120 continue
c
      d(1) = a(1,1)
      e(1) = 0.0
      a(1,1) = 1.0
c
c     accumulate the transformation matrix
c
      DO 180 i = 2, n
      l = i - 1
      IF (d(i) .eq. 0.0) go to 160
c
      DO 150 j = 1, l
      gg = 0.0d0
      DO 130 k = 1, l
  130 gg = gg + a(i,k)*a(k,j)
      g = gg
      DO 140 k = 1, l
  140 a(k,j) = a(k,j) - g*a(k,i)
  150 continue
c
  160 continue
      d(i) = a(i,i)
      a(i,i) = 1.0
      DO 170 j = 1, l
      a(i,j) = 0.0
  170 a(j,i) = 0.0
  180 continue
c
      return
      END
      subroutine tql2 (nn, mm, d, e, z, conv)
      dimension z(mm,nn), d(nn), e(nn)
      logical conv
c
c     Real symmetric tridiagonal ql eigenvalue and vector program
c
c     Num. Math. 11, 293-306(1968)
c
c     The principal virtue of this routine is absolutely
c     orthogonal eigenvectors even for multiple or pathologically
c     close eigenvalues.  Of course, IF the matrix has very large
c     (relatively) elements in the lower right corner they won't
c     be correct; in such a case the matrix should be turned
c     END-for-END.
c
c     d(1:n) contains the diagonal elements of the matrix on input
c     and the eigenvalues sorted in increasing absolute value on
c     output.
c     e(2:n) contains the subdiagonal elements of the matrix on
c     input and is destroyed.
c     z(1:n,1:n) contains the orthogonal matrix which reduces the
c     original matrix to tridiagonal form, e.g. the output matrix
c     from tred2, or the identity matrix if the original matrix
c     is tridiagonal.  on output z contains the eigenvectors by
c     columns.
c     Eps is the smallest positive number such that
c     (1.0 + eps .ne. 1.0) is .true.
c     conv is set to .true. UNLESS the algorithm fails to converge
c     in thirty iterations.
c
c     Eps for the VAX 11/780
c
      data eps / '00003500'x /
c
      n = nn
      IF (n .gt. 1) go to 10
      conv = .false.
      IF (n .lt. 1) return
      conv = .true.
      z(1,1) = 1.0
      return
c
   10 continue
      conv = .true.
      DO 20 i = 2, n
   20 e(i-1) = e(i)
      e(n) = 0.0
      b = 0.0
      f = 0.0
      DO 180 l = 1, n
      j = 0
      h = eps*(abs(d(l)) + abs(e(l)))
      IF (b .lt. h) b = h
c
c     search for small subdiagonal element
c
      DO 30 m = l, n
      IF (abs(e(m)) .le. b) go to 40
   30 continue
   40 IF (m .eq. l) go to 120
c
   50 continue
      IF (j .lt. 30) go to 60
      conv = .false.
      return
c
   60 continue
      j = j + 1
c
c     calculate origin shift for this iteration
c
      p = (d(l+1) - d(l))/(2.0*e(l))
      r = sqrt(p*p + 1.0)
      h = r + abs(p)
      IF (p .lt. 0.0) h = -h
      h = d(l) - e(l)/h
      DO 70 i = l, n
   70 d(i) = d(i) - h
      f = f + h
c
c     ql transformation
c
      p = d(m)
      c = 1.0
      s = 0.0
      DO 110 ii = l, m
      IF (ii .eq. m) go to 110
      i = m + l - ii - 1
      g = c*e(i)
      h = c*p
      IF (abs(p) .ge. abs(e(i))) go to 80
c
      c = p/e(i)
      r = sqrt(c*c + 1.0)
      e(i+1) = s*e(i)*r
      s = 1.0/r
      c = c/r
      go to 90
c
   80 continue
      c = e(i)/p
      r = sqrt(c*c + 1.0)
      e(i+1) = s*p*r
      s = c/r
      c = 1.0/r
c
   90 continue
      p = c*d(i) - s*g
      d(i+1) = h + s*(c*g + s*d(i))
c
c     form vector
c
      DO 100 k = 1, n
      h = z(k,i+1)
      z(k,i+1) = s*z(k,i) + c*h
      z(k,i) = c*z(k,i) - s*h
  100 continue
c
  110 continue
      e(l) = s*p
      d(l) = c*p
      IF (abs(e(l)) .gt. b) go to 50
c
c     found an eigenvalue
c
  120 continue
      d(l) = d(l) + f
c
c     order eigenvalues and eigenvectors
c
      IF (l .eq. 1) go to 180
      p = d(l)
      r = abs(p)
      DO 130 j = 1, l
      IF (abs(d(j)) .gt. r) go to 140
  130 continue
      go to 180
c
  140 continue
      k = l
  150 continue
      d(k) = d(k-1)
      k = k - 1
      IF (k .gt. j) go to 150
      d(j) = p
c
      DO 170 i = 1, n
      p = z(i,l)
      k = l
  160 z(i,k) = z(i,k-1)
      k = k - 1
      IF (k .gt. j) go to 160
  170 z(i,j) = p
c
  180 continue
      return
c
      END
