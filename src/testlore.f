      PROGRAM MULTILORE
C
C LODR stands for Locally-overlayed dipeptide residual and is a tool for 
C     analysing the differences in conformation in various structures at the
C     single residue phi,psi angle level.  The LODR score naturally combines 
C     differences in phi and psi into a single number that relfects the level
C     of difference in chain path.
C
C The program reads 1 multi-model coordinate file and writes 2 output files 
C     1) a log file - LODR.log  
C     2) a user specifed output file with the results of the LODR analysis
C        this file has the colums
C
C The user must also specify two integer values (N,M) defining how to compare 
C the coordinate sets.  N defines how many coordinate sets are in the first
C group and M is the number of coordinate sets in the second group.  Forseen 
C uses are as follows:
C    N   M    situation
C    0,  #    one ensemble of # equivalent coordinate sets is being analyzed 
C             to assess their local agreement along the chain.  This would 
C             normally be used for the analysis of an NMR ensemble.
C    1,  1    Two single structures are being compared.
C    1,  #    a single structure is being compared to an enzemble of # 
C             coordinate sets.  This would normally be used for comparing one
C             X-ray structure with an NMR ensemble of the same molecule 
C    #n, #m   Two ensembles are being compared.  This could be either a pair of 
C             NMR determinations of the same molecule, or any two sets of
C             structures that represent different states of a protein 
C             (i.e. oxidized/reduced or bound/free)
C
C The output file has the information
C
C  Col 1: residue number  
C  Col 2: average spread first (N) group of coordinate sets
C  Col 3: rms spread first (N) group of coordinate sets
C  Col 4: average spread second (M) group of coordinate sets
C  Col 5: rms spread second (M) group of coordinate sets
C  Col 6: average spread between first (N) and second (M) sets
C  Col 7: rms spread between first (N) and second (M) sets
C  Col 8: closest approach between first (N) and second (M) sets 
C
      CHARACTER*20 BLANK
      CHARACTER*80 INFILE,LINE,OUTFILE
      CHARACTER*5 AT1,ATREF
      CHARACTER*4 RES1,ATYP,RESNAM,LABEL,RESREF,RES
      REAL RMS1,RMS2,RMS12,AVE1,AVE2,AVE12
      DIMENSION XIN(3,5001,31),XPEP(3,13,31),XRA(3,5001)
      DIMENSION X(3,5001),Y(3,5001),WT(5000),ISITIN(5000)
      DIMENSION ROT(3,3),DIST(5000),RMSPEP(31,31), ATREF(5000)
      DIMENSION RESREF(5000),RES(2000)
      DIMENSION RES1(5000),AT1(5000),B(5000),P(5000),RESNAM(5000)
C
      BLANK='                    '
      WRITE(*,*) 'Program Multicore version 1.0 - Dec 2004)'
      WRITE(*,*) '**********************************'
      WRITE(*,*) 'All coordinate sets must have the exact same' 
      WRITE(*,*) 'number of atoms. Each model must start with a '
      WRITE(*,*) 'MODEL line and end with a TER line'
      WRITE(*,*) ' '
      WRITE(*,*) 'Name of the coordinate file?'
      READ(5,1000) INFILE
      OPEN(1,FILE=INFILE,STATUS='OLD',READONLY) 
      WRITE(*,*) 'How many coordinates files in each of the two sets?'
      WRITE(*,*) 'Options are 0,# OR 1,1 OR 1,# OR #n,#m'
      READ(*,*) NN, MM
      WRITE(*,*) 'Name for the results file?'
      READ(5,1000) OUTFILE
      OPEN(21,FILE=OUTFILE,STATUS='NEW') 
      I = 0
 5    READ(1,1000,END=25) LINE
      WRITE(*,*) LINE
      IF (LINE(1:5).EQ.'MODEL') GOTO 7
      GOTO 5
 7    I=I+1
      IPT = 1
 10   READ (1,1010,END=20,ERR=888) LABEL,IAT,AT1(IPT),RESNAM(IPT),
     $ RES1(IPT),XIN(1,IPT,I),XIN(2,IPT,I),XIN(3,IPT,I),P(IPT),B(IPT)
      IF (LABEL.EQ.'TER '.OR.LABEL.EQ.'END ') GOTO 20
C
C Select backbone atoms only
C
      WRITE(22,*) AT1(IPT),RESNAM(IPT),RES1(IPT)
      WRITE(22,*) BLANK, RESREF(IPT), ATREF(IPT)
C     IF (AT1(IPT).NE.' N   '.AND.AT1(IPT).NE.' CA  '.AND.
C    $    AT1(IPT).NE.' C   '.AND.AT1(IPT).NE.' O   '.AND.
C    $    AT1(IPT).NE.' OT1 ') GOTO 10
C
C check to make sure the atoms from each model are the same ones
C
      IF (I.EQ.1) THEN
        ATREF(IPT) = AT1(IPT)
        RESREF(IPT) = RES1(IPT)
      ELSE
        IF (AT1(IPT).NE.ATREF(IPT).OR.RES1(IPT).NE.RESREF(IPT)) THEN
          WRITE(*,*) IPT,' mismatch',RES1(IPT),AT1(IPT),' is not ',
     $                               RESREF(IPT),ATREF(IPT)
          WRITE(*,*) I, RESNAM(IPT), RES1(IPT)
          PAUSE 'To continue anyway enter "C"'
        ENDIF
      ENDIF
C
C Since the atom matched, increment the atom counter
C
      IPT = IPT + 1
      GOTO 10
C
C All done with this coordinate set
C
 20   NATOMS = IPT - 1
      WRITE(*,*) NATOMS,' backbone atoms read for set',I
      IF (I.EQ.1) IPT1 = NATOMS
c
      IF (NATOMS.NE.IPT1) THEN
        WRITE(*,*) 'wrong number of atoms in file',I
        PAUSE 'To continue anyway enter "C"'
        IPT1 = MIN(IPT1,NATOMS)
        WRITE(*,*) 'Fitting based on the first',IPT1,' coordinates.'
        WRITE(*,1060) IPT1
      ENDIF
      GOTO 5
 25   CONTINUE
c
 30   CONTINUE
      NFILE = I
      WRITE(*,*) nfile,' models read in'
      WRITE(*,*) 'Coordinates to be grouped into two groups of'
      STOP
C    $             ,nn,' and ',mm,' sets'
      ITOT= NN + MM
      IF (NFILE.GT.ITOT) then
        NFILE = ITOT
        WRITE(*,*) '!!!!Only first',NFILE,' sets used in analysis!!!!'
      ELSEIF (NFILE.LT.ITOT) then
        stop 'not enough coordinate sets to break up as specified'
      ENDIF
C
      DMAX = 50.   ! number has no effect on the results as no cycling is done
      NCYC = 1
C
      NRES = IPT1/4
      WRITE(*,*) NRES,' residues in the coordinate set'
      WRITE(21,1008)  
      DO 300 IPEP = 2 , NRES-1
       IBEG = IPEP*4 - 6
       IEND = IPEP*4 + 2
       RES(IPEP) = RESREF(IBEG+5)  ! CA-atom of central residue
       NPTS = IEND - IBEG + 1
C
       KK = 0
       DO 39 K = IBEG , IEND      ! fill the reference array from the 
         KK = KK+1                ! first coordinate set
         DO 39 L = 1 , 3         
           XPEP(L,KK,1) = XIN(L,K,1)
           X(L,KK) = XIN(L,K,1)
 39     CONTINUE
C
      DO 200 IFILE = 2 , NFILE
        KK = 0
        DO 40 K = IBEG , IEND   ! fill the coordinate arrays with the 
        KK = KK + 1             ! set for comparison
        WT(KK) = 0.0  
        IF (KK.LE.5) WT(KK) = 1.0  ! overlay based on first peptide bond
        DO 40 L = 1 , 3         
          Y(L,KK) = XIN(L,K,IFILE)
 40     CONTINUE
C
C calculate and carry out the rotation
C
        CALL ITERATE(NCYC,DMAX,NPTS,Y,X,WT,ROT,RMS)
C
C transfer the rotated coordinates back into the storage array for averaging
C
        if (RMS.GT.0.1) then  ! this indicates a fundamental difference in
	  WRITE(*,*)'RMS',RMS,' FOR FILE,PEP',IFILE,IPEP ! peptide geometry
          do 1998 kkkk=1,5
            write(*,*) kkkk,(y(ll,kkkk),x(ll,kkkk), ll=1,3)
 1998     continue
        endif
        DO 38 K = 1 , NPTS
        DO 38 L = 1 , 3  
 38       XPEP(L,K,IFILE) = Y(L,K)
C

 200  CONTINUE
C
C Calculate the rms deviations for all coordinate pairs 
C
      AVE1 = 0.
      AVE12 = 0.
      AVE2 = 0.
      RMS1 = 0.
      RMS12 = 0.
      RMS2 = 0.
      NUM1 = 0
      NUM12 = 0
      NUM2 = 0
      dmin = 10.
      imin = 0
      jmin = 0
      DO 283 JFILE = 1 , NFILE
      DO 290 IFILE = JFILE+1 , NFILE
      SUMD2 = 0.0
      NUSED = NPTS - 5
      If (nused.ne.4) stop 'nused.ne.4!! nused='
      DO 285 I = 6, NPTS   ! LEAVE OUT FIRST PEPTIDE BOND
        D2 = (XPEP(1,I,IFILE) - XPEP(1,I,JFILE))**2 + 
     $       (XPEP(2,I,IFILE) - XPEP(2,I,JFILE))**2 +
     $       (XPEP(3,I,IFILE) - XPEP(3,I,JFILE))**2 
        SUMD2 = SUMD2 + D2
 285  CONTINUE
      RMSPEP(IFILE,JFILE) = SQRT(SUMD2/NUSED)
C
C compile the rms, average and minimal differences
C
      IF (jfile.le.nn.and.ifile.le.nn) then     ! intra first group
        AVE1 = AVE1 + RMSPEP(IFILE,JFILE)
        RMS1 = RMS1 + RMSPEP(IFILE,JFILE)**2
        NUM1 = NUM1 + 1
      ELSEIF (jfile.le.nn.and.ifile.gt.nn) then ! between groups
        AVE12 = AVE12 + RMSPEP(IFILE,JFILE)
        RMS12= RMS12 + RMSPEP(IFILE,JFILE)**2
        NUM12 = NUM12 + 1
        if (dmin.gt.rmspep(ifile,jfile)) then
          dmin = rmspep(ifile,jfile)
          imin = ifile
          jmin = jfile
        endif
      ELSE ! (jfile.gt.nn.and.ifile.gt.nn)        intra second group
        AVE2 = AVE2 + RMSPEP(IFILE,JFILE)
        RMS2 = RMS2 + RMSPEP(IFILE,JFILE)**2
        NUM2 = NUM2 + 1
      ENDIF
 290  CONTINUE
 283  CONTINUE
      IF (NUM1.NE.0) then
        AVE1 = AVE1/NUM1
        RMS1 = SQRT(RMS1/NUM1)
      ENDIF
      IF (NUM12.NE.0) then
        AVE12 = AVE12/NUM12
        RMS12 = SQRT(RMS12/NUM12)
      ENDIF
      IF (NUM2.NE.0) then
        AVE2 = AVE2/NUM2
        RMS2 = SQRT(RMS2/NUM2)
      ENDIF
      if (dmin.eq.10.) dmin = 0.
      WRITE(21,1011) RES(IPEP),AVE1,RMS1,AVE2,RMS2,AVE12,RMS12,
     $                         dmin,jmin,imin 
 300  CONTINUE
      STOP
 888  stop 'format error reading coordinate file'
 1000 FORMAT(A)
 1010 FORMAT(A4,2X,I5,1X,A5,A4,1X,A4,4X,3F8.3,2F6.2)
 1011 FORMAT(A4,3(2F6.2,8x),f6.2,2I4)
 1008 FORMAT('RES,AVE1,RMS1,AVE2,RMS2,AVE12,RMS12,jmin,imin,dmin')
 1009 FORMAT(I4,F6.2)
 1012 FORMAT(<NFILE-1>F4.1)
 1013 FORMAT('TER   ',I5,6X,A4,1X,A4)
 1014 FORMAT('END     ')
 1015 FORMAT('ATOM',10X,I5)
 1020 FORMAT(10X,' atom-1    atom-2     dist')
 1030 FORMAT(10X,2(2A4,2X),F7.2)
 1040 FORMAT(/,I6,' atomic coordinates read from model',I5) 
 1050 FORMAT(/,I6,' atomic coordinates read from ',A80) 
 1060 FORMAT(' !!!!!Fit based only on first',I5,' common coordinates',/)
 1070 FORMAT(/,I6,' atomic coordinates written to ',A80)
 1080 FORMAT(/' Best rotation of',I5,' equivalent atoms based on',
     $ I5,' atoms with',/,' final deviation less than',F6.2,' A.')
 1090 FORMAT(/,' The center of mass of set 1 before superposition is:',
     $ /,/,20X,3F10.4)
 1100 FORMAT(/,' The center of mass of set 2 is:',/,/,20X,3F10.4)
 1110 FORMAT(3(3F10.6/))
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
C       TYPE *,NUSE,' ATOMS TO BE USED FOR THE NEXT TRY.'
        IF (NSAME.EQ.NPT) THEN
C          WRITE(6,*) 'PROCESS CONVERGED IN',ICYC,'CYCLES.'
          RETURN
        ENDIF
C
 100  CONTINUE
C
C All done
C
C     WRITE(6,*) 'COMPLETED',NCYC,' ITERATIONS WITHOUT CONVERGENCE'
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
C      WRITE (8,*) 'NUMBER OF POINTS',NPT
C      WRITE (8,*) 'CENTER OF MASS SET 1'
C      WRITE (8,1000) CENX
C      WRITE (8,*) 'CENTER OF MASS SET 2'
C      WRITE (8,1000) CENY
C      WRITE (8,*) 'RMS DEVIATION',RMS
C      WRITE (8,*) 'ROTATION MATRIX'
C      WRITE (8,1000) ROT
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
