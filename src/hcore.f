      PROGRAM homologcore

C
C Reads the two PROTEIN DATA BANK files for the homologous proteins
C and loops through each equivalent secondary structural element, finding 
C the best rotation based on that element and then using a user defined
C cutoff evaluates the extent of the structural equivalence.
C Then using the full set of "structurally equivalent" residues, a final
C overlay is calculated and used to evaluate the structural divergence 
C of the two proteins.
C
C The second file should be the complete coordinate set that one desires
C to have rotated, even though only CA atoms will be used for the rotation
C
C Requires the two PDB files and one input file called "hcore.in"
C having the following information:
C	line 1: y or n (y=only use ca, n=use all atoms)
C	line 2: filename for PDB file 1
C	line 3: filename for PDB file 2 to be rotated onto pdb 1
C	line 4: cutoff in Angstrom for structural equivalence
C	lines 5-n: 1 line per structural element. Format(2(A5),1X,I4,1X)
C            residue name for first residue in element in protein 1
C            residue name for first residue in element in protein 2
C            number of residues in element (inclusive of ends)
C            example: 'x 203x 169x   6x
C
      INTEGER POINTER(20000,100),LENGTH(100),JPOINT(20000),NCA(2)
      INTEGER ICA,IPT
      CHARACTER*80 FILNAM(2),NAME
      CHARACTER*4 ATOMID,ATYP,CSEGID
      CHARACTER*4 RESNUM,NUMCA,RESBEG(2),OUTRES(0:10)
      CHARACTER*3 RESNAM
      CHARACTER*3 NAMCA
      CHARACTER*1 ONLYCA, achain(20000,2)
      REAL  XIN(3,20001,2),WT(20000),OUTDIST(10),XOUT(3),Q(20000,2)
      DIMENSION X(3,20001),Y(3,20001),ISITIN(20000),IBEG(2,100),
     $    NATOMS(2)
      DIMENSION ROT(3,3),DIST(20000,100),CENOLD(3),CENNEW(3),XREL(3)
      DIMENSION RESNUM(20000,2),ATOMID(20000,2),B(20000,2),
     $  P(20000,2),ATOMNUM(20000,2),
     $  RESNAM(20000,2),CSEGID(20000,2)
      DIMENSION NUMCA(20000,2),NAMCA(20000,2),XCA(3,20000,2)

      write(*,*)'Make sure you have the file hcore.in'
      write(*,*) ' '
      write(*,*) ' '
      write(*,*)'1) y on n (to or not to use only ca)'
      write(*,*)'2) PDB file 1'
      write(*,*)'3) PDB file 2 to be rotated onto pdb 1'
      write(*,*)'4) cutoff in Angstrom for structural equivalence'
      write(*,*)'5-n) 1 line per structural element. '
      write(*,*)'                        Format(2(1x,A4),1X,I4,1X)'
      write(*,*)' residue name for first residue in element in pdb 1'
      write(*,*)' residue name for first residue in element in pdb 2'
      write(*,*)' number of residues in element (inclusive of ends)'
      write(*,*) ' '
C
      OPEN(UNIT=5,FILE='hcore.in',STATUS='OLD')
      OPEN(UNIT=80,FILE='hcore.out')
C See if only want to use CA  y=only ca  n=use all atoms
        READ(5,1005) onlyca   
C
C Read in the coordinates for the two structures
C
      DO 30 I = 1 , 2
        READ(5,1000) FILNAM(I)
        OPEN(1,FILE=FILNAM(I),STATUS='OLD',READONLY)
        IPT = 1
        ICA = 0
 5      READ(1,1000,END=20) NAME
        IF (NAME(1:4).EQ.'ATOM') GOTO 7
        GOTO 5
 7      BACKSPACE(1)
 10     READ (1,1299,END=20,ERR=20) ATOMNUM(IPT,I),ATOMID(IPT,I),
     $     RESNAM(IPT,I), achain(IPT,I), RESNUM(IPT,I),
     $     XIN(1,IPT,I),XIN(2,IPT,I),XIN(3,IPT,I),B(IPT,I),
     $     Q(IPT,I),CSEGID(IPT,I)
 1299   FORMAT(6x,i5,1X,A4,1X,A3,1X,a1,A4,4x,3F8.3,2F6.2,6x,a4)

       IF (onlyca .eq. 'n' ) then
          ICA = ICA + 1
          XCA(1,ICA,I) = XIN(1,IPT,I)
          XCA(2,ICA,I) = XIN(2,IPT,I)
          XCA(3,ICA,I) = XIN(3,IPT,I)
          NUMCA(ICA,I) = RESNUM(IPT,I)
          NAMCA(ICA,I) = RESNAM(IPT,I)
       ENDif

       IF (onlyca .eq. 'y' ) then
        IF (ATOMID(IPT,I).EQ.' CA ') THEN
          ICA = ICA + 1
          XCA(1,ICA,I) = XIN(1,IPT,I)
          XCA(2,ICA,I) = XIN(2,IPT,I)
          XCA(3,ICA,I) = XIN(3,IPT,I)
          NUMCA(ICA,I) = RESNUM(IPT,I)
          NAMCA(ICA,I) = RESNAM(IPT,I)
        ENDIF
       ENDif

        IPT = IPT + 1
        GOTO 10
 20     NATOMS(I) = IPT - 1
        NCA(I) = ICA
        WRITE(80,*) NATOMS(I),' total  atoms read from  ',FILNAM(I)
        WRITE(80,*) NCA(I),   ' C-alpha atoms read from ',FILNAM(I)
        CLOSE(1)
 30   CONTINUE
C
C Get the cutoff in A for assessing structural equivalence
C
      READ(5,*) CUTOFF
      WRITE(80,*) CUTOFF,' A cutoff to be used'
C
C Do loop over each structural element listed in the input file
C
      ISEG = 0
 33   READ(5,1015,END=900,ERR=999) RESBEG(1),RESBEG(2),ILEN
      ISEG = ISEG + 1
      LENGTH(ISEG) = ILEN
C
      DO 35 I = 1 , 2
      DO 34 J = 1 , NCA(I)
        IF (RESBEG(I).EQ.NUMCA(J,I)) THEN
          IBEG(I,ISEG) = J
          GOTO 35
        ENDIF
 34   CONTINUE
      WRITE(*,*) 'RESIDUE NOT FOUND',RESBEG(I)
 35   CONTINUE
C
      DMAX = 50.
      NCYC = 1
C
C Overlay the whole structures based on the equivalent elements and
C not allowing any gaps
C
      IBEFOR = MIN(IBEG(1,ISEG),IBEG(2,ISEG)) - 1
      IAFTER = MIN(NCA(1)-IBEG(1,ISEG),NCA(2)-IBEG(2,ISEG))
      ITOTAL = IBEFOR + IAFTER + 1
      ICNT = 0
      DO 40 K = -IBEFOR , IAFTER, 1
        ICNT = ICNT + 1
      DO 40 L = 1 , 3
        X(L,ICNT) = XCA(L,IBEG(2,ISEG)+K,2)
        Y(L,ICNT) = XCA(L,IBEG(1,ISEG)+K,1)
        WT(ICNT) = 0.0
        IF (K.GE.0.AND.K.LT.LENGTH(ISEG)) WT(ICNT) = 1.0
 40   CONTINUE
      IF (ICNT.NE.ITOTAL) STOP 'Help, ICOUNT.ne.ITOTAL !' 
      CALL ITERATE(NCYC,DMAX,ITOTAL,X,Y,WT,ROT,RMS)
C
C write out locally rotated coordinates for coordinate file 2
C
C      OPEN(1,FILE='TEMP.PDB',STATUS='NEW')
C      ICNT = 0
C      DO 45 IPT = IBEG(2,ISEG)-IBEFOR , IBEG(2,ISEG)+IAFTER
C        ICNT = ICNT + 1
C        WRITE(1,1299) ATOMNUM(IPT,2),ATOMID(IPT,2),RESNAM(IPT,2),
C     $  RESNUM(IPT,2),X(1,ICNT),X(2,ICNT),X(3,ICNT),B(IPT,2),
C     $  Q(IPT,2),CSEGID(IPT,2)
C 45   CONTINUE
C      CLOSE(1)
C      
C use the rotated coordinates get separations of all aligned atoms
C
      DO 50 I = 1 , ITOTAL
        D2 = (X(1,I) - Y(1,I))**2 + 
     $       (X(2,I) - Y(2,I))**2 +
     $       (X(3,I) - Y(3,I))**2 
        ITEMP1 = IBEG(1,ISEG) - (IBEFOR+1) + I
        ITEMP2 = IBEG(2,ISEG) - (IBEFOR+1) + I
        POINTER(ITEMP1,ISEG) = ITEMP2
        DIST(ITEMP1,ISEG) = SQRT(D2)
 50   CONTINUE
C
C On to the next segment
C
      GOTO 33
C
C Done with the individual overlays. Summarize the results.
C
 900  CONTINUE
      WRITE(80,*) ISEG,' total elements compared'
C
C This is the place to write out all of the distances if desired
C      
C
C
C Now extend each element until an atom is found separated more than the cutoff
C
      DO 950 I = 1 , ISEG
        MIDPT = LENGTH(I)/2
        IMULT = 1.0
        DO 910 J = IBEG(1,I)+MIDPT , 1 , -1
          IF (DIST(J,I).GT.CUTOFF) IMULT = 0
          POINTER(J,I) = POINTER(J,I)*IMULT
 910    CONTINUE
        IMULT = 1.0
        DO 920 J = IBEG(1,I)+MIDPT , NCA(1)
          IF (DIST(J,I).GT.CUTOFF) IMULT = 0
          POINTER(J,I) = POINTER(J,I)*IMULT
 920    CONTINUE
 950  CONTINUE
C
C Synthesis of the complete equivalence list and output the summarized results
C
      WRITE(80,1232)
      DO 960 I = 1 , NCA(1)
        IEQUIV = 0
        POINTER(I,100) = 0
        OUTRES(0) = NUMCA(I,1)
        DO 955 J = 1 , ISEG
          IF (POINTER(I,J).NE.0) THEN
            IEQUIV = IEQUIV + 1
            OUTRES(IEQUIV) = NUMCA(POINTER(I,J),2)
            OUTDIST(IEQUIV) = DIST(I,J)
            POINTER(I,100) = POINTER(I,J)   ! Save the last one to be used
          ENDIF
 955    CONTINUE
        if (iequiv.eq.0) THEN
          WRITE(80,1236) OUTRES(0)
        ELSE
          WRITE(80,1233) OUTRES(0),(OUTRES(K),OUTDIST(K),K=1,IEQUIV)
        ENDIF
 960  CONTINUE
 1232 FORMAT(' Residue     equiv #1      equiv #2      equiv #3')
 1233 FORMAT(3X,A4,<IEQUIV>(5X,A4,F5.1))
 1236 FORMAT(3X,A4)
C
C Here one should now apply the final overlay and calculations
C
      ICNT = 0
      DO 450 K = 1 , NCA(1)
        IF (POINTER(K,100).EQ.0) GOTO 450
        ICNT = ICNT + 1
        DO 440 L = 1 , 3
          X(L,ICNT) = XCA(L,POINTER(K,100),2)
          Y(L,ICNT) = XCA(L,K,1)
          WT(ICNT) = 1.0
          JPOINT(ICNT) = K
 440    CONTINUE
 450  CONTINUE
      ITOTAL = ICNT  
      NCYC = 1
      DMAX = 50
C
      CALL ITERATE(NCYC,DMAX,ITOTAL,X,Y,WT,ROT,RMS)
C
      WRITE(80,*) ITOTAL,' core atoms were superimposed.'
      WRITE(80,1090)  X(1,ICNT+1),X(2,ICNT+1),X(3,ICNT+1)
      WRITE(80,1100)  Y(1,ICNT+1),Y(2,ICNT+1),Y(3,ICNT+1)
      WRITE(80,1110)  ROT
      WRITE(80,1120)  ITOTAL,RMS
      WRITE(80,1020) 
      DO 645 I = 1 , 3
        CENOLD(I) = X(I,ICNT+1)
 645    CENNEW(I) = Y(I,ICNT+1)
      DO 650 I = 1 , ICNT
         J = JPOINT(I)
         D2 = (X(1,I) - Y(1,I))**2 + 
     $        (X(2,I) - Y(2,I))**2 +
     $        (X(3,I) - Y(3,I))**2 
         XDIST = SQRT(D2)
         WRITE(80,1030) NUMCA(J,1),NUMCA(POINTER(J,100),2),XDIST
 650   CONTINUE
C
 1155 FORMAT(21I5)
 1160 FORMAT(21F5.2)
C
C Apply final rotation matrix to the complete coordinate set of file 2
C and write them out
C
      OPEN(1,FILE='rotated.pdb',STATUS='NEW')
      ICNT = 0
      DO 700 IPT = 1 , NATOMS(2)
        ICNT = ICNT + 1
        DO 705 I = 1 , 3
 705      XREL(I) = XIN(I,IPT,2) - CENOLD(I)
        DO 720 I = 1,3
        XOUT(I) = CENNEW(I)
        DO 710 J = 1,3
 710      XOUT(I) = XOUT(I) + ROT(J,I)*XREL(J)
 720    CONTINUE
        WRITE(1,1301) ATOMNUM(IPT,2),ATOMID(IPT,2),RESNAM(IPT,2),
     $    achain(IPT,2), RESNUM(IPT,2),
     $    XOUT,B(IPT,2),Q(IPT,2),CSEGID(IPT,2)
 700  CONTINUE
      WRITE(6,*) NATOMS(2),' atoms written to file ROTATED.PDB'
      CLOSE(1)
c
 999  STOP
 1301 FORMAT('ATOM',2X,i5,1x,A4,1x,A3,1X,a1,A4,4X,3F8.3,2F6.2,6x,a4)
 1000 FORMAT(A)
 1005 FORMAT (A1)
 1010 FORMAT(13X,A4,4X,A4,4X,3F8.3)
 1015 FORMAT(2(1x,A4),1X,I4)
 1020 FORMAT(10X,' atom-1    atom-2     dist')
 1030 FORMAT(10X,2(2X,A4,4X),F7.2)
 1040 FORMAT(/,I6,' atomic coordinates read from ',I5,' files') 
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
      DIMENSION XINIT(3,20000)
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
      DIMENSION XREL(3,20000),YREL(3,20000)
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
