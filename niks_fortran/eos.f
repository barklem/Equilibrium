C=========================================================================
C EQCOUNT: Counts the list of species for solving the equation of state by
C          merging the default list and species present in the line list.
C
C We assume that only neutral molecules can appear in the line list.
C For atoms, all the ions present in the table of partition functions
C are added to the list. Atomic names are case sensitive, that is the first
C character must be uppercase and for 2 character names the second character
C must be lower case.
C
C Inputs:
C   ELEMEN - the names of chemical elements in the periodic table
C   SPNAME - the names of the species present in the line lists + continuous
C            absorbers
C   ION    - ionization stage (1 -neutral, 2 - first ion etc.)
C   NLINES - the length of the line list, also dimensions of arrays SPNAME,
C            ION, SPINDX
C   NLIST  - if >0 on input, indicates that the default list of species have
C            been loaded, otherwise EQLIST loads the default list to SPLIST.
C   ELESIZ - Size of all arrays related to atomic list.
C
C   Return code  0: OK
C                1: illegal species name
C               >1: SPLSIZ is too small
C
      integer function eqcount(elemen,spname,ion,nlines,nlist,ELESIZ)
      INCLUDE 'SIZES.EOS'

      integer nlines,nlist,ELESIZ
      character*(3) elemen(ELESIZ)
      character*2 tmp
      character*(SPCHAR) spname(nlines)
      character*(SPCHAR) tmplist(SPLSIZ),chname
      integer ion(nlines),ionmax
      real a(IONSIZ)
      double precision b(IONSIZ)
      INCLUDE 'DEFAULT.EOS'
C
      eqcount=0
      ionmax=0
      do 1 ispec=1,NDEF
        tmplist(ispec)=default(ispec)
   1  continue
      ncount=NDEF
C
C Associate each species in SPNAME with an entry in SPLIST. If SPNAME
C  contains a new species not in SPLIST, then add that new species at
C  the end of SPLIST.
C
      if(nlines.gt.0) then
        do 5 ilin=1,nlines
          call mbuild(spname(ilin),ion(ilin)-1,chname)
c          write(*,*) ncount,ilin,ionmax,spname(ilin),chname
          do 2 ispec=1,ncount
            if(tmplist(ispec).eq.chname) goto 5
   2      continue
c          write(*,*) ncount,ilin,chname
C
C Look for atomic species. Negative ions (e.g. H-) are treated as molecules
C
          if((spname(ilin)(2:2).EQ.' '.OR.
     *       (spname(ilin)(3:3).EQ.' '.AND.
     *        spname(ilin)(2:2).GE.'a'.AND.
     *        spname(ilin)(2:2).LE.'z')).AND.
     *        ion(ilin).GT.0) then
            iel=0
            tmp=spname(ilin)(1:2)
            do 3 i=1,ELESIZ
              if(tmp.eq.elemen(i)(1:2)) iel=i
   3        continue
            if(iel.lt.1) then
              eqcount=1
              return
c              write(*,*) 'eqlist: Wrong species: ',spname(ilin)
c              stop
            end if
            call XSAHA(iel,1.,1.,1.,ionmax,a,b,IONSIZ,5)
            if(ncount+ionmax.gt.SPLSIZ) then
              eqcount=ncount+ionmax
              return
            endif
            tmplist(ncount+1)=elemen(iel)(1:2)
            do 4 i=2,ionmax
              ncount=ncount+1
              i1=index(tmplist(ncount),' ')
              tmplist(ncount+1)=tmplist(ncount)(1:i1-1)//'+'
   4        continue
            ncount=ncount+1
          else
C
C Molecules are counted here
C
            if(ncount.ge.SPLSIZ) then
              eqcount=ncount+1
              return
            endif
            tmplist(ncount+1)=chname
            ncount=ncount+1
          end if
   5    continue
      endif
C
C All lines have been processed, add free electrons and return
C
      nlist=ncount+1
      eqcount=0
C
      return
      end

C=========================================================================
C EQLIST: Creates the list of species for solving the equation of state by
C         merging the default list and species present in the line list.
C
C We assume that only neutral molecules can appear in the line list.
C For atoms, all the ions present in the table of partition functions
C are added to the list. Atomic names are case sensitive, that is the first
C character must be uppercase and for 2 character names the second character
C must be lower case.
C
C Inputs:
C   ELEMEN - the names of chemical elements in the periodic table
C   SPNAME - the names of the species present in the line lists + continuous
C            absorbers
C   ION    - ionization stage (1 -neutral, 2 - first ion etc.)
C   NLINES - the length of the line list, also dimensions of arrays SPNAME,
C            ION, SPINDX
C   NLIST  - if >0 on input, indicates that the default list of species have
C            been loaded, otherwise EQLIST loads the default list to SPLIST.
C   SPLDIM - maximum length of the compiled lists of species SPLIST (must
C            be smaller than SPLSIZ).
C   ELESIZ - Size of all arrays related to atomic list.
C
C Outputs:
C   SPINDX - index array of size NLINES which upon return holds pointers to
C            the complete list of species SPLIST: line L is produced by
C            species SPLIST(SPINDEX(L))
C   SPLIST - upon return contains the compiled list of all species (default
C            list + species in the line list + continuous absorbers)
C   NLIST  - the size of the compiled list of species SPLIST
C
C   Return code  0: OK
C                1: illegal species name
C                2: SPLDIM is too small)
C                3: Missing ionization stage
C                4: e- is not the last item in the list
C                5: Unreasonable abundances
C
C  2006.12.27 - converted eqlist to a function for compatibility with the SME
C
C
      integer*4 function eqlist(abund,elemen,spname,ion,spindx,splist,
     &                  nlines,nlist,SPLDIM,ELESIZ)
      INCLUDE 'SIZES.EOS'

      integer nlines,nlist,SPLDIM,ELESIZ
      character*(SPCHAR) spname(nlines),splist(SPLDIM)
      character*(3) elemen(ELESIZ)
      character*2 tmp
      integer ion(nlines),spindx(nlines)
      dimension abund(ELESIZ)
      real a(IONSIZ)
      double precision b(IONSIZ)
C
C SPLIST should contain all the major contributors to the electron pressure,
C and all the molecules which significantly affect atomic partial pressures.
C For each call to EQSTAT, the base set of species at the beginning of SPLIST
C are supplemented by any new species that appear in SPNAME. It is common
C for some of the species in the base set (at the beginning of SPNAME) to be
C duplicated in SPNAME. This allows one to get ZETA for these species and is
C not a problem.
C
      integer splmax
      character*(SPCHAR) chname
      INCLUDE 'DEFAULT.EOS'
C
C Determine maximum allowed number of species, based on sizes of arrays
C  defined locally (using SPLSIZ) and passed by argument (using spldim).
C
      splmax=min(SPLSIZ,SPLDIM)
C
C Load base set of species (SPLIST) with default set of species (DEFAULT),
C  if passed value of NLIST is 0. Be sure to include "e-" at the end of
C  SPLIST.
C
      idef=0
      ionmax=0
      if(nlist.eq.0) then
C
C  nlines set to -1 indicates that we need to get partial pressures for all atoms
C  This mode is meant for use within VALD
C
        if(nlines.eq.-1) then
c
c Add all atoms first (the call to XSAHA is dummy,
C just to get the number of ions available in the table)
c
          do 20 iel=1,ELESIZ
            call XSAHA(iel,1.,1.,1.,ionmax,a,b,IONSIZ,5)
            idef=idef+1
            if(idef.gt.splmax) goto 900
            splist(idef)=elemen(iel)(1:2)
            do 10 i=2,ionmax
              idef=idef+1
              if(idef.gt.splmax) goto 900
              splist(idef)=splist(idef-1)
              isp=index(splist(idef),' ')
              if(isp.le.0) then
c                write(*,*) 'eqlist: Insufficient length of splist ',
c     *                     'elements to store ion',elemen(iel)(1:2),i
c                stop
                 eqlist=2
                 return
              endif
              splist(idef)(isp:isp)='+'
  10        continue
  20      continue
        endif
C
C  Copy the default list
C
        if(NDEF+idef.gt.splmax) goto 900
        do 30 jdef=1,NDEF
          splist(idef+jdef)=default(jdef)
  30    continue
        nlist=NDEF+idef
c        nlist=idef
      endif
C
C Check that abundances are sensible.
C
      absum=0.0
      do 50 ielem=1,ELESIZ
        if(abund(ielem).lt.0.0.or.abund(ielem).gt.1.0) then
          write(*,40) ielem,abund(ielem)
  40      format('eqlist: bad abundance for element',i3,':',1pe13.4)
          write(*,*) (abund(ispec),ispec=1,99)
c          stop
          eqlist=5
          return
        endif
        absum=absum+abund(ielem)
  50  continue
c      do 60 ielem=1,ELESIZ
c        abund(ielem)=abund(ielem)/absum
c  60  continue
c      if(abs(absum-1.0).gt.1.0e-3) then
c        write(*,70) absum
c  70    format('eqlist: warning! abundances are not normalized:'
c     &           ,1pe13.5)
c      endif

C
C Associate each species in SPNAME with an entry in SPLIST. If SPNAME
C  contains a new species not in SPLIST, then add that new species at
C  the end of SPLIST.
C
      do 100 ispec=nlist+1,splmax
        splist(ispec)='        '
 100  continue
      inew=nlist+1
      if(nlines.gt.0) then
        do 150 ilin=1,nlines
          call mbuild(spname(ilin),ion(ilin)-1,chname)
          do 110 ispec=1,nlist
            if(splist(ispec).eq.chname) then
              spindx(ilin)=ispec
              goto 150
            endif
 110      continue
C
C Look for atomic species. Negative ions (e.g. H-) are treated as molecules
C
c          write(*,*) spname(ilin),ion(ilin),chname
          if((spname(ilin)(2:2).EQ.' '.OR.
     *       (spname(ilin)(3:3).EQ.' '.AND.
     *        spname(ilin)(2:2).GE.'a'.AND.
     *        spname(ilin)(2:2).LE.'z')).AND.
     *        ion(ilin).GT.0) then
            iel=0
            tmp=spname(ilin)(1:2)
            do 120 i=1,ELESIZ
              if(tmp.eq.elemen(i)(1:2)) iel=i
 120        continue
            if(iel.lt.1) then
c              write(*,*) 'eqlist: Wrong species: "'//spname(ilin)//'"'
c              stop
              eqlist=1
              return
            end if
            call XSAHA(iel,1.,1.,1.,ionmax,a,b,IONSIZ,5)
C
C  Make sure that neutral atoms are included as well as all
C  the intermediate ions
C
            do 130 ii=0,ionmax-1
            if(inew.gt.splmax) goto 900
            call mbuild(spname(ilin),ii,chname)
            splist(inew)=chname
            if(ii.eq.ion(ilin)-1) spindx(ilin)=inew
 130        inew=inew+1
          else
c       write(*,*) 'Molecule: '//chname,inew
            if(inew.gt.splmax) goto 900
            splist(inew)=chname
            spindx(ilin)=inew
            inew=inew+1
          end if
          nlist=inew-1
 150    continue
      endif
C
C Make sure free electrons are the last species in the list.
C
      do 170 ispec=1,nlist-1
      if(splist(ispec).eq.'e-') then
c        write(*,*) 'eqlist: "e-" may only occur at the end of the'
c     &          // ' species list (SPLIST).'
c        stop
        eqlist=4
        return
      endif
 170  continue
      if(splist(nlist).ne.'e-') then
        nlist=nlist+1
        if(nlist.gt.splmax) goto 900
        splist(nlist)='e-'
      endif
C
C Make sure neutral hydrogen and neutral helium are in SPLIST. These
C  species are needed for H1FRCT and HE1FRCT. Remember the locations
C  of these species in SPLIST for later use. Code is optimized for
C  the case where H and He both occur early in SPLIST list.
C
c      ih1=-1
c      do 200 ispec=1,nlist
c        if(splist(ispec).eq.'H') then
c          ih1=ispec
c          goto 210
c        endif
c 200  continue
c      write(*,*) 'eqlist: "H" must be in species list (SPLIST)'
c      stop
c 210  ihe1=-1
c      do 220 ispec=1,nlist
c        if(splist(ispec).eq.'He') then
c          ihe1=ispec
c          goto 230
c        endif
c 220  continue
c      write(*,*) 'eqlist: "He" must be in species list (SPLIST)'
c      stop
c 230  continue
C
C Sort the list
C
      call sort2(nlist,splist,nlines,spindx,elemen,ELESIZ)
c      do 250 ispec=1,nlist
c 250  write(*,*) ispec,' "',splist(ispec),'"'
C
      eqlist=0
      return
C
C Error handlers.
C
 900  continue
c      write(*,905) splmax
c 905  format('eqlist: species list (SPLIST) not long enough:',i5,a)
c      stop
      eqlist=2
      return
      end

c
C=========================================================================
C EQSTAT: Determine thermodynamic quantities required for spectroscopy.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   PTOTAL [real] Total gas pressure (in dyne/cm^2), given by NTOTAL*K*T,
C     which is to be used in calculating chemical and ionization equilibrium,
C     and partial pressures.
C   PELEC [real] Electron pressure (in dyne/cm^2), given by NELEC*K*T,
C     which is to be used in calculating ionization equilibrium.
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPNAME [character*(*) array(NLINES)] Case-sensitive species name of atom
C     or molecule. The first letter of each atom name must be uppercase. The
C     second letter of each atom name, if present, must be lowercase. Each
C     atom name may optionally be followed by a multiplicity number between
C     1 and 4. If no multiplicity number is given for a particular atom, then
C     its multiplicity is assumed to be 1. All atomic and molecular species
C     in SPNAME must be neutral, with the charge state specified separately
C     in the ION input argument.
C   ION [integer array(NLINES)] Charge state for each of the atomic and
C     molecular species specified in SPNAME. ION=-1 for negative ions (e.g.
C     H minus), ION=0 for neutrals, ION=1 for singly ionized species, etc.
C   NLINES [integer] Number of valid entries in SPNAME and ION. From an
C     external perspective, each entry in SPNAME and ION will correspond to
C     a single spectral line, so some specie/charge combinations may appear
C     more than once, while others may not appear at all.
C   SPLDIM [integer] Array sizes for the arguments SPLIST and XFRACT, which
C     contain information for each species. The maximum allowed number of
C     species is SPLMAX=MIN(SPLSIZ,SPLDIM), where SPLSIZ is a parameter
C     defined in the file SIZES.SYN and used to dimension the local arrays
C     XNPF, PFUNC, and POTION. SPLMAX must be large enough to handle the
C     base set of species used when computing the molecular equilibrium and
C     also any additional species that appear only in the line list. Ideally,
C     the calling routine will <1> Include SIZES.SYN, <2> Use SPLSIZ to
C     dimension SPLIST and XFRACT, and <3> Pass SPLSIZ in place of SPLDIM.
C     However, SPLDIM is passed separately to allow for error checking in
C     the cases when this is not done (e.g. when called from IDL).
C   MODE [integer] Determines the content of the output:
C      1    - number densities
C      2    - partition functions
C      3    - partial pressures
C  other<10 - number densities/partition functions
C       +10 - the same as above but electron density is assumed to be known
C             precisely and not re-determined in the process
C
C Input/Output:
C   SPLIST [character*(*) array(SPLDIM)] If NLIST is nonzero upon entry,
C     then SPLIST must contain the base set of species that must be included
C     in the molecular equilibrium calculation, regardless of which species
C     are represented by lines in SPNAME. Until the code is cleaned up, the
C     species list in SPLIST must include "e-" after the NLIST element.
C     If NLIST is zero upon entry, then SPLIST is loaded with the base set
C     of species coded into EQSTAT below (in the variable DEFAULT). Again,
C     an "e-" is appended after the base set.
C     Regardless of the whether SPLIST is valid upon entry or needs to be
C     loaded with the defaults, species that are in the lines list SPNAME,
C     but are not in the base set of species will be inserted into SPLIST
C     after the "e-" entry. Currently, the extended list is not used, but
C     in the future, we may solve for the equilibrium of all species in the
C     extended SPLIST.
C   NLIST [integer] If nonzero upon entry, NLIST is the number of species
C     in the base set of species passed in SPLIST (including the mandatory
C     "e-" at the beginning of the list). If NLIST is zero upon entry, this
C     indicates that the default base set of species coded in EQSTAT should
C     be used. Upon exit, NLIST is set to the number of species in SPLIST,
C     which contains the base set plus any additional species that occur
C     in the line list.
C
C Outputs:
C   SPINDX [integer array(NLINES)] Species index assigned to each line in
C     the input line list (specified by the input arguments SPNAME and ION).
C     The species index is used to reconstruct the species name (in SPLIST)
C     or "zeta" value (in XFRACT) computed for each line in the input line
C     list. For example, ZETA(SPINDX(370)) contains the zeta value for the
C     line corresponding to SPNAME(370) and ION(370).
C   XFRACT [real array(SPLDIM)] Zeta (in cm^-3) for the atomic or molecular
C     species in the corresponding entry of SPNAME and the charge state in
C     corresponding entry of ION. Zeta is the number density divided by the
C     partition function, and is required for spectrum synthesis.
C   POTI   [real array(SPLDIM)] ionization potential in eV for the
C     corresponding species.
C   ATWGHT [real array(SPLDIM-1)] molecular weights in AMU for the
C     corresponding species.
C   H1FRCT [real] Number density (in cm^-3) of neutral atomic hydgrogen,
C     used in computing damping constants (and continuous opacities?).
C   HE1FRCT [real] Number density (in cm^-3) of neutral atomic helium,
C     used in computing damping constants (and continuous opacities?).
C
      subroutine eqstat(mode,temp,Pg,Pe,abund,elemen,amass,
     &                  ELESIZ,spindx,splist,xfract,poti,atwght,
     &                  h1frct,he1frct,nlines,nlist,xne,xna,rho,niter)
c      subroutine eqstat(mode,temp,xnatom,xnelec,abund,elemen,amass,
c     &                  ELESIZ,spindx,splist,xfract,poti,atwght,
c     &                  h1frct,he1frct,nlines,nlist,xne,xna,rho,niter)
cc     &                  FAILED)
      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'

      integer mode,ELESIZ,niter
      integer nlines,nlist
      real temp,Tk,Pg,Pe,Pgas,Pelec,xna,xne,rho,
     *     h1frct,he1frct
c      real xnatom,xnelec,xne_old,xna_old
      real Pg_old,Pe_old
      character*(SPCHAR) splist(nlist)
      character*(3) elemen(ELESIZ)
      integer spindx(nlines)
      real xfract(nlist),poti(nlist),atwght(nlist)
      real abund(ELESIZ),amass(ELESIZ)
      logical FAILED

      integer Anum(4),Natm(4),maxion,nelm,nchg
      real xnpf(SPLSIZ),pfunc(SPLSIZ),tol,tol1,xtotal
      real potion(IONSIZ),wtmol
      double precision awt(SPLSIZ-1),fract(IONSIZ)
      integer icharge,iter,ispec,IH1,IHe1,mmode

      INTEGER MAXITER
      REAL kBol
      DOUBLE PRECISION PSI
      PARAMETER (kBol=1.38065E-16,MAXITER=5000)
C
C Call equation of state solver.
C
c      open(87,file='dumpb.dat',form='unformatted',status='old')
c      read(87) temp,Pgas,Pelec,abund,elemen,amass,
c     &         mmode,spindx(nlines),splist,nlines,nlist
c      close(87)
      TOL=1.E-5
      TOL1=1.E-3
      Pgas=Pg
      Pelec=Pe
      PSI=2.d0/(1.d0+SQRT(5.d0))
      DO ISPEC=1,NLIST
        IF(SPLIST(ISPEC).EQ.'H  ') IH1 =ISPEC
        IF(SPLIST(ISPEC).EQ.'He ') IHE1=ISPEC
        XNPF(ISPEC)=-1.
      END DO
      Tk=temp*kBol
c      Pgas=(xnatom+xnelec)*Tk
      mmode=mod(mode,10)

      if(temp.gt.12000.) then
C
C Hot gas: assume no molecules and use Saha equation
C
        niter=1
        if(mode.ge.10) then
          call Nelect(temp,Pgas,abund,amass,ELESIZ,
     *                xna,xne,h1frct,he1frct,wtmol)
          Pelec=xne*Tk
        else
          xne=Pelec/Tk
        endif
c        xne=xnelec
c        xna=xnatom
        xna=(Pgas-Pelec)/Tk
        call xsaha(1,temp,xne,xna,2,potion,fract,IONSIZ,2)
        h1frct=fract(1)
        call xsaha(2,temp,xne,xna,3,potion,fract,IONSIZ,2)
        he1frct=fract(1)

        rho=xna*wtmol
        do 2 ispec=1,nlist-1
        CALL MPARSE(elemen,splist(ispec),Nelm,Nchg,Anum,Natm,ELESIZ)
        icharge=Nchg+1
        if(Nelm.eq.1.and.Natm(1).eq.1.and.Nchg.ge.0) then
C
C Get the number of ionization stages available in XSAHA
C
          call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,
     *               IONSIZ,5)
C
C Atom. Parser returns atomic number in Anum(1)
C
          if(mmode.eq.1) then
C
C  MODE=1, Return number densities
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,
     *                 IONSIZ,2)
            xfract(ispec)=fract(icharge)*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          else if(mmode.eq.2) then
C
C  MODE=2, Return partition functions
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,
     *                IONSIZ,3)
            xfract(ispec)=fract(icharge)
            poti(ispec)=potion(icharge)
          else if(mmode.eq.3) then
C
C  MODE=3, Return partial pressures
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,
     *                 IONSIZ,2)
            xfract(ispec)=fract(icharge)*kBol*temp*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          else
C
C  Any other MODE: Return number densities / partition functions
C
            call xsaha(Anum(1),temp,xne,xna,maxion,potion,fract,
     *                 IONSIZ,1)
            xfract(ispec)=fract(icharge)*xna*abund(Anum(1))
            poti(ispec)=potion(icharge)
          endif
          atwght(ispec)=amass(Anum(1))
        else
C
C Ignore molecules
C
          poti(ispec)  =1.
          atwght(ispec)=1.
          xfract(ispec)=0.
        endif
  2     continue
C
C Electrons
C
        if(mmode.eq.1) then
          xfract(nlist)=xne
        else if(mmode.eq.2) then
          xfract(nlist)=1.
        else if(mmode.eq.3) then
          xfract(nlist)=xne*Tk
        else
          xfract(nlist)=xne
        endif
      else
C
C Cold gas: solve molecular equilibrium
C
        niter=0
c        write(*,*) NLINES,NLIST,temp,Pgas,Pelec,mmode
c        write(*,'(10f8.3)') log10(abund)
        IF(mode.ge.10) then
          if(temp.gt.4000.) then
            Pe_old=Pgas*0.1
          else if(temp.gt.2000.) then
            Pe_old=Pgas*0.01
          else
            Pe_old=Pgas*0.001
          endif
        else
          Pe_old=Pelec
        endif
        Pg_old=Pg
c        IF(mode.ge.10) then
c          if(temp.gt.4000.) then
c            xne_old=xnatom*0.1
c          else if(temp.gt.2000.) then
c            xne_old=xnatom*0.01
c          else
c            xne_old=xnatom*0.001
c          endif
c        else
c          xne_old=xnelec
c        endif
  3     continue
        if(Pelec/Pgas.lt.0.01) then
c        if(xnelec/xnatom.lt.0.01) then
          call GAS(temp,Pg_old,Pe_old,abund,elemen,amass,
     *             ELESIZ,tol,splist,nlist,
     *             xne,xna,rho,Pgas,xnpf,pfunc,poti,xtotal,
     *             awt,iter,FAILED)
c          call lnGAS(temp,xne_old,xna_old,abund,elemen,amass,
c     *             ELESIZ,tol,
c     *             splist,nlist,xne,xna,rho,xnpf,pfunc,poti,xtotal,
c     *             awt,iter,FAILED)
        else
          call GAS(temp,Pg_old,Pe_old,abund,elemen,amass,
     *             ELESIZ,tol,splist,nlist,
     *             xne,xna,rho,Pgas,xnpf,pfunc,poti,xtotal,
     *             awt,iter,FAILED)
c          call GAS(temp,xne_old,xna_old,abund,elemen,amass,
c     *             ELESIZ,tol,
c     *             splist,nlist,xne,xna,rho,xnpf,pfunc,poti,xtotal,
c     *             awt,iter,FAILED)
        endif
        niter=niter+iter
        IF(niter.ge.MAXITER) THEN
          Pelec=xne*Tk
          WRITE(*,*) 'T,Pgas,Pnew,Pelec,Pe_in,Pe_out,NITER=',
     *                Temp,Pgas,Pg,Pelec,Pe_old,Pe,niter,FAILED
c          WRITE(*,*) 'T,Pgas,Pnew,XNA_in,XNA_out,XNE_in,XNE_out=',
c     *                Temp,Pgas,Pnew,xna_old,xna,xne_old,xne,niter,
c     *                FAILED
          IF(niter.gt.MAXITER*20) STOP
        END IF
        IF(mode.lt.10.and.
     *    (abs(Pgas -Pg_old)/max(1.E-20,Pgas ).gt.tol1.or.
     *     abs(Pelec-Pe_old)/max(1.E-20,Pelec).gt.tol1)) THEN
          Pe_old=Pelec
          Pg_old=Pg
c     *    (abs(xne-xne_old)/max(1.E-20,xne).gt.tol1.or.
c     *     abs(xna-xna_old)/max(1.E-20,xna).gt.tol1)) THEN
c          xne_old=xne
c          xna_old=xnatom
          GOTO 3
        END IF
c        write(*,*) 'T, P', Temp, Pg
c        do ispec=1,nlist-1
c          write(*,*) ispec,splist(ispec),xnpf(ispec)
c        enddo
c      write(*,'(F10.1,13E11.4)') Temp,xnpf(1),
c     &                                xnpf(2),
c     &                                xnpf(3),
c     &                                xnpf(4),
c     &                                xnpf(5),
c     &                                xnpf(6),
c     & (Pgas-Pelec)/Tk,xna,Pelec/Tk,xne,rho
C
C Fill return arrays.
C
        do 4 ispec=1,nlist-1
  4     atwght(ispec)=awt(ispec)
C
        if(mmode.eq.1) then
C
C  MODE=1, Return number densities
C
          do 5 ispec=1,nlist-1
  5       xfract(ispec)=xnpf(ispec)*pfunc(ispec)
          xfract(nlist)=xne
        else if(mmode.eq.2) then
C
C  MODE=2, Return partition functions
C
          do 6 ispec=1,nlist
  6       xfract(ispec)=pfunc(ispec)
          xfract(nlist)=1.
        else if(mmode.eq.3) then
C
C  MODE=3, Return partial pressures
C
          do 7 ispec=1,nlist-1
  7       xfract(ispec)=xnpf(ispec)*pfunc(ispec)*Tk
          xfract(nlist)=xne*Tk
        else
C
C  Any other MODE: Return number densities / partition functions
C
          do 8 ispec=1,nlist-1
  8       xfract(ispec)=xnpf(ispec)
          xfract(nlist)=xne
        endif
C
        h1frct=xnpf(ih1)*pfunc(ih1)
        he1frct=xnpf(ihe1)*pfunc(ihe1)
      endif
C
      return
      end

      function llength(name,elemen,ELESIZ)
C
C  Returns an almost unique integer for molecule "name" which
C  is assumed to include up to 4 different types of atoms.
C  For molecule A1_n1 A2_n2 A3_n3 A4_n4 Ch
C  llength = (n1 + n2 + n3 + n4)*10000 + (Z1 + Z2 + Z3 + Z4)*10 + charge
C  Charge of -1 corresponds to 9. Positive charge is limited to +8.
C
      integer iel(4),nat(4),charge,ELESIZ
      character*(*) name
      character*(3) elemen(ELESIZ)
C
      call mparse(elemen,name,nel,charge,iel,nat,ELESIZ)
      llength=0
      do 1 i=1,nel
   1  llength=llength+iel(i)*10+10000*nat(i)
      if(charge.gt.0) then
        llength=llength+charge
      else if(charge.lt.0) then
        llength=llength+9
      end if
C
      return
      end

C=========================================================================
C NELECT: Finds consistent electron number density.
C
C Inputs:
C   T [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   P [real] Total gas pressure (in dyne/cm^2), given by NTOTAL*K*T,
C     which is to be used in calculating chemical and ionization equilibrium,
C     and partial pressures.
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   AMASS [real array(ELESIZ)] atomic weights in AMU.
C Outputs:
C   XNA    [real] Atomic number density
C   XNE    [real] Electron number density
C   H1FRC  [real] Number density (in cm^-3) of neutral atomic hydgrogen,
C      used in computing damping constants.
C   HE1FRC [real] Number density (in cm^-3) of neutral atomic helium,
C      used in computing damping constants.
C   WTMOLE [real] Mean molecular weight in AMU.
C
      SUBROUTINE NELECT(T,P,ABUND,AMASS,ELESIZ,
     *                  XNA,XNE,H1FRC,HE1FRC,WTMOLE)
C
C
C  AUTHOR: N.Piskunov
C
C  LAST UPDATE: 29 January 1993
C
      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      INTEGER ELESIZ
      REAL T,P,XNE,XNA,H1FRC,HE1FRC,WTMOLE
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      DOUBLE PRECISION kBol
      PARAMETER (kBol=1.38065D-16)

      DOUBLE PRECISION FRACT(IONSIZ)
      DOUBLE PRECISION TK,XNTOT,XNENEW,X,XA,XE,ERROR
      REAL POTI(IONSIZ)
      INTEGER L,IEL,ION,MAXION
C
      TK=kBol*T
      XNTOT=P/TK
      XE=XNTOT*0.5D0
      XA=XE
      DO 4 L=1,200
        XNENEW=0.D0
        DO 2 IEL=1,ELESIZ
          X=0.D0
          XNE=XE
          XNA=XA
C
C  Get the number of known ions
C
          CALL XSAHA(IEL,T,XNE,XNA,MAXION,POTI,FRACT,IONSIZ,5)
C
C  Get the number of electrons contributed by all ions of atom IEL
C
          CALL XSAHA(IEL,T,XNE,XNA,MAXION,POTI,FRACT,IONSIZ,2)
          IF(IEL.EQ.1) H1FRC =FRACT(1)
          IF(IEL.EQ.2) HE1FRC=FRACT(1)
          DO 1 ION=1,MIN(MAXION,IEL+1)
            X=X+FRACT(ION)*(ION-1)
   1      CONTINUE
          XNENEW=XNENEW+X*XA*ABUND(IEL)
   2    CONTINUE
        XNENEW=(XNENEW+XE)*0.5D0
        ERROR=ABS((XE-XNENEW)/XNENEW)
        XE=XNENEW
        XA=XNTOT-XE
c        write(*,'('' T,XNE,XNA,ERROR='',F8.1,3E14.6)') T,XNE,XNA,ERROR
        IF(ERROR.LT.1.D-5) THEN
          X=0.D0
          DO 3 IEL=1,99
            X=X+ABUND(IEL)*AMASS(IEL)
   3      CONTINUE
          WTMOLE=X*1.660E-24
          RETURN
        END IF
   4  CONTINUE
      WRITE(*,*) 'Can''t converge calculating electron density'
C
      STOP
      END

C=========================================================================
C SORT2: sorts two arrays in atomic element order of the first (character) array.
C Hydrogen first, Helium next etc. All atoms/ions must end up before molecules
C that contain this atoms.
C
      subroutine sort2(nlist,list1,nlines,list2,elemen,ELESIZ)
      include 'SIZES.EOS'
c
      integer nlist,nlines,ELESIZ
      character*(*) list1(nlist)
      character*(3) elemen(ELESIZ)
      character*(SPCHAR) name,name1,name2
      integer list2(nlines)
c
c Go through the list (except the last item which is e-)
c
      i=0
   1  if(i.lt.nlist-2) then
c
c Set the first entry as the minimum rank in the remaining part of the list
c
        i=i+1
        imin=i
        name2=list1(imin)
        l2=llength(name2,elemen,ELESIZ)
c
c Go through other entries. Look for smaller or identical ranks.
c
        j=i
   2    if(j.lt.nlist-1) then
          j=j+1
          name1=list1(j)
          l1=llength(name1,elemen,ELESIZ)
          if(l1.lt.l2.or.(l1.eq.l2.and.name1.lt.name2)) then
c
c Found smaller rank. Store the location of the new winner.
c
            imin=j
            name2=list1(imin)
            l2=llength(name2,elemen,ELESIZ)
c            if(list1(list2(4)).eq.'e-') write(*,*) 'A',name1,name2,
c     *      imin,list1(imin),(list2(k),k=1,nlines)
          else if(name1.eq.name2) then
c
c Found more than one candidate: kill the latter and update the index vector
c
            do 3 k=j,nlist-1
   3        list1(k)=list1(k+1)
            nlist=nlist-1
            if(nlines.gt.0) then
              do 4 k=1,nlines
              if(list2(k).eq.j) list2(k)=imin
              if(list2(k).gt.j) list2(k)=list2(k)-1
   4          continue
            endif
          end if
          go to 2
        end if
c
c Put entries in the correct order and update the index vector
c
        name=list1(i)
c        if(list1(list2(4)).eq.'e-') write(*,*) 'C',name,
c     *    list1(imin),imin,list1(imin),(list2(k),k=1,nlines)
        list1(i)=list1(imin)
        list1(imin)=name
        if(nlines.gt.0) then
          do 5 k=1,nlines
          l=list2(k)
          if(l.eq.i)    list2(k)=imin
          if(l.eq.imin) list2(k)=i
   5      continue
        endif
        go to 1
      end if
c
      return
      end

C=========================================================================
C MBUILD: Build complete name from charge value and neutral species name.
C
C Inputs:
C   SPNAME [character] Name of neutral atom or molecule,
C   ICHARGE [integer] Desired charge value (-1, 0, 1 - 4) for output
C   atomic or molecular species. The charge value is interpreted as follows:
C       -1:  negative ion
C        0:  neutral species
C       +1:  singly ionized species
C       +2:  doubly ionized species, etc.
C
C     All other charge values are invalid and generate fatal errors.
C
C Outputs:
C   CHNAME [character] Complete name of species constructed from input charge
C     value and neutral species name.
C
C 96-Jun-01 Valenti  Wrote.
C 96-Dec-12 Piskunov Expanded to IONSIZ ionization stage
C
      subroutine mbuild(spname,icharge,chname)
      INCLUDE 'SIZES.EOS'

      character*(*) spname,chname
C
C Generate a fatal error if the neutral species begins with a space.
C
      if(spname(1:1).eq.' ') then
        write(*,*) 'mbuild: species name is blank'
        stop
      endif
C
C Check that requested charge value is allowed.
C
      if(icharge.lt.-1 .or. icharge.gt.IONSIZ-1) then
        write(*,200) spname,icharge
 200    format('mbuild: invalid charge value for ',a,':',i4)
        stop
      endif
C
C Initialize the output string with spaces.
C
      chname=' '
C
C Handle the simple case where a neutral charge state was requested.
C Just copy the input neutral species name up to the first space or
C   until SPCHAR characters have been copied.
C
      if(icharge.eq.0) then
        chname=spname
        return
      endif
C
C Find location of the first space, which is where the charge string will go.
C A fatal error occurs if the output requires more than SPCHAR characters.
C
      ispace=index(spname,' ')
      if(ispace.le.0.or.ispace+abs(icharge)-1.gt.len(chname)) then
        write(*,201) spname,icharge
 201    format('mbuild: no room in string "',a,'" for charge:',i4)
        stop
      end if
C
C Copy neutral species name.
C
      chname=spname
C
C Insert charge string beginning at first space.
C
      if(icharge.lt.0) then
        chname(ispace:ispace)='-'
      else if(icharge.gt.0.and.icharge.lt.IONSIZ) then
        chname(ispace:ispace+icharge-1)='++++++++++++++++++++++++++++++'
      else
        write(*,*) 'The charge is too large. Must be less than',IONSIZ,
     *             spname,icharge
        stop
      endif
C
c      write(*,*) icharge,'"',chname,'"'
      return
      end

C=========================================================================
C MPARSE: Parse molecular name. Get number and type of atomic constituents.
C
C Inputs:
C   SPNAME [character array(*)] Case-sensitive species name of molecule.
C     First letter of each atom name must be uppercase. The second letter
C     of each atom name, if present, must be lowercase. Each atom name may
C     optionally be followed by a multiplicity number between 1 and 4. If
C     no multiplicity number is given for a particular atom, then its
C     multiplicity is assumed to be 1. Finally, a non-neutral charge state
C     for the molecule may be specified with a trailing "-", "+", or "++".
C     In the absence of such a charge indicator, the molecule is assumed
C     to be neutral.
C   ELEMEN [character array(*)] Case-sensitive list of atoms participating
C     in molecule formation (periodic table).
C
C Outputs:
C   NEL [integer] Number of elements comprising molecule. Also gives the
C     maximum valid index for IEL and NAT.
C   CHARGE [integer] Charge state of the molecule (-1, 0, +1,...,+(IONSIZ-1)).
C   IEL [integer array(4)] atomic number(s) of the atomic types comprising
C     the molecule in SPNAME.
C   NAT [integer array(4)] multiplicity (up to 4) for each of the atomic
C     types in IEL.
C
      SUBROUTINE MPARSE(ELEMEN,SPNAME,NEL,CHARGE,IEL,NAT,ELESIZ)
      INCLUDE 'SIZES.EOS'
C
      INTEGER IEL(4),NAT(4),NEL,CHARGE,ELESIZ
      CHARACTER SPNAME*(SPCHAR),TMP*2
      CHARACTER*(3) ELEMEN(ELESIZ)
C
C  Set pointer I1 to beginning of first atom name.
C
c      write(*,*) LEN(ELEMEN(1))
      CHARGE=0
      I1=1
C
C  Loop through (up to four) different atoms in a molecule.
C
      DO 4 J=1,4
C
C  Set pointer I2 to the end of the next atom's name.
C
      I2=I1
      IF(ICHAR(SPNAME(I1+1:I1+1)).GE.ICHAR('a').AND.
     *   ICHAR(SPNAME(I1+1:I1+1)).LE.ICHAR('z')) I2=I1+1
C
C  Update number of atomic species in molecule.
C
      NEL=J
C
C  Find atomic the atomic number of current atom.
C
      TMP='  '
      TMP=SPNAME(I1:I2)
      DO 1 I=1,ELESIZ
      IF(TMP.EQ.ELEMEN(I)(1:2)) GO TO 2
   1  CONTINUE
C
C  Fall through to here if atom name was not in ELEMEN list.
C
c      WRITE(*,*) 'Unknown element: ',SPNAME,i1,i2,' ',SPNAME(i1:i2)
      WRITE(*,*) 'Unknown element: ',SPNAME(I1:I2),' "',SPNAME(1:I2),'"'
      STOP
C
C  Save atomic number of current atom.
C
   2  IEL(NEL)=I
C
C  Check for optional atomic multiplicity. Default is 1; maximum is 5.
C
      I1=I2+1
      NAT(NEL)=1
      IF(SPNAME(I1:I1).EQ.'1') THEN
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'2') THEN
        NAT(NEL)=2
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'3') THEN
        NAT(NEL)=3
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'4') THEN
        NAT(NEL)=4
        I1=I1+1
      ELSE IF(SPNAME(I1:I1).EQ.'5') THEN
        NAT(NEL)=5
        I1=I1+1
      END IF
C
C   Check for optional charge on molecule. Default is neutral; "-", "+",
C   "++", etc. up to IONSIZ are allowed.
C
      IF(I1.GT.SPCHAR) RETURN
      IF(SPNAME(I1:I1).EQ.' ') RETURN
      IF(SPNAME(I1:I1).EQ.'-') THEN
        CHARGE=-1
        RETURN
      ENDIF
      IF(SPNAME(I1:I1).EQ.'+') THEN
        CHARGE=1
        DO 3 IONN=1,IONSIZ-1
        IF(SPNAME(I1+IONN:I1+IONN).NE.'+') RETURN
   3    CHARGE=CHARGE+1
      END IF
C
C  Fall through if we didn't just find a charge state and return. Loop
C  back and interpret character pointed at by I1 as beginning of atom.
C
   4  CONTINUE
C
C  There were 4 different atomic types, but presumably we are done.
C
      RETURN
      END

C=========================================================================
C GAS: Determines the equilibrium ionization and partial pressure for every
C      atom and molecule in the species list, assuming no other species are
C      present. Temperature, total pressure, and elemental abundances must
C      be specified, but all atomic and molecular constants are determined
C      internally.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   XNELEC [real] Estimated electron number density (in 1/cm^3)
C   XNATOM [real] Number density (in 1/cm^3) of all particles other than
C     electrons (i.e. atoms or molecules), used to calculate total pressure?
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPLIST [character*(*) array(NLIST)] List of species to consider in
C     solving for the molecular equilibrium, including both the base set,
C     plus any additional species that occur in the line list.
C   NLIST [integer] the number of valid species in SPLIST.
C   TOL [real] iterative solution of the population equations terminates
C     successfully when the largest fractional change in pressure for all
C     species (including electrons) drops below this "tolerance".
C Outputs:
C   XNE [real] electron number density (in 1/cm^3) determined in GAS routine.
C   AWT [real*8] atomic weights of each species
C
      SUBROUTINE GAS(TEMP,Pgas,Pelec,ABUND,ELEMEN,AMASS,ELESIZ,
     *               TOL,SPLIST,NLIST,XNE,XNA,RHO,Pgnew,
     *               XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
     *               FAILED)
c      SUBROUTINE GAS(TEMP,XNELEC,XNATOM,ABUND,ELEMEN,AMASS,ELESIZ,
c     *               TOL,SPLIST,NLIST,
c     *               XNE,XNA,RHO,XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
c     *               FAILED)

      IMPLICIT NONE
      INCLUDE 'SIZES.EOS'
C
      CHARACTER ENAME*(SPCHAR),BLANK*1
      INTEGER MAXIT,MAXREF
      DOUBLE PRECISION KBOL,HMASS,AMULOG
      PARAMETER (BLANK=' ',ENAME='e-',KBOL=1.38065D-16,MAXIT=1000,
     *           HMASS=1.66053D-24,AMULOG=-23.779751D0,MAXREF=10)
      LOGICAL PRINT,FAILED

      INTEGER NLIST,ELESIZ
      CHARACTER*(SPCHAR) SPLIST(NLIST)
      CHARACTER*(3) ELEMEN(ELESIZ)
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      CHARACTER NAMEMX*(SPCHAR),NAMET*(SPCHAR)
      INTEGER JATOM, TYPE(SPLSIZ-1),NCH(SPLSIZ-1),IATOM(ELEDIM),
     *  INDSP(ELEDIM),NAT(4,SPLSIZ-1),ZAT(4,SPLSIZ-1),NTOT(SPLSIZ-1),
     *  NEL(SPLSIZ-1),IAT(SPLSIZ-1),INDZAT(99)
      REAL T,TEMP,XNELEC,XNATOM,TOL,XNE,XNA,RHO,Pgas,Pelec,Pgnew,
     *  POTI(IONSIZ),XNPF(*),PFUNC(*),POTION(*),XTOTAL
      DOUBLE PRECISION FRACT(IONSIZ),IT(SPLSIZ-1),KT(SPLSIZ-1),
     *  AWT(SPLSIZ-1)

      DOUBLE PRECISION A(ELEDIM+1,ELEDIM+1),RHS(ELEDIM+1),
     *  AA(ELEDIM+1,ELEDIM+1),
     *  B(ELEDIM+1),BB(ELEDIM+1),
     *  P(ELEDIM+1),PP(SPLSIZ-1),PP0(SPLSIZ-1),PART(SPLSIZ-1),ND

      DOUBLE PRECISION PE,PG,PF,PNEW,PENEW,DP,DPE,PION,PENQ
      DOUBLE PRECISION RNF(ELEDIM),AL(ELEDIM+1)
      INTEGER NELM,NCHG,ANUM(4),NATM(4),IPIV(ELEDIM+1),IWORK(ELEDIM+1),
     *  INFO,REPEAT,ISPEC,NSP1,NELT,NQ,K,KK,IDIR,KMAX,I,J,NEQ,IELM,NP,
     *  IIH2,IICO,IIH2O,IDAMAX,NGIT
      DOUBLE PRECISION RATIOM,QPRD,RHSTOT,SCALE,FACTOR,PNOTE,PDTOT,PU,
     *  PD,GMU,PTOT,DELP,DELPE,PQ,RCOND,DASUM,DELMAX,PE0,
     *  PTOTH,PHyd,PTOTC,PTOTO,WATCOR,AQUAD,BQUAD,CQUAD,DPQ,DPTOT
c      DOUBLE PRECISION PZS,COMPZ

      DOUBLE PRECISION RSCL(ELEDIM+1),CSCL(ELEDIM+1)
      DOUBLE PRECISION FERR(1),BERR(1),WORK(5*(ELEDIM+1))
      CHARACTER*1 EQUED

cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      real ttt(101)
c      real*8 Kttt(101)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C
C Initialize the Reciprocal Neutral Fraction (RNF). The RNF is used to
C adjust the initial neutral atomic partial pressures used in the linear
C solver. Originally, atomic species were assumed to be predominantly
C neutral, but at low electron pressures, this is a poor assumption for
C species with low ionization potentials.
C
      DO 1 I=1,ELEDIM
   1  RNF(I)=1.0D0
C
C Total gas and electron pressure
C
c      T=MAX(1200.,TEMP)
      T=TEMP
      PG=Pgas
      PE=Pelec
      XNELEC=PE/(KBOL*TEMP)
      XNATOM=PG/(KBOL*TEMP)
C
C Avoid unpleasant surprises
C
      IF(PG.GT.PE) THEN
        XNATOM=XNATOM-XNELEC
      ELSE
        XNELEC=XNATOM*0.01
      END IF
c      PG=(XNATOM+XNELEC)*KBOL*TEMP
c      PE=XNELEC*KBOL*TEMP
C
C  Calculate equilibrium constants for each species in list (except 'e-').
C
c      PRINT=.TRUE.
      PRINT=.FALSE.
      PION=0
      IIH2=0
      IICO=0
      IIH2O=0
      JATOM=0
      NP=0
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      open(13,file='KT_eos.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')
c      write(13) NLIST,LEN(SPLIST(1))
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO 4 ISPEC=1,NLIST-1
      PP0(ISPEC)=0.D0
      CALL MPARSE(ELEMEN,SPLIST(ISPEC),NELM,NCHG,ANUM,NATM,ELESIZ)
      IF(NCHG.EQ.0) NP=ISPEC
      IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.EQ.0) THEN
C
C  Neutral atom
C
        TYPE(ISPEC)=1
        KT(ISPEC)=1.0
        IT(ISPEC)=1.0
        JATOM=JATOM+1
        IF(JATOM.GT.ELEDIM) THEN
          write(*,'(a,2i4)') 'gas: too many element types,' //
     *      ' increase ELEDIM:',ELEDIM,JATOM
          stop
        END IF
        IATOM(JATOM)=ANUM(1)
        INDSP(JATOM)=ISPEC
        IAT(ISPEC)=JATOM
        AWT(ISPEC)=AMASS(ANUM(1))
        INDZAT(ANUM(1))=JATOM
        NTOT(ISPEC)=1
        CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,
     *             IONSIZ,3)
        PART(ISPEC)=FRACT(1)
        POTION(ISPEC)=POTI(1)
      ELSE IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.NE.0) THEN
C
C  Ionized atom
C
        TYPE(ISPEC)=3
        IF(NCHG.GT.0) THEN
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,
     *               IONSIZ,2)
          IT(ISPEC)=FRACT(NCHG+1)/FRACT(1)*PE**NCHG
          RNF(ANUM(1))=RNF(ANUM(1))+FRACT(NCHG+1)/FRACT(1)
c          if(ANUM(1).eq.26) write(*,*) SPLIST(ISPEC),NCHG,
c     *                      (FRACT(I),I=1,IONSIZ)
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,
     *               IONSIZ,3)
          PART(ISPEC)=FRACT(NCHG+1)
c          if(ANUM(1).eq.62) write(*,*) 'pf: ',SPLIST(ISPEC),NCHG,FRACT
          POTION(ISPEC)=POTI(NCHG+1)
C
C Negative ions
C
        ELSE IF(NCHG.LT.0) THEN
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,
     *               IONSIZ,3)
          QPRD=2.D0*FRACT(1)
          CALL NEGION(ANUM(1),T,QPRD,
     *                IT(ISPEC),PART(ISPEC),POTION(ISPEC))
          KT(ISPEC)=1.
        END IF
C
        KT(ISPEC)=1.0
        AWT(ISPEC)=AMASS(ANUM(1))
        NTOT(ISPEC)=1
      ELSE IF(NELM.GT.1.OR.NATM(1).GT.1) THEN
C
C  Neutral or ionized molecule
C
        TYPE(ISPEC)=2
C
C  Calculate mass ratio (RATIOM) and partition function product (QPRD)
C  needed by MOLCON. See MOLCON header for decription of these quantities.
C  While we are at it, calculate the atomic weight (AWT) of the molecule
C  and the total number of atoms (NTOT) of any type in the molecule.
C
        NTOT(ISPEC)=0
        AWT(ISPEC)=0.0D0
        RATIOM=0.0D0
C
C  Fixed the partition function ratio for ionized molecules.
C  Now we start with a product of partition functions for free
C  electrons in denominator. NP 29-12-2006.
C        QPRD=0.0D0
        QPRD=-NCHG*LOG10(2.0)
        DO 2 IELM=1,NELM
        NTOT(ISPEC)=NTOT(ISPEC)+NATM(IELM)
        AWT(ISPEC)=AWT(ISPEC)+NATM(IELM)*AMASS(ANUM(IELM))
        RATIOM=RATIOM+NATM(IELM)*LOG10(AMASS(ANUM(IELM)))
        CALL XSAHA(ANUM(IELM),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,
     *             IONSIZ,3)
        IF(SPLIST(ISPEC).EQ.'H2')  IIH2=ISPEC
        IF(SPLIST(ISPEC).EQ.'CO')  IICO=ISPEC
        IF(SPLIST(ISPEC).EQ.'H2O') IIH2O=ISPEC
c       if(splist(ispec).eq.'N2')write(*,*)
c     &    anum(ielm),(fract(i),i=1,2)
   2    QPRD=QPRD+NATM(IELM)*LOG10(FRACT(1))
        RATIOM=RATIOM-LOG10(AWT(ISPEC))+(NTOT(ISPEC)-1)*AMULOG
C
C  Now get the molecular constants from MOLCON.
C
        CALL MOLCON(SPLIST(ISPEC),T,NTOT(ISPEC),
     &              RATIOM,QPRD,KT(ISPEC),PART(ISPEC),PION)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        do ittt=0,100
c          ttt(ittt+1)=20.*ittt+1000.
c          CALL MOLCON(SPLIST(ISPEC),ttt(ittt+1),NTOT(ISPEC),
c     &                RATIOM,QPRD,Kttt(ittt+1),PART(ISPEC),PION)
c        enddo
c        write(13) SPLIST(ispec),ttt,Kttt
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  Finally, record the charge state of the molecule.
C
        IF(NCHG.EQ.0) THEN
          IT(ISPEC)=1
        ELSE
          write(*,*) T,ISPEC,SPLIST(ISPEC),KT(ISPEC),IT(ISPEC)
C
C  The first option was used with Sauval & Tatum constants.
C  JV fits to NextGen pressures needed IT(ISPEC)=1.0 for positive
C  molecular ions.
C
          IF(SPLIST(ISPEC).EQ.'H2+'.OR.SPLIST(ISPEC).EQ.'NO+') THEN
            K=1
            DO IELM=2,NELM
              IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
     *          K=IELM
            ENDDO
            IT(ISPEC)=IT(INDSP(ANUM(K))+1)
            KT(ISPEC)=KT(ISPEC)/IT(ISPEC)
          ENDIF
          IT(ISPEC)=1.0
        END IF
C
C  Store ionization potential (needed e.g. for broadening calculations)
C
        IF(PION.GT.0.D0) THEN
          POTION(ISPEC)=PION
        ELSE
c
c  If ionization potential is not available use the one for TiO!
c
          POTION(ISPEC)=6.4
        ENDIF
      ELSE
C
C  Fall through to here when the molecular formula doesn't make sense.
C
        WRITE(*,*) 'Wrong formula for the species: ',splist(ISPEC)
        STOP
      END IF
C
C  Now save results of MPARSE into arrays.
C
      NEL(ISPEC)=NELM
      NCH(ISPEC)=NCHG
      DO 3 IELM=1,NELM
      ZAT(IELM,ISPEC)=ANUM(IELM)
   3  NAT(IELM,ISPEC)=NATM(IELM)
C
C  Go back for next species.
C
   4  CONTINUE
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      close(13)
c      stop
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      NEQ=JATOM+1
C==================================
C== End of species list parsing. ==
C==================================
C
C Print diagnostic: neutral fractions.
C
c     write(*,*) 'Reciprocal Neutral Fractions'
c     do 850 i=1,JATOM/7
c       write(*,860) (jeff(iatom(j)),j=7*i-6,7*i)
c850  continue
c860  format(1p7e10.3,a)
c     if(JATOM.gt.7*(JATOM/7)) write(*,860)
c    *  (jeff(iatom(j)),j=7*(JATOM/7)+1,JATOM)
c      do 52 i=1,nlist-1
c  52  write(*,'(I4,1P2E12.4,3I3,A6,0Pf8.2,8I4)')
c     *  i,IT(i),KT(i),NCH(i),NTOT(i),NEL(i),SPLIST(i),AWT(i),
c     *  (ZAT(j,i),NAT(j,i),j=1,NEL(i))
C================================================================
C== UPDATE MAIN ARRAYS                                         ==
C================================================================
c
c Make the initial estimate of the partial pressures for neutral atoms. These
c pressures are used as input to the linear solver. When only abundances are
c considered, the largest errors occur for low ionization elements, which can
c be highly ionized at low electron pressures. Thus, we apply a correction
c to recover the neutral fraction for each atom. The neutral fraction only
c corrects for losses into ionization states included in the species list.
c When the ionization correction is included, the largest error in the inital
c guess for carbon, which has unaccounted for losses into CO. Late in the
c convergence process, nitrogen becomes the dominant source of error.
c
      DO 5 J=1,JATOM
      P(J)=PG*ABUND(IATOM(J))/RNF(IATOM(J))
      ISPEC=INDSP(J)
      PP0(ISPEC)=P(J)
   5  CONTINUE
c
c Make an initial guess at the balance between H and H2.
c Assumes pressures of species other than H, H2, He, and Ne are negligible.
c Constraints:
c   KT(IIH2)*PP(IIH2)=P(1)**2           <-- chemical equilibrium
c   P(1)+2*PP(IIH2)=ABUND(1)*(PG-PE)    <-- H particle conservation
c
      IF(IIH2.GT.0) THEN
        PHyd=0.5*(-KT(IIH2)+SQRT(KT(IIH2)**2
     &        +4.0*KT(IIH2)*(PG-PE-P(2)-P(10))))
      ELSE
        PHyd=(PG-PE)*ABUND(1)
      ENDIF
c      IF(PHyd.GT.0.) P(1)=PHyd
c
c Make an initial guess at the balance between C, O, CO, and H2O.
c Constraints:
c   KT(IICO)*PP(IICO)=P(6)*P(8)         <-- chemical equilibrium
c   KT(IIH2O)*PP(IIH2O)=P(1)**2*P(8)    <-- chemical equilibrium
c   PTOTH=P(1)+2*PP(IIH2)       <-- defines density of H nuclei
c   PTOTC=P(6)+PP(IICO)                 <-- defines density of C nuclei
c   PTOTO=P(8)+PP(IICO)+PP(IIH2O)       <-- defines density of O nuclei
c   PTOTC=PTOTH*ABUND(6)/ABUND(1)       <-- abundance constraint
c   PTOTO=PTOTH*ABUND(8)/ABUND(1)       <-- abundance constraint
c
      PTOTH=P(1)
      IF(IIH2.GT.0) PTOTH=PTOTH+2.0*P(1)**2/KT(IIH2)
      PTOTC=PTOTH*ABUND(6)/ABUND(1)
      PTOTO=PTOTH*ABUND(8)/ABUND(1)
      IF(IIH2O.GT.0) THEN
        WATCOR=1.0+P(1)**2/KT(IIH2O)
        AQUAD=1.0/WATCOR
        IF(IICO.GT.0) THEN
          BQUAD=KT(IICO)+(PTOTO-PTOTC)/WATCOR
          CQUAD=-KT(IICO)*PTOTC
c          P(6)=(-BQUAD+SQRT(BQUAD**2-4.0*AQUAD*CQUAD))/(2.0*AQUAD)
c          P(8)=(P(6)+PTOTO-PTOTC)/WATCOR
        ELSE
c          P(6)=PTOTC
c          P(8)=PTOTO
        ENDIF
      ELSE
c        P(6)=PTOTC
c        P(8)=PTOTO
      ENDIF
c      IF(P(6).LE.0.) P(6)=PTOTC
c      IF(P(8).LE.0.) P(8)=PTOTO
      PE0=PE
      NAMEMX=BLANK
      DELMAX=0.0D0
c      COMPZ=0.0D0
c      PZS=0.0D0
c      write(*,*) SPLIST(1),P(1),SPLIST(IIH2),P(IIH2),
c     *           SPLIST(IIH2+1),P(IIH2+1),
c     *           SPLIST(IIH2+2),P(IIH2+2)
c      DO 6 J=1,JATOM
c      NN=INDSP(J)
c      IF(IPR(NN).NE.2) GOTO 3
c      NNP=INDX(3,ITAB(ZAT(1,NN)),1,1,1)
c      COMPZ=COMPZ+ABUND(IATOM(J))
c      IF(PE.EQ.0.0D0) PZS= PZS + P(J)
c      IF(PE.GT.0.0D0) PZS= PZS + (1.0D0+IT(NNP)/PE)*P(J)
c   6  CONTINUE
c      do J=1,JATOM
c        write(*,*) J,P(J),ABUND(IATOM(J)),SPLIST(INDSP(J))
c      enddo
c      write(*,*) JATOM+1,PE,'e-'
c      stop
C================================================================
C== MAIN LOOP: FILL LINEARIZED COEFFICIENT MATRIX AND RHS VECTOR,
C== AND SOLVE SYSTEM FOR PARTIAL PRESSURE CORRECTIONS.         ==
C== ISOLV=1: LINEARIZE ONLY THE PARTIAL PRESSURES OF THE NEUTRAL=
C== ATOMS FOR WHICH IPR(J)=1 (MAJOR SPECIES). THE ELECTRON     ==
C== PRESSURE PE IS ASSUMED TO BE GIVEN IN THIS CASE, AND SO IS ==
C== NOT INCLUDED IN THE LINEARIZATION. THIS IS NECESSARY SINCE ==
C== MOST OF THESE ELECTRONS (AT COOL TEMPS.) ORIGINATE FROM    ==
C== ELEMENTS NOT CONSIDERED IN THE LINEARIZATION. IN ORDER TO  ==
C== OBTAIN A GOOD VALUE FOR PE IN THE FIRST PLACE, IT IS       ==
C== NECESSARY TO CALL GAS WITH ISOLV=2.                        ==
C== ISOLV=2: THIS LINEARIZES THE PARTIAL PRESSURES OF THE NEUTRAL
C== ATOMS FOR WHICH IPR(J)=1 OR 2. THIS LIST OF ELEMENTS SHOULD==
C== INCLUDE ALL THE SIGNIFICANT CONTRIBUTORS TO THE TOTAL      ==
C== PRESSURE PG, AS WELL AS THE ELECTON PRESSURE PE. ANY ELEMENT=
C== (IPR(J)=3) NOT INCLUDED IS ASSUMED TO HAVE A NEGLIGIBLE    ==
C== EFFECT ON BOTH P AND PE.                                   ==
C== IN BOTH CASES, THE PARTIAL PRESSURES OF THE NEUTRAL ATOMS  ==
C== FOR ELEMENTS NOT INCLUDED IN THE LINEARIZATION ARE         ==
C== CALCULATED DIRECTLY FROM THE NOW DETERMINED PRESSURES OF   ==
C== THE LINEARIZED ELEMENTS.                                   ==
C================================================================
      NGIT=0
      RHSTOT=1.D99
C
C Top of loop in which linearized equations are solved recursively.
C
      REPEAT=0
      KMAX=1
   7  IF(NGIT.GE.MAXIT) THEN
        WRITE(6,208)
 208    FORMAT('*** ERROR: TOO MANY ITERATIONS IN ROUTINE "GAS"')
        WRITE(6,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),RHSTOT
        write(*,*) TEMP,PG,P(1),XNATOM,XNELEC
        STOP
      END IF
      NGIT=NGIT+1
      P(NEQ)=PE

      SCALE=10.D0
      IDIR=0
   9  CALL EOSFCN(NEQ,P,B,A,1,PG,NCH,NLIST,
     *  IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)

c      write(*,*) 'Pe,SCALE,B(1),Ptot,Pg=',PE,SCALE,B(1),PTOT,PG

      IF(B(1).GT.1.D2) THEN
        IF(IDIR.NE.-1) THEN
          SCALE=SQRT(SCALE)
          IDIR=-1
        ENDIF
C
C Neutral atomic pressures are too high. Scale them down until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)/SCALE
        ENDDO
        GOTO 9
      ELSE IF(B(1).LT.-1.D2) THEN
        IF(IDIR.NE.1) THEN
          SCALE=SQRT(SCALE)
          IDIR=1
        ENDIF
C
C Neutral atomic pressures are too low. Scale them up until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)*SCALE
        ENDDO
        GOTO 9
      ENDIF

      CALL EOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,
     *  IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
C================================================================
C== NOW SOLVE THE LINEARIZED EQUATIONS (USING ROUTINE "LINEQ") ==
C================================================================
      IF(PRINT) THEN
        WRITE(*,200) NGIT
 200    FORMAT('LOG OF COEFFICIENT MATRIX AT ITERATION #',I5//)
        KK=MIN(30,NEQ-1)
        WRITE(*,201) (SPLIST(INDSP(K)),K=1,KK-1),'e-','RHS'
 201    FORMAT(4x,31(1x,a3,2x))
        DO 21 I=1,KK-1
        DO 20 J=1,KK-1
  20    AL(J)=LOG10(ABS(A(J,I))+1.0D-50)
        AL(KK)=LOG10(ABS(A(NEQ,I))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(I))+1.0D-50)
        NAMET=SPLIST(INDSP(I))
        WRITE(*,202) NAMET,(AL(J),J=1,KK+1)
  21    CONTINUE
        DO 22 J=1,KK-1
        AL(J)=LOG10(ABS(A(J,NEQ))+1.0D-50)
  22    CONTINUE
        AL(KK)=LOG10(ABS(A(NEQ,NEQ))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(NEQ))+1.0D-50)
        NAMET='e-'
        WRITE(*,202) NAMET,(AL(J),J=1,KK+1)
 202    FORMAT(A2,31F6.1)
        WRITE(*,'(/)')
      END IF
C
C  Save a copy of the RHS for future step refinement
C
      DO 23 I=1,NEQ
  23  RHS(I)=B(I)
      RHSTOT=DASUM(NEQ,RHS,1)
C
C  Solve linear system for corrections
C  In order not to solve for Pelect, one should use NEQ-1 as the first
C  argument. NEQ solves the whole system including electron pressure
C
c
c  Using LAPACK routine
c
c        open(unit=4,file='dump.bin',form='UNFORMATTED')
c        write(4) NEQ
c        write(4) ((A(i,j),i=1,NEQ),j=1,NEQ)
c        write(4) (B(i),i=1,NEQ)
      CALL DGESVX('E','N',NEQ,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,EQUED,
     *            RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,FERR,BERR,
     *            WORK,IWORK,INFO)
      CALL DCOPY(NEQ,BB,1,B,1)
c
c  The same thing using LINEQ2 or LINEQ and BLAS 2/3
c      CALL LINEQ(NEQ,1,A,ELEDIM+1,IPIV,B,ELEDIM+1,INFO)

      IF(INFO.NE.0) THEN
        IF(REPEAT.LT.2) THEN
          DO J=1,NEQ
           P(J)=P(J)*0.999D0
          END DO
          REPEAT=REPEAT+1
          GO TO 9
        ELSE IF(REPEAT.LT.4) THEN
          DO J=1,NEQ
           P(J)=P(J)*1.001D0
          END DO
          REPEAT=REPEAT+1
          GO TO 9
        ELSE
c          WRITE(*,*) 'EOS: LINEQ failed to solved for corrections to'
c          WRITE(*,*) '     the partial pressures. Matrix is degenerate'
c          WRITE(*,*) '     Temp=',TEMP,', Natom=',XNATOM,', Ne=',XNELEC
c          WRITE(*,*) '     INFO=',INFO,' Iter=',NGIT,' EQUED=',EQUED
cc          open(unit=4,file='dump.bin',form='UNFORMATTED')
cc          write(4) NEQ,((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
cc          close(4)
cc          write(1) 0
cc          close(1)
c          STOP
          CALL DGESVX('E','N',NEQ-1,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,EQUED,
     *                RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,FERR,BERR,
     *                WORK,IWORK,INFO)
          CALL DCOPY(NEQ,BB,1,B,1)
          PTOT=0.D0
          DO J=1,NEQ-1
            PTOT=PTOT+P(J)
          END DO
          PE=MAX(PG-PTOT,1.D-20)
        END IF
      END IF
      REPEAT=0

c
C=================================================================
C== FINALLY, UPDATE THE PARTIAL PRESSURES FOR THE MAJOR SPECIES ==
C== BY ADDING THE PRESSURE CORRECTIONS OBTAINED FOR EACH ATOM   ==
C== FROM THE LINEARIZATION PROCEDURE.                           ==
C=================================================================
      DELMAX=-1.0D0
      KMAX=1
      DO 31 K=1,NEQ
      ISPEC=INDSP(K)
C
C Compute the maximum correction in order to computer the under-relaxation factor
C
      DP=B(K)
      IF(P(K).GT.1.0D-20.AND.ABS(P(K)/PG).GE.1.0D-15) THEN
        DELP=ABS(DP/P(K))
        IF(DELP.GT.DELMAX) THEN
          DELMAX=DELP
        END IF
      END IF
  31  CONTINUE
C
C  Under-relaxation factor
C
      FACTOR=1.D0/(DELMAX+1.D0)
C
C Apply corrections
C
      DELMAX=-1.0D0
      KMAX=1
      DO 32 K=1,JATOM
      ISPEC=INDSP(K)
C
C  Restrict the correction to avoid getting negative pressures
C
      PNEW=P(K)-B(K)*FACTOR
c      write(*,*) K,P(K),PNEW,B(K)
      IF(PNEW.LT.0.D0) PNEW=MIN(MIN(P(K),ABS(PNEW)),PG)
c      IF(PNEW.LT.0.D0) PNEW=ABS(PNEW)
      DP=PNEW-P(K)
      IF(ABS(DP).GT.1.D-15) DP=DP*MIN(1.D0,0.4D0*P(K)/ABS(DP))
      P(K)=PNEW
      IF(P(K).GT.1.0D-20.AND.ABS(P(K)/PG).GE.1.0D-15) THEN
        DELP=ABS(DP/P(K))
        IF(DELP.GT.DELMAX) THEN
          NAMEMX=SPLIST(ISPEC)
          DELMAX=DELP
          KMAX=K
        END IF
      END IF
  32  CONTINUE

c      PENEW=BBB(NEQ)
      PENEW=PE-B(NEQ)*FACTOR
c      write(*,*) NEQ,PE,PENEW,B(NEQ)
      IF(PENEW.LT.0.D0) PENEW=MIN(PE,ABS(PENEW))
c      IF(PENEW.LT.0.D0) PENEW=ABS(PENEW)
      DPE=PENEW-PE
      IF(ABS(DPE).GT.1.D-15) DPE=DPE*MIN(1.D0,0.4D0*PE/ABS(DPE))
      PE=PENEW
      IF(ABS(PE/PG).GE.1.0D-15) THEN
        DELPE=ABS(DPE/PE)
        IF(DELPE.GT.DELMAX) NAMEMX=ENAME
        IF(DELPE.GT.DELMAX) DELMAX=DELPE
      END IF
C================================================================
C== PRINT OUT SUMMARY LINE FOR EACH ITERATION                  ==
C================================================================
      PTOT=PE
      PQ=0.0D0
      DO ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=LOG(MAX(IT(ISPEC),1.D-115))-LOG(KT(ISPEC))-
     -     LOG(MAX(PE,1.D-115))*NQ
        DO I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF+LOG(MAX(P(J),1.D-115))*NAT(I,ISPEC)
        ENDDO
c        PENQ=1.0D0
c        IF(PE.GT.0.0D0.AND.NQ.NE.0) PENQ=PE**NQ
c        PP(ISPEC)=IT(ISPEC)/(KT(ISPEC)*PENQ)*PF
        PP(ISPEC)=EXP(PF)
        PTOT=PTOT+PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
c        write(*,*) ISPEC,SPLIST(ISPEC),PP(ISPEC),PTOT,PG
      ENDDO
c      stop
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(PE-PQ)/PG
c      write(*,*) PG,PTOT,DELMAX,DPTOT,DPQ,FACTOR
      IF(PRINT) THEN
        WRITE(6,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),
     *               PTOT/TEMP/KBOL,DPTOT,PE/TEMP/KBOL,DPQ
 203    FORMAT(I10,2X,A8,1P9E11.3)
      END IF
      IF((DPTOT.GT.TOL.OR.DPQ.GT.TOL.OR.DELMAX.GT.TOL)
     *   .AND.NGIT.LT.MAXIT) GOTO 7
C
C Bottom of the loop in which linearized equations are solved recursively.
C
C================================================================
C== CALCULATE FINAL PARTIAL PRESSURES AFTER CONVERGENCE OBTAINED=
C================================================================
      PTOT=PE
      PD=0.0D0
      PU=0.0D0
      PQ=0.0D0
      DO 34 ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=1.0D0
        DO 33 I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF*P(J)**NAT(I,ISPEC)
  33    CONTINUE
        PENQ=1.0D0
        IF(PE.GT.0.0D0) PENQ=PE**NQ
        PP(ISPEC)=IT(ISPEC)/(KT(ISPEC)*PENQ)*PF
        PTOT=PTOT+PP(ISPEC)
        PD=PD+NTOT(ISPEC)*PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
        PU=PU+AWT(ISPEC)*PP(ISPEC)
  34  CONTINUE
      PP(NLIST)=PE
      PDTOT=PD+PE
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(PQ-PE)/PG
      GMU=PU/PTOT
      ND=PTOT/(TEMP*KBOL)
      RHO=ND*GMU*HMASS
      XNE=PE/(TEMP*KBOL)
C================================================================
C== WRITE OUT FINAL PARTIAL PRESSURES                          ==
C================================================================
      IF(PRINT) THEN
c      IF(DASUM(NLIST-1,PP,1)+PE.GT.PG*1.01D0) THEN
        write(*,'(''AFTER '',I3,'' iterations.   Max change of:'',G10.3,
     #      ''  in element:'',A)') NGIT,DELMAX,NAMEMX
        WRITE(6,'(''AFTER '',I3,'' ITERATIONS WITH ''/
     #            ''T='',1PE10.3,''   P='',E10.3)') NGIT,TEMP,
     #    DASUM(NLIST-1,PP,1)+PE
        WRITE(6,'(''PDTOT='',1PE10.3,''   DPTOT='',E10.3,
     #            ''  DPQ='',E10.3,''  Nelectron='',E10.3,'' cm^3''/
     #    '' Nparticle='',1PE10.3,'' cm^3   Mean At.Wt.='',
     #    0PF7.3,''   Density='',1PE10.3,'' g/cm^3''//
     #    '' # Species   Abundance   Initial P   Final P'',
     #    ''      IT         KT         pf''/)')
     #    PDTOT,DPTOT,DPQ,XNE,ND-XNE,GMU,RHO
        NSP1=NLIST
        DO 35 ISPEC=1,NLIST-1
        IF(TYPE(ISPEC).NE.1) THEN
          WRITE(*,206) ISPEC,SPLIST(ISPEC),PP0(ISPEC),PP(ISPEC),
     #                 IT(ISPEC),KT(ISPEC),PART(ISPEC)
 206      FORMAT(I3,1X,A8,11X,1P5E11.3)
        ELSE
          J=IAT(ISPEC)
          WRITE(*,207) ISPEC,splist(ISPEC),ABUND(IATOM(J)),PP0(ISPEC),
     #                 PP(ISPEC),IT(ISPEC),KT(ISPEC),PART(ISPEC)
 207      FORMAT(I3,1X,A8,1P6E11.3)
        END IF
  35    CONTINUE
        WRITE(*,206) NSP1,ENAME,PE0,PE
        WRITE(*,*) IDAMAX(NLIST-1,PP,1),SPLIST(IDAMAX(NLIST-1,PP,1))
c        stop
      END IF
C
C Fill up the output array and set up flags
C PNOTE is the partial pressure due to everything except electrons.
C XNA is the number density of everything except electrons.
C
      PNOTE=0.D0
      DO 36 ISPEC=1,NLIST-1
      IF(PART(ISPEC).GT.0.) THEN
        IF(PP(ISPEC)/(KBOL*TEMP*PART(ISPEC)).GE.1.D-20) THEN
          XNPF(ISPEC)=PP(ISPEC)/(KBOL*TEMP*PART(ISPEC))
        ELSE
          XNPF(ISPEC)=0.0
        END IF
        PFUNC(ISPEC)=PART(ISPEC)
      ELSE
        XNPF(ISPEC)=0.
        PFUNC(ISPEC)=1.
      END IF
      PNOTE=PNOTE+PP(ISPEC)
c      write(*,*) ISPEC,PNOTE,PP(ISPEC),SPLIST(ISPEC)
  36  CONTINUE
      XNPF(NLIST)=XNE
      PFUNC(NLIST)=1.0
      XTOTAL=PD/(KBOL*TEMP)
      XNA=PNOTE/(KBOL*TEMP)
      Pgnew=PTOT
C
      RETURN
      END

C=========================================================================
C LOGARITHMIC version: the solution is found for the logs of ficticious
C                      partial pressures.
C GAS: Determines the equilibrium ionization and partial pressure for every
C      atom and molecule in the species list, assuming no other species are
C      present. Temperature, total pressure, and elemental abundances must
C      be specified, but all atomic and molecular constants are determined
C      internally.
C
C Inputs:
C   TEMP [real] Temperature (in K) which is to be used in calculating the
C     equilibrium constants and partition functions.
C   XNELEC [real] Estimated electron number density (in 1/cm^3)
C   XNATOM [real] Number density (in 1/cm^3) of all particles other than
C     electrons (i.e. atoms or molecules), used to calculate total pressure?
C   ABUND [real array(ELESIZ)] The fraction of all atomic species with respect
C     to the total number of atomic nuclei in any form. Thus, hydrogen has
C     an abundance slightly less than 1, molecules contain multiple atomic
C     nuclei each of which contributes separately to the "total number of
C     atomic nuclei", and ionization state and electrons are irrelevant.
C     All abundances should be greater than or equal to 0 and less than or
C     equal to 1. Value outside this range generate a fatal error. A warning
C     is issued if the sum of ABUND is not equal to 1. Atomic number is used
C     to index a particular element, e.g. ABUND(26) corresponds to iron.
C   SPLIST [character*(*) array(NLIST)] List of species to consider in
C     solving for the molecular equilibrium, including both the base set,
C     plus any additional species that occur in the line list.
C   NLIST [integer] the number of valid species in SPLIST.
C   TOL [real] iterative solution of the population equations terminates
C     successfully when the largest fractional change in pressure for all
C     species (including electrons) drops below this "tolerance".
C Outputs:
C   XNE [real] electron number density (in 1/cm^3) determined in GAS routine.
C   AWT [real*8] atomic weights of each species
C
      SUBROUTINE lnGAS(TEMP,Pgas,Pelec,ABUND,ELEMEN,AMASS,ELESIZ,
     *                 TOL,SPLIST,NLIST,XNE,XNA,RHO,Pgnew,XNPF,
     *                 PFUNC,POTION,XTOTAL,AWT,NGIT,
     *                 FAILED)
c      SUBROUTINE lnGAS(TEMP,XNELEC,XNATOM,ABUND,ELEMEN,AMASS,ELESIZ,
c     *                 TOL,SPLIST,NLIST,
c     *                 XNE,XNA,RHO,XNPF,PFUNC,POTION,XTOTAL,AWT,NGIT,
c     *                 FAILED)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INCLUDE 'SIZES.EOS'
C
      CHARACTER ENAME*(SPCHAR),BLANK*1
      INTEGER MAXIT,MAXREF
      DOUBLE PRECISION KBOL,HMASS,AMULOG
      PARAMETER (BLANK=' ',ENAME='e-',KBOL=1.38065D-16,MAXIT=1000,
     *           HMASS=1.66053D-24,AMULOG=-23.779751D0,MAXREF=10)

      LOGICAL PRINT,FAILED

      INTEGER NLIST,ELESIZ
      CHARACTER*(SPCHAR) SPLIST(NLIST)
      CHARACTER*(3) ELEMEN(ELESIZ)
      REAL ABUND(ELESIZ),AMASS(ELESIZ)

      CHARACTER NAMEMX*(SPCHAR),NAMET*(SPCHAR)
      INTEGER JATOM, TYPE(SPLSIZ-1),NCH(SPLSIZ-1),IATOM(ELEDIM),
     *  INDSP(ELEDIM),NAT(4,SPLSIZ-1),ZAT(4,SPLSIZ-1),NTOT(SPLSIZ-1),
     *  NEL(SPLSIZ-1),IAT(SPLSIZ-1),INDZAT(99)
      REAL T,TEMP,XNELEC,XNATOM,TOL,Pgas,Pelec,Pgnew,XNE,XNA,RHO,
     *  POTI(IONSIZ),XNPF(*),PFUNC(*),POTION(*),XTOTAL
      DOUBLE PRECISION FRACT(IONSIZ),IT(SPLSIZ-1),KT(SPLSIZ-1),
     *  AWT(SPLSIZ-1)

      DOUBLE PRECISION A(ELEDIM+1,ELEDIM+1),RHS(ELEDIM+1),
     *  AA(ELEDIM+1,ELEDIM+1),
     *  B(ELEDIM+1),BB(ELEDIM+1),
     *  P(ELEDIM+1),PP(SPLSIZ-1),PP0(SPLSIZ-1),PART(SPLSIZ-1),ND

      DOUBLE PRECISION PE,PG,PF,PNEW,PENEW,
     *  DP,DPE,PION
c      DOUBLE PRECISION AT,BT,PN,DPF(4),CRATIO,BBB(ELEDIM+1),
c     *  PENQ,DPP,DPPE
      DOUBLE PRECISION RNF(ELEDIM),AL(ELEDIM+1)
      INTEGER NELM,NCHG,ANUM(4),NATM(4),IPIV(ELEDIM+1),IWORK(ELEDIM+1),
     *  INFO
      DOUBLE PRECISION RATIOM,QPRD,RHSTOT,SCALE
c      DOUBLE PRECISION DUMMY,SCOLD,RHS0,RHS1,RHS2

c      DOUBLE PRECISION BOLD(ELEDIM+1),S(ELEDIM+1),GAMMA,BNORM,BOLDN
      DOUBLE PRECISION RSCL(ELEDIM+1),CSCL(ELEDIM+1)
c      DOUBLE PRECISION ROWCND,COLCND,AMX
      DOUBLE PRECISION FERR(1),BERR(1),WORK(5*(ELEDIM+1))
      CHARACTER*1 EQUED

      PARAMETER (NFIELDS=40)
      CHARACTER*(*) FORMAT201,FORMAT202
c      CHARACTER*(*) AFIELDS
c      PARAMETER (AFIELDS=CHAR(NFIELDS/10+ICHAR('0'))//
c     *                   CHAR(MOD(NFIELDS,10)+ICHAR('0')))
c      PARAMETER (FORMAT201='(4x,'//AFIELDS//'(1X,A3,2X))')
c      PARAMETER (FORMAT202='(A2,'//AFIELDS//'F6.1)')
      PARAMETER (FORMAT201='(4x,40(1X,A3,2X))')
      PARAMETER (FORMAT202='(A2,40F6.1)')

cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      real ttt(101)
c      real*8 Kttt(101)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

C
C Initialize the Reciprocal Neutral Fraction (RNF). The RNF is used to
C adjust the initial neutral atomic partial pressures used in the linear
C solver. Originally, atomic species were assumed to be predominantly
C neutral, but at low electron pressures, this is a poor assumption for
C species with low ionization potentials.
C
      DO 1 I=1,ELEDIM
   1  RNF(I)=1.0D0
C
C Total gas and electron pressure
C
      T=MAX(1200.,TEMP)
c      T=TEMP
      PG=Pgas
      PE=Pelec
      XNELEC=PE/(KBOL*TEMP)
      XNATOM=PG/(KBOL*TEMP)
C
C Avoid unpleasant surprises
C
      if(PG.GT.PE) THEN
        XNATOM=XNATOM-XNELEC
      ELSE
        XNELEC=XNATOM*0.01
      END IF
c      PG=(XNATOM+XNELEC)*KBOL*TEMP
c      PE=XNELEC*KBOL*TEMP
C
C  Calculate equilibrium constants for each species in list (except 'e-').
C
c      PRINT=.TRUE.
      PRINT=.FALSE.
      PION=0
      IIH2=0
      IICO=0
      IIH2O=0
      JATOM=0
      NP=0
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      open(13,file='KT_eos.dat',FORM='UNFORMATTED',STATUS='UNKNOWN')
c      write(13) NLIST,LEN(SPLIST(1))
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO 4 ISPEC=1,NLIST-1
      PP0(ISPEC)=0.D0
      CALL MPARSE(ELEMEN,SPLIST(ISPEC),NELM,NCHG,ANUM,NATM,ELESIZ)
c      write(*,*) ISPEC,'"'//SPLIST(ISPEC)//'"',NELM,NCHG,
c     *           ANUM,NATM,ELESIZ
      IF(NCHG.EQ.0) NP=ISPEC
      IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.EQ.0) THEN
C
C  Neutral atom
C
        TYPE(ISPEC)=1
        KT(ISPEC)=1.0
        IT(ISPEC)=1.0
        JATOM=JATOM+1
        IF(JATOM.GT.ELEDIM) THEN
          write(*,'(a,2i4)') 'gas: too many element types,' //
     *      ' increase ELEDIM:',ELEDIM,JATOM
          stop
        END IF
        IATOM(JATOM)=ANUM(1)
        INDSP(JATOM)=ISPEC
        IAT(ISPEC)=JATOM
        AWT(ISPEC)=AMASS(ANUM(1))
        INDZAT(ANUM(1))=JATOM
        NTOT(ISPEC)=1
        CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,
     *             IONSIZ,3)
        PART(ISPEC)=FRACT(1)
        POTION(ISPEC)=POTI(1)
      ELSE IF(NELM.EQ.1.AND.NATM(1).EQ.1.AND.NCHG.NE.0) THEN
C
C  Ionized atom
C
        TYPE(ISPEC)=3
        IF(NCHG.GT.0) THEN
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,
     *               IONSIZ,2)
          IT(ISPEC)=FRACT(NCHG+1)/FRACT(1)*PE**NCHG
          RNF(ANUM(1))=RNF(ANUM(1))+FRACT(NCHG+1)/FRACT(1)
c          if(ANUM(1).eq.26) write(*,*) SPLIST(ISPEC),NCHG,
c     *                      (FRACT(I),I=1,IONSIZ)
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,
     *               IONSIZ,3)
          PART(ISPEC)=FRACT(NCHG+1)
c          if(ANUM(1).eq.62) write(*,*) 'pf: ',SPLIST(ISPEC),NCHG,FRACT
          POTION(ISPEC)=POTI(NCHG+1)
C
C Negative ions
C
        ELSE IF(NCHG.LT.0) THEN
          CALL XSAHA(ANUM(1),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,
     *               IONSIZ,3)
          QPRD=2.D0*FRACT(1)
          CALL NEGION(ANUM(1),T,QPRD,
     *                IT(ISPEC),PART(ISPEC),POTION(ISPEC))
          KT(ISPEC)=1.
        END IF
C
        KT(ISPEC)=1.0
        AWT(ISPEC)=AMASS(ANUM(1))
        NTOT(ISPEC)=1
      ELSE IF(NELM.GT.1.OR.NATM(1).GT.1) THEN
C
C  Neutral or ionized molecule
C
        TYPE(ISPEC)=2
C
C  Calculate mass ratio (RATIOM) and partition function product (QPRD)
C  needed by MOLCON. See MOLCON header for decription of these quantities.
C  While we are at it, calculate the atomic weight (AWT) of the molecule
C  and the total number of atoms (NTOT) of any type in the molecule.
C
        NTOT(ISPEC)=0
        AWT(ISPEC)=0.0D0
        RATIOM=0.0D0
C
C  Fixed the partition function ratio for ionized molecules.
C  Now we start with a product of partition functions for free
C  electrons in denominator. NP 29-12-2006.
C       QPRD=0.0D0
        QPRD=-NCHG*LOG10(2.0)
        DO 2 IELM=1,NELM
        NTOT(ISPEC)=NTOT(ISPEC)+NATM(IELM)
        AWT(ISPEC)=AWT(ISPEC)+NATM(IELM)*AMASS(ANUM(IELM))
        RATIOM=RATIOM+NATM(IELM)*LOG10(AMASS(ANUM(IELM)))
        CALL XSAHA(ANUM(IELM),T,XNELEC,XNATOM,IONSIZ,POTI,FRACT,
     *             IONSIZ,3)
        IF(SPLIST(ISPEC).EQ.'H2')  IIH2=ISPEC
        IF(SPLIST(ISPEC).EQ.'CO')  IICO=ISPEC
        IF(SPLIST(ISPEC).EQ.'H2O') IIH2O=ISPEC
c       if(splist(ispec).eq.'N2')write(*,*)
c     *    anum(ielm),(fract(i),i=1,2)
   2    QPRD=QPRD+NATM(IELM)*LOG10(FRACT(1))
        RATIOM=RATIOM-LOG10(AWT(ISPEC))+(NTOT(ISPEC)-1)*AMULOG
C
C  Now get the molecular constants from MOLCON.
C
        CALL MOLCON(SPLIST(ISPEC),T,NTOT(ISPEC),
     *              RATIOM,QPRD,KT(ISPEC),PART(ISPEC),PION)
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        do ittt=0,100
c          ttt(ittt+1)=20.*ittt+1000.
c          CALL MOLCON(SPLIST(ISPEC),ttt(ittt+1),NTOT(ISPEC),
c     *                RATIOM,QPRD,Kttt(ittt+1),PART(ISPEC),PION)
c        END DO
c        write(13) SPLIST(ispec),ttt,Kttt
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C
C  Finally, record the charge state of the molecule.
C
        IF(NCHG.EQ.0) THEN
          IT(ISPEC)=1
        ELSE
C
C  The first option was used with Sauval & Tatum constants.
C  JV fits to NextGen pressures needed IT(ISPEC)=1.0 for positive
C  molecular ions.
C
          IF(SPLIST(ISPEC).EQ.'H2+'.OR.SPLIST(ISPEC).EQ.'NO+') THEN
            K=1
            DO IELM=2,NELM
              IF(POTION(INDSP(ANUM(IELM))).LT.POTION(INDSP(ANUM(K))))
     *          K=IELM
            ENDDO
            IT(ISPEC)=IT(INDSP(ANUM(K))+1)
            KT(ISPEC)=KT(ISPEC)/IT(ISPEC)
          endif
          IT(ISPEC)=1.0
        END IF
C
C  Store ionization potential (needed e.g. for broadening calculations)
C
        IF(PION.GT.0.D0) THEN
          POTION(ISPEC)=PION
        ELSE
c
c  If ionization potential is not available use the one for TiO!
c
          POTION(ISPEC)=6.4
        END IF
      ELSE
C
C  Fall through to here when the molecular formula doesn't make sense.
C
        WRITE(*,*) 'Wrong formula for the species: ',splist(ISPEC)
        STOP
      END IF
C
C  Now save results of MPARSE into arrays.
C
      NEL(ISPEC)=NELM
      NCH(ISPEC)=NCHG
      DO 3 IELM=1,NELM
      ZAT(IELM,ISPEC)=ANUM(IELM)
c      if(ANUM(IELM).eq.6.or.ANUM(IELM).eq.8) then
c        write(*,*) ISPEC,SPLIST(ISPEC),IT(ISPEC),KT(ISPEC)
c      endif
   3  NAT(IELM,ISPEC)=NATM(IELM)
C
C  Go back for next species.
C
      IT(ISPEC)=MAX(1.D-50,IT(ISPEC))
      KT(ISPEC)=MAX(1.D-50,KT(ISPEC))
   4  CONTINUE
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c      close(13)
c      stop
cC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      NEQ=JATOM+1
C==================================
C== End of species list parsing. ==
C==================================
C
C Print diagnostic: neutral fractions.
C
c     write(*,*) 'Reciprocal Neutral Fractions'
c     do 850 i=1,JATOM/7
c       write(*,860) (jeff(iatom(j)),j=7*i-6,7*i)
c850  continue
c860  format(1p,7e10.3,a)
c     if(JATOM.gt.7*(JATOM/7)) write(*,860)
c    *  (jeff(iatom(j)),j=7*(JATOM/7)+1,JATOM)
c      do 52 i=1,nlist-1
c  52  write(*,'(I4,1P2E12.4,3I3,A6,0Pf8.2,8I4)')
c     *  i,IT(i),KT(i),NCH(i),NTOT(i),NEL(i),SPLIST(i),AWT(i),
c     *  (ZAT(j,i),NAT(j,i),j=1,NEL(i))
C================================================================
C== UPDATE MAIN ARRAYS                                         ==
C================================================================
c
c Make the initial estimate of the partial pressures for neutral atoms. These
c pressures are used as input to the linear solver. When only abundances are
c considered, the largest errors occur for low ionization elements, which can
c be highly ionized at low electron pressures. Thus, we apply a correction
c to recover the neutral fraction for each atom. The neutral fraction only
c corrects for losses into ionization states included in the species list.
c When the ionization correction is included, the largest error in the inital
c guess for carbon, which has unaccounted for losses into CO. Late in the
c convergence process, nitrogen becomes the dominant source of error.
c
      DO 5 J=1,JATOM
      P(J)=PG*ABUND(IATOM(J))/RNF(IATOM(J))
      ISPEC=INDSP(J)
      PP0(ISPEC)=P(J)
   5  CONTINUE
c
c Make an initial guess at the balance between H and H2.
c Assumes pressures of species other than H, H2, He, and Ne are negligible.
c Constraints:
c   KT(IIH2)*PP(IIH2)=P(1)**2           <-- chemical equilibrium
c   P(1)+2*PP(IIH2)=ABUND(1)*(PG-PE)    <-- H particle conservation
c
      IF(IIH2.GT.0) THEN
        PHyd=0.5*(-KT(IIH2)+SQRT(KT(IIH2)**2
     *        +4.0*KT(IIH2)*(PG-PE-P(2)-P(10))))
      ELSE
        PHyd=(PG-PE)*ABUND(1)
      END IF
      IF(PHyd.GT.0.) P(1)=PHyd
c
c Make an initial guess at the balance between C, O, CO, and H2O.
c Constraints:
c   KT(IICO)*PP(IICO)=P(6)*P(8)         <-- chemical equilibrium
c   KT(IIH2O)*PP(IIH2O)=P(1)**2*P(8)    <-- chemical equilibrium
c   PTOTH=P(1)+2*PP(IIH2)               <-- defines density of H nuclei
c   PTOTC=P(6)+PP(IICO)                 <-- defines density of C nuclei
c   PTOTO=P(8)+PP(IICO)+PP(IIH2O)       <-- defines density of O nuclei
c   PTOTC=PTOTH*ABUND(6)/ABUND(1)       <-- abundance constraint
c   PTOTO=PTOTH*ABUND(8)/ABUND(1)       <-- abundance constraint
c
      PTOTH=P(1)
      IF(IIH2.GT.0) PTOTH=PTOTH+2.0*P(1)**2/KT(IIH2)
      PTOTC=PTOTH*ABUND(6)/ABUND(1)
      PTOTO=PTOTH*ABUND(8)/ABUND(1)
      IF(IIH2O.GT.0) THEN
        WATCOR=1.0+P(1)**2/KT(IIH2O)
        AQUAD=1.0/WATCOR
        IF(IICO.GT.0) THEN
          BQUAD=KT(IICO)+(PTOTO-PTOTC)/WATCOR
          CQUAD=-KT(IICO)*PTOTC
          P(6)=(-BQUAD+SQRT(BQUAD**2-4.0*AQUAD*CQUAD))/(2.0*AQUAD)
          P(8)=(P(6)+PTOTO-PTOTC)/WATCOR
        ELSE
          P(6)=PTOTC
          P(8)=PTOTO
        END IF
      ELSE
        P(6)=PTOTC
        P(8)=PTOTO
      END IF
      IF(P(6).LE.0.) P(6)=PTOTC
      IF(P(8).LE.0.) P(8)=PTOTO
      PE0=PE
      NAMEMX=BLANK
      DELMAX=0.0D0
c      COMPZ=0.0D0
c      PZS=0.0D0
c      DO 6 J=1,JATOM
c      NN=INDSP(J)
c      IF(IPR(NN).NE.2) GOTO 3
c      NNP=INDX(3,ITAB(ZAT(1,NN)),1,1,1)
c      COMPZ=COMPZ+ABUND(IATOM(J))
c      IF(PE.EQ.0.0D0) PZS= PZS + P(J)
c      IF(PE.GT.0.0D0) PZS= PZS + (1.0D0+IT(NNP)/PE)*P(J)
c   6  CONTINUE
c      do J=1,JATOM
c        write(*,*) J,P(J),ABUND(IATOM(J)),SPLIST(INDSP(J))
c      END DO
c      write(*,*) JATOM+1,PE,'e-'
c      stop
C================================================================
C== MAIN LOOP: FILL LINEARIZED COEFFICIENT MATRIX AND RHS VECTOR,
C== AND SOLVE SYSTEM FOR PARTIAL PRESSURE CORRECTIONS.         ==
C== ISOLV=1: LINEARIZE ONLY THE PARTIAL PRESSURES OF THE NEUTRAL=
C== ATOMS FOR WHICH IPR(J)=1 (MAJOR SPECIES). THE ELECTRON     ==
C== PRESSURE PE IS ASSUMED TO BE GIVEN IN THIS CASE, AND SO IS ==
C== NOT INCLUDED IN THE LINEARIZATION. THIS IS NECESSARY SINCE ==
C== MOST OF THESE ELECTRONS (AT COOL TEMPS.) ORIGINATE FROM    ==
C== ELEMENTS NOT CONSIDERED IN THE LINEARIZATION. IN ORDER TO  ==
C== OBTAIN A GOOD VALUE FOR PE IN THE FIRST PLACE, IT IS       ==
C== NECESSARY TO CALL GAS WITH ISOLV=2.                        ==
C== ISOLV=2: THIS LINEARIZES THE PARTIAL PRESSURES OF THE NEUTRAL
C== ATOMS FOR WHICH IPR(J)=1 OR 2. THIS LIST OF ELEMENTS SHOULD==
C== INCLUDE ALL THE SIGNIFICANT CONTRIBUTORS TO THE TOTAL      ==
C== PRESSURE PG, AS WELL AS THE ELECTON PRESSURE PE. ANY ELEMENT=
C== (IPR(J)=3) NOT INCLUDED IS ASSUMED TO HAVE A NEGLIGIBLE    ==
C== EFFECT ON BOTH P AND PE.                                   ==
C== IN BOTH CASES, THE PARTIAL PRESSURES OF THE NEUTRAL ATOMS  ==
C== FOR ELEMENTS NOT INCLUDED IN THE LINEARIZATION ARE         ==
C== CALCULATED DIRECTLY FROM THE NOW DETERMINED PRESSURES OF   ==
C== THE LINEARIZED ELEMENTS.                                   ==
C================================================================
      NGIT=0
      RHSTOT=1.D99
C
C Top of loop in which linearized equations are solved recursively.
C
      KMAX=1
      DO J=1,NEQ-1
        P(J)=LOG(MAX(P(J),1.D-50))
      END DO
      PE=LOG(MAX(PE,1.D-50))
      IDIR=0
c      open(unit=4,file='dump.bin',form='UNFORMATTED')
c      write(4) NEQ
   7  IF(NGIT.GE.MAXIT) THEN
        WRITE(6,208)
 208    FORMAT('*** ERROR: TOO MANY ITERATIONS IN ROUTINE "GAS"')
        WRITE(6,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),RHSTOT
        write(*,*) TEMP,PG,P(1),XNATOM,XNELEC
        STOP
      END IF
      NGIT=NGIT+1
      P(NEQ)=PE

c      do J=1,NEQ
c        p(J)=exp(p(j))
c      enddo
c      CALL EOSFCN(NEQ,P,B,A,1,PG,NCH,NLIST,
c     *     IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
c      CALL EOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,
c     *     IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
c      do j=1,NEQ
c        SCALE=P(J)
c        P(J)=P(J)+0.001d0
c        CALL EOSFCN(NEQ,P,BB,A,1,PG,NCH,NLIST,
c     *    IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
c        write(*,*) J,SCALE
c        write(*,'(40e10.3)')(a(i,j)-(bb(i)-b(i))/0.001d0
c     *                            ,i=1,40)
c        write(*,'(40e10.3)')(a(i,j),i=1,40)
c        write(*,'(40e10.3)')((bb(i)-b(i))/(0.001d0),i=1,40)
c        write(*,'(40e10.3)')(bb(i),i=1,40)
c        P(J)=SCALE
c      enddo
c      stop

      SCALE=10.D0
   9  CALL lnEOSFCN(NEQ,P,B,A,1,PG,NCH,NLIST,
     *              IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
      IF(B(1).GT.0.05D0*PG) THEN
        IF(IDIR.NE.-1) THEN
          SCALE=ABS(SCALE*0.1D0)
          IDIR=-1
        END IF
C
C Neutral atomic pressures are too high. Scale them down until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)-SCALE
        END DO
        GOTO 9
      ELSE IF(B(1).LT.-0.05D0*PG) THEN
        IF(IDIR.NE.1) THEN
          SCALE=ABS(SCALE*0.1D0)
          IDIR=1
        END IF
C
C Neutral atomic pressures are too low. Scale them up until
C partical conservation equation will become negative
C
        DO J=1,NEQ-1
          P(J)=P(J)+SCALE
        END DO
        GOTO 9
      END IF
      CALL lnEOSFCN(NEQ,P,B,A,2,PG,NCH,NLIST,
     *              IATOM,INDSP,NAT,ZAT,NTOT,NEL,IAT,INDZAT,ABUND,KT,IT)
C
C================================================================
C== NOW SOLVE THE LINEARIZED EQUATIONS (USING ROUTINE "LINEQ") ==
C================================================================
      IF(PRINT) THEN
        WRITE(*,200) NGIT
 200    FORMAT('LOG OF COEFFICIENT MATRIX AT ITERATION #',I5//)
        KK=MIN(NFIELDS,NEQ-1)
        WRITE(*,FORMAT201) (SPLIST(INDSP(K)),K=1,KK-1),'e-','RHS'
        DO 21 I=1,KK-1
        DO 20 J=1,KK-1
  20    AL(J)=LOG10(ABS(A(J,I))+1.0D-50)
        AL(KK)=LOG10(ABS(A(NEQ,I))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(I))+1.0D-50)
        NAMET=SPLIST(INDSP(I))
        WRITE(*,FORMAT202) NAMET,(AL(J),J=1,KK+1)
  21    CONTINUE
        DO 22 J=1,KK-1
        AL(J)=LOG10(ABS(A(J,NEQ))+1.0D-50)
  22    CONTINUE
        AL(KK)=LOG10(ABS(A(NEQ,NEQ))+1.0D-50)
        AL(KK+1)=LOG10(ABS(B(NEQ))+1.0D-50)
        NAMET='e-'
        WRITE(*,FORMAT202) NAMET,(AL(J),J=1,KK+1)
        WRITE(*,'(/)')
      END IF
C
C  Save a copy of the RHS for future step refinement
C
      DO 23 I=1,NEQ
  23  RHS(I)=B(I)
      RHSTOT=DASUM(NEQ,RHS,1)
C
C  Solve linear system for corrections
C  In order not to solve for Pelect, one should use NEQ-1 as the first
C  argument. NEQ solves the whole system including electron pressure
C
c
c  Using LAPACK routine
c
c        open(unit=4,file='dump.bin',form='UNFORMATTED')
c        write(4) NEQ
c        write(4) ((A(i,j),i=1,NEQ),j=1,NEQ)
c        write(4) (B(i),i=1,NEQ)
c      write(4) ((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
      CALL DGESVX('E','N',NEQ,1,A,ELEDIM+1,AA,ELEDIM+1,IPIV,EQUED,
     *            RSCL,CSCL,B,ELEDIM+1,BB,ELEDIM+1,RCOND,FERR,BERR,
     *            WORK,IWORK,INFO)
      CALL DCOPY(NEQ,BB,1,B,1)
c      write(4) ((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
c
c  The same thing using LINEQ2 or LINEQ and BLAS 2/3
c      CALL LINEQ(NEQ,1,A,ELEDIM+1,IPIV,B,ELEDIM+1,INFO)
      IF(INFO.NE.0) THEN
        WRITE(*,*) 'EOS: LINEQ failed to solved for corrections to'
        WRITE(*,*) '     the partial pressures. Matrix is degenerate'
        WRITE(*,*) '     Temp=',TEMP,', Natom=',XNATOM,', Nelec=',XNELEC
        WRITE(*,*) '     Pg=',PG,', INFO=',INFO,
     *             ', Element: ',SPLIST(INDSP(INFO)),
     *             ', Iter=',NGIT,' EQUED=',EQUED
c        open(unit=4,file='dump.bin',form='UNFORMATTED')
c        write(4) NEQ,((A(i,j),i=1,NEQ),j=1,NEQ),(B(i),i=1,NEQ)
c        close(4)
c        write(1) 0
c        close(1)
        IF(PRINT) THEN
c          close(4)
          STOP
        END IF
c        DO J=1,NEQ
c          P(J)=MAX(P(J)+0.1D0,-115.d0)
c          write(*,*) J,P(J),B(J),B(J)*FACTOR
c        END DO
        write(*,*) P(INFO),B(INFO),B(INFO)*FACTOR
        P(INFO)=MAX(P(INFO)+0.1D0,-115.d0)
        PRINT=.TRUE.
        GO TO 9
      END IF
c
C=================================================================
C== FINALLY, UPDATE THE PARTIAL PRESSURES FOR THE MAJOR SPECIES ==
C== BY ADDING THE PRESSURE CORRECTIONS OBTAINED FOR EACH ATOM   ==
C== FROM THE LINEARIZATION PROCEDURE.                           ==
C=================================================================
      DELMAX=-200.0D0
      KMAX=1
      DO K=1,JATOM
c        write(*,*) K,P(K),B(K)
        ISPEC=INDSP(K)
        DP=ABS(P(K))
        DELP=ABS(B(K))
        IF(DP.GT.1.D-10) DELP=DELP/DP
        IF(DELP.GT.DELMAX) THEN
          NAMEMX=SPLIST(ISPEC)
          DELMAX=DELP
          KMAX=K
        END IF
      END DO
      DPE=ABS(P(NEQ))
      DELPE=ABS(B(NEQ))
      IF(DPE.GT.1.D-10) DELPE=DELPE/DPE
      IF(DELPE.GT.DELMAX) THEN
        NAMEMX=ENAME
        DELMAX=DELPE
        KMAX=NEQ
      END IF
C
C  Under-relaxation factor
C
      FACTOR=0.2D0/(DELMAX+0.2D0)
      DO K=1,JATOM
C
C  Apply corrections
C
        DP=B(K)*FACTOR
c        DP=10.D0*DP/MAX(10.D0,ABS(DP))
        PNEW=P(K)-DP
        P(K)=MAX(PNEW,-115.D0)
      END DO
      DP=B(NEQ)*FACTOR
c      DP=10.D0*DP/MAX(10.D0,ABS(DP))
      PENEW=PE-DP
      PE=MAX(PENEW,-115.D0)
C================================================================
C== PRINT OUT SUMMARY LINE FOR EACH ITERATION                  ==
C================================================================
      PTOT=EXP(PE)
      PQ=0.0D0
      DO ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=-PE*NQ+LOG(IT(ISPEC))-LOG(KT(ISPEC))
        DO I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF+P(J)*NAT(I,ISPEC)
        END DO
        PP(ISPEC)=EXP(PF)
        PTOT=PTOT+PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
c        write(*,*) ISPEC,SPLIST(ISPEC),PP(ISPEC),PTOT,PG,NQ,PQ,EXP(PE)
      END DO
c      stop
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(EXP(PE)-PQ)/PG
c      write(*,*) DELMAX,DPTOT,DPQ
      IF(PRINT) THEN
        WRITE(6,203) NGIT,NAMEMX,DELMAX,PE,B(KMAX),P(KMAX),
     *               PTOT/TEMP/KBOL,DPTOT,EXP(PE)/TEMP/KBOL,DPQ,FACTOR
 203    FORMAT(I10,2X,A8,1P,9E11.3)
      END IF
      write(*,*) NGIT,TOL,RHSTOT,DPTOT,DELMAX,PTOT,PG,DPQ,FACTOR
      IF((RHSTOT.GT.TOL.OR.DPTOT.GT.TOL.OR.DELMAX.GT.TOL)
     *   .AND.NGIT.LT.MAXIT) GO TO 7
C
C Bottom of the loop in which linearized equations are solved recursively.
C
C================================================================
C== CALCULATE FINAL PARTIAL PRESSURES AFTER CONVERGENCE OBTAINED=
C================================================================
c      write(*,*) RHSTOT,DELMAX,DPTOT,DPQ,TOL
      PTOT=EXP(PE)
      PD=0.0D0
      PU=0.0D0
      PQ=0.0D0
      DO ISPEC=1,NLIST-1
        NELT=NEL(ISPEC)
        NQ=NCH(ISPEC)
        PF=-PE*NQ+LOG(IT(ISPEC))-LOG(KT(ISPEC))
        DO I=1,NELT
          J=INDZAT(ZAT(I,ISPEC))
          PF=PF+P(J)*NAT(I,ISPEC)
        END DO
        PP(ISPEC)=EXP(PF)
        PTOT=PTOT+PP(ISPEC)
        PD=PD+NTOT(ISPEC)*PP(ISPEC)
        PQ=PQ+NQ*PP(ISPEC)
        PU=PU+AWT(ISPEC)*PP(ISPEC)
c        write(*,*) ISPEC,SPLIST(ISPEC),PP(ISPEC),PTOT,PG,NQ,PQ,EXP(PE)
      END DO
      PE=EXP(PE)
      DO J=1,JATOM
        P(J)=EXP(P(J))
      END DO
      PP(NLIST)=PE
      PDTOT=PD+PE
      DPTOT=DABS(PTOT-PG)/PG
      DPQ=DABS(PQ-PE)/PG
      GMU=PU/PTOT
      ND=PTOT/(TEMP*KBOL)
      RHO=ND*GMU*HMASS
      XNE=PE/(TEMP*KBOL)
C================================================================
C== WRITE OUT FINAL PARTIAL PRESSURES                          ==
C================================================================
      IF(PRINT) THEN
        write(*,'(''AFTER '',I3,'' iterations.   Max change of:'',G10.3,
     #      ''  in element:'',A)') NGIT,DELMAX,NAMEMX
        WRITE(6,'(''AFTER '',I3,'' ITERATIONS WITH ''/
     #            ''T='',1PE10.3,''   P='',E10.3)') NGIT,TEMP,PG
        WRITE(6,'(''PDTOT='',1PE10.3,''   DPTOT='',E10.3,
     #            ''  DPQ='',E10.3,''  Nelectron='',E10.3,'' cm^3''/
     #    '' Nparticle='',1PE10.3,'' cm^3   Mean At.Wt.='',
     #    0PF7.3,''   Density='',1PE10.3,'' g/cm^3''/
     #    '' # Species   Abundance   Initial P   Final P'',
     #    ''      IT         KT         pf''//)')
     #   PDTOT,DPTOT,DPQ,XNE,ND-XNE,GMU,RHO
        NSP1=NLIST
        DO 35 ISPEC=1,NLIST-1
        IF(TYPE(ISPEC).NE.1) THEN
          WRITE(*,206) ISPEC,SPLIST(ISPEC),PP0(ISPEC),PP(ISPEC),
     #                 IT(ISPEC),KT(ISPEC),PART(ISPEC)
 206      FORMAT(I3,1X,A8,11X,1P,5E11.3)
        ELSE
          J=IAT(ISPEC)
          WRITE(*,207) ISPEC,splist(ISPEC),ABUND(IATOM(J)),PP0(ISPEC),
     #                 PP(ISPEC),IT(ISPEC),KT(ISPEC),PART(ISPEC)
 207      FORMAT(I3,1X,A8,1P,6E11.3)
        END IF
  35    CONTINUE
        WRITE(*,206) NSP1,ENAME,PE0,EXP(PE)
      END IF
C
C Fill up the output array and set up flags
C PNOTE is the partial pressure due to everything except electrons.
C XNA is the number density of everything except electrons.
C
      PNOTE=0.0
      DO 36 ISPEC=1,NLIST-1
      IF(PART(ISPEC).GT.0.) THEN
        IF(PP(ISPEC)/KBOL/TEMP/PART(ISPEC).GE.1.D-20) THEN
          XNPF(ISPEC)=PP(ISPEC)/(KBOL*TEMP*PART(ISPEC))
        ELSE
          XNPF(ISPEC)=0.0
        END IF
        PFUNC(ISPEC)=PART(ISPEC)
      ELSE
        XNPF(ISPEC)=0.
        PFUNC(ISPEC)=1.
      END IF
      PNOTE=PNOTE+PP(ISPEC)
c      write(*,*) ISPEC,PNOTE,PP(ISPEC),SPLIST(ISPEC)
  36  CONTINUE
      XNPF(NLIST)=XNE
      PFUNC(NLIST)=1.0
      XTOTAL=PD/(TEMP*KBOL)
      XNA=PNOTE/(TEMP*KBOL)
c      write(*,*) 'Pg,PD,PNOTE,PE,PNOTE+PE',Pg,PD,PTOT,PE,PNOTE+PE
      Pgnew=Ptot
C
      RETURN
      END


C=========================================================================
C MOLCON: Returns equilibrium constant and partition function for a given
C   molecule and temperature.
C
C Inputs:
C   SPNAME [character(*)] Name of molecule, chosen from SPLIST below.
C   T [real] Temperature (in K) at which EQK and PART are to be found.
C   NTOT [real] Total number of atoms in the molecule.
C   RATIOM [real] Logarithm (base 10) of mass ratio (in g^(natoms-1)):
C     ratiom = Sum{log10(Atomic Masses)} - log10(Sum{Atomic Masses})
C   QPRD [double] Logarithm of product of atomic partition functions:
C     qprd = Sum{log10(Atomic Partition Functions)}
C
C Outputs:
C   EQK [real] Equilibrium constant (in dynes/cm/cm) at temperature T,
C     calculated from dissociation energy and partition function.
C   PART [real] Partition function at temperature T, calculated from
C     expressions in the references cited below.
C
C References:
C   For diatomic molecules: Sauval & Tatum (1984, ApJS, 56, 193).
C
      SUBROUTINE MOLCON(SPNAME,T,NTOT,RATIOM,QPRD,EQK,PART,PION)
C
      INTEGER MSPEC,NTOT
      DOUBLE PRECISION KERG,KEV
      DOUBLE PRECISION RATIOM,QPRD,PION
      PARAMETER (KERG=1.38065D-16,KEV=KERG/1.60219D-12)
      PARAMETER (CONST=25947.256)
C
      REAL T
      DOUBLE PRECISION TH,LOGTH,EQK_ST,EQK,PART,Qm_spln,Kp_spln
C
C Combine equilibrium constant coefficients into one large array.
C
      PARAMETER (MSPEC=196)
      PARAMETER (NEQCOE=7)
      DOUBLE PRECISION COEF(NEQCOE,MSPEC)
      DOUBLE PRECISION C01(NEQCOE,50),C02(NEQCOE,50),
     *                 C03(NEQCOE,50),C04(NEQCOE,46)
      EQUIVALENCE (C01(1,1),COEF(1,  1)),(C02(1,1),COEF(1, 51))
      EQUIVALENCE (C03(1,1),COEF(1,101)),(C04(1,1),COEF(1,151))
C
C Combine partition function coefficients into one large array.
C
      PARAMETER (NPCOEF=11)
      DOUBLE PRECISION PCOEF(NPCOEF,MSPEC)
      DOUBLE PRECISION P01(NPCOEF,50),P02(NPCOEF,50),
     *                 P03(NPCOEF,50),P04(NPCOEF,46)
      EQUIVALENCE (P01(1,1),PCOEF(1,  1)),(P02(1,1),PCOEF(1, 51))
      EQUIVALENCE (P03(1,1),PCOEF(1,101)),(P04(1,1),PCOEF(1,151))
C
      CHARACTER SPNAME*(*),SPLIST(MSPEC)*8
      SAVE
C
C Molecular species list from NextGen models (Allard & Hauschildt).
C See old/eos.4.f for molecular species list from Sauval & Tatum (1984).
C
      DATA SPLIST/
     * 'H2      ','CO      ','H2O     ','OH      ','N2      ',
     * 'SiO     ','HS      ','H2S     ','NH      ','SiH     ',
     * 'CH      ','H2+     ','NO      ','MgH     ','HCl     ',
     * 'SiS     ','AlOH    ','NH2     ','AlH     ','CN      ',
     * 'CO2     ','SO      ','TiO     ','S2      ','FeH     ',
     * 'NH3     ','HCN     ','HCO     ','O2      ','CH2     ',
     * 'HF      ','H3+     ','CaH     ','Al2O    ','AlO     ',
     * 'CH3     ','SiH2    ','MgO     ','C2      ','TiO2    ',
     * 'VO2     ','NaH     ','AlCl    ','AlF     ','VO      ',
     * 'CS      ','MgOH    ','PO2     ','CaOH    ','PH2     ',
     * 'C2H     ','ScO     ','AlO2H   ','AlS     ','FeO     ',
     * 'CrO     ','CH4     ','NS      ','SO2     ','SiN     ',
     * 'OH-     ','ZrO     ','NO+     ','ZrO2    ','BO      ',
     * 'SiO2    ','HBO     ','SiC     ','YO2     ','TiS     ',
     * 'HBO2    ','C2H2    ','OCS     ','ZrO+    ','NaOH    ',
     * 'CaCl    ','AlOF    ','YO      ','NaCl    ','C2O     ',
     * 'CHP     ','HS-     ','H2-     ','TiH     ','PH3     ',
     * 'MgS     ','TiO+    ','LaO2    ','Si2     ','SiH4    ',
     * 'BH2     ','AlOCl   ','LaO     ','C2N     ','AlBO2   ',
     * 'KCl     ','SiH-    ','CaF     ','CaO2H2  ','KOH     ',
     * 'CN-     ','Al2O2   ','BaOH    ','SrOH    ','BO2     ',
     * 'SiF     ','CH-     ','C3      ','C2-     ','MgO2H2  ',
     * 'BeOH    ','HBS     ','SiC2    ','FeO2H2  ','CrO2    ',
     * 'BeH2O2  ','BH3     ','NaCN    ','BeH2    ','Si2N    ',
     * 'CaCl2   ','NaBO2   ','C3H     ','OBF     ','CS2     ',
     * 'LiOH    ','Al2     ','LiCl    ','TiOCl   ','C2H4    ',
     * 'CHCl    ','TiCl    ','AlOF2   ','KBO2    ','Si2C    ',
     * 'CHF     ','BO-     ','AlO2    ','BaO2H2  ','OTiF    ',
     * 'CS-     ','C2N2    ','SrO2H2  ','ClCN    ','AlClF   ',
     * 'KCN     ','AlCl2   ','BaCl2   ','AlF2    ','MgCl2   ',
     * 'FeO-    ','BO2H2   ','SiH3Cl  ','FeCl2   ','Si3     ',
     * 'SiH3F   ','CH3Cl   ','SrCl2   ','CaF2    ','TiF2    ',
     * 'LiBO2   ','MgClF   ','BeBO2   ','C2HCl   ','TiCl2   ',
     * 'C4      ','H3BO3   ','MgF2    ','BaClF   ','BeF2    ',
     * 'C2HF    ','BeCl2   ','TiOCl2  ','ZrCl2   ','BaF2    ',
     * 'BeC2    ','Be2O    ','SrF2    ','ZrF2    ','FeF2    ',
     * 'P4      ','SiH2F2  ','H3O+    ','C5      ','TiF3    ',
     * 'TiCl3   ','ZrCl3   ','Na2Cl2  ','Na2O2H2 ','Be3O3   ',
     * 'K2Cl2   ','K2O2H2  ','ZrCl4   ','Na2C2N2 ','ZrF4    ',
     * 'Li2O2H2 '/
C
C Dissociation energy (first column, in eV) and equilibrium constant
C   coefficients. See the file "atomiz.notes" for the information on the
C   origin of the dissociation energies. The polynomial fit coefficients
C   for the equilibrium constants were determined with "ng_kfit.pro" and
C   are meant to reproduce the constants used in constructing the NextGen
C   models. The NextGen equilibrium constants were fit over the temperature
C   range 1600 < T < 7730 K. The fits are likely to diverge rapidly from
C   the truth outside this temperature range.
C Equilibrium constants may be constructed from the coefficients using:
C
C     log10(Kp) = Sum{i=2,7}{COEF(i)*log10(THETA)**(i-2)} - COEF(1)*THETA
C
      DATA C01/
     *   4.4781, 12.1354, -0.7752, -0.7821,  0.1464,  0.1603, -0.0626,  H2
     *  11.0920, 13.2368, -0.8342, -0.0477, -0.2923, -0.4557,  0.6108,  CO
     *   9.6221, 24.7774, -2.3428,  1.6868, -1.2845, -2.9925,  3.6555,  H2O
     *   4.3920, 11.8016, -0.8507, -0.5193,  0.0502, -0.3409,  0.4836,  OH
     *   9.7594, 12.8868, -0.8813,  0.2639, -1.5912,  1.5866, -0.5407,  N2
     *   8.2600, 12.9252, -0.7608, -0.3541,  1.5620, -3.5952,  2.5962,  SiO
     *   3.5500, 11.4382, -0.7816, -0.4659,  0.4314, -1.2144,  0.9648,  HS
     *   7.5946, 23.8543, -0.9525, -0.8118,  0.2051, -1.0299,  1.1555,  H2S
     *   3.4700, 11.4658, -0.7258, -0.6418, -0.0442,  0.2836, -0.1618,  NH
     *   3.0600, 11.2595, -0.6962, -0.6435,  0.6663, -0.3357, -0.4151,  SiH
     *   3.4650, 11.5333, -0.5255, -0.7105,  0.2264, -0.9271,  0.9577,  CH
     *   2.6508, 15.8052, 33.7578, 34.5956, 27.3455, 16.6214,  9.9717,  H2+
     *   6.4968, 11.9347, -0.7596,  0.0953, -0.9731,  0.8265, -0.2151,  NO
     *   1.3400, 10.2911, -0.3698, -0.0655, -2.9771,  6.1325, -4.3869,  MgH
     *   4.4336, 11.9041, -0.8281, -0.6163,  0.1580, -0.5068,  0.5164,  HCl
     *   6.4200, 12.6363, -0.7355,  0.0488,  0.8442, -2.0131,  1.3603,  SiS
     *  10.1252, 25.2575, -0.6810, -0.3051, -1.5765,  2.7536, -1.8355,  AlOH
     *   7.4400, 23.7389, -1.0179, -0.9947, -1.4353,  3.2530, -1.9224,  NH2
     *   3.0600, 11.4907, -0.4322, -0.6561, -0.5978,  2.4923, -2.4038,  AlH
     *   7.7600, 12.4438, -0.4756, -0.4909, -1.4623,  2.6823, -1.5396,  CN
     *  16.5382, 26.9571, -0.7464, -0.4921, -0.8506, -0.1365,  0.2358,  CO2
     *   5.3590, 12.3380, -0.4956, -0.2251, -0.1907, -0.2038,  0.2579,  SO
     *   6.8700, 11.9229, -1.4044,  0.7899, -0.7317, -0.0193, -0.4994,  TiO
     *   4.3693, 12.3190, -0.5050, -0.0290, -0.0266, -0.6002,  0.4572,  S2
     *   2.4100, 12.1214,  0.9438,  2.2756, -0.1086,  4.1281, -1.9952,  FeH
     *  12.1388, 36.6661, -1.4062, -0.9258, -1.6969,  0.6005,  1.2302,  NH3
     *  13.2363, 25.1318, -0.5532, -0.0850, -0.9817,  0.6676,  0.3054,  HCN
     *  11.8560, 24.6414, -0.9415, -0.1856, -0.2948, -0.1630,  0.5836,  HCO
     *   5.1156, 12.8758, -0.4856, -0.5054, -0.0776, -0.0713,  0.2369,  O2
     *   7.9400, 23.8609, -1.0762, -0.4928, -0.4092,  0.0031,  0.3761,  CH2
     *   5.8690, 12.2896, -0.9180, -0.6238,  0.1243, -0.3525,  0.4767,  HF
     *   0.0000, 18.8343, 12.4131, 11.9991,  6.8079,  8.4071,  2.6202,  H3+
     *   1.7000, 10.1982, -0.9309,  1.8315, -5.6059,  6.9571, -3.5023,  CaH
     *  10.9653, 24.8807, -0.0033,  0.4796, -1.6979,  3.5631, -2.5414,  Al2O
     *   5.2700, 12.2132, -0.5246, -0.1918, -0.6810,  1.7287, -1.5839,  AlO
     *  12.6885, 36.6540, -1.3373, -1.0064, -0.5880, -0.2362,  0.8764,  CH3
     *   0.0000, 17.8513,-15.5361,-17.6144,-13.1604, -6.4819, -5.6361,  SiH2
     *   3.5300, 10.7940,  0.0122,  1.1189, -1.8758,  2.9976, -2.7758,  MgO
     *   6.2100, 12.4672, -0.4452, -0.0100, -0.1868, -0.3860,  0.6230,  C2
     *  13.2915, 25.9340, -1.4243,  1.6519, -0.7240, -0.7271,  0.7518,  TiO2
     *  12.9619, 25.9238, -1.2927,  1.3710, -2.4073,  2.2875, -0.5486,  VO2
     *   1.8800, 10.7184, -0.3642,  0.7843, -6.5309, 13.2912, -9.9502,  NaH
     *   5.1200, 11.8277, -0.3468, -1.0735,  1.8038, -1.7748,  0.4333,  AlCl
     *   6.8900, 12.2422, -0.4905, -0.4198,  0.0242,  0.3868, -0.5765,  AlF
     *   6.4100, 12.8108, -0.5811, -0.7895, -2.6766,  8.5158, -6.9993,  VO
     *   7.3550, 12.8487, -0.7627, -0.2538,  1.5240, -4.0119,  3.0234,  CS
     *   8.0735, 23.3256, -0.5884,  0.3637, -2.4401,  3.3936, -1.7121,  MgOH
     *  11.7451, 25.2051, -0.9105,  1.0031, -0.7207, -1.1064,  1.6239,  PO2
     *   8.7035, 23.1900, -1.0964,  2.5340, -5.9823,  5.3416, -1.1946,  CaOH
     *   6.4895, 23.0863, -1.3781,  0.2539, -0.6746, -1.2341,  1.5623/  PH2
      DATA C02/
     *  12.2087, 24.9752, -0.3204, -0.5640, -0.8997,  1.6927, -0.7771,  C2H
     *   6.9600, 12.5225, -1.2695,  1.7628, -2.0543, -1.2215,  2.3706,  ScO
     *  15.6364, 37.7022, -0.5885, -0.0823, -1.7283,  3.0502, -2.0176,  AlO2H
     *   3.8400, 11.9140, -0.5187, -0.1193, -0.3886,  1.1704, -1.2299,  AlS
     *   4.2000, 12.5326, -1.0657,  1.0360, -1.5641,  0.9560, -0.3218,  FeO
     *   4.4000, 11.0587, -1.3926,  1.4461, -2.1552,  3.3409, -3.1078,  CrO
     *  17.2173, 49.9426, -0.9720, -2.4957, -0.0017, -2.3299,  3.1042,  CH4
     *   4.8000, 11.9223, -0.6951,  0.1870, -0.7158,  0.4121,  0.0296,  NS
     *  11.1405, 25.9246, -0.5809,  0.0734, -0.3333,  0.1699,  0.0529,  SO2
     *   6.6880, 14.0972,  4.2904,  4.9608,  2.9390,  3.9789,  0.8908,  SiN
     *   4.7600, 19.9888, -6.7088, -4.3846, -2.8142, -2.3004, -0.3157,  OH-
     *   7.8500, 12.4674, -1.1280,  0.0368,  0.2221,  1.1043, -1.8804,  ZrO
     *  10.8500, 17.5169, 33.0097, 36.2110, 26.7396, 15.2392, 11.4130,  NO+
     *  14.4650, 25.6324, -1.5339,  1.1586, -0.9355,  1.6114, -1.2154,  ZrO2
     *   8.2800, 12.6246, -0.6966, -0.3874,  0.2531, -0.7582,  0.5307,  BO
     *  13.0355, 26.5610, -0.2891,  0.3006, -0.4009,  0.5864, -0.4006,  SiO2
     *  12.7425, 25.2283, -0.4780, -0.3611, -0.2189, -0.2108,  0.5883,  HBO
     *   4.6400, 11.8909, -0.8762,  0.1138,  0.0665, -0.5226,  0.3331,  SiC
     *  15.2000, 25.8617, -1.4050, -0.3896,  1.0805,  2.9269, -3.7531,  YO2
     *   4.7500, 11.6628, -1.4463,  1.3742, -0.8127, -0.4623,  0.2288,  TiS
     *  19.0991, 38.4541, -0.7808, -0.4220, -0.9239,  1.0793, -0.2304,  HBO2
     *  16.9704, 37.7481, -0.2529, -1.0622, -0.1485, -0.7058,  1.1910,  C2H2
     *  14.3762, 26.3815, -0.1712,  0.1197,  0.0059, -0.9891,  1.1946,  OCS
     *   0.0000,  2.5576, -0.5567, -4.5109, -4.3690, -0.1528, -3.1319,  ZrO+
     *   8.0150, 23.3420, -0.6139,  1.4091, -6.8466, 13.0407, -9.2977,  NaOH
     *   4.0900, 10.6268, -1.1367,  2.5278, -5.6022,  4.8741, -1.1616,  CaCl
     *  12.9003, 25.5751, -0.0730,  0.2808, -1.1757,  2.3733, -1.6726,  AlOF
     *   7.2900, 12.4422, -1.3547,  1.3087,  0.1688, -5.4106,  5.1158,  YO
     *   4.2300, 11.0864, -0.4463,  1.1926, -7.5820, 15.2552,-11.1116,  NaCl
     *  14.5371, 25.6134, -0.0508,  0.3710, -0.6246, -0.7682,  0.5868,  C2O
     *  11.4442, 24.7107, -0.5678, -0.0389,  1.0076, -4.6514,  4.3893,  CHP
     *   3.7900, 19.0227, -8.0668, -5.9821, -3.8685, -3.1838, -1.0364,  HS-
     *   0.7300, 19.7162, -5.0018, -2.7680, -1.2845, -0.9859, -0.3380,  H2-
     *   2.1200, 12.4717,  0.1601,  1.4596, -0.2012,  5.0788, -4.5487,  TiH
     *   9.7800, 35.8044, -1.3937, -0.2650, -0.6732, -2.5437,  2.9710,  PH3
     *   2.4000, 11.3146, -0.5595,  0.3619, -2.0065,  3.8766, -2.9900,  MgS
     *   0.0000,  4.5751,  3.4421,  0.7560, -1.7011,  1.4510, -1.3922,  TiO+
     *  21.1510, 31.0805, 10.7070, 12.8687, 10.5799,  6.4414,  3.6171,  LaO2
     *   3.2100, 12.1817, -0.7102, -0.2403,  1.1042, -1.3644,  0.3198,  Si2
     *  13.2716, 48.6914, -1.0602, -1.2802, -0.8603,  0.1159, -0.0701,  SiH4
     *   8.2349, 24.0157, -0.6514, -0.6064, -0.6542,  0.9096, -0.5839,  BH2
     *  10.9011, 25.1839, -0.1060,  0.2530, -1.1850,  2.3355, -1.6111,  AlOCl
     *   8.2300, 12.1920,  0.1751, -0.7678, -1.3836,  1.7704, -0.0265,  LaO
     *  14.0629, 25.1475, -0.2270,  0.7024, -0.8499,  0.4583,  0.1889,  C2N
     *  20.0747, 38.6719, -0.2664,  0.2782, -1.2642,  1.6020, -0.5248,  AlBO2
     *   4.3400, 10.9561, -0.8720,  3.4218,-12.2306, 18.7863,-11.1011,  KCl
     *   3.2300, 19.3359, -5.7570, -3.5853, -1.3882, -2.3313, -0.4930,  SiH-
     *   5.4800, 11.0459, -0.8574,  2.3137, -4.6777,  4.4532, -1.1716,  CaF
     *  17.8875, 47.4921, -1.1390,  2.7534, -7.2248,  6.3242, -1.1381,  CaO2H2
     *   8.1892, 23.3129, -1.0581,  3.5131,-11.3115, 16.9078, -9.8867/  KOH
      DATA C03/
     *  10.3100, 21.7682, -5.8992, -3.8627, -4.0284,  1.2924, -2.5856,  CN-
     *  16.1405, 37.9519, -0.0230,  0.6639, -2.4910,  5.5385, -4.2945,  Al2O2
     *   9.0621, 23.3478, -2.1422,  1.7058, -1.6807, 10.3429,-14.0183,  BaOH
     *   8.6837, 23.1042, -1.2656,  3.2436, -7.2017,  6.5067, -1.7129,  SrOH
     *  13.9839, 25.6721, -0.0784,  0.0544, -0.2755,  0.6140, -0.3673,  BO2
     *   5.5700, 12.0158, -0.5187, -0.1216,  0.6738, -0.6377,  0.1588,  SiF
     *   0.0000, 16.4621,-13.8562,-13.1896, -9.2577, -6.3354, -2.5704,  CH-
     *  13.8610, 26.3081, -1.3134,  0.1185, -0.0461, -0.4056,  0.8088,  C3
     *   8.4800, 21.1413, -5.8697, -3.3745, -2.7491, -1.8902, -0.2441,  C2-
     *  17.1545, 48.1845, -0.5683,  0.1125, -3.0973,  4.3727, -2.1978,  MgO2H2
     *   9.3961, 23.7967, -0.6500,  0.2061, -1.9381,  2.1259, -0.6451,  BeOH
     *  10.4305, 24.8357, -0.4930, -0.4550,  0.8862, -2.7257,  2.4025,  HBS
     *  13.1966, 25.7392,  0.0961, -0.7979, -0.1515,  4.2750, -4.6336,  SiC2
     *  17.4231, 48.8561, -0.4831,  0.9575, -1.9798, -0.0476,  1.2346,  FeO2H2
     *  10.0930, 25.0689, -1.5784,  2.2605, -3.1152,  3.7375, -2.5596,  CrO2
     *  20.0817, 49.3051, -0.2203,  0.6123, -1.9159,  3.0362, -0.6588,  BeH2O2
     *  11.4541, 36.8342, -1.3068, -1.2283, -0.7130, -0.1039,  0.8121,  BH3
     *  12.5346, 24.2744, -0.4230,  2.1003, -7.6565, 14.5171,-10.4377,  NaCN
     *   6.5483, 23.5736, -0.7830, -0.0881, -2.2398,  2.7050, -1.5244,  BeH2
     *  10.1248, 24.8268, -0.3784,  0.5561, -0.7324,  1.7508, -1.6977,  Si2N
     *   9.3132, 22.5681, -0.7730,  3.2979, -6.3686,  5.5210, -0.9987,  CaCl2
     *  18.8913, 37.0212, -0.3881,  1.7934, -7.5472, 14.9782,-11.0505,  NaBO2
     *   0.0000, 19.8338,-46.6804,-50.9308,-35.9059,-13.5611,-23.8103,  C3H
     *  15.5315, 26.0301, -0.1824,  0.0109, -0.3944,  0.5184, -0.0882,  OBF
     *  11.9993, 26.2368, -0.1708,  0.2491,  0.4220, -2.2962,  2.2409,  CS2
     *   8.9381, 23.5703, -0.6263,  1.0060, -4.3983,  7.4665, -4.8955,  LiOH
     *   1.5500, 11.3681, -0.1946, -0.0669, -2.3347,  5.3477, -4.0343,  Al2
     *   4.8400, 11.3090, -0.5602,  0.5886, -3.9705,  7.3873, -5.2571,  LiCl
     *  11.3225, 25.4462, -1.0487,  1.8142, -1.5110,  0.4282, -0.0240,  TiOCl
     *  23.3326, 62.7915, -1.3095, -1.6903, -0.9624, -1.6171,  2.5521,  C2H4
     *   7.4689, 23.8059, -0.5629,  0.0019, -0.3896, -0.7781,  0.3890,  CHCl
     *   6.6900, 14.8883,  5.3193,  8.9551,  3.7271,  5.1452,  1.0391,  TiCl
     *  19.2284, 37.1933,  0.1308, -0.0614, -0.9981,  2.9770, -2.1833,  AlOF2
     *  18.9713, 36.8674, -0.8338,  3.8816,-11.3916, 16.8414, -9.6911,  KBO2
     *  11.2271, 25.9412,  0.1074, -0.8813, -0.2594,  4.4112, -4.4861,  Si2C
     *   9.2183, 24.5270, -0.6453, -1.0757, -0.7155,  2.2944, -1.4513,  CHF
     *   0.0000, 11.8175,-29.4442,-30.6402,-22.9279,-13.1209, -8.8023,  BO-
     *  10.9760, 27.6834,  5.5082,  6.6402,  5.5692,  2.7324,  1.9375,  AlO2
     *  18.0802, 47.0050, -2.3587,  2.3466, -2.2753,  8.4432,-11.3032,  BaO2H2
     *  12.8526, 25.8889, -1.0260,  1.8361, -1.5017,  0.3478,  0.0486,  OTiF
     *   6.5000, 20.6745, -7.9942, -5.7057, -2.6759, -6.1649,  1.2656,  CS-
     *  21.5636, 39.0495, -0.1190,  0.7088, -1.5184,  0.4914,  0.9277,  C2N2
     *  17.5958, 46.9386, -1.3295,  3.5725, -8.4710,  7.5694, -1.8456,  SrO2H2
     *  12.2076, 25.3442, -0.0379, -0.1189, -0.8276,  1.3188, -0.6986,  ClCN
     *  10.6135, 23.6489, -0.5207,  0.0519, -0.6538,  1.9149, -1.5058,  AlClF
     *  12.5010, 24.1386, -0.8692,  4.1888,-11.7377, 17.1662, -9.8522,  KCN
     *   8.8688, 23.5425, -0.5528,  0.0031, -0.7346,  2.3344, -1.9878,  AlCl2
     *   9.6070, 22.2204, -2.5275,  2.8555, -1.4987,  7.7865,-11.3039,  BaCl2
     *  12.3143, 24.3964, -0.4940,  0.0699, -0.5475,  1.6261, -1.2695,  AlF2
     *   8.1536, 22.9187, -0.1815,  0.6847, -2.4792,  4.3296, -2.7691/  MgCl2
      DATA C04/
     *   0.0000, 17.5598,-16.6727,-14.0707,-13.0780, -5.4193, -4.7856,  FeO-
     *  20.4537, 49.9913, -0.5362, -0.7176, -1.2169,  1.1206, -0.3773,  BO2H2
     *  14.1133, 48.5194, -0.8436, -1.0629, -0.7362,  0.3080, -0.3403,  SiH3Cl
     *   8.3239, 23.6272, -0.2108,  1.1105, -2.1105,  1.5380, -0.1684,  FeCl2
     *   7.3840, 24.8600, -0.1499, -0.1631,  0.1378,  1.6604, -1.9986,  Si3
     *  16.1268, 48.9782, -0.8260, -1.0380, -0.6452, -0.1029,  0.1199,  SiH3F
     *  16.2992, 49.7196, -1.2716, -1.4752, -1.1626,  0.6516, -0.0837,  CH3Cl
     *   9.1791, 22.1133, -1.4891,  4.1050, -7.6534,  6.6694, -1.5355,  SrCl2
     *  11.6845, 23.2600, -1.2039,  3.3661, -6.2828,  5.1661, -0.6547,  CaF2
     *  13.7563, 25.2856, -0.4137,  1.0746, -1.1248,  0.2935,  0.3807,  TiF2
     *  19.4163, 36.9346, -0.3977,  1.3814, -4.7577,  8.2956, -5.5779,  LiBO2
     *   9.5422, 23.6489, -0.6541,  0.7042, -2.5258,  4.5411, -3.0359,  MgClF
     *  19.3953, 37.4967, -0.4103,  0.6249, -2.5737,  3.7334, -2.0769,  BeBO2
     *  16.1988, 37.8077, -0.3545, -0.2428, -0.1731, -1.4896,  1.9844,  C2HCl
     *   9.9277, 24.6274, -0.5062,  0.9860, -1.3100,  0.8075, -0.0931,  TiCl2
     *  19.7168, 40.3256, -0.2533,  0.3731, -0.5863, -0.6939,  0.9337,  C4
     *  30.6562, 75.8041, -1.6269, -1.1205, -1.8109,  2.1354, -0.8357,  H3BO3
     *  10.7510, 23.8686, -0.6130,  0.7434, -2.6657,  5.0507, -3.5509,  MgF2
     *   0.0000, 13.8534,-28.5088,-27.6557,-25.0420, -4.2145,-21.0916,  BaClF
     *  13.3200, 24.6323, -0.2099,  0.5174, -1.9085,  2.9836, -1.7351,  BeF2
     *  16.6788, 38.1093, -0.3632, -0.2642, -0.4287, -0.5573,  0.9863,  C2HF
     *   9.6498, 23.7877, -0.2606,  0.4816, -1.7048,  2.1226, -0.8176,  BeCl2
     *  15.7352, 37.1910, -1.0480,  1.8371, -1.1420, -0.7526,  1.2880,  TiOCl2
     *  10.7683, 24.3508, -0.5859,  0.0972, -0.3635,  0.9082, -0.3338,  ZrCl2
     *  11.9101, 22.9073, -2.4413,  2.9420, -1.3655,  7.3312,-10.8692,  BaF2
     *  12.4073, 25.2586, -0.5256,  0.7548, -2.0655,  2.2598, -0.9944,  BeC2
     *   9.9676, 24.0020, -0.4765,  1.0925, -3.6131,  4.2582, -1.8225,  Be2O
     *  11.3542, 22.8132, -1.4157,  4.1790, -7.3508,  5.5696, -0.4507,  SrF2
     *  13.7587, 24.7160, -1.0103,  0.2376, -0.4664, -0.9114,  6.9672,  ZrF2
     *  13.0910, 27.6502,  6.5468,  8.2502,  7.3334,  4.1191,  1.2402,  FeF2
     *  12.5389, 37.9053, -1.3490,  3.1985, -1.1165, -6.7253,  7.3584,  P4
     *  19.0240, 49.7099, -0.5565, -0.7375, -0.2251, -1.1324,  1.2457,  SiH2F2
     *   3.2806, 41.7329, 32.0127, 34.5233, 27.1981, 13.3168, 13.4808,  H3O+
     *  27.0859, 54.0398,  0.0077,  0.4169, -0.9261, -0.3135,  0.6322,  C5
     *  19.7864, 37.9176, -0.7063,  1.7895, -1.5401,  0.9448, -0.6313,  TiF3
     *  14.3199, 37.3165, -0.8450,  1.6603, -1.6009,  0.8934, -0.5070,  TiCl3
     *  15.5540, 36.5254, -0.7361,  0.8503, -0.3688,  0.0324,  0.0881,  ZrCl3
     *  10.6603, 34.6664, -0.4567,  3.2641,-13.6211, 27.6173,-20.7914,  Na2Cl2
     *  18.1954, 60.7438, -0.7643,  2.2577,-14.4187, 28.3225,-20.4866,  (NaOH)2
     *  28.8149, 64.3940, -0.2174,  1.3367, -6.6368,  8.6309, -4.6284,  Be3O3
     *  10.8345, 33.9871, -1.3140,  7.4840,-21.9583, 33.6428,-20.3143,  K2Cl2
     *  18.3196, 60.4179, -1.6298,  6.4524,-22.9230, 33.8810,-20.0092,  (KOH)2
     *  20.4364, 49.7173, -0.6667,  0.8064, -0.1308, -0.4433,  0.8970,  ZrCl4
     *  27.1266, 62.7471, -0.3813,  3.6624,-15.0927, 27.0694,-18.7738,  (NaCN)2
     *  27.0557, 51.2712, -0.5271,  0.8930, -0.5666,  1.5292, -1.3568,  ZrF4
     *  20.3442, 61.3686, -0.8410,  1.3617, -9.5297, 16.1158,-11.1739/  (LiOH)2
C
C Coefficients for constructing partition functions (and then equilibrium
C   constants, perhaps). For diatomic molecules other than H2 and CO, the
C   data are from Sauval & Tatum (1984, ApJS, 56, 193). For H2 and CO, the
C   data are from Irwin (1987, A&A, 182, 348). For polyatomic molecules,
C   the coefficients are from Irwin (1988, A&AS, 74,145).
C Coefficients used to construct the partition function, as follows:
C
C     log10(Q) = Sum{i=0,9}{PCOEF(i+1)*log10(THETA)**i}
C                                                           Ioniz. pot.
      DATA P01/
     *   1.69179,      -1.72270,       0.798033,     -0.157089,         H2
     *  -0.535313,      1.75818,      -2.63895,       1.35708,          H2
     *   0.0,           0.0,                                 15.42593,  H2
     *   3.615300,     -1.773848,      0.3516181,     0.08620792,       CO
     *   0.2911791,    -1.141469,      2.513133,     -2.886502,         CO
     *   1.238932,      0.0,                                 14.01400,  CO
     *   4.344711818,  -3.6343233,     1.415963,      0.01594,          H2O
     *   0.56542,      -1.2583,        0.53796,       3*0.0, 12.62100,  H2O
     *   3.0929, -1.6778,  0.6743, -0.1874,  0.0000,  5*0.0, 13.01700,  OH
     *   3.2643, -1.7303,  0.4192,  0.0000,  0.0000,  5*0.0, 15.58100,  N2
     *   4.2275, -1.9144,  0.7201, -1.3099,  1.1657,  5*0.0, 11.49000,  SiO
     *  1.0, 9*0.,                                           10.42200,  HS
     *   5.117210341,  -3.94844146,    1.23193,       0.076156,         H2S
     *   0.42163,      -0.453534,      0.0,           3*0.0, 10.45700,  H2S
     *   3.0735, -1.8501,  0.9607, -0.3935,  0.0000,  5*0.0, 13.49000,  NH
     *   3.6908, -1.9801,  0.7704, -0.2247,  0.0000,  5*0.0,  7.91000,  SiH
     *   3.3586, -2.0656,  0.9624, -0.2239,  0.0000,  5*0.0, 10.64000,  CH
     *   2.5410, -2.4336,  1.4979,  0.0192, -0.7483,  5*0.0, -1.00000,  H2+
     *   4.3073, -1.8255,  0.3765,  0.0000,  0.0000,  5*0.0,  9.26420,  NO
     *   3.6704, -2.2682,  0.9354, -0.2597,  0.0000,  5*0.0,  7.20000,  MgH
     *   2.8005, -1.7476,  0.5310,  0.0000,  0.0000,  5*0.0, 12.74400,  HCl
     *   4.8026, -1.9753,  0.2600,  0.0000,  0.0000,  5*0.0, 10.53000,  SiS
     *   6.103792598,  -4.3938712,     0.662588,      0.3751,           AlOH
     *   0.38386,      -0.2147,        0.0,           3*0.0, -1.00000,  AlOH
     *   4.819621858,  -3.84200734,    1.5386462,     0.784399,         NH2
     *  -2.34404,       2.50803,      -1.13304,       3*0.0, 11.14000,  NH2
     *   3.3209, -2.5909,  1.7415, -0.7636,  0.0000,  5*0.0,  5.50000,  AlH
     *   4.0078, -2.1514,  0.9226, -0.1671,  0.0000,  5*0.0, 13.59800,  CN
     *   6.01081285,   -4.438833,      0.840462,      0.2945,           CO2
     *   0.3694,       -0.273,         0.0,           3*0.0, 13.77700,  CO2
     *   4.7963, -2.1308,  0.5224,  0.0000,  0.0000,  5*0.0, 10.29400,  SO
c    *   5.7765, -2.3739,  0.8940, -0.3641,  0.0000,  5*0.0,  6.40000,  TiO
     *   5.3051, -2.3739,  0.8940, -0.3641,  0.0000,  5*0.0,  6.40000,  TiO
     *   5.0796, -2.1967,  0.4101,  0.0000,  0.0000,  5*0.0,  9.35600,  S2
     *   4.6265980,    -2.5625800,     0.38885943,    0.40219820,       FeH
     *  -0.21386399,    0.027845045,   0.0,           3*0.0,  7.37000,  FeH
     *   5.884176216,  -5.8364867,     1.608417,      1.50876,          NH3
     *  -0.59607,      -0.58961,       0.2459,        3*0.0, -1.00000,  NH3
     *   5.434042379,  -4.2409874,     0.988745,      0.49464,          HCN
     *   0.03719,      -0.22924,       0.0,           3*0.0, 13.60000,  HCN
     *   6.298781639,  -3.85672804,    0.8551678,     0.321901,         HCO
     *   0.020274,      0.15254,      -0.25298,       3*0.0,  8.12000,  HCO
     *   4.0636, -2.0779,  0.7660, -0.2111,  0.0000,  5*0.0, 12.06970,  O2
     *  1.0, 9*0.,                                           10.39600,  CH2
     *   2.4164, -1.6132,  0.6357, -0.1767,  0.0000,  5*0.0, 16.03000,  HF
     *  1.0, 9*0.,                                           -1.00000,  H3+
     *   3.8411, -2.3891,  1.3578, -0.6893,  0.0000,  5*0.0,  5.86000,  CaH
     *  1.0, 9*0.,                                           -1.00000,  Al2O
     *   4.9191, -2.6291,  0.5831,  0.3163,  0.0000,  5*0.0,  9.46000,  AlO
     *  1.0, 9*0.,                                            9.84000,  CH3
     *  1.0, 9*0.,                                            8.80000,  SiH2
     *   5.3182, -2.6502, -0.2781, -0.7823,  1.3107,  5*0.0,  8.76000,  MgO
     *   4.3091, -2.2406,  0.4865, -0.2049,  0.0000,  5*0.0, 11.40000,  C2
     *  1.0, 9*0.,                                            9.50000,  TiO2
     *   8.457240767,  -4.1987868,     0.334575,      0.20744,          VO2
     *   0.18226,      -0.053465,      0.0,           3*0.0, -1.00000,  VO2
     *   3.5453, -2.3457,  0.8557, -0.1685,  0.0000,  5*0.0,  4.70000,  NaH
     *   5.1115, -2.2303,  0.8001, -0.5192,  0.0000,  5*0.0,  9.40000,  AlCl
     *   4.5405, -2.1033,  0.6208, -0.2930,  0.0000,  5*0.0, -1.00000,  AlF
     *   5.0687, -2.2186,  0.9545, -0.4592,  0.0000,  5*0.0,  7.23860,  VO
     *   4.1646, -1.9348,  0.8034, -1.3669,  1.1561,  5*0.0, 11.33000,  CS
     *   6.8401894714, -4.338616427,   0.71600166,    0.128126,         MgOH
     *   0.5978087,    -0.8658369,     0.385049,      3*0.0,  7.50000,  MgOH
     *  1.0, 9*0.,                                           11.90000,  PO2
     *   7.1623971155, -4.471282563,   1.1221899,    -0.558812,         CaOH
     *   0.2294,        1.78658,      -2.95118,       1.41591,          CaOH
     *   2*0.0,                                               5.80000,  CaOH
     *  1.0, 9*0.,                                            9.82400/  PH2
      DATA P02/
     *  1.0, 9*0.,                                           11.61000,  C2H
     *   4.8065, -2.2129,  0.9991, -0.5414,  0.0000,  5*0.0, -1.00000,  ScO
     *  1.0, 9*0.,                                           -1.00000,  AlO2H
     *   5.2461, -2.1319,  0.5340, -0.2309,  0.0000,  5*0.0, -1.00000,  AlS
     *   5.5642, -2.1947,  0.5065,  0.0000,  0.0000,  5*0.0,  8.90000,  FeO
     *   5.5270, -2.1311,  0.6523, -0.2533,  0.0000,  5*0.0,  7.85000,  CrO
     *  1.0, 9*0.,                                           12.61000,  CH4
     *   4.8052, -1.9619,  0.3140,  0.0000,  0.0000,  5*0.0,  8.87000,  NS
     *  1.0, 9*0.,                                           12.34900,  SO2
     *   4.6570, -2.3587,  0.8819, -0.1642,  0.0000,  5*0.0, -1.00000,  SiN
     *  1.0, 9*0.,                                           -1.00000,  OH-
     *   5.3279, -2.4694,  0.2164, -0.2313,  0.0000,  5*0.0,  6.00000,  ZrO
     *   3.5649, -1.7328,  0.4241,  0.0000,  0.0000,  5*0.0, -1.00000,  NO+
     *   8.72011985,   -4.247295,      0.2758,        0.20738,          ZrO2
     *   0.09406,       0.0,           0.0,           3*0.0, -1.00000,  ZrO2
     *   3.9953, -1.8665,  0.5965, -0.1617,  0.0000,  5*0.0, 13.30000,  BO
     *  1.0, 9*0.,                                           -1.00000,  SiO2
     *  1.0, 9*0.,                                           -1.00000,  HBO
     *   5.1477, -1.8671,  0.2404,  0.0000,  0.0000,  5*0.0,  9.20000,  SiC
     *  1.0, 9*0.,                                           -1.00000,  YO2
     *   5.8948, -2.2183,  0.5928, -0.3106,  0.0000,  5*0.0,  7.10000,  TiS
     *  1.0, 9*0.,                                           -1.00000,  HBO2
     *   7.1220464309, -6.966653604,   1.9668235,     0.362597,         C2H2
     *   0.608996,     -0.920435,      0.271892,      3*0.0, 11.40000,  C2H2
     *  1.0, 9*0.,                                           11.18500,  OCS
     *  1.0, 9*0.,                                           -1.00000,  ZrO+
     *  1.0, 9*0.,                                           -1.00000,  NaOH
     *   5.7494, -2.3340,  0.8685, -0.5306,  0.0000,  5*0.0,  5.86000,  CaCl
     *  1.0, 9*0.,                                           -1.00000,  AlOF
     *   4.9515, -2.0866,  0.6565, -0.3082,  0.0000,  5*0.0,  6.00000,  YO
     *   5.3364, -2.2844,  0.2820,  0.1185,  0.0000,  5*0.0, -1.00000,  NaCl
     *  1.0, 9*0.,                                           -1.00000,  C2O
     *  1.0, 9*0.,                                           10.79000,  CHP
     *  1.0, 9*0.,                                           -1.00000,  HS-
     *  1.0, 9*0.,                                           -1.00000,  H2-
     *  1.0, 9*0.,                                            6.00000,  TiH
     *  1.0, 9*0.,                                            9.86900,  PH3
     *   5.0367, -2.1625,  0.4859, -0.1780,  0.0000,  5*0.0, -1.00000,  MgS
     *  1.0, 9*0.,                                           -1.00000,  TiO+
     *  1.0, 9*0.,                                           -1.00000,  LaO2
     *   5.2617, -2.1485,  0.5647, -0.2985,  0.0000,  5*0.0, -1.00000,  Si2
     *  1.0, 9*0.,                                           -1.00000,  SiH4
     *  1.0, 9*0.,                                            9.80000,  BH2
     *  1.0, 9*0.,                                           -1.00000,  AlOCl
     *   5.1147, -2.5016,  1.0445, -0.3135,  0.0000,  5*0.0,  4.95000,  LaO
     *  1.0, 9*0.,                                           12.00000,  C2N
     *  1.0, 9*0.,                                           -1.00000,  AlBO2
     *   5.6860, -2.3016,  0.2086,  0.1763,  0.0000,  5*0.0, -1.00000,  KCl
     *  1.0, 9*0.,                                           -1.00000,  SiH-
     *   5.2010, -2.2653,  0.8941, -0.5384,  0.0000,  5*0.0, -1.00000,  CaF
     *  1.0, 9*0.,                                           -1.00000,  CaO2H2
     *  1.0, 9*0.,                                            7.50000/  KOH
      DATA P03/
     *  1.0, 9*0.,                                           -1.00000,  CN-
     *  1.0, 9*0.,                                           -1.00000,  Al2O2
     *  1.0, 9*0.,                                           -1.00000,  BaOH
     *  1.0, 9*0.,                                           -1.00000,  SrOH
     *  1.0, 9*0.,                                           -1.00000,  BO2
     *   5.0871, -2.0375,  0.4478, -0.1243,  0.0000,  5*0.0,  7.54000,  SiF
     *  1.0, 9*0.,                                           -1.00000,  CH-
     *   6.618407932,  -3.576399,      0.883642,      0.087548,         C3
     *   0.04817,      -0.16471,       0.0,           3*0.0, -1.00000,  C3
     *  1.0, 9*0.,                                           -1.00000,  C2-
     *  1.0, 9*0.,                                           -1.00000,  MgO2H2
     *  1.0, 9*0.,                                           -1.00000,  BeOH
     *  1.0, 9*0.,                                           -1.00000,  HBS
     *   7.54651307623,-5.075563869,   1.82960795,    0.0983258,        SiC2
     *  -6.335157,     14.33103,     -13.01689,       4.428233,         SiC2
     *   2*0.0,                                              10.20000,  SiC2
     *  1.0, 9*0.,                                           -1.00000,  FeO2H2
     *  1.0, 9*0.,                                           -1.00000,  CrO2
     *  1.0, 9*0.,                                           -1.00000,  BeH2O2
     *  1.0, 9*0.,                                           -1.00000,  BH3
     *  1.0, 9*0.,                                           -1.00000,  NaCN
     *  1.0, 9*0.,                                           -1.00000,  BeH2
     *  1.0, 9*0.,                                           -1.00000,  Si2N
     *  1.0, 9*0.,                                           -1.00000,  CaCl2
     *  1.0, 9*0.,                                           -1.00000,  NaBO2
     *  1.0, 9*0.,                                           -1.00000,  C3H
     *  1.0, 9*0.,                                           -1.00000,  OBF
     *  1.0, 9*0.,                                           10.07300,  CS2
     *  1.0, 9*0.,                                           -1.00000,  LiOH
     *   5.5538, -2.3365,  0.5754, -0.2119,  0.0000,  5*0.0,  5.40000,  Al2
     *   4.5605, -2.2216,  0.5760, -0.1706,  0.0000,  5*0.0,  9.57000,  LiCl
     *  1.0, 9*0.,                                           -1.00000,  TiOCl
     *  1.0, 9*0.,                                           -1.00000,  C2H4
     *  1.0, 9*0.,                                           -1.00000,  CHCl
     *  1.0, 9*0.,                                           -1.00000,  TiCl
     *  1.0, 9*0.,                                           -1.00000,  AlOF2
     *  1.0, 9*0.,                                           -1.00000,  KBO2
     *  1.0, 9*0.,                                           -1.00000,  Si2C
     *  1.0, 9*0.,                                           10.06000,  CHF
     *  1.0, 9*0.,                                           -1.00000,  BO-
     *  1.0, 9*0.,                                           -1.00000,  AlO2
     *  1.0, 9*0.,                                           -1.00000,  BaO2H2
     *  1.0, 9*0.,                                           -1.00000,  OTiF
     *  1.0, 9*0.,                                           -1.00000,  CS-
     *  1.0, 9*0.,                                           -1.00000,  C2N2
     *  1.0, 9*0.,                                           -1.00000,  SrO2H2
     *  1.0, 9*0.,                                           12.36000,  ClCN
     *  1.0, 9*0.,                                           -1.00000,  AlClF
     *  1.0, 9*0.,                                           -1.00000,  KCN
     *  1.0, 9*0.,                                           -1.00000,  AlCl2
     *  1.0, 9*0.,                                           -1.00000,  BaCl2
     *  1.0, 9*0.,                                           -1.00000,  AlF2
     *  1.0, 9*0.,                                           -1.00000/  MgCl2
      DATA P04/
     *  1.0, 9*0.,                                           -1.00000,  FeO-
     *  1.0, 9*0.,                                           -1.00000,  BO2H2
     *  1.0, 9*0.,                                           -1.00000,  SiH3Cl
     *  1.0, 9*0.,                                           -1.00000,  FeCl2
     *  1.0, 9*0.,                                           -1.00000,  Si3
     *  1.0, 9*0.,                                           -1.00000,  SiH3F
     *  1.0, 9*0.,                                           -1.00000,  CH3Cl
     *  1.0, 9*0.,                                           -1.00000,  SrCl2
     *  1.0, 9*0.,                                           -1.00000,  CaF2
     *  1.0, 9*0.,                                           -1.00000,  TiF2
     *  1.0, 9*0.,                                           -1.00000,  LiBO2
     *  1.0, 9*0.,                                           -1.00000,  MgClF
     *  1.0, 9*0.,                                           -1.00000,  BeBO2
     *  1.0, 9*0.,                                           -1.00000,  C2HCl
     *  1.0, 9*0.,                                           -1.00000,  TiCl2
     *  1.0, 9*0.,                                           -1.00000,  C4
     *  1.0, 9*0.,                                           -1.00000,  H3BO3
     *  1.0, 9*0.,                                           -1.00000,  MgF2
     *  1.0, 9*0.,                                           -1.00000,  BaClF
     *  1.0, 9*0.,                                           -1.00000,  BeF2
     *  1.0, 9*0.,                                           -1.00000,  C2HF
     *  1.0, 9*0.,                                           -1.00000,  BeCl2
     *  1.0, 9*0.,                                           -1.00000,  TiOCl2
     *  1.0, 9*0.,                                           -1.00000,  ZrCl2
     *  1.0, 9*0.,                                           -1.00000,  BaF2
     *  1.0, 9*0.,                                           -1.00000,  BeC2
     *  1.0, 9*0.,                                           -1.00000,  Be2O
     *  1.0, 9*0.,                                           -1.00000,  SrF2
     *  1.0, 9*0.,                                           -1.00000,  ZrF2
     *  1.0, 9*0.,                                           -1.00000,  FeF2
     *  1.0, 9*0.,                                           -1.00000,  P4
     *  1.0, 9*0.,                                           -1.00000,  SiH2F2
     *  1.0, 9*0.,                                           -1.00000,  H3O+
     *  1.0, 9*0.,                                           -1.00000,  C5
     *  1.0, 9*0.,                                           -1.00000,  TiF3
     *  1.0, 9*0.,                                           -1.00000,  TiCl3
     *  1.0, 9*0.,                                           -1.00000,  ZrCl3
     *  1.0, 9*0.,                                           -1.00000,  Na2Cl2
     *  1.0, 9*0.,                                           -1.00000,  Na2O2H2
     *  1.0, 9*0.,                                           -1.00000,  Be3O3
     *  1.0, 9*0.,                                           -1.00000,  K2Cl2
     *  1.0, 9*0.,                                           -1.00000,  K2O2H2
     *  1.0, 9*0.,                                           -1.00000,  ZrCl4
     *  1.0, 9*0.,                                           -1.00000,  Na2C2N2
     *  1.0, 9*0.,                                           -1.00000,  ZrF4
     *  1.0, 9*0.,                                           -1.00000/  Li2O2H2
C
C Try to find the input speicies name (SPNAME) in the list (SPLIST) of
C species for which we have equilibrium constant coefficients. Note that
C the index is stored in a new variable J, rather than using the loop
C variable I, because some optimizers don't save the loop variable after
C normal termination of the loop.
C
      DO 1 I=1,MSPEC
      J=I
      IF(SPLIST(J).EQ.SPNAME) GO TO 2
   1  CONTINUE
C
C Fall through to here, if requested molecule was not in SPLIST.
C Print a warning, but return anyway.
C
      WRITE(*,*) 'MOLCON: Don''t have dissociation constant for ',
     *           'molecule: "', SPNAME, '"'
      EQK =1.D20
      PART=1.D0
      RETURN
C
C Calculate independent variable for polynomial expansions.
C Note that the polynomial expansions in Sauval & Tatum (1984) and Irwin
C (1987,1988) are in terms of log10(5040/T), not log10(5039.7475/T), but
C the more accuate value of 5039.7475 should be used in converting the
C partition function into an equilibrium constant.
C
   2  TH=5040.D0/T
      LOGTH=LOG10(TH)
C
C Construct equilibrium constant from polynomial coefficients and
C dissociation constant. A "+1" term at the end would convert from
C pascals (i.e. N/m/m as in Sauval) to dynes/cm/cm.
C
c      if (t.lt.1600) logth=log10(5040.0/1600.0)
c      if (t.gt.7730) logth=log10(5040.0/7730.0)
      EQK=COEF(2,J)+LOGTH*(COEF(3,J)+LOGTH*(COEF(4,J)+
     &              LOGTH*(COEF(5,J)+LOGTH*(COEF(6,J)+
     &              LOGTH*(COEF(7,J))))))
     &             -TH*COEF(1,J)
C    &             +1.0D0
      EQK=10.D0**EQK
C
C Just for the reference, the relation between partition functions
C and equilibrium constant:
C
C            P(A)*P(B)*...      N(A)*N(B)*...
C K(AB...) = ------------- = kT-------------- =
C              P(AB...)           N(AB...)
C
C             2*pi*kT 3/2    M(A)*M(B)*... 3/2   Q(A)*Q(B)*...
C       = kT*(-------)    * (-------------)    * ------------- * exp(-D(AB)/kT)
C               h^2            M(AB...)           Q(AB...)
C
C where, K - equilibrium constant, Q - partition functions, M - masses
C        P - partial pressures, N - number densities, T - temperature,
C        D - complete dissociation energy, h - plank constant. Remember
C        to use masses in grams (1 amu = 1.660540E-24 g) and energy in
C        ergs (1 eV = 1.60219E-12 ergs). Also, k = 1.38065E-16 erg/K,
C        h = 6.626076E-27 erg s, and pi = 3.1415926536.
C
C Construct partition function from polynomial coefficients.
C
      PART=PCOEF(NPCOEF-1,J)
      DO 3 I=NPCOEF-2,1,-1
    3 PART=LOGTH*PART+PCOEF(I,J)
C
C Copy ionization potential
C
      PION=PCOEF(NPCOEF,J)
C
C Calculate equilibrium constant (EQK) from partition function, dissociation
C constant, and other information passed into subroutine. The constants used
C are:  79.733501 = 1.5*log10(2*pi/h/h)  [in cgs units] and
C      -15.859914 = alog10(k)            [in cgs units].
C       5039.7475 = alog10(e)*k*(eV/erg)
C
      EQK_ST=(NTOT-1)*(79.733501D0+2.5D0*(LOG10(T)-15.859914D0))+
     &       1.5D0*RATIOM+QPRD-PART-COEF(1,J)*5039.7475D0/T
C
C Convert equilibrium constant and partition function from logarithms.
C
      EQK_ST=10.D0**EQK_ST
      PART=10.D0**PART
C
C Check if there is relevant data in Paul Barklem's tables
C
      CALL KP_Q_SPLN(SPNAME,T,Qm_spln,Kp_spln)
c      write(*,'(F10.1,A9,5G13.6)') T,SPNAME,EQK,EQK_ST,Kp_spln,
c     &                             PART,Qm_spln
      IF(Kp_spln.GE.0.d0) THEN
        EQK =Kp_spln
        PART=Qm_spln
      ENDIF

c      if(spname.eq.'H2+') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'NO') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'C3') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'SiN') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'SiS') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'TiO') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'N2') write(*,'(a,f10.2,1p3e14.6,i3,1p2e14.6)')
c     &   spname,t , eqk, eqk_st, part,NTOT,QPRD,RATIOM
c
c Don't use EQK_ST based on partition function - use direct fit to EQK.
c
c      EQK=EQK_ST
C
C Done.
C
      RETURN
      END


C=========================================================================
C MOLCON: Returns equilibrium constant and partition function for a given
C   molecule and temperature.
C
C Inputs:
C   SPNAME [character(*)] Name of molecule, chosen from SPLIST below.
C   T [real] Temperature (in K) at which EQK and PART are to be found.
C   NTOT [real] Total number of atoms in the molecule.
C   RATIOM [real] Logarithm (base 10) of mass ratio (in g^(natoms-1)):
C     ratiom = Sum{log10(Atomic Masses)} - log10(Sum{Atomic Masses})
C   QPRD [double] Logarithm of product of atomic partition functions:
C     qprd = Sum{log10(Atomic Partition Functions)}
C
C Outputs:
C   EQK [real] Equilibrium constant (in dynes/cm/cm) at temperature T,
C     calculated from dissociation energy and partition function.
C   PART [real] Partition function at temperature T, calculated from
C     expressions in the references cited below.
C
C References:
C   For diatomic molecules: Sauval & Tatum (1984, ApJS, 56, 193).
C
      SUBROUTINE MOLCON2(SPNAME,T,NTOT,RATIOM,QPRD,EQK,PART,PION)
C
      INTEGER MSPEC,NTOT
      DOUBLE PRECISION KERG,KEV
      DOUBLE PRECISION RATIOM,QPRD,PION
      PARAMETER (KERG=1.38065D-16,KEV=KERG/1.60219D-12)
      PARAMETER (CONST=25947.256)
C
      REAL T
      DOUBLE PRECISION TH,LOGTH,EQK_ST,EQK,PART,Qm_spln,Kp_spln
C
C Combine equilibrium constant coefficients into one large array.
C
      PARAMETER (MSPEC=196)
      PARAMETER (NEQCOE=7)
      DOUBLE PRECISION COEF(NEQCOE,MSPEC)
      DOUBLE PRECISION C01(NEQCOE,50),C02(NEQCOE,50),
     *                 C03(NEQCOE,50),C04(NEQCOE,46)
      EQUIVALENCE (C01(1,1),COEF(1,  1)),(C02(1,1),COEF(1, 51))
      EQUIVALENCE (C03(1,1),COEF(1,101)),(C04(1,1),COEF(1,151))
C
C Combine partition function coefficients into one large array.
C
      PARAMETER (NPCOEF=11)
      DOUBLE PRECISION PCOEF(NPCOEF,MSPEC)
      DOUBLE PRECISION P01(NPCOEF,50),P02(NPCOEF,50),
     *                 P03(NPCOEF,50),P04(NPCOEF,46)
      EQUIVALENCE (P01(1,1),PCOEF(1,  1)),(P02(1,1),PCOEF(1, 51))
      EQUIVALENCE (P03(1,1),PCOEF(1,101)),(P04(1,1),PCOEF(1,151))
C
      CHARACTER SPNAME*(*),SPLIST(MSPEC)*8
      SAVE
C
C Molecular species list from NextGen models (Allard & Hauschildt).
C See old/eos.4.f for molecular species list from Sauval & Tatum (1984).
C
      DATA SPLIST/
     * 'H2      ','CO      ','H2O     ','OH      ','N2      ',
     * 'SiO     ','HS      ','H2S     ','NH      ','SiH     ',
     * 'CH      ','H2+     ','NO      ','MgH     ','HCl     ',
     * 'SiS     ','AlOH    ','NH2     ','AlH     ','CN      ',
     * 'CO2     ','SO      ','TiO     ','S2      ','FeH     ',
     * 'NH3     ','HCN     ','HCO     ','O2      ','CH2     ',
     * 'HF      ','H3+     ','CaH     ','Al2O    ','AlO     ',
     * 'CH3     ','SiH2    ','MgO     ','C2      ','TiO2    ',
     * 'VO2     ','NaH     ','AlCl    ','AlF     ','VO      ',
     * 'CS      ','MgOH    ','PO2     ','CaOH    ','PH2     ',
     * 'C2H     ','ScO     ','AlO2H   ','AlS     ','FeO     ',
     * 'CrO     ','CH4     ','NS      ','SO2     ','SiN     ',
     * 'OH-     ','ZrO     ','NO+     ','ZrO2    ','BO      ',
     * 'SiO2    ','HBO     ','SiC     ','YO2     ','TiS     ',
     * 'HBO2    ','C2H2    ','OCS     ','ZrO+    ','NaOH    ',
     * 'CaCl    ','AlOF    ','YO      ','NaCl    ','C2O     ',
     * 'CHP     ','HS-     ','H2-     ','TiH     ','PH3     ',
     * 'MgS     ','TiO+    ','LaO2    ','Si2     ','SiH4    ',
     * 'BH2     ','AlOCl   ','LaO     ','C2N     ','AlBO2   ',
     * 'KCl     ','SiH-    ','CaF     ','CaO2H2  ','KOH     ',
     * 'CN-     ','Al2O2   ','BaOH    ','SrOH    ','BO2     ',
     * 'SiF     ','CH-     ','C3      ','C2-     ','MgO2H2  ',
     * 'BeOH    ','HBS     ','SiC2    ','FeO2H2  ','CrO2    ',
     * 'BeH2O2  ','BH3     ','NaCN    ','BeH2    ','Si2N    ',
     * 'CaCl2   ','NaBO2   ','C3H     ','OBF     ','CS2     ',
     * 'LiOH    ','Al2     ','LiCl    ','TiOCl   ','C2H4    ',
     * 'CHCl    ','TiCl    ','AlOF2   ','KBO2    ','Si2C    ',
     * 'CHF     ','BO-     ','AlO2    ','BaO2H2  ','OTiF    ',
     * 'CS-     ','C2N2    ','SrO2H2  ','ClCN    ','AlClF   ',
     * 'KCN     ','AlCl2   ','BaCl2   ','AlF2    ','MgCl2   ',
     * 'FeO-    ','BO2H2   ','SiH3Cl  ','FeCl2   ','Si3     ',
     * 'SiH3F   ','CH3Cl   ','SrCl2   ','CaF2    ','TiF2    ',
     * 'LiBO2   ','MgClF   ','BeBO2   ','C2HCl   ','TiCl2   ',
     * 'C4      ','H3BO3   ','MgF2    ','BaClF   ','BeF2    ',
     * 'C2HF    ','BeCl2   ','TiOCl2  ','ZrCl2   ','BaF2    ',
     * 'BeC2    ','Be2O    ','SrF2    ','ZrF2    ','FeF2    ',
     * 'P4      ','SiH2F2  ','H3O+    ','C5      ','TiF3    ',
     * 'TiCl3   ','ZrCl3   ','Na2Cl2  ','Na2O2H2 ','Be3O3   ',
     * 'K2Cl2   ','K2O2H2  ','ZrCl4   ','Na2C2N2 ','ZrF4    ',
     * 'Li2O2H2 '/
C
C Dissociation energy (first column, in eV) and equilibrium constant
C   coefficients. See the file "atomiz.notes" for the information on the
C   origin of the dissociation energies. The polynomial fit coefficients
C   for the equilibrium constants were determined with "ng_kfit.pro" and
C   are meant to reproduce the constants used in constructing the NextGen
C   models. The NextGen equilibrium constants were fit over the temperature
C   range 1600 < T < 7730 K. The fits are likely to diverge rapidly from
C   the truth outside this temperature range.
C Equilibrium constants may be constructed from the coefficients using:
C
C     log10(Kp) = Sum{i=2,7}{COEF(i)*log10(THETA)**(i-2)} - COEF(1)*THETA
C
      DATA C01/
     *   4.4781, 12.1354, -0.7752, -0.7821,  0.1464,  0.1603, -0.0626,  H2
     *  11.0920, 13.2368, -0.8342, -0.0477, -0.2923, -0.4557,  0.6108,  CO
     *   9.6221, 24.7774, -2.3428,  1.6868, -1.2845, -2.9925,  3.6555,  H2O
     *   4.3920, 11.8016, -0.8507, -0.5193,  0.0502, -0.3409,  0.4836,  OH
     *   9.7594, 12.8868, -0.8813,  0.2639, -1.5912,  1.5866, -0.5407,  N2
     *   8.2600, 12.9252, -0.7608, -0.3541,  1.5620, -3.5952,  2.5962,  SiO
     *   3.5500, 11.4382, -0.7816, -0.4659,  0.4314, -1.2144,  0.9648,  HS
     *   7.5946, 23.8543, -0.9525, -0.8118,  0.2051, -1.0299,  1.1555,  H2S
     *   3.4700, 11.4658, -0.7258, -0.6418, -0.0442,  0.2836, -0.1618,  NH
     *   3.0600, 11.2595, -0.6962, -0.6435,  0.6663, -0.3357, -0.4151,  SiH
     *   3.4650, 11.5333, -0.5255, -0.7105,  0.2264, -0.9271,  0.9577,  CH
     *   2.6508, 15.8052, 33.7578, 34.5956, 27.3455, 16.6214,  9.9717,  H2+
     *   6.4968, 11.9347, -0.7596,  0.0953, -0.9731,  0.8265, -0.2151,  NO
     *   1.3400, 10.2911, -0.3698, -0.0655, -2.9771,  6.1325, -4.3869,  MgH
     *   4.4336, 11.9041, -0.8281, -0.6163,  0.1580, -0.5068,  0.5164,  HCl
     *   6.4200, 12.6363, -0.7355,  0.0488,  0.8442, -2.0131,  1.3603,  SiS
     *  10.1252, 25.2575, -0.6810, -0.3051, -1.5765,  2.7536, -1.8355,  AlOH
     *   7.4400, 23.7389, -1.0179, -0.9947, -1.4353,  3.2530, -1.9224,  NH2
     *   3.0600, 11.4907, -0.4322, -0.6561, -0.5978,  2.4923, -2.4038,  AlH
     *   7.7600, 12.4438, -0.4756, -0.4909, -1.4623,  2.6823, -1.5396,  CN
     *  16.5382, 26.9571, -0.7464, -0.4921, -0.8506, -0.1365,  0.2358,  CO2
     *   5.3590, 12.3380, -0.4956, -0.2251, -0.1907, -0.2038,  0.2579,  SO
     *   6.8700, 11.9229, -1.4044,  0.7899, -0.7317, -0.0193, -0.4994,  TiO
     *   4.3693, 12.3190, -0.5050, -0.0290, -0.0266, -0.6002,  0.4572,  S2
     *   2.4100, 12.1214,  0.9438,  2.2756, -0.1086,  4.1281, -1.9952,  FeH
     *  12.1388, 36.6661, -1.4062, -0.9258, -1.6969,  0.6005,  1.2302,  NH3
     *  13.2363, 25.1318, -0.5532, -0.0850, -0.9817,  0.6676,  0.3054,  HCN
     *  11.8560, 24.6414, -0.9415, -0.1856, -0.2948, -0.1630,  0.5836,  HCO
     *   5.1156, 12.8758, -0.4856, -0.5054, -0.0776, -0.0713,  0.2369,  O2
     *   7.9400, 23.8609, -1.0762, -0.4928, -0.4092,  0.0031,  0.3761,  CH2
     *   5.8690, 12.2896, -0.9180, -0.6238,  0.1243, -0.3525,  0.4767,  HF
     *   0.0000, 18.8343, 12.4131, 11.9991,  6.8079,  8.4071,  2.6202,  H3+
     *   1.7000, 10.1982, -0.9309,  1.8315, -5.6059,  6.9571, -3.5023,  CaH
     *  10.9653, 24.8807, -0.0033,  0.4796, -1.6979,  3.5631, -2.5414,  Al2O
     *   5.2700, 12.2132, -0.5246, -0.1918, -0.6810,  1.7287, -1.5839,  AlO
     *  12.6885, 36.6540, -1.3373, -1.0064, -0.5880, -0.2362,  0.8764,  CH3
     *   0.0000, 17.8513,-15.5361,-17.6144,-13.1604, -6.4819, -5.6361,  SiH2
     *   3.5300, 10.7940,  0.0122,  1.1189, -1.8758,  2.9976, -2.7758,  MgO
     *   6.2100, 12.4672, -0.4452, -0.0100, -0.1868, -0.3860,  0.6230,  C2
     *  13.2915, 25.9340, -1.4243,  1.6519, -0.7240, -0.7271,  0.7518,  TiO2
     *  12.9619, 25.9238, -1.2927,  1.3710, -2.4073,  2.2875, -0.5486,  VO2
     *   1.8800, 10.7184, -0.3642,  0.7843, -6.5309, 13.2912, -9.9502,  NaH
     *   5.1200, 11.8277, -0.3468, -1.0735,  1.8038, -1.7748,  0.4333,  AlCl
     *   6.8900, 12.2422, -0.4905, -0.4198,  0.0242,  0.3868, -0.5765,  AlF
     *   6.4100, 12.8108, -0.5811, -0.7895, -2.6766,  8.5158, -6.9993,  VO
     *   7.3550, 12.8487, -0.7627, -0.2538,  1.5240, -4.0119,  3.0234,  CS
     *   8.0735, 23.3256, -0.5884,  0.3637, -2.4401,  3.3936, -1.7121,  MgOH
     *  11.7451, 25.2051, -0.9105,  1.0031, -0.7207, -1.1064,  1.6239,  PO2
     *   8.7035, 23.1900, -1.0964,  2.5340, -5.9823,  5.3416, -1.1946,  CaOH
     *   6.4895, 23.0863, -1.3781,  0.2539, -0.6746, -1.2341,  1.5623/  PH2
      DATA C02/
     *  12.2087, 24.9752, -0.3204, -0.5640, -0.8997,  1.6927, -0.7771,  C2H
     *   6.9600, 12.5225, -1.2695,  1.7628, -2.0543, -1.2215,  2.3706,  ScO
     *  15.6364, 37.7022, -0.5885, -0.0823, -1.7283,  3.0502, -2.0176,  AlO2H
     *   3.8400, 11.9140, -0.5187, -0.1193, -0.3886,  1.1704, -1.2299,  AlS
     *   4.2000, 12.5326, -1.0657,  1.0360, -1.5641,  0.9560, -0.3218,  FeO
     *   4.4000, 11.0587, -1.3926,  1.4461, -2.1552,  3.3409, -3.1078,  CrO
     *  17.2173, 49.9426, -0.9720, -2.4957, -0.0017, -2.3299,  3.1042,  CH4
     *   4.8000, 11.9223, -0.6951,  0.1870, -0.7158,  0.4121,  0.0296,  NS
     *  11.1405, 25.9246, -0.5809,  0.0734, -0.3333,  0.1699,  0.0529,  SO2
     *   6.6880, 14.0972,  4.2904,  4.9608,  2.9390,  3.9789,  0.8908,  SiN
     *   4.7600, 19.9888, -6.7088, -4.3846, -2.8142, -2.3004, -0.3157,  OH-
     *   7.8500, 12.4674, -1.1280,  0.0368,  0.2221,  1.1043, -1.8804,  ZrO
     *  10.8500, 17.5169, 33.0097, 36.2110, 26.7396, 15.2392, 11.4130,  NO+
     *  14.4650, 25.6324, -1.5339,  1.1586, -0.9355,  1.6114, -1.2154,  ZrO2
     *   8.2800, 12.6246, -0.6966, -0.3874,  0.2531, -0.7582,  0.5307,  BO
     *  13.0355, 26.5610, -0.2891,  0.3006, -0.4009,  0.5864, -0.4006,  SiO2
     *  12.7425, 25.2283, -0.4780, -0.3611, -0.2189, -0.2108,  0.5883,  HBO
     *   4.6400, 11.8909, -0.8762,  0.1138,  0.0665, -0.5226,  0.3331,  SiC
     *  15.2000, 25.8617, -1.4050, -0.3896,  1.0805,  2.9269, -3.7531,  YO2
     *   4.7500, 11.6628, -1.4463,  1.3742, -0.8127, -0.4623,  0.2288,  TiS
     *  19.0991, 38.4541, -0.7808, -0.4220, -0.9239,  1.0793, -0.2304,  HBO2
     *  16.9704, 37.7481, -0.2529, -1.0622, -0.1485, -0.7058,  1.1910,  C2H2
     *  14.3762, 26.3815, -0.1712,  0.1197,  0.0059, -0.9891,  1.1946,  OCS
     *   0.0000,  2.5576, -0.5567, -4.5109, -4.3690, -0.1528, -3.1319,  ZrO+
     *   8.0150, 23.3420, -0.6139,  1.4091, -6.8466, 13.0407, -9.2977,  NaOH
     *   4.0900, 10.6268, -1.1367,  2.5278, -5.6022,  4.8741, -1.1616,  CaCl
     *  12.9003, 25.5751, -0.0730,  0.2808, -1.1757,  2.3733, -1.6726,  AlOF
     *   7.2900, 12.4422, -1.3547,  1.3087,  0.1688, -5.4106,  5.1158,  YO
     *   4.2300, 11.0864, -0.4463,  1.1926, -7.5820, 15.2552,-11.1116,  NaCl
     *  14.5371, 25.6134, -0.0508,  0.3710, -0.6246, -0.7682,  0.5868,  C2O
     *  11.4442, 24.7107, -0.5678, -0.0389,  1.0076, -4.6514,  4.3893,  CHP
     *   3.7900, 19.0227, -8.0668, -5.9821, -3.8685, -3.1838, -1.0364,  HS-
     *   0.7300, 19.7162, -5.0018, -2.7680, -1.2845, -0.9859, -0.3380,  H2-
     *   2.1200, 12.4717,  0.1601,  1.4596, -0.2012,  5.0788, -4.5487,  TiH
     *   9.7800, 35.8044, -1.3937, -0.2650, -0.6732, -2.5437,  2.9710,  PH3
     *   2.4000, 11.3146, -0.5595,  0.3619, -2.0065,  3.8766, -2.9900,  MgS
     *   0.0000,  4.5751,  3.4421,  0.7560, -1.7011,  1.4510, -1.3922,  TiO+
     *  21.1510, 31.0805, 10.7070, 12.8687, 10.5799,  6.4414,  3.6171,  LaO2
     *   3.2100, 12.1817, -0.7102, -0.2403,  1.1042, -1.3644,  0.3198,  Si2
     *  13.2716, 48.6914, -1.0602, -1.2802, -0.8603,  0.1159, -0.0701,  SiH4
     *   8.2349, 24.0157, -0.6514, -0.6064, -0.6542,  0.9096, -0.5839,  BH2
     *  10.9011, 25.1839, -0.1060,  0.2530, -1.1850,  2.3355, -1.6111,  AlOCl
     *   8.2300, 12.1920,  0.1751, -0.7678, -1.3836,  1.7704, -0.0265,  LaO
     *  14.0629, 25.1475, -0.2270,  0.7024, -0.8499,  0.4583,  0.1889,  C2N
     *  20.0747, 38.6719, -0.2664,  0.2782, -1.2642,  1.6020, -0.5248,  AlBO2
     *   4.3400, 10.9561, -0.8720,  3.4218,-12.2306, 18.7863,-11.1011,  KCl
     *   3.2300, 19.3359, -5.7570, -3.5853, -1.3882, -2.3313, -0.4930,  SiH-
     *   5.4800, 11.0459, -0.8574,  2.3137, -4.6777,  4.4532, -1.1716,  CaF
     *  17.8875, 47.4921, -1.1390,  2.7534, -7.2248,  6.3242, -1.1381,  CaO2H2
     *   8.1892, 23.3129, -1.0581,  3.5131,-11.3115, 16.9078, -9.8867/  KOH
      DATA C03/
     *  10.3100, 21.7682, -5.8992, -3.8627, -4.0284,  1.2924, -2.5856,  CN-
     *  16.1405, 37.9519, -0.0230,  0.6639, -2.4910,  5.5385, -4.2945,  Al2O2
     *   9.0621, 23.3478, -2.1422,  1.7058, -1.6807, 10.3429,-14.0183,  BaOH
     *   8.6837, 23.1042, -1.2656,  3.2436, -7.2017,  6.5067, -1.7129,  SrOH
     *  13.9839, 25.6721, -0.0784,  0.0544, -0.2755,  0.6140, -0.3673,  BO2
     *   5.5700, 12.0158, -0.5187, -0.1216,  0.6738, -0.6377,  0.1588,  SiF
     *   0.0000, 16.4621,-13.8562,-13.1896, -9.2577, -6.3354, -2.5704,  CH-
     *  13.8610, 26.3081, -1.3134,  0.1185, -0.0461, -0.4056,  0.8088,  C3
     *   8.4800, 21.1413, -5.8697, -3.3745, -2.7491, -1.8902, -0.2441,  C2-
     *  17.1545, 48.1845, -0.5683,  0.1125, -3.0973,  4.3727, -2.1978,  MgO2H2
     *   9.3961, 23.7967, -0.6500,  0.2061, -1.9381,  2.1259, -0.6451,  BeOH
     *  10.4305, 24.8357, -0.4930, -0.4550,  0.8862, -2.7257,  2.4025,  HBS
     *  13.1966, 25.7392,  0.0961, -0.7979, -0.1515,  4.2750, -4.6336,  SiC2
     *  17.4231, 48.8561, -0.4831,  0.9575, -1.9798, -0.0476,  1.2346,  FeO2H2
     *  10.0930, 25.0689, -1.5784,  2.2605, -3.1152,  3.7375, -2.5596,  CrO2
     *  20.0817, 49.3051, -0.2203,  0.6123, -1.9159,  3.0362, -0.6588,  BeH2O2
     *  11.4541, 36.8342, -1.3068, -1.2283, -0.7130, -0.1039,  0.8121,  BH3
     *  12.5346, 24.2744, -0.4230,  2.1003, -7.6565, 14.5171,-10.4377,  NaCN
     *   6.5483, 23.5736, -0.7830, -0.0881, -2.2398,  2.7050, -1.5244,  BeH2
     *  10.1248, 24.8268, -0.3784,  0.5561, -0.7324,  1.7508, -1.6977,  Si2N
     *   9.3132, 22.5681, -0.7730,  3.2979, -6.3686,  5.5210, -0.9987,  CaCl2
     *  18.8913, 37.0212, -0.3881,  1.7934, -7.5472, 14.9782,-11.0505,  NaBO2
     *   0.0000, 19.8338,-46.6804,-50.9308,-35.9059,-13.5611,-23.8103,  C3H
     *  15.5315, 26.0301, -0.1824,  0.0109, -0.3944,  0.5184, -0.0882,  OBF
     *  11.9993, 26.2368, -0.1708,  0.2491,  0.4220, -2.2962,  2.2409,  CS2
     *   8.9381, 23.5703, -0.6263,  1.0060, -4.3983,  7.4665, -4.8955,  LiOH
     *   1.5500, 11.3681, -0.1946, -0.0669, -2.3347,  5.3477, -4.0343,  Al2
     *   4.8400, 11.3090, -0.5602,  0.5886, -3.9705,  7.3873, -5.2571,  LiCl
     *  11.3225, 25.4462, -1.0487,  1.8142, -1.5110,  0.4282, -0.0240,  TiOCl
     *  23.3326, 62.7915, -1.3095, -1.6903, -0.9624, -1.6171,  2.5521,  C2H4
     *   7.4689, 23.8059, -0.5629,  0.0019, -0.3896, -0.7781,  0.3890,  CHCl
     *   6.6900, 14.8883,  5.3193,  8.9551,  3.7271,  5.1452,  1.0391,  TiCl
     *  19.2284, 37.1933,  0.1308, -0.0614, -0.9981,  2.9770, -2.1833,  AlOF2
     *  18.9713, 36.8674, -0.8338,  3.8816,-11.3916, 16.8414, -9.6911,  KBO2
     *  11.2271, 25.9412,  0.1074, -0.8813, -0.2594,  4.4112, -4.4861,  Si2C
     *   9.2183, 24.5270, -0.6453, -1.0757, -0.7155,  2.2944, -1.4513,  CHF
     *   0.0000, 11.8175,-29.4442,-30.6402,-22.9279,-13.1209, -8.8023,  BO-
     *  10.9760, 27.6834,  5.5082,  6.6402,  5.5692,  2.7324,  1.9375,  AlO2
     *  18.0802, 47.0050, -2.3587,  2.3466, -2.2753,  8.4432,-11.3032,  BaO2H2
     *  12.8526, 25.8889, -1.0260,  1.8361, -1.5017,  0.3478,  0.0486,  OTiF
     *   6.5000, 20.6745, -7.9942, -5.7057, -2.6759, -6.1649,  1.2656,  CS-
     *  21.5636, 39.0495, -0.1190,  0.7088, -1.5184,  0.4914,  0.9277,  C2N2
     *  17.5958, 46.9386, -1.3295,  3.5725, -8.4710,  7.5694, -1.8456,  SrO2H2
     *  12.2076, 25.3442, -0.0379, -0.1189, -0.8276,  1.3188, -0.6986,  ClCN
     *  10.6135, 23.6489, -0.5207,  0.0519, -0.6538,  1.9149, -1.5058,  AlClF
     *  12.5010, 24.1386, -0.8692,  4.1888,-11.7377, 17.1662, -9.8522,  KCN
     *   8.8688, 23.5425, -0.5528,  0.0031, -0.7346,  2.3344, -1.9878,  AlCl2
     *   9.6070, 22.2204, -2.5275,  2.8555, -1.4987,  7.7865,-11.3039,  BaCl2
     *  12.3143, 24.3964, -0.4940,  0.0699, -0.5475,  1.6261, -1.2695,  AlF2
     *   8.1536, 22.9187, -0.1815,  0.6847, -2.4792,  4.3296, -2.7691/  MgCl2
      DATA C04/
     *   0.0000, 17.5598,-16.6727,-14.0707,-13.0780, -5.4193, -4.7856,  FeO-
     *  20.4537, 49.9913, -0.5362, -0.7176, -1.2169,  1.1206, -0.3773,  BO2H2
     *  14.1133, 48.5194, -0.8436, -1.0629, -0.7362,  0.3080, -0.3403,  SiH3Cl
     *   8.3239, 23.6272, -0.2108,  1.1105, -2.1105,  1.5380, -0.1684,  FeCl2
     *   7.3840, 24.8600, -0.1499, -0.1631,  0.1378,  1.6604, -1.9986,  Si3
     *  16.1268, 48.9782, -0.8260, -1.0380, -0.6452, -0.1029,  0.1199,  SiH3F
     *  16.2992, 49.7196, -1.2716, -1.4752, -1.1626,  0.6516, -0.0837,  CH3Cl
     *   9.1791, 22.1133, -1.4891,  4.1050, -7.6534,  6.6694, -1.5355,  SrCl2
     *  11.6845, 23.2600, -1.2039,  3.3661, -6.2828,  5.1661, -0.6547,  CaF2
     *  13.7563, 25.2856, -0.4137,  1.0746, -1.1248,  0.2935,  0.3807,  TiF2
     *  19.4163, 36.9346, -0.3977,  1.3814, -4.7577,  8.2956, -5.5779,  LiBO2
     *   9.5422, 23.6489, -0.6541,  0.7042, -2.5258,  4.5411, -3.0359,  MgClF
     *  19.3953, 37.4967, -0.4103,  0.6249, -2.5737,  3.7334, -2.0769,  BeBO2
     *  16.1988, 37.8077, -0.3545, -0.2428, -0.1731, -1.4896,  1.9844,  C2HCl
     *   9.9277, 24.6274, -0.5062,  0.9860, -1.3100,  0.8075, -0.0931,  TiCl2
     *  19.7168, 40.3256, -0.2533,  0.3731, -0.5863, -0.6939,  0.9337,  C4
     *  30.6562, 75.8041, -1.6269, -1.1205, -1.8109,  2.1354, -0.8357,  H3BO3
     *  10.7510, 23.8686, -0.6130,  0.7434, -2.6657,  5.0507, -3.5509,  MgF2
     *   0.0000, 13.8534,-28.5088,-27.6557,-25.0420, -4.2145,-21.0916,  BaClF
     *  13.3200, 24.6323, -0.2099,  0.5174, -1.9085,  2.9836, -1.7351,  BeF2
     *  16.6788, 38.1093, -0.3632, -0.2642, -0.4287, -0.5573,  0.9863,  C2HF
     *   9.6498, 23.7877, -0.2606,  0.4816, -1.7048,  2.1226, -0.8176,  BeCl2
     *  15.7352, 37.1910, -1.0480,  1.8371, -1.1420, -0.7526,  1.2880,  TiOCl2
     *  10.7683, 24.3508, -0.5859,  0.0972, -0.3635,  0.9082, -0.3338,  ZrCl2
     *  11.9101, 22.9073, -2.4413,  2.9420, -1.3655,  7.3312,-10.8692,  BaF2
     *  12.4073, 25.2586, -0.5256,  0.7548, -2.0655,  2.2598, -0.9944,  BeC2
     *   9.9676, 24.0020, -0.4765,  1.0925, -3.6131,  4.2582, -1.8225,  Be2O
     *  11.3542, 22.8132, -1.4157,  4.1790, -7.3508,  5.5696, -0.4507,  SrF2
     *  13.7587, 24.7160, -1.0103,  0.2376, -0.4664, -0.9114,  6.9672,  ZrF2
     *  13.0910, 27.6502,  6.5468,  8.2502,  7.3334,  4.1191,  1.2402,  FeF2
     *  12.5389, 37.9053, -1.3490,  3.1985, -1.1165, -6.7253,  7.3584,  P4
     *  19.0240, 49.7099, -0.5565, -0.7375, -0.2251, -1.1324,  1.2457,  SiH2F2
     *   3.2806, 41.7329, 32.0127, 34.5233, 27.1981, 13.3168, 13.4808,  H3O+
     *  27.0859, 54.0398,  0.0077,  0.4169, -0.9261, -0.3135,  0.6322,  C5
     *  19.7864, 37.9176, -0.7063,  1.7895, -1.5401,  0.9448, -0.6313,  TiF3
     *  14.3199, 37.3165, -0.8450,  1.6603, -1.6009,  0.8934, -0.5070,  TiCl3
     *  15.5540, 36.5254, -0.7361,  0.8503, -0.3688,  0.0324,  0.0881,  ZrCl3
     *  10.6603, 34.6664, -0.4567,  3.2641,-13.6211, 27.6173,-20.7914,  Na2Cl2
     *  18.1954, 60.7438, -0.7643,  2.2577,-14.4187, 28.3225,-20.4866,  (NaOH)2
     *  28.8149, 64.3940, -0.2174,  1.3367, -6.6368,  8.6309, -4.6284,  Be3O3
     *  10.8345, 33.9871, -1.3140,  7.4840,-21.9583, 33.6428,-20.3143,  K2Cl2
     *  18.3196, 60.4179, -1.6298,  6.4524,-22.9230, 33.8810,-20.0092,  (KOH)2
     *  20.4364, 49.7173, -0.6667,  0.8064, -0.1308, -0.4433,  0.8970,  ZrCl4
     *  27.1266, 62.7471, -0.3813,  3.6624,-15.0927, 27.0694,-18.7738,  (NaCN)2
     *  27.0557, 51.2712, -0.5271,  0.8930, -0.5666,  1.5292, -1.3568,  ZrF4
     *  20.3442, 61.3686, -0.8410,  1.3617, -9.5297, 16.1158,-11.1739/  (LiOH)2
C
C Coefficients for constructing partition functions (and then equilibrium
C   constants, perhaps). For diatomic molecules other than H2 and CO, the
C   data are from Sauval & Tatum (1984, ApJS, 56, 193). For H2 and CO, the
C   data are from Irwin (1987, A&A, 182, 348). For polyatomic molecules,
C   the coefficients are from Irwin (1988, A&AS, 74,145).
C Coefficients used to construct the partition function, as follows:
C
C     log10(Q) = Sum{i=0,9}{PCOEF(i+1)*log10(THETA)**i}
C                                                           Ioniz. pot.
      DATA P01/
     *   1.69179,      -1.72270,       0.798033,     -0.157089,         H2
     *  -0.535313,      1.75818,      -2.63895,       1.35708,          H2
     *   0.0,           0.0,                                 15.42593,  H2
     *   3.615300,     -1.773848,      0.3516181,     0.08620792,       CO
     *   0.2911791,    -1.141469,      2.513133,     -2.886502,         CO
     *   1.238932,      0.0,                                 14.01400,  CO
     *   4.344711818,  -3.6343233,     1.415963,      0.01594,          H2O
     *   0.56542,      -1.2583,        0.53796,       3*0.0, 12.62100,  H2O
     *   3.0929, -1.6778,  0.6743, -0.1874,  0.0000,  5*0.0, 13.01700,  OH
     *   3.2643, -1.7303,  0.4192,  0.0000,  0.0000,  5*0.0, 15.58100,  N2
     *   4.2275, -1.9144,  0.7201, -1.3099,  1.1657,  5*0.0, 11.49000,  SiO
     *  1.0, 9*0.,                                           10.42200,  HS
     *   5.117210341,  -3.94844146,    1.23193,       0.076156,         H2S
     *   0.42163,      -0.453534,      0.0,           3*0.0, 10.45700,  H2S
     *   3.0735, -1.8501,  0.9607, -0.3935,  0.0000,  5*0.0, 13.49000,  NH
     *   3.6908, -1.9801,  0.7704, -0.2247,  0.0000,  5*0.0,  7.91000,  SiH
     *   3.3586, -2.0656,  0.9624, -0.2239,  0.0000,  5*0.0, 10.64000,  CH
     *   2.5410, -2.4336,  1.4979,  0.0192, -0.7483,  5*0.0, -1.00000,  H2+
     *   4.3073, -1.8255,  0.3765,  0.0000,  0.0000,  5*0.0,  9.26420,  NO
     *   3.6704, -2.2682,  0.9354, -0.2597,  0.0000,  5*0.0,  7.20000,  MgH
     *   2.8005, -1.7476,  0.5310,  0.0000,  0.0000,  5*0.0, 12.74400,  HCl
     *   4.8026, -1.9753,  0.2600,  0.0000,  0.0000,  5*0.0, 10.53000,  SiS
     *   6.103792598,  -4.3938712,     0.662588,      0.3751,           AlOH
     *   0.38386,      -0.2147,        0.0,           3*0.0, -1.00000,  AlOH
     *   4.819621858,  -3.84200734,    1.5386462,     0.784399,         NH2
     *  -2.34404,       2.50803,      -1.13304,       3*0.0, 11.14000,  NH2
     *   3.3209, -2.5909,  1.7415, -0.7636,  0.0000,  5*0.0,  5.50000,  AlH
     *   4.0078, -2.1514,  0.9226, -0.1671,  0.0000,  5*0.0, 13.59800,  CN
     *   6.01081285,   -4.438833,      0.840462,      0.2945,           CO2
     *   0.3694,       -0.273,         0.0,           3*0.0, 13.77700,  CO2
     *   4.7963, -2.1308,  0.5224,  0.0000,  0.0000,  5*0.0, 10.29400,  SO
c    *   5.7765, -2.3739,  0.8940, -0.3641,  0.0000,  5*0.0,  6.40000,  TiO
     *   5.3051, -2.3739,  0.8940, -0.3641,  0.0000,  5*0.0,  6.40000,  TiO
     *   5.0796, -2.1967,  0.4101,  0.0000,  0.0000,  5*0.0,  9.35600,  S2
     *   4.6265980,    -2.5625800,     0.38885943,    0.40219820,       FeH
     *  -0.21386399,    0.027845045,   0.0,           3*0.0,  7.37000,  FeH
     *   5.884176216,  -5.8364867,     1.608417,      1.50876,          NH3
     *  -0.59607,      -0.58961,       0.2459,        3*0.0, -1.00000,  NH3
     *   5.434042379,  -4.2409874,     0.988745,      0.49464,          HCN
     *   0.03719,      -0.22924,       0.0,           3*0.0, 13.60000,  HCN
     *   6.298781639,  -3.85672804,    0.8551678,     0.321901,         HCO
     *   0.020274,      0.15254,      -0.25298,       3*0.0,  8.12000,  HCO
     *   4.0636, -2.0779,  0.7660, -0.2111,  0.0000,  5*0.0, 12.06970,  O2
     *  1.0, 9*0.,                                           10.39600,  CH2
     *   2.4164, -1.6132,  0.6357, -0.1767,  0.0000,  5*0.0, 16.03000,  HF
     *  1.0, 9*0.,                                           -1.00000,  H3+
     *   3.8411, -2.3891,  1.3578, -0.6893,  0.0000,  5*0.0,  5.86000,  CaH
     *  1.0, 9*0.,                                           -1.00000,  Al2O
     *   4.9191, -2.6291,  0.5831,  0.3163,  0.0000,  5*0.0,  9.46000,  AlO
     *  1.0, 9*0.,                                            9.84000,  CH3
     *  1.0, 9*0.,                                            8.80000,  SiH2
     *   5.3182, -2.6502, -0.2781, -0.7823,  1.3107,  5*0.0,  8.76000,  MgO
     *   4.3091, -2.2406,  0.4865, -0.2049,  0.0000,  5*0.0, 11.40000,  C2
     *  1.0, 9*0.,                                            9.50000,  TiO2
     *   8.457240767,  -4.1987868,     0.334575,      0.20744,          VO2
     *   0.18226,      -0.053465,      0.0,           3*0.0, -1.00000,  VO2
     *   3.5453, -2.3457,  0.8557, -0.1685,  0.0000,  5*0.0,  4.70000,  NaH
     *   5.1115, -2.2303,  0.8001, -0.5192,  0.0000,  5*0.0,  9.40000,  AlCl
     *   4.5405, -2.1033,  0.6208, -0.2930,  0.0000,  5*0.0, -1.00000,  AlF
     *   5.0687, -2.2186,  0.9545, -0.4592,  0.0000,  5*0.0,  7.23860,  VO
     *   4.1646, -1.9348,  0.8034, -1.3669,  1.1561,  5*0.0, 11.33000,  CS
     *   6.8401894714, -4.338616427,   0.71600166,    0.128126,         MgOH
     *   0.5978087,    -0.8658369,     0.385049,      3*0.0,  7.50000,  MgOH
     *  1.0, 9*0.,                                           11.90000,  PO2
     *   7.1623971155, -4.471282563,   1.1221899,    -0.558812,         CaOH
     *   0.2294,        1.78658,      -2.95118,       1.41591,          CaOH
     *   2*0.0,                                               5.80000,  CaOH
     *  1.0, 9*0.,                                            9.82400/  PH2
      DATA P02/
     *  1.0, 9*0.,                                           11.61000,  C2H
     *   4.8065, -2.2129,  0.9991, -0.5414,  0.0000,  5*0.0, -1.00000,  ScO
     *  1.0, 9*0.,                                           -1.00000,  AlO2H
     *   5.2461, -2.1319,  0.5340, -0.2309,  0.0000,  5*0.0, -1.00000,  AlS
     *   5.5642, -2.1947,  0.5065,  0.0000,  0.0000,  5*0.0,  8.90000,  FeO
     *   5.5270, -2.1311,  0.6523, -0.2533,  0.0000,  5*0.0,  7.85000,  CrO
     *  1.0, 9*0.,                                           12.61000,  CH4
     *   4.8052, -1.9619,  0.3140,  0.0000,  0.0000,  5*0.0,  8.87000,  NS
     *  1.0, 9*0.,                                           12.34900,  SO2
     *   4.6570, -2.3587,  0.8819, -0.1642,  0.0000,  5*0.0, -1.00000,  SiN
     *  1.0, 9*0.,                                           -1.00000,  OH-
     *   5.3279, -2.4694,  0.2164, -0.2313,  0.0000,  5*0.0,  6.00000,  ZrO
     *   3.5649, -1.7328,  0.4241,  0.0000,  0.0000,  5*0.0, -1.00000,  NO+
     *   8.72011985,   -4.247295,      0.2758,        0.20738,          ZrO2
     *   0.09406,       0.0,           0.0,           3*0.0, -1.00000,  ZrO2
     *   3.9953, -1.8665,  0.5965, -0.1617,  0.0000,  5*0.0, 13.30000,  BO
     *  1.0, 9*0.,                                           -1.00000,  SiO2
     *  1.0, 9*0.,                                           -1.00000,  HBO
     *   5.1477, -1.8671,  0.2404,  0.0000,  0.0000,  5*0.0,  9.20000,  SiC
     *  1.0, 9*0.,                                           -1.00000,  YO2
     *   5.8948, -2.2183,  0.5928, -0.3106,  0.0000,  5*0.0,  7.10000,  TiS
     *  1.0, 9*0.,                                           -1.00000,  HBO2
     *   7.1220464309, -6.966653604,   1.9668235,     0.362597,         C2H2
     *   0.608996,     -0.920435,      0.271892,      3*0.0, 11.40000,  C2H2
     *  1.0, 9*0.,                                           11.18500,  OCS
     *  1.0, 9*0.,                                           -1.00000,  ZrO+
     *  1.0, 9*0.,                                           -1.00000,  NaOH
     *   5.7494, -2.3340,  0.8685, -0.5306,  0.0000,  5*0.0,  5.86000,  CaCl
     *  1.0, 9*0.,                                           -1.00000,  AlOF
     *   4.9515, -2.0866,  0.6565, -0.3082,  0.0000,  5*0.0,  6.00000,  YO
     *   5.3364, -2.2844,  0.2820,  0.1185,  0.0000,  5*0.0, -1.00000,  NaCl
     *  1.0, 9*0.,                                           -1.00000,  C2O
     *  1.0, 9*0.,                                           10.79000,  CHP
     *  1.0, 9*0.,                                           -1.00000,  HS-
     *  1.0, 9*0.,                                           -1.00000,  H2-
     *  1.0, 9*0.,                                            6.00000,  TiH
     *  1.0, 9*0.,                                            9.86900,  PH3
     *   5.0367, -2.1625,  0.4859, -0.1780,  0.0000,  5*0.0, -1.00000,  MgS
     *  1.0, 9*0.,                                           -1.00000,  TiO+
     *  1.0, 9*0.,                                           -1.00000,  LaO2
     *   5.2617, -2.1485,  0.5647, -0.2985,  0.0000,  5*0.0, -1.00000,  Si2
     *  1.0, 9*0.,                                           -1.00000,  SiH4
     *  1.0, 9*0.,                                            9.80000,  BH2
     *  1.0, 9*0.,                                           -1.00000,  AlOCl
     *   5.1147, -2.5016,  1.0445, -0.3135,  0.0000,  5*0.0,  4.95000,  LaO
     *  1.0, 9*0.,                                           12.00000,  C2N
     *  1.0, 9*0.,                                           -1.00000,  AlBO2
     *   5.6860, -2.3016,  0.2086,  0.1763,  0.0000,  5*0.0, -1.00000,  KCl
     *  1.0, 9*0.,                                           -1.00000,  SiH-
     *   5.2010, -2.2653,  0.8941, -0.5384,  0.0000,  5*0.0, -1.00000,  CaF
     *  1.0, 9*0.,                                           -1.00000,  CaO2H2
     *  1.0, 9*0.,                                            7.50000/  KOH
      DATA P03/
     *  1.0, 9*0.,                                           -1.00000,  CN-
     *  1.0, 9*0.,                                           -1.00000,  Al2O2
     *  1.0, 9*0.,                                           -1.00000,  BaOH
     *  1.0, 9*0.,                                           -1.00000,  SrOH
     *  1.0, 9*0.,                                           -1.00000,  BO2
     *   5.0871, -2.0375,  0.4478, -0.1243,  0.0000,  5*0.0,  7.54000,  SiF
     *  1.0, 9*0.,                                           -1.00000,  CH-
     *   6.618407932,  -3.576399,      0.883642,      0.087548,         C3
     *   0.04817,      -0.16471,       0.0,           3*0.0, -1.00000,  C3
     *  1.0, 9*0.,                                           -1.00000,  C2-
     *  1.0, 9*0.,                                           -1.00000,  MgO2H2
     *  1.0, 9*0.,                                           -1.00000,  BeOH
     *  1.0, 9*0.,                                           -1.00000,  HBS
     *   7.54651307623,-5.075563869,   1.82960795,    0.0983258,        SiC2
     *  -6.335157,     14.33103,     -13.01689,       4.428233,         SiC2
     *   2*0.0,                                              10.20000,  SiC2
     *  1.0, 9*0.,                                           -1.00000,  FeO2H2
     *  1.0, 9*0.,                                           -1.00000,  CrO2
     *  1.0, 9*0.,                                           -1.00000,  BeH2O2
     *  1.0, 9*0.,                                           -1.00000,  BH3
     *  1.0, 9*0.,                                           -1.00000,  NaCN
     *  1.0, 9*0.,                                           -1.00000,  BeH2
     *  1.0, 9*0.,                                           -1.00000,  Si2N
     *  1.0, 9*0.,                                           -1.00000,  CaCl2
     *  1.0, 9*0.,                                           -1.00000,  NaBO2
     *  1.0, 9*0.,                                           -1.00000,  C3H
     *  1.0, 9*0.,                                           -1.00000,  OBF
     *  1.0, 9*0.,                                           10.07300,  CS2
     *  1.0, 9*0.,                                           -1.00000,  LiOH
     *   5.5538, -2.3365,  0.5754, -0.2119,  0.0000,  5*0.0,  5.40000,  Al2
     *   4.5605, -2.2216,  0.5760, -0.1706,  0.0000,  5*0.0,  9.57000,  LiCl
     *  1.0, 9*0.,                                           -1.00000,  TiOCl
     *  1.0, 9*0.,                                           -1.00000,  C2H4
     *  1.0, 9*0.,                                           -1.00000,  CHCl
     *  1.0, 9*0.,                                           -1.00000,  TiCl
     *  1.0, 9*0.,                                           -1.00000,  AlOF2
     *  1.0, 9*0.,                                           -1.00000,  KBO2
     *  1.0, 9*0.,                                           -1.00000,  Si2C
     *  1.0, 9*0.,                                           10.06000,  CHF
     *  1.0, 9*0.,                                           -1.00000,  BO-
     *  1.0, 9*0.,                                           -1.00000,  AlO2
     *  1.0, 9*0.,                                           -1.00000,  BaO2H2
     *  1.0, 9*0.,                                           -1.00000,  OTiF
     *  1.0, 9*0.,                                           -1.00000,  CS-
     *  1.0, 9*0.,                                           -1.00000,  C2N2
     *  1.0, 9*0.,                                           -1.00000,  SrO2H2
     *  1.0, 9*0.,                                           12.36000,  ClCN
     *  1.0, 9*0.,                                           -1.00000,  AlClF
     *  1.0, 9*0.,                                           -1.00000,  KCN
     *  1.0, 9*0.,                                           -1.00000,  AlCl2
     *  1.0, 9*0.,                                           -1.00000,  BaCl2
     *  1.0, 9*0.,                                           -1.00000,  AlF2
     *  1.0, 9*0.,                                           -1.00000/  MgCl2
      DATA P04/
     *  1.0, 9*0.,                                           -1.00000,  FeO-
     *  1.0, 9*0.,                                           -1.00000,  BO2H2
     *  1.0, 9*0.,                                           -1.00000,  SiH3Cl
     *  1.0, 9*0.,                                           -1.00000,  FeCl2
     *  1.0, 9*0.,                                           -1.00000,  Si3
     *  1.0, 9*0.,                                           -1.00000,  SiH3F
     *  1.0, 9*0.,                                           -1.00000,  CH3Cl
     *  1.0, 9*0.,                                           -1.00000,  SrCl2
     *  1.0, 9*0.,                                           -1.00000,  CaF2
     *  1.0, 9*0.,                                           -1.00000,  TiF2
     *  1.0, 9*0.,                                           -1.00000,  LiBO2
     *  1.0, 9*0.,                                           -1.00000,  MgClF
     *  1.0, 9*0.,                                           -1.00000,  BeBO2
     *  1.0, 9*0.,                                           -1.00000,  C2HCl
     *  1.0, 9*0.,                                           -1.00000,  TiCl2
     *  1.0, 9*0.,                                           -1.00000,  C4
     *  1.0, 9*0.,                                           -1.00000,  H3BO3
     *  1.0, 9*0.,                                           -1.00000,  MgF2
     *  1.0, 9*0.,                                           -1.00000,  BaClF
     *  1.0, 9*0.,                                           -1.00000,  BeF2
     *  1.0, 9*0.,                                           -1.00000,  C2HF
     *  1.0, 9*0.,                                           -1.00000,  BeCl2
     *  1.0, 9*0.,                                           -1.00000,  TiOCl2
     *  1.0, 9*0.,                                           -1.00000,  ZrCl2
     *  1.0, 9*0.,                                           -1.00000,  BaF2
     *  1.0, 9*0.,                                           -1.00000,  BeC2
     *  1.0, 9*0.,                                           -1.00000,  Be2O
     *  1.0, 9*0.,                                           -1.00000,  SrF2
     *  1.0, 9*0.,                                           -1.00000,  ZrF2
     *  1.0, 9*0.,                                           -1.00000,  FeF2
     *  1.0, 9*0.,                                           -1.00000,  P4
     *  1.0, 9*0.,                                           -1.00000,  SiH2F2
     *  1.0, 9*0.,                                           -1.00000,  H3O+
     *  1.0, 9*0.,                                           -1.00000,  C5
     *  1.0, 9*0.,                                           -1.00000,  TiF3
     *  1.0, 9*0.,                                           -1.00000,  TiCl3
     *  1.0, 9*0.,                                           -1.00000,  ZrCl3
     *  1.0, 9*0.,                                           -1.00000,  Na2Cl2
     *  1.0, 9*0.,                                           -1.00000,  Na2O2H2
     *  1.0, 9*0.,                                           -1.00000,  Be3O3
     *  1.0, 9*0.,                                           -1.00000,  K2Cl2
     *  1.0, 9*0.,                                           -1.00000,  K2O2H2
     *  1.0, 9*0.,                                           -1.00000,  ZrCl4
     *  1.0, 9*0.,                                           -1.00000,  Na2C2N2
     *  1.0, 9*0.,                                           -1.00000,  ZrF4
     *  1.0, 9*0.,                                           -1.00000/  Li2O2H2
C
C Try to find the input speicies name (SPNAME) in the list (SPLIST) of
C species for which we have equilibrium constant coefficients. Note that
C the index is stored in a new variable J, rather than using the loop
C variable I, because some optimizers don't save the loop variable after
C normal termination of the loop.
C
      DO 1 I=1,MSPEC
      J=I
      IF(SPLIST(J).EQ.SPNAME) GO TO 2
   1  CONTINUE
C
C Fall through to here, if requested molecule was not in SPLIST.
C Print a warning, but return anyway.
C
      WRITE(*,*) 'MOLCON: Don''t have dissociation constant for ',
     *           'molecule: "', SPNAME, '"'
      EQK =1.D20
      PART=1.D0
      RETURN
C
C Calculate independent variable for polynomial expansions.
C Note that the polynomial expansions in Sauval & Tatum (1984) and Irwin
C (1987,1988) are in terms of log10(5040/T), not log10(5039.7475/T), but
C the more accuate value of 5039.7475 should be used in converting the
C partition function into an equilibrium constant.
C
   2  TH=5040.D0/T
      LOGTH=LOG10(TH)
C
C Construct equilibrium constant from polynomial coefficients and
C dissociation constant. A "+1" term at the end would convert from
C pascals (i.e. N/m/m as in Sauval) to dynes/cm/cm.
C
c      if (t.lt.1600) logth=log10(5040.0/1600.0)
c      if (t.gt.7730) logth=log10(5040.0/7730.0)
      EQK=COEF(2,J)+LOGTH*(COEF(3,J)+LOGTH*(COEF(4,J)+
     &              LOGTH*(COEF(5,J)+LOGTH*(COEF(6,J)+
     &              LOGTH*(COEF(7,J))))))
     &             -TH*COEF(1,J)
C    &             +1.0D0
      EQK=10.D0**EQK
C
C Just for the reference, the relation between partition functions
C and equilibrium constant:
C
C            P(A)*P(B)*...      N(A)*N(B)*...
C K(AB...) = ------------- = kT-------------- =
C              P(AB...)           N(AB...)
C
C             2*pi*kT 3/2    M(A)*M(B)*... 3/2   Q(A)*Q(B)*...
C       = kT*(-------)    * (-------------)    * ------------- * exp(-D(AB)/kT)
C               h^2            M(AB...)           Q(AB...)
C
C where, K - equilibrium constant, Q - partition functions, M - masses
C        P - partial pressures, N - number densities, T - temperature,
C        D - complete dissociation energy, h - plank constant. Remember
C        to use masses in grams (1 amu = 1.660540E-24 g) and energy in
C        ergs (1 eV = 1.60219E-12 ergs). Also, k = 1.38065E-16 erg/K,
C        h = 6.626076E-27 erg s, and pi = 3.1415926536.
C
C Construct partition function from polynomial coefficients.
C
      PART=PCOEF(NPCOEF-1,J)
      DO 3 I=NPCOEF-2,1,-1
    3 PART=LOGTH*PART+PCOEF(I,J)
C
C Copy ionization potential
C
      PION=PCOEF(NPCOEF,J)
C
C Calculate equilibrium constant (EQK) from partition function, dissociation
C constant, and other information passed into subroutine. The constants used
C are:  79.733501 = 1.5*log10(2*pi/h/h)  [in cgs units] and
C      -15.859914 = alog10(k)            [in cgs units].
C       5039.7475 = alog10(e)*k*(eV/erg)
C
      EQK_ST=(NTOT-1)*(79.733501D0+2.5D0*(LOG10(T)-15.859914D0))+
     &       1.5D0*RATIOM+QPRD-PART-COEF(1,J)*5039.7475D0/T
C
C Convert equilibrium constant and partition function from logarithms.
C
      EQK_ST=10.D0**EQK_ST
      PART=10.D0**PART
C
C Check if there is relevant data in Paul Barklem's tables
C
c      CALL KP_Q_SPLN(SPNAME,T,Qm_spln,Kp_spln)
c      write(*,'(F10.1,A9,5G13.6)') T,SPNAME,EQK,EQK_ST,Kp_spln,
c     &                             PART,Qm_spln
c      IF(Kp_spln.GE.0.d0) THEN
c        EQK =Kp_spln
c        PART=Qm_spln
c      ENDIF

c      if(spname.eq.'H2+') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'NO') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'C3') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'SiN') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'SiS') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'TiO') write(*,'(a,f10.2,1p3e16.8)')
c     &   spname,t , eqk, eqk_st, part
c      if(spname.eq.'N2') write(*,'(a,f10.2,1p3e14.6,i3,1p2e14.6)')
c     &   spname,t , eqk, eqk_st, part,NTOT,QPRD,RATIOM
c
c Don't use EQK_ST based on partition function - use direct fit to EQK.
c
c      EQK=EQK_ST
C
C Done.
C
      RETURN
      END


C=========================================================================
C Kp_spln: Returns equilibrium constant for a given molecule and temperature.
C
C Inputs:
C   ISPEC [integer] Species number according to the table below.
C   THETA=5040/T [real] where T is the temperature (in K) at which Kp is
C         interpolated.
C History:
C  28-jun-2007: First version written by N. Piskunov including 57 species.
C               Molecular equilibium tabulated by P. Barklem, resampled
C               for optimal spline interpolation and converted to Fortran
C               DATA statements by J. Valenti
C
C Outputs:
C   Kp_spln [real*8] Equilibrium constant (in dynes/cm/cm) at temperature T,
C     interpolated from Paul Berklem's tables.
C
C Reference:
C   Paul Barklem 2007, in preparaton.
C
      SUBROUTINE KP_Q_SPLN(SPNAME,TEMP,Qm_spln,Kp_spln)
C
      IMPLICIT NONE
      INTEGER ISPEC
      CHARACTER SPNAME*(*)
      REAL TEMP
      REAL*8 Qm_spln,Kp_spln
C
C  Local variables
C
      LOGICAL FIRST
      INTEGER MSPEC,NTH,KLO,KHI,I
      PARAMETER(MSPEC=58,NTH=36)
      REAL TEMP0
      REAL*8 THETA,A,U(NTH),SPL_INTERP
C
C Molecular equilibrium constants, alog10(Kp)
C  as a function of theta, 5040/temperature
C
C D2 = SPL_INIT(THETA,KP_<species>)
C KP2 = SPL_INTERP(THETA,Kp_<species>,D2,theta2)
C
      CHARACTER SPLIST(MSPEC)*8
      REAL*8 Qm(NTH,MSPEC),Qm2(NTH,MSPEC),Kp(NTH,MSPEC),Kp2(NTH,MSPEC)
      REAL*8 TH(NTH),
     * Qm_H2p (NTH),Qm_H2  (NTH),Qm_H2m (NTH),Qm_CH  (NTH),Qm_CHm (NTH),
     * Qm_C2  (NTH),Qm_C2m (NTH),Qm_CN  (NTH),Qm_CNm (NTH),Qm_NH  (NTH),
     * Qm_N2  (NTH),Qm_OH  (NTH),Qm_OHm (NTH),Qm_BO  (NTH),Qm_CO  (NTH),
     * Qm_NOp (NTH),Qm_NO  (NTH),Qm_O2  (NTH),Qm_HF  (NTH),Qm_NaH (NTH),
     * Qm_MgH (NTH),Qm_MgO (NTH),Qm_AlH (NTH),Qm_AlO (NTH),Qm_AlF (NTH),
     * Qm_Al2 (NTH),Qm_SiH (NTH),Qm_SiHm(NTH),Qm_SiC (NTH),Qm_SiN (NTH),
     * Qm_SiO (NTH),Qm_SiF (NTH),Qm_Si2 (NTH),Qm_SH  (NTH),Qm_SHm (NTH),
     * Qm_CS  (NTH),Qm_NS  (NTH),Qm_SO  (NTH),Qm_MgS (NTH),Qm_AlS (NTH),
     * Qm_SiS (NTH),Qm_S2  (NTH),Qm_HCl (NTH),Qm_LiCl(NTH),Qm_NaCl(NTH),
     * Qm_AlCl(NTH),Qm_CaH (NTH),Qm_CaF (NTH),Qm_CaCl(NTH),Qm_ScO (NTH),
     * Qm_TiO (NTH),Qm_TiS (NTH),Qm_VO  (NTH),Qm_CrO (NTH),Qm_FeO (NTH),
     * Qm_YO  (NTH),Qm_ZrO (NTH),Qm_LaO (NTH),
     * Kp_H2p (NTH),Kp_H2  (NTH),Kp_H2m (NTH),Kp_CH  (NTH),Kp_CHm (NTH),
     * Kp_C2  (NTH),Kp_C2m (NTH),Kp_CN  (NTH),Kp_CNm (NTH),Kp_NH  (NTH),
     * Kp_N2  (NTH),Kp_OH  (NTH),Kp_OHm (NTH),Kp_BO  (NTH),Kp_CO  (NTH),
     * Kp_NOp (NTH),Kp_NO  (NTH),Kp_O2  (NTH),Kp_HF  (NTH),Kp_NaH (NTH),
     * Kp_MgH (NTH),Kp_MgO (NTH),Kp_AlH (NTH),Kp_AlO (NTH),Kp_AlF (NTH),
     * Kp_Al2 (NTH),Kp_SiH (NTH),Kp_SiHm(NTH),Kp_SiC (NTH),Kp_SiN (NTH),
     * Kp_SiO (NTH),Kp_SiF (NTH),Kp_Si2 (NTH),Kp_SH  (NTH),Kp_SHm (NTH),
     * Kp_CS  (NTH),Kp_NS  (NTH),Kp_SO  (NTH),Kp_MgS (NTH),Kp_AlS (NTH),
     * Kp_SiS (NTH),Kp_S2  (NTH),Kp_HCl (NTH),Kp_LiCl(NTH),Kp_NaCl(NTH),
     * Kp_AlCl(NTH),Kp_CaH (NTH),Kp_CaF (NTH),Kp_CaCl(NTH),Kp_ScO (NTH),
     * Kp_TiO (NTH),Kp_TiS (NTH),Kp_VO  (NTH),Kp_CrO (NTH),Kp_FeO (NTH),
     * Kp_YO  (NTH),Kp_ZrO (NTH),Kp_LaO (NTH)
C
      EQUIVALENCE (Qm(1, 1),Qm_H2p ),(Qm(1, 2),Qm_H2  )
      EQUIVALENCE (Qm(1, 3),Qm_H2m ),(Qm(1, 4),Qm_CH  )
      EQUIVALENCE (Qm(1, 5),Qm_CHm ),(Qm(1, 6),Qm_C2  )
      EQUIVALENCE (Qm(1, 7),Qm_C2m ),(Qm(1, 8),Qm_CN  )
      EQUIVALENCE (Qm(1, 9),Qm_CNm ),(Qm(1,10),Qm_NH  )
      EQUIVALENCE (Qm(1,11),Qm_N2  ),(Qm(1,12),Qm_OH  )
      EQUIVALENCE (Qm(1,13),Qm_OHm ),(Qm(1,14),Qm_BO  )
      EQUIVALENCE (Qm(1,15),Qm_CO  ),(Qm(1,16),Qm_NOp )
      EQUIVALENCE (Qm(1,17),Qm_NO  ),(Qm(1,18),Qm_O2  )
      EQUIVALENCE (Qm(1,19),Qm_HF  ),(Qm(1,20),Qm_NaH )
      EQUIVALENCE (Qm(1,21),Qm_MgH ),(Qm(1,22),Qm_MgO )
      EQUIVALENCE (Qm(1,23),Qm_AlH ),(Qm(1,24),Qm_AlO )
      EQUIVALENCE (Qm(1,25),Qm_AlF ),(Qm(1,26),Qm_Al2 )
      EQUIVALENCE (Qm(1,27),Qm_SiH ),(Qm(1,28),Qm_SiHm)
      EQUIVALENCE (Qm(1,29),Qm_SiC ),(Qm(1,30),Qm_SiN )
      EQUIVALENCE (Qm(1,31),Qm_SiO ),(Qm(1,32),Qm_SiF )
      EQUIVALENCE (Qm(1,33),Qm_Si2 ),(Qm(1,34),Qm_SH  )
      EQUIVALENCE (Qm(1,35),Qm_SHm ),(Qm(1,36),Qm_CS  )
      EQUIVALENCE (Qm(1,37),Qm_NS  ),(Qm(1,38),Qm_SO  )
      EQUIVALENCE (Qm(1,39),Qm_MgS ),(Qm(1,40),Qm_AlS )
      EQUIVALENCE (Qm(1,41),Qm_SiS ),(Qm(1,42),Qm_S2  )
      EQUIVALENCE (Qm(1,43),Qm_HCl ),(Qm(1,44),Qm_LiCl)
      EQUIVALENCE (Qm(1,45),Qm_NaCl),(Qm(1,46),Qm_AlCl)
      EQUIVALENCE (Qm(1,47),Qm_CaH ),(Qm(1,48),Qm_CaF )
      EQUIVALENCE (Qm(1,49),Qm_CaCl),(Qm(1,50),Qm_ScO )
      EQUIVALENCE (Qm(1,51),Qm_TiO ),(Qm(1,52),Qm_TiS )
      EQUIVALENCE (Qm(1,53),Qm_VO  ),(Qm(1,54),Qm_CrO )
      EQUIVALENCE (Qm(1,55),Qm_FeO ),(Qm(1,56),Qm_YO  )
      EQUIVALENCE (Qm(1,57),Qm_ZrO ),(Qm(1,58),Qm_LaO )
      EQUIVALENCE (Kp(1, 1),Kp_H2p ),(Kp(1, 2),Kp_H2  )
      EQUIVALENCE (Kp(1, 3),Kp_H2m ),(Kp(1, 4),Kp_CH  )
      EQUIVALENCE (Kp(1, 5),Kp_CHm ),(Kp(1, 6),Kp_C2  )
      EQUIVALENCE (Kp(1, 7),Kp_C2m ),(Kp(1, 8),Kp_CN  )
      EQUIVALENCE (Kp(1, 9),Kp_CNm ),(Kp(1,10),Kp_NH  )
      EQUIVALENCE (Kp(1,11),Kp_N2  ),(Kp(1,12),Kp_OH  )
      EQUIVALENCE (Kp(1,13),Kp_OHm ),(Kp(1,14),Kp_BO  )
      EQUIVALENCE (Kp(1,15),Kp_CO  ),(Kp(1,16),Kp_NOp )
      EQUIVALENCE (Kp(1,17),Kp_NO  ),(Kp(1,18),Kp_O2  )
      EQUIVALENCE (Kp(1,19),Kp_HF  ),(Kp(1,20),Kp_NaH )
      EQUIVALENCE (Kp(1,21),Kp_MgH ),(Kp(1,22),Kp_MgO )
      EQUIVALENCE (Kp(1,23),Kp_AlH ),(Kp(1,24),Kp_AlO )
      EQUIVALENCE (Kp(1,25),Kp_AlF ),(Kp(1,26),Kp_Al2 )
      EQUIVALENCE (Kp(1,27),Kp_SiH ),(Kp(1,28),Kp_SiHm)
      EQUIVALENCE (Kp(1,29),Kp_SiC ),(Kp(1,30),Kp_SiN )
      EQUIVALENCE (Kp(1,31),Kp_SiO ),(Kp(1,32),Kp_SiF )
      EQUIVALENCE (Kp(1,33),Kp_Si2 ),(Kp(1,34),Kp_SH  )
      EQUIVALENCE (Kp(1,35),Kp_SHm ),(Kp(1,36),Kp_CS  )
      EQUIVALENCE (Kp(1,37),Kp_NS  ),(Kp(1,38),Kp_SO  )
      EQUIVALENCE (Kp(1,39),Kp_MgS ),(Kp(1,40),Kp_AlS )
      EQUIVALENCE (Kp(1,41),Kp_SiS ),(Kp(1,42),Kp_S2  )
      EQUIVALENCE (Kp(1,43),Kp_HCl ),(Kp(1,44),Kp_LiCl)
      EQUIVALENCE (Kp(1,45),Kp_NaCl),(Kp(1,46),Kp_AlCl)
      EQUIVALENCE (Kp(1,47),Kp_CaH ),(Kp(1,48),Kp_CaF )
      EQUIVALENCE (Kp(1,49),Kp_CaCl),(Kp(1,50),Kp_ScO )
      EQUIVALENCE (Kp(1,51),Kp_TiO ),(Kp(1,52),Kp_TiS )
      EQUIVALENCE (Kp(1,53),Kp_VO  ),(Kp(1,54),Kp_CrO )
      EQUIVALENCE (Kp(1,55),Kp_FeO ),(Kp(1,56),Kp_YO  )
      EQUIVALENCE (Kp(1,57),Kp_ZrO ),(Kp(1,58),Kp_LaO )
      SAVE
C
      DATA SPLIST/
     * 'H2+     ','H2      ','H2-     ','CH      ','CH-     ',
     * 'C2      ','C2-     ','CN      ','CN-     ','NH      ',
     * 'N2      ','OH      ','OH-     ','BO      ','CO      ',
     * 'NO+     ','NO      ','O2      ','HF      ','NaH     ',
     * 'MgH     ','MgO     ','AlH     ','AlO     ','AlF     ',
     * 'Al2     ','SiH     ','SiH-    ','SiC     ','SiN     ',
     * 'SiO     ','SiF     ','Si2     ','SH      ','SH-     ',
     * 'CS      ','NS      ','SO      ','MgS     ','AlS     ',
     * 'SiS     ','S2      ','HCl     ','LiCl    ','NaCl    ',
     * 'AlCl    ','CaH     ','CaF     ','CaCl    ','ScO     ',
     * 'TiO     ','TiS     ','VO      ','CrO     ','FeO     ',
     * 'YO      ','ZrO     ','LaO     '/
C
      DATA TH/
     1  5.04000000d-01, 5.32816666d-01, 5.89833336d-01, 6.72826987d-01, 5040/T
     2  7.86336892d-01, 9.38437080d-01, 1.14109393d+00, 1.41138056d+00, 5040/T
     3  1.77345739d+00, 2.26155251d+00, 2.92438188d+00, 3.83170538d+00, 5040/T
     4  5.08409278d+00, 6.82755387d+00, 9.27559118d+00, 1.27426522d+01, 5040/T
     5  1.76951993d+01, 2.48301764d+01, 3.51963422d+01, 5.03830850d+01, 5040/T
     6  7.28161181d+01, 1.06223484d+02, 1.56374576d+02, 2.32259457d+02, 5040/T
     7  3.47982501d+02, 5.25821797d+02, 8.01202354d+02, 1.23082959d+03, 5040/T
     8  1.90607197d+03, 2.97511316d+03, 4.67983950d+03, 7.41762986d+03, 5040/T
     9  1.18454701d+04, 1.90564126d+04, 3.08803914d+04, 5.03999992d+04/ 5040/T
C
C Molecular partition functions
C
      DATA Qm_H2p/
     1  3.16267828d+00, 3.11731572d+00, 3.03301461d+00, 2.92135062d+00, H2+
     2  2.78592091d+00, 2.62927368d+00, 2.45493684d+00, 2.26897367d+00, H2+
     3  2.08004824d+00, 1.89656355d+00, 1.72321037d+00, 1.56106148d+00, H2+
     4  1.40953276d+00, 1.26663134d+00, 1.12870288d+00, 9.91797882d-01, H2+
     5  8.53629822d-01, 7.14250795d-01, 5.72467338d-01, 4.37375094d-01, H2+
     6  3.05453704d-01, 1.84598818d-01, 8.06881265d-02, 2.38057269d-02, H2+
     7  4.39343320d-03, 3.17882761d-04, 2.78602393d-04,-1.35412088d-04, H2+
     8  9.67906370d-05, 2.83810105d-04, 6.21340947d-05,-1.48294987d-04, H2+
     9  1.19806307d-04, 4.49942650d-04, 5.42173714d-04, 2.63960122d-19/ H2+
      DATA Qm_H2/
     1  2.29424176d+00, 2.24556879d+00, 2.15660080d+00, 2.04203230d+00, H2
     2  1.90870863d+00, 1.76299415d+00, 1.61165102d+00, 1.46082610d+00, H2
     3  1.31459049d+00, 1.17465271d+00, 1.04101026d+00, 9.12250950d-01, H2
     4  7.85762722d-01, 6.58821989d-01, 5.29994245d-01, 3.99441359d-01, H2
     5  2.68276046d-01, 1.38511027d-01, 1.07111745d-02,-1.03221077d-01, H2
     6 -2.01061903d-01,-2.67229131d-01,-3.06685987d-01,-3.05106300d-01, H2
     7 -2.91672611d-01,-3.02510144d-01,-3.06339103d-01,-2.98255093d-01, H2
     8 -3.03013192d-01,-3.06845137d-01,-3.02303095d-01,-2.97991498d-01, H2
     9 -3.03484773d-01,-3.10249119d-01,-3.12138892d-01,-3.01029996d-01/ H2
      DATA Qm_H2m/
     1  2.59526213d+00, 2.54659436d+00, 2.45762984d+00, 2.34306219d+00, H2-
     2  2.20973862d+00, 2.06402415d+00, 1.91268102d+00, 1.76185609d+00, H2-
     3  1.61562048d+00, 1.47568270d+00, 1.34204025d+00, 1.21328095d+00, H2-
     4  1.08679272d+00, 9.59851985d-01, 8.31024240d-01, 7.00471354d-01, H2-
     5  5.69306042d-01, 4.39541023d-01, 3.11741170d-01, 1.97808919d-01, H2-
     6  9.99680931d-02, 3.38008644d-02,-5.65599168d-03,-4.07630432d-03, H2-
     7  9.35738509d-03,-1.48014811d-03,-5.30910765d-03, 2.77490263d-03, H2-
     8 -1.98319670d-03,-5.81514166d-03,-1.27309971d-03, 3.03849771d-03, H2-
     9 -2.45477746d-03,-9.21912291d-03,-1.11088960d-02,-5.40842441d-18/ H2-
      DATA Qm_CH/
     1  4.09454845d+00, 4.03220440d+00, 3.91951148d+00, 3.77721044d+00, CH
     2  3.61570249d+00, 3.44308749d+00, 3.26553980d+00, 3.08767363d+00, CH
     3  2.91356046d+00, 2.74694132d+00, 2.59019141d+00, 2.44339380d+00, CH
     4  2.30442717d+00, 2.16970558d+00, 2.03530309d+00, 1.89847102d+00, CH
     5  1.75825842d+00, 1.61502940d+00, 1.46708991d+00, 1.32239482d+00, CH
     6  1.17569105d+00, 1.03142213d+00, 8.87607004d-01, 7.68285866d-01, CH
     7  6.77576841d-01, 6.19286828d-01, 5.98162143d-01, 6.05344058d-01, CH
     8  5.99794356d-01, 5.95413239d-01, 6.00604749d-01, 6.05533212d-01, CH
     9  5.99254004d-01, 5.91521872d-01, 5.89361726d-01, 6.02059991d-01/ CH
      DATA Qm_CHm/
     1  4.09454845d+00, 4.03220440d+00, 3.91951148d+00, 3.77721044d+00, CH-
     2  3.61570249d+00, 3.44308749d+00, 3.26553980d+00, 3.08767363d+00, CH-
     3  2.91356046d+00, 2.74694132d+00, 2.59019141d+00, 2.44339380d+00, CH-
     4  2.30442717d+00, 2.16970558d+00, 2.03530309d+00, 1.89847102d+00, CH-
     5  1.75825842d+00, 1.61502940d+00, 1.46708991d+00, 1.32239482d+00, CH-
     6  1.17569105d+00, 1.03142213d+00, 8.87607004d-01, 7.68285866d-01, CH-
     7  6.77576841d-01, 6.19286828d-01, 5.98162143d-01, 6.05344058d-01, CH-
     8  5.99794356d-01, 5.95413239d-01, 6.00604749d-01, 6.05533212d-01, CH-
     9  5.99254004d-01, 5.91521872d-01, 5.89361726d-01, 6.02059991d-01/ CH-
      DATA Qm_C2/
     1  5.04115178d+00, 4.97854091d+00, 4.86608447d+00, 4.72408284d+00, C2
     2  4.56073360d+00, 4.38122886d+00, 4.18950267d+00, 3.98882965d+00, C2
     3  3.78180125d+00, 3.56991527d+00, 3.35313664d+00, 3.12987577d+00, C2
     4  2.89740336d+00, 2.65231907d+00, 2.39194918d+00, 2.11991938d+00, C2
     5  1.85568913d+00, 1.62969215d+00, 1.43817729d+00, 1.28755301d+00, C2
     6  1.12865144d+00, 9.67371508d-01, 7.98455286d-01, 6.35106239d-01, C2
     7  4.72003797d-01, 3.04045464d-01, 1.39473875d-01,-8.65658395d-03, C2
     8 -1.48227559d-01,-2.63386656d-01,-2.96336836d-01,-2.90010350d-01, C2
     9 -3.09027672d-01,-3.31175062d-01,-3.37371872d-01,-3.01029996d-01/ C2
      DATA Qm_C2m/
     1  4.45648843d+00, 4.38979169d+00, 4.27093699d+00, 4.12446863d+00, C2-
     2  3.96367459d+00, 3.79784756d+00, 3.63128021d+00, 3.46398784d+00, C2-
     3  3.29503238d+00, 3.12512084d+00, 2.95666634d+00, 2.79264407d+00, C2-
     4  2.63548539d+00, 2.48599990d+00, 2.34250202d+00, 2.20118238d+00, C2-
     5  2.05810169d+00, 1.91136422d+00, 1.75824111d+00, 1.60573317d+00, C2-
     6  1.44733943d+00, 1.28540016d+00, 1.11457265d+00, 9.52576772d-01, C2-
     7  7.90940804d-01, 6.21008511d-01, 4.55062833d-01, 3.06567630d-01, C2-
     8  1.63956233d-01, 4.44591200d-02, 6.94832691d-03, 1.06481798d-02, C2-
     9 -7.47103092d-03,-2.82236622d-02,-3.40356963d-02,-1.65737082d-17/ C2-
      DATA Qm_CN/
     1  4.73930577d+00, 4.67609227d+00, 4.56237616d+00, 4.41860097d+00, CN
     2  4.25381340d+00, 4.07539867d+00, 3.89073259d+00, 3.70677075d+00, CN
     3  3.52844994d+00, 3.35758358d+00, 3.19387161d+00, 3.03702382d+00, CN
     4  2.88730979d+00, 2.74410326d+00, 2.60469507d+00, 2.46518981d+00, CN
     5  2.32259529d+00, 2.17597053d+00, 2.02291928d+00, 1.87052795d+00, CN
     6  1.71228406d+00, 1.55056287d+00, 1.38006048d+00, 1.21855074d+00, CN
     7  1.05766058d+00, 8.88890514d-01, 7.24775740d-01, 5.79311060d-01, CN
     8  4.41799935d-01, 3.30270404d-01, 3.03291682d-01, 3.13384156d-01, CN
     9  2.91758320d-01, 2.66140616d-01, 2.58977780d-01, 3.01029996d-01/ CN
      DATA Qm_CNm/
     1  4.72601398d+00, 4.66373624d+00, 4.55150458d+00, 4.40932544d+00, CN-
     2  4.24608199d+00, 4.06907886d+00, 3.88565040d+00, 3.70273812d+00, CN-
     3  3.52528581d+00, 3.35512671d+00, 3.19198430d+00, 3.03559021d+00, CN-
     4  2.88623303d+00, 2.74330350d+00, 2.60410754d+00, 2.46476279d+00, CN-
     5  2.32228834d+00, 2.17575213d+00, 2.02276850d+00, 1.87042088d+00, CN-
     6  1.71221034d+00, 1.55051260d+00, 1.38002829d+00, 1.21852894d+00, CN-
     7  1.05764484d+00, 8.88881270d-01, 7.24770584d-01, 5.79307141d-01, CN-
     8  4.41798167d-01, 3.30270387d-01, 3.03291625d-01, 3.13383655d-01, CN-
     9  2.91758698d-01, 2.66142040d-01, 2.58979496d-01, 3.01029996d-01/ CN-
      DATA Qm_NH/
     1  3.77089299d+00, 3.70850776d+00, 3.59744910d+00, 3.46037343d+00, NH
     2  3.30887589d+00, 3.15109496d+00, 2.99206916d+00, 2.83445638d+00, NH
     3  2.67987579d+00, 2.52975595d+00, 2.38532657d+00, 2.24696307d+00, NH
     4  2.11337650d+00, 1.98162616d+00, 1.84852289d+00, 1.71227558d+00, NH
     5  1.57258714d+00, 1.43005318d+00, 1.28309934d+00, 1.13983393d+00, NH
     6  9.95257463d-01, 8.54229399d-01, 7.15659424d-01, 6.04918159d-01, NH
     7  5.27619579d-01, 4.86036812d-01, 4.73622924d-01, 4.79501747d-01, NH
     8  4.75442735d-01, 4.72199066d-01, 4.76043638d-01, 4.79693195d-01, NH
     9  4.75043405d-01, 4.69317716d-01, 4.67718114d-01, 4.77121255d-01/ NH
      DATA Qm_N2/
     1  3.82533357d+00, 3.77704498d+00, 3.69097394d+00, 3.58260865d+00, N2
     2  3.45768555d+00, 3.31992809d+00, 3.17252167d+00, 3.01853841d+00, N2
     3  2.86108775d+00, 2.70333858d+00, 2.54833923d+00, 2.39852586d+00, N2
     4  2.25490848d+00, 2.11629216d+00, 1.97945556d+00, 1.84087564d+00, N2
     5  1.69848038d+00, 1.55191345d+00, 1.39891505d+00, 1.24659595d+00, N2
     6  1.08845001d+00, 9.26870251d-01, 7.56577798d-01, 5.95383392d-01, N2
     7  4.34976861d-01, 2.66958810d-01, 1.04030161d-01,-3.94444033d-02, N2
     8 -1.73616385d-01,-2.80206194d-01,-3.01059584d-01,-2.87483161d-01, N2
     9 -3.11447151d-01,-3.40191537d-01,-3.48225307d-01,-3.01029996d-01/ N2
      DATA Qm_OH/
     1  3.70528328d+00, 3.65138235d+00, 3.55539924d+00, 3.43646432d+00, OH
     2  3.30357820d+00, 3.16255675d+00, 3.01717650d+00, 2.87025768d+00, OH
     3  2.72439070d+00, 2.58190037d+00, 2.44429498d+00, 2.31156479d+00, OH
     4  2.18184249d+00, 2.05210767d+00, 1.91983190d+00, 1.78404186d+00, OH
     5  1.64485834d+00, 1.50301265d+00, 1.35703876d+00, 1.21520068d+00, OH
     6  1.07275441d+00, 9.34997971d-01, 8.01912652d-01, 7.00092802d-01, OH
     7  6.35201135d-01, 6.06794258d-01, 5.99929861d-01, 6.03436627d-01, OH
     8  6.01082602d-01, 5.99194031d-01, 6.01432550d-01, 6.03557501d-01, OH
     9  6.00850165d-01, 5.97516388d-01, 5.96585022d-01, 6.02059991d-01/ OH
      DATA Qm_OHm/
     1  2.95023047d+00, 2.90796763d+00, 2.83150173d+00, 2.73415374d+00, OH-
     2  2.62155552d+00, 2.49774312d+00, 2.36621000d+00, 2.23024907d+00, OH-
     3  2.09298432d+00, 1.95713134d+00, 1.82447857d+00, 1.69523589d+00, OH-
     4  1.56776546d+00, 1.43937059d+00, 1.30791004d+00, 1.17269357d+00, OH-
     5  1.03396479d+00, 8.92504970d-01, 7.46898833d-01, 6.05426844d-01, OH-
     6  4.63428741d-01, 3.26277155d-01, 1.94155423d-01, 9.38364309d-02, OH-
     7  3.08673244d-02, 4.28699547d-03,-1.89094292d-03, 1.22280590d-03, OH-
     8 -8.68650427d-04,-2.54709963d-03,-5.57633334d-04, 1.33089943d-03, OH-
     9 -1.07522277d-03,-4.03808941d-03,-4.86583331d-03,-2.36895652d-18/ OH-
      DATA Qm_BO/
     1  4.63386799d+00, 4.57378147d+00, 4.46913601d+00, 4.34243673d+00, BO
     2  4.20268790d+00, 4.05406189d+00, 3.89829926d+00, 3.73658927d+00, BO
     3  3.57074746d+00, 3.40342587d+00, 3.23772886d+00, 3.07667403d+00, BO
     4  2.92246882d+00, 2.77552645d+00, 2.63375399d+00, 2.49324386d+00, BO
     5  2.35039126d+00, 2.20369911d+00, 2.05059409d+00, 1.89811862d+00, BO
     6  1.73976147d+00, 1.57787411d+00, 1.40712319d+00, 1.24523963d+00, BO
     7  1.08377547d+00, 9.14111854d-01, 7.48588444d-01, 6.00787118d-01, BO
     8  4.59345720d-01, 3.41721176d-01, 3.06752397d-01, 3.12041955d-01, BO
     9  2.93148325d-01, 2.71294504d-01, 2.65177453d-01, 3.01029996d-01/ BO
      DATA Qm_CO/
     1  4.17035582d+00, 4.12074763d+00, 4.03322600d+00, 3.92389753d+00, CO
     2  3.79820183d+00, 3.65948444d+00, 3.51071048d+00, 3.35483953d+00, CO
     3  3.19491751d+00, 3.03410521d+00, 2.87555064d+00, 2.72197200d+00, CO
     4  2.57488432d+00, 2.43369631d+00, 2.29551138d+00, 2.15653859d+00, CO
     5  2.01413327d+00, 1.86760435d+00, 1.71462672d+00, 1.56229613d+00, CO
     6  1.40411180d+00, 1.24245572d+00, 1.07203582d+00, 9.10635529d-01, CO
     7  7.49904848d-01, 5.81380503d-01, 4.17647890d-01, 2.72817240d-01, CO
     8  1.36370879d-01, 2.64389359d-02, 1.47533892d-03, 1.27326996d-02, CO
     9 -9.64185307d-03,-3.62682048d-02,-4.37118437d-02,-2.12824680d-17/ CO
      DATA Qm_NOp/
     1  4.12793386d+00, 4.07909441d+00, 3.99237755d+00, 3.88355685d+00, NO+
     2  3.75832826d+00, 3.62035670d+00, 3.47281734d+00, 3.31878813d+00, NO+
     3  3.16138106d+00, 3.00376099d+00, 2.84896188d+00, 2.69938979d+00, NO+
     4  2.55600736d+00, 2.41756797d+00, 2.28082498d+00, 2.14227550d+00, NO+
     5  1.99988594d+00, 1.85332015d+00, 1.70032203d+00, 1.54800252d+00, NO+
     6  1.38985547d+00, 1.22827365d+00, 1.05797782d+00, 8.96778005d-01, NO+
     7  7.36362989d-01, 5.68331708d-01, 4.05382059d-01, 2.61871961d-01, NO+
     8  1.27640430d-01, 2.09644034d-02, 6.93950025d-06, 1.35256156d-02, NO+
     9 -1.03972479d-02,-3.90871705d-02,-4.71057597d-02,-2.29344608d-17/ NO+
      DATA Qm_NO/
     1  4.89800724d+00, 4.84728380d+00, 4.75667892d+00, 4.64261676d+00, NO
     2  4.51126882d+00, 4.36638748d+00, 4.21092606d+00, 4.04761425d+00, NO
     3  3.87918485d+00, 3.70846377d+00, 3.53830656d+00, 3.37127397d+00, NO
     4  3.20893273d+00, 3.05084091d+00, 2.89383833d+00, 2.73281923d+00, NO
     5  2.56305218d+00, 2.38245899d+00, 2.18949156d+00, 1.99526895d+00, NO
     6  1.80096692d+00, 1.61596868d+00, 1.43167966d+00, 1.27050899d+00, NO
     7  1.11347083d+00, 9.38999569d-01, 7.70753764d-01, 6.23758316d-01, NO
     8  4.76920044d-01, 3.51957876d-01, 3.10472539d-01, 3.12116272d-01, NO
     9  2.93504590d-01, 2.72512633d-01, 2.66625999d-01, 3.01029996d-01/ NO
      DATA Qm_O2/
     1  4.78101842d+00, 4.71464304d+00, 4.59844054d+00, 4.45685785d+00, O2
     2  4.29972545d+00, 4.13207248d+00, 3.95719531d+00, 3.77824431d+00, O2
     3  3.59820965d+00, 3.41944724d+00, 3.24379760d+00, 3.07313948d+00, O2
     4  2.90950872d+00, 2.75421080d+00, 2.60644772d+00, 2.46287070d+00, O2
     5  2.31902640d+00, 2.17208710d+00, 2.01883909d+00, 1.86611621d+00, O2
     6  1.70743247d+00, 1.54506163d+00, 1.37358461d+00, 1.21060998d+00, O2
     7  1.04747134d+00, 8.75199555d-01, 7.05570125d-01, 5.51105608d-01, O2
     8  3.98399270d-01, 2.61712541d-01, 1.99242516d-01, 1.85560979d-01, O2
     9  1.71502753d-01, 1.57880075d-01, 1.53989653d-01, 1.76091259d-01/ O2
      DATA Qm_HF/
     1  2.99954398d+00, 2.94822515d+00, 2.85693539d+00, 2.74379061d+00, HF
     2  2.61702682d+00, 2.48191172d+00, 2.34212616d+00, 2.20068522d+00, HF
     3  2.06030371d+00, 1.92317691d+00, 1.79041743d+00, 1.66147438d+00, HF
     4  1.53411680d+00, 1.40549700d+00, 1.27370016d+00, 1.13827218d+00, HF
     5  9.99542577d-01, 8.58326353d-01, 7.13254335d-01, 5.72737081d-01, HF
     6  4.32269913d-01, 2.97578963d-01, 1.69765121d-01, 7.64283732d-02, HF
     7  2.20194794d-02, 2.75782968d-03,-8.70595991d-04, 5.93707491d-04, HF
     8 -4.22282782d-04,-1.23822787d-03,-2.71083434d-04, 6.46992849d-04, HF
     9 -5.22700233d-04,-1.96304462d-03,-2.36543745d-03,-1.15162565d-18/ HF
      DATA Qm_NaH/
     1  4.33405348d+00, 4.27041387d+00, 4.15549726d+00, 4.00963974d+00, NaH
     2  3.84114155d+00, 3.65521401d+00, 3.45639501d+00, 3.25035396d+00, NaH
     3  3.04337919d+00, 2.83964069d+00, 2.64058683d+00, 2.44732840d+00, NaH
     4  2.26190116d+00, 2.08655631d+00, 1.92243482d+00, 1.76823905d+00, NaH
     5  1.61974816d+00, 1.47206207d+00, 1.31971479d+00, 1.16882051d+00, NaH
     6  1.01309234d+00, 8.55316977d-01, 6.90896957d-01, 5.38806489d-01, NaH
     7  3.92550503d-01, 2.47116160d-01, 1.22062714d-01, 4.32111212d-02, NaH
     8  3.50152377d-03,-1.64862059d-02,-4.22607185d-03, 1.01266106d-02, NaH
     9 -8.18106426d-03,-3.07246739d-02,-3.70227419d-02,-1.80247164d-17/ NaH
      DATA Qm_MgH/
     1  4.46189498d+00, 4.39564667d+00, 4.27580789d+00, 4.12403578d+00, MgH
     2  3.95062071d+00, 3.76347982d+00, 3.56932573d+00, 3.37356092d+00, MgH
     3  3.17947067d+00, 2.98847655d+00, 2.80181367d+00, 2.62145627d+00, MgH
     4  2.44966532d+00, 2.28792529d+00, 2.13570885d+00, 1.98984357d+00, MgH
     5  1.84554644d+00, 1.69932781d+00, 1.54762716d+00, 1.39740643d+00, MgH
     6  1.24259503d+00, 1.08613507d+00, 9.23674403d-01, 7.74562523d-01, MgH
     7  6.32920067d-01, 4.95221296d-01, 3.84131960d-01, 3.24488174d-01, MgH
     8  3.01892478d-01, 2.92173518d-01, 2.98828230d-01, 3.06300794d-01, MgH
     9  2.96771772d-01, 2.85037878d-01, 2.81759748d-01, 3.01029996d-01/ MgH
      DATA Qm_MgO/
     1  6.06134597d+00, 5.99916727d+00, 5.88582789d+00, 5.73958811d+00, MgO
     2  5.56658224d+00, 5.36924187d+00, 5.14768897d+00, 4.90068070d+00, MgO
     3  4.62673551d+00, 4.32581257d+00, 4.00249991d+00, 3.67145987d+00, MgO
     4  3.36074957d+00, 3.09873865d+00, 2.88842028d+00, 2.71081674d+00, MgO
     5  2.54921945d+00, 2.39507779d+00, 2.23934587d+00, 2.08632798d+00, MgO
     6  1.92680691d+00, 1.76326642d+00, 1.59010178d+00, 1.42430644d+00, MgO
     7  1.25689659d+00, 1.07829069d+00, 8.98632295d-01, 7.28215808d-01, MgO
     8  5.49185408d-01, 3.65572362d-01, 2.25660262d-01, 1.16815728d-01, MgO
     9  7.99498261d-03,-9.96339338d-02,-1.40913384d-01, 9.31473670d-08/ MgO
      DATA Qm_AlH/
     1  4.19411233d+00, 4.12133771d+00, 3.98823767d+00, 3.81718484d+00, AlH
     2  3.61968136d+00, 3.40749475d+00, 3.19327774d+00, 2.98698685d+00, AlH
     3  2.79196131d+00, 2.60613554d+00, 2.42706694d+00, 2.25479414d+00, AlH
     4  2.09088649d+00, 1.93621576d+00, 1.78939420d+00, 1.64666577d+00, AlH
     5  1.50363115d+00, 1.35784037d+00, 1.20642724d+00, 1.05662168d+00, AlH
     6  9.02368892d-01, 7.46727059d-01, 5.85497631d-01, 4.38263366d-01, AlH
     7  2.99577129d-01, 1.66900729d-01, 6.47140265d-02, 1.59610884d-02, AlH
     8  3.21431435d-04,-5.69362707d-03,-1.39971940d-03, 3.34983874d-03, AlH
     9 -2.70630524d-03,-1.01637566d-02,-1.22471646d-02,-5.96259643d-18/ AlH
      DATA Qm_AlO/
     1  5.73420675d+00, 5.66722573d+00, 5.54534631d+00, 5.38876934d+00, AlO
     2  5.20534148d+00, 5.00017567d+00, 4.77776773d+00, 4.54372167d+00, AlO
     3  4.30606959d+00, 4.07445249d+00, 3.85622162d+00, 3.65283260d+00, AlO
     4  3.46154785d+00, 3.28037992d+00, 3.10963286d+00, 2.94934674d+00, AlO
     5  2.79670669d+00, 2.64680904d+00, 2.49271040d+00, 2.33953690d+00, AlO
     6  2.18011062d+00, 2.01664820d+00, 1.84354510d+00, 1.67804523d+00, AlO
     7  1.51105105d+00, 1.33287298d+00, 1.15396283d+00, 9.84814889d-01, AlO
     8  8.07695079d-01, 6.27234177d-01, 4.92798315d-01, 3.93888313d-01, AlO
     9  3.00591807d-01, 2.09669821d-01, 1.76511635d-01, 3.01030009d-01/ AlO
      DATA Qm_AlF/
     1  5.26257440d+00, 5.19561062d+00, 5.07848727d+00, 4.93607656d+00, AlF
     2  4.77865265d+00, 4.61111117d+00, 4.43519539d+00, 4.25149082d+00, AlF
     3  4.06085687d+00, 3.86484678d+00, 3.66560444d+00, 3.46579261d+00, AlF
     4  3.26857309d+00, 3.07743473d+00, 2.89565651d+00, 2.72523213d+00, AlF
     5  2.56524467d+00, 2.41169863d+00, 2.25621332d+00, 2.10315790d+00, AlF
     6  1.94362352d+00, 1.78005174d+00, 1.60683090d+00, 1.44098022d+00, AlF
     7  1.27347995d+00, 1.09470232d+00, 9.14788709d-01, 7.43983696d-01, AlF
     8  5.64317113d-01, 3.79677081d-01, 2.38017994d-01, 1.25925396d-01, AlF
     9  1.18923373d-02,-1.01431253d-01,-1.45722582d-01, 1.74483121d-07/ AlF
      DATA Qm_Al2/
     1  6.42074270d+00, 6.36394465d+00, 6.25787458d+00, 6.11655420d+00, Al2
     2  5.94428834d+00, 5.74448342d+00, 5.52219867d+00, 5.28565340d+00, Al2
     3  5.04535805d+00, 4.80941760d+00, 4.57903135d+00, 4.35075649d+00, Al2
     4  4.12233928d+00, 3.89461787d+00, 3.67030972d+00, 3.45308543d+00, Al2
     5  3.24681213d+00, 3.05484516d+00, 2.87292598d+00, 2.71080811d+00, Al2
     6  2.54815797d+00, 2.38386955d+00, 2.21025982d+00, 2.04295332d+00, Al2
     7  1.87338820d+00, 1.69235110d+00, 1.50857127d+00, 1.33130334d+00, Al2
     8  1.14191284d+00, 9.41648380d-01, 7.73626630d-01, 6.14286597d-01, Al2
     9  4.24040321d-01, 2.15611363d-01, 8.11371015d-02, 1.79670158d-01/ Al2
      DATA Qm_SiH/
     1  4.39169203d+00, 4.33331223d+00, 4.22713564d+00, 4.09218492d+00, SiH
     2  3.93837302d+00, 3.77365144d+00, 3.60356062d+00, 3.43107193d+00, SiH
     3  3.25799578d+00, 3.08637324d+00, 2.91870708d+00, 2.75745961d+00, SiH
     4  2.60421623d+00, 2.45867656d+00, 2.31813989d+00, 2.17854135d+00, SiH
     5  2.03662359d+00, 1.89134483d+00, 1.74045476d+00, 1.59140584d+00, SiH
     6  1.43821641d+00, 1.28414607d+00, 1.12528736d+00, 9.81712482d-01, SiH
     7  8.48939693d-01, 7.26290850d-01, 6.40796324d-01, 6.10070788d-01, SiH
     8  6.01821565d-01, 5.99027837d-01, 6.01342294d-01, 6.03776110d-01, SiH
     9  6.00673553d-01, 5.96853106d-01, 5.95785777d-01, 6.02059991d-01/ SiH
      DATA Qm_SiHm/
     1  4.01495803d+00, 3.97034959d+00, 3.88930594d+00, 3.78544401d+00, SiH-
     2  3.66422129d+00, 3.52936112d+00, 3.38399534d+00, 3.23111266d+00, SiH-
     3  3.07377414d+00, 2.91517336d+00, 2.75850082d+00, 2.60651380d+00, SiH-
     4  2.46075423d+00, 2.32065629d+00, 2.18339521d+00, 2.04535129d+00, SiH-
     5  1.90412138d+00, 1.75925310d+00, 1.60868403d+00, 1.45988608d+00, SiH-
     6  1.30693592d+00, 1.15311654d+00, 9.94567437d-01, 8.51402695d-01, SiH-
     7  7.19262554d-01, 5.97666058d-01, 5.13793501d-01, 4.84654034d-01, SiH-
     8  4.76829050d-01, 4.74165305d-01, 4.76425642d-01, 4.78784345d-01, SiH-
     9  4.75777658d-01, 4.72075264d-01, 4.71040916d-01, 4.77121255d-01/ SiH-
      DATA Qm_SiC/
     1  5.87867313d+00, 5.82587254d+00, 5.72923298d+00, 5.60408126d+00, SiC
     2  5.45629324d+00, 5.29003422d+00, 5.10928175d+00, 4.91841155d+00, SiC
     3  4.72216109d+00, 4.52494934d+00, 4.32996430d+00, 4.13898899d+00, SiC
     4  3.95340262d+00, 3.77518425d+00, 3.60641215d+00, 3.44750596d+00, SiC
     5  3.29569665d+00, 3.14617760d+00, 2.99223611d+00, 2.83909336d+00, SiC
     6  2.67968635d+00, 2.51621207d+00, 2.34306368d+00, 2.17748071d+00, SiC
     7  2.01034363d+00, 1.83193043d+00, 1.65264939d+00, 1.48291429d+00, SiC
     8  1.30484920d+00, 1.12283297d+00, 9.85728058d-01, 8.82008234d-01, SiC
     9  7.81221049d-01, 6.82277969d-01, 6.45349570d-01, 7.78151285d-01/ SiC
      DATA Qm_SiN/
     1  5.45648830d+00, 5.38625221d+00, 5.26094368d+00, 5.10430651d+00, SiN
     2  4.92625169d+00, 4.73331168d+00, 4.53133905d+00, 4.32653437d+00, SiN
     3  4.12445296d+00, 3.92805020d+00, 3.73737251d+00, 3.55178282d+00, SiN
     4  3.37216975d+00, 3.20070771d+00, 3.03891456d+00, 2.88586117d+00, SiN
     5  2.73759706d+00, 2.58936723d+00, 2.43572560d+00, 2.28254018d+00, SiN
     6  2.12320640d+00, 1.95986499d+00, 1.78692737d+00, 1.62173522d+00, SiN
     7  1.45520098d+00, 1.27767205d+00, 1.09979719d+00, 9.32304879d-01, SiN
     8  7.57818710d-01, 5.81773378d-01, 4.54927005d-01, 3.69036534d-01, SiN
     9  2.95143451d-01, 2.24667133d-01, 2.00286837d-01, 3.01029997d-01/ SiN
      DATA Qm_SiO/
     1  4.93333326d+00, 4.85762066d+00, 4.73261894d+00, 4.59149427d+00, SiO
     2  4.44402885d+00, 4.29026034d+00, 4.12820759d+00, 3.95772107d+00, SiO
     3  3.78048485d+00, 3.59904042d+00, 3.41634779d+00, 3.23567513d+00, SiO
     4  3.06036106d+00, 2.89319056d+00, 2.73531321d+00, 2.58515075d+00, SiO
     5  2.43840212d+00, 2.29063894d+00, 2.13708000d+00, 1.98389514d+00, SiO
     6  1.82456526d+00, 1.66122234d+00, 1.48827693d+00, 1.32307735d+00, SiO
     7  1.15652846d+00, 9.78968146d-01, 8.01046277d-01, 6.33482383d-01, SiO
     8  4.58875859d-01, 2.82626968d-01, 1.55434285d-01, 6.89720993d-02, SiO
     9 -5.74510431d-03,-7.70808512d-02,-1.01809258d-01, 1.15935281d-09/ SiO
      DATA Qm_SiF/
     1  5.76452930d+00, 5.70451434d+00, 5.59800380d+00, 5.46586964d+00, SiF
     2  5.31666132d+00, 5.15488960d+00, 4.98285177d+00, 4.80199160d+00, SiF
     3  4.61380944d+00, 4.42020223d+00, 4.22351338d+00, 4.02655369d+00, SiF
     4  3.83258409d+00, 3.64509569d+00, 3.46718517d+00, 3.30039616d+00, SiF
     5  3.14317427d+00, 2.99110768d+00, 2.83623030d+00, 2.68312865d+00, SiF
     6  2.52363108d+00, 2.36009251d+00, 2.18690281d+00, 2.02117245d+00, SiF
     7  1.85384391d+00, 1.67525527d+00, 1.49566541d+00, 1.32540267d+00, SiF
     8  1.14655828d+00, 9.63260019d-01, 8.23923415d-01, 7.16128241d-01, SiF
     9  6.08938008d-01, 5.03078442d-01, 4.62712408d-01, 6.02060068d-01/ SiF
      DATA Qm_Si2/
     1  6.35928555d+00, 6.29767441d+00, 6.18669331d+00, 6.04671130d+00, Si2
     2  5.88662202d+00, 5.71155248d+00, 5.52373886d+00, 5.32360303d+00, Si2
     3  5.11104478d+00, 4.88611882d+00, 4.64916736d+00, 4.40100634d+00, Si2
     4  4.14350248d+00, 3.88068225d+00, 3.62026469d+00, 3.37428812d+00, Si2
     5  3.15518020d+00, 2.96734274d+00, 2.79575598d+00, 2.64212773d+00, Si2
     6  2.48193048d+00, 2.31811117d+00, 2.14485723d+00, 1.97746945d+00, Si2
     7  1.80783276d+00, 1.62724739d+00, 1.44392754d+00, 1.26711941d+00, Si2
     8  1.07876591d+00, 8.80082687d-01, 7.14432073d-01, 5.59121001d-01, Si2
     9  3.75656000d-01, 1.77430463d-01, 5.76900576d-02, 1.77458636d-01/ Si2
      DATA Qm_SH/
     1  4.14906525d+00, 4.09195586d+00, 3.98935375d+00, 3.86102335d+00, SH
     2  3.71680164d+00, 3.56349216d+00, 3.40534192d+00, 3.24495364d+00, SH
     3  3.08455768d+00, 2.92659022d+00, 2.77343919d+00, 2.62673523d+00, SH
     4  2.48647740d+00, 2.35052172d+00, 2.21530853d+00, 2.07780460d+00, SH
     5  1.93669190d+00, 1.79209854d+00, 1.64210968d+00, 1.49433963d+00, SH
     6  1.34301598d+00, 1.19172645d+00, 1.03706402d+00, 9.00100633d-01, SH
     7  7.78590572d-01, 6.74570448d-01, 6.14771497d-01, 6.06095424d-01, SH
     8  6.00717587d-01, 5.97758900d-01, 6.01109935d-01, 6.04327984d-01, SH
     9  6.00227698d-01, 5.95178661d-01, 5.93768098d-01, 6.02059991d-01/ SH
      DATA Qm_SHm/
     1  3.54193175d+00, 3.48584399d+00, 3.38472823d+00, 3.25767670d+00, SH-
     2  3.11425856d+00, 2.96130725d+00, 2.80326198d+00, 2.64289194d+00, SH-
     3  2.48249763d+00, 2.32453022d+00, 2.17137920d+00, 2.02467524d+00, SH-
     4  1.88441741d+00, 1.74846173d+00, 1.61324854d+00, 1.47574461d+00, SH-
     5  1.33463191d+00, 1.19003855d+00, 1.04004969d+00, 8.92279641d-01, SH-
     6  7.40955993d-01, 5.89666461d-01, 4.35004033d-01, 2.98040642d-01, SH-
     7  1.76530581d-01, 7.25104567d-02, 1.27115056d-02, 4.03543298d-03, SH-
     8 -1.34240450d-03,-4.30109138d-03,-9.50056781d-04, 2.26799306d-03, SH-
     9 -1.83229305d-03,-6.88133043d-03,-8.29189338d-03,-4.03695187d-18/ SH-
      DATA Qm_CS/
     1  4.88712667d+00, 4.81090713d+00, 4.68327897d+00, 4.53742749d+00, CS
     2  4.38514027d+00, 4.22830350d+00, 4.06513858d+00, 3.89469477d+00, CS
     3  3.71800605d+00, 3.53740535d+00, 3.35583874d+00, 3.17657604d+00, CS
     4  3.00290898d+00, 2.83748208d+00, 2.68118849d+00, 2.53215415d+00, CS
     5  2.38596058d+00, 2.23836123d+00, 2.08484806d+00, 1.93171607d+00, CS
     6  1.77246735d+00, 1.60924924d+00, 1.43649384d+00, 1.27158594d+00, CS
     7  1.10548446d+00, 9.28614610d-01, 7.51772194d-01, 5.85908326d-01, CS
     8  4.14068043d-01, 2.42581030d-01, 1.23516775d-01, 4.98145303d-02, CS
     9 -7.43209267d-03,-6.07363554d-02,-7.84662453d-02, 8.01415524d-11/ CS
      DATA Qm_NS/
     1  5.44544822d+00, 5.38784554d+00, 5.28667840d+00, 5.16201509d+00, NS
     2  5.02092256d+00, 4.86655571d+00, 4.70073189d+00, 4.52515467d+00, NS
     3  4.34175891d+00, 4.15271131d+00, 3.96038928d+00, 3.76732701d+00, NS
     4  3.57597940d+00, 3.38814672d+00, 3.20404864d+00, 3.02158075d+00, NS
     5  2.83701036d+00, 2.64870564d+00, 2.45739074d+00, 2.27790151d+00, NS
     6  2.10542793d+00, 1.93834495d+00, 1.76432407d+00, 1.59969272d+00, NS
     7  1.43405825d+00, 1.25625403d+00, 1.07862409d+00, 9.12250224d-01, NS
     8  7.38649780d-01, 5.64340599d-01, 4.41038692d-01, 3.60704663d-01, NS
     9  2.94023817d-01, 2.30974742d-01, 2.09578366d-01, 3.01029996d-01/ NS
      DATA Qm_SO/
     1  5.50262626d+00, 5.44098046d+00, 5.33145345d+00, 5.19465597d+00, SO
     2  5.03829857d+00, 4.86688958d+00, 4.68430490d+00, 4.49466566d+00, SO
     3  4.30202981d+00, 4.10968780d+00, 3.91985782d+00, 3.73416413d+00, SO
     4  3.55447214d+00, 3.38300621d+00, 3.22120911d+00, 3.06814612d+00, SO
     5  2.91987529d+00, 2.77164263d+00, 2.61799844d+00, 2.46480812d+00, SO
     6  2.30546592d+00, 2.14211104d+00, 1.96915261d+00, 1.80392823d+00, SO
     7  1.63734439d+00, 1.45973884d+00, 1.28174421d+00, 1.11406313d+00, SO
     8  9.39271546d-01, 7.62708512d-01, 6.34975737d-01, 5.47619343d-01, SO
     9  4.71613264d-01, 3.98932366d-01, 3.73655359d-01, 4.77121256d-01/ SO
      DATA Qm_MgS/
     1  5.75604548d+00, 5.69476131d+00, 5.58406199d+00, 5.44400743d+00, MgS
     2  5.28379141d+00, 5.10970643d+00, 4.92574198d+00, 4.73386568d+00, MgS
     3  4.53493451d+00, 4.32977684d+00, 4.11973766d+00, 3.90675536d+00, MgS
     4  3.69336976d+00, 3.48273702d+00, 3.27852021d+00, 3.08444748d+00, MgS
     5  2.90305835d+00, 2.73443192d+00, 2.57054535d+00, 2.41651786d+00, MgS
     6  2.25647132d+00, 2.09257732d+00, 1.91902125d+00, 1.75204725d+00, MgS
     7  1.58292082d+00, 1.40225567d+00, 1.21916134d+00, 1.04308783d+00, MgS
     8  8.55420292d-01, 6.57930890d-01, 4.94520426d-01, 3.42918704d-01, MgS
     9  1.65113618d-01,-2.49762251d-02,-1.33715294d-01, 5.98134506d-04/ MgS
      DATA Qm_AlS/
     1  5.94686136d+00, 5.88585994d+00, 5.77645108d+00, 5.63916759d+00, AlS
     2  5.48312096d+00, 5.31400979d+00, 5.13504377d+00, 4.94769826d+00, AlS
     3  4.75286173d+00, 4.55171689d+00, 4.34599479d+00, 4.13796846d+00, AlS
     4  3.93047306d+00, 3.72688617d+00, 3.53089571d+00, 3.34583221d+00, AlS
     5  3.17316628d+00, 3.01135446d+00, 2.85165161d+00, 2.69839545d+00, AlS
     6  2.53854802d+00, 2.37466024d+00, 2.20104341d+00, 2.03419951d+00, AlS
     7  1.86523325d+00, 1.68458206d+00, 1.50159898d+00, 1.32579845d+00, AlS
     8  1.13843565d+00, 9.41457503d-01, 7.78960190d-01, 6.28897461d-01, AlS
     9  4.53485990d-01, 2.66723052d-01, 1.62175545d-01, 3.01454659d-01/ AlS
      DATA Qm_SiS/
     1  5.48323443d+00, 5.41715542d+00, 5.30430126d+00, 5.17050859d+00, SiS
     2  5.02406852d+00, 4.86667168d+00, 4.69826981d+00, 4.51954689d+00, SiS
     3  4.33207827d+00, 4.13788727d+00, 3.93931987d+00, 3.73912806d+00, SiS
     4  3.54053362d+00, 3.34712523d+00, 3.16241661d+00, 2.98887218d+00, SiS
     5  2.82626007d+00, 2.67114232d+00, 2.51490546d+00, 2.36178080d+00, SiS
     6  2.20203331d+00, 2.03814724d+00, 1.86445957d+00, 1.69780820d+00, SiS
     7  1.52908938d+00, 1.34851371d+00, 1.16576843d+00, 9.90472875d-01, SiS
     8  8.03741021d-01, 6.07802521d-01, 4.47114066d-01, 3.00143714d-01, SiS
     9  1.29598458d-01,-5.05790210d-02,-1.47321463d-01, 2.14205530d-04/ SiS
      DATA Qm_S2/
     1  5.86889648d+00, 5.79614754d+00, 5.66806995d+00, 5.51182236d+00, S2
     2  5.33943349d+00, 5.15670309d+00, 4.96544499d+00, 4.76636439d+00, S2
     3  4.56118078d+00, 4.35288180d+00, 4.14465632d+00, 3.93889481d+00, S2
     4  3.73728945d+00, 3.54180724d+00, 3.35514140d+00, 3.17969456d+00, S2
     5  3.01551215d+00, 2.85943918d+00, 2.70275484d+00, 2.54964499d+00, S2
     6  2.38988443d+00, 2.22599356d+00, 2.05230895d+00, 1.88561424d+00, S2
     7  1.71683599d+00, 1.53621652d+00, 1.35338400d+00, 1.17793063d+00, S2
     8  9.90975630d-01, 7.94676238d-01, 6.33378300d-01, 4.85355726d-01, S2
     9  3.13133410d-01, 1.30712335d-01, 3.13889833d-02, 1.76361725d-01/ S2
      DATA Qm_HCl/
     1  3.40870792d+00, 3.35562881d+00, 3.26080519d+00, 3.14252139d+00, HCl
     2  3.00896539d+00, 2.86536243d+00, 2.71533390d+00, 2.56186082d+00, HCl
     3  2.40782426d+00, 2.25603181d+00, 2.10885443d+00, 1.96755549d+00, HCl
     4  1.83154475d+00, 1.69823251d+00, 1.56421721d+00, 1.42716121d+00, HCl
     5  1.28634408d+00, 1.14211504d+00, 9.92633904d-01, 8.45590155d-01, HCl
     6  6.95337685d-01, 5.45653544d-01, 3.93420180d-01, 2.60350543d-01, HCl
     7  1.45727412d-01, 5.26613886d-02, 4.80512121d-03, 3.75311913d-03, HCl
     8 -1.91409628d-03,-5.73792999d-03,-1.25909165d-03, 3.00523647d-03, HCl
     9 -2.42790598d-03,-9.11820480d-03,-1.09872913d-02,-5.34922051d-18/ HCl
      DATA Qm_LiCl/
     1  5.28461641d+00, 5.22784859d+00, 5.12275990d+00, 4.98547155d+00, LiCl
     2  4.82367198d+00, 4.64447296d+00, 4.45427041d+00, 4.25729714d+00, LiCl
     3  4.05541169d+00, 3.84947985d+00, 3.64066587d+00, 3.43085649d+00, LiCl
     4  3.22266404d+00, 3.01931490d+00, 2.82433823d+00, 2.64081581d+00, LiCl
     5  2.46981581d+00, 2.30931786d+00, 2.15043409d+00, 1.99741677d+00, LiCl
     6  1.83790624d+00, 1.67454295d+00, 1.50174029d+00, 1.33622744d+00, LiCl
     7  1.16930913d+00, 9.91762390d-01, 8.13645666d-01, 6.45536224d-01, LiCl
     8  4.70367620d-01, 2.93099204d-01, 1.63969906d-01, 7.43326041d-02, LiCl
     9 -4.79803769d-03,-8.07314005d-02,-1.07367430d-01, 2.16536226d-09/ LiCl
      DATA Qm_NaCl/
     1  6.05484882d+00, 6.00044430d+00, 5.89891560d+00, 5.76441400d+00, NaCl
     2  5.60273344d+00, 5.41936595d+00, 5.22039179d+00, 5.01160297d+00, NaCl
     3  4.79709126d+00, 4.57864134d+00, 4.35664700d+00, 4.13162587d+00, NaCl
     4  3.90495348d+00, 3.67888944d+00, 3.45650696d+00, 3.24157886d+00, NaCl
     5  3.03793582d+00, 2.84869150d+00, 2.66914903d+00, 2.50844762d+00, NaCl
     6  2.34634955d+00, 2.18217376d+00, 2.00860414d+00, 1.84134090d+00, NaCl
     7  1.67183570d+00, 1.49089476d+00, 1.30725966d+00, 1.13021393d+00, NaCl
     8  9.41176977d-01, 7.41474174d-01, 5.74360421d-01, 4.16535057d-01, NaCl
     9  2.28762767d-01, 2.41344825d-02,-1.04634091d-01, 2.50427115d-03/ NaCl
      DATA Qm_AlCl/
     1  5.87385640d+00, 5.80471391d+00, 5.68222026d+00, 5.53156670d+00, AlCl
     2  5.36438975d+00, 5.18716180d+00, 5.00236756d+00, 4.81026078d+00, AlCl
     3  4.61080639d+00, 4.40463365d+00, 4.19306147d+00, 3.97795903d+00, AlCl
     4  3.76175513d+00, 3.54749441d+00, 3.33879714d+00, 3.13954184d+00, AlCl
     5  2.95278357d+00, 2.77953550d+00, 2.61249319d+00, 2.45753898d+00, AlCl
     6  2.29725264d+00, 2.13332298d+00, 1.95976103d+00, 1.79266651d+00, AlCl
     7  1.62337444d+00, 1.44256974d+00, 1.25921379d+00, 1.08268099d+00, AlCl
     8  8.94351070d-01, 6.95793789d-01, 5.30602168d-01, 3.75987151d-01, AlCl
     9  1.93356369d-01,-3.62148156d-03,-1.21468206d-01, 1.19123479d-03/ AlCl
      DATA Qm_CaH/
     1  4.73726728d+00, 4.66251358d+00, 4.52672222d+00, 4.35438645d+00, CaH
     2  4.15878795d+00, 3.95212925d+00, 3.74483446d+00, 3.54265546d+00, CaH
     3  3.34588531d+00, 3.15264928d+00, 2.96255670d+00, 2.77728828d+00, CaH
     4  2.59942258d+00, 2.43118671d+00, 2.27318555d+00, 2.12333917d+00, CaH
     5  1.97706889d+00, 1.83001625d+00, 1.67765318d+00, 1.52646901d+00, CaH
     6  1.37024920d+00, 1.21167555d+00, 1.04600309d+00, 8.91965879d-01, CaH
     7  7.42674159d-01, 5.92222257d-01, 4.58228152d-01, 3.65483132d-01, CaH
     8  3.09563861d-01, 2.78760963d-01, 2.95068337d-01, 3.15344213d-01, CaH
     9  2.89466465d-01, 2.57602183d-01, 2.48700175d-01, 3.01029996d-01/ CaH
      DATA Qm_CaF/
     1  5.99137474d+00, 5.91955936d+00, 5.79167613d+00, 5.63374902d+00, CaF
     2  5.45879210d+00, 5.27519172d+00, 5.08689719d+00, 4.89434395d+00, CaF
     3  4.69662252d+00, 4.49343181d+00, 4.28574827d+00, 4.07556495d+00, CaF
     4  3.86558414d+00, 3.65909981d+00, 3.45980289d+00, 3.27120392d+00, CaF
     5  3.09521756d+00, 2.93092608d+00, 2.76977985d+00, 2.61634526d+00, CaF
     6  2.45648528d+00, 2.29267402d+00, 2.11920512d+00, 1.95250280d+00, CaF
     7  1.78377205d+00, 1.60357793d+00, 1.42127090d+00, 1.24649305d+00, CaF
     8  1.06079647d+00, 8.66461519d-01, 7.08287629d-01, 5.65775875d-01, CaF
     9  4.02573234d-01, 2.31944440d-01, 1.45347586d-01, 3.01109593d-01/ CaF
      DATA Qm_CaCl/
     1  6.57733019d+00, 6.50103170d+00, 6.36486194d+00, 6.19650504d+00, CaCl
     2  6.01060541d+00, 5.81741464d+00, 5.62190359d+00, 5.42390965d+00, CaCl
     3  5.22096548d+00, 5.01144344d+00, 4.79556499d+00, 4.57478143d+00, CaCl
     4  4.35120173d+00, 4.12751144d+00, 3.90704451d+00, 3.69373812d+00, CaCl
     5  3.49150903d+00, 3.30346530d+00, 3.12485062d+00, 2.96461503d+00, CaCl
     6  2.80265895d+00, 2.63844224d+00, 2.46475481d+00, 2.29730119d+00, CaCl
     7  2.12749273d+00, 1.94607399d+00, 1.76169449d+00, 1.58348214d+00, CaCl
     8  1.39260048d+00, 1.18995970d+00, 1.01809659d+00, 8.52424277d-01, CaCl
     9  6.51740703d-01, 4.26344792d-01, 2.64185602d-01, 3.17225685d-01/ CaCl
      DATA Qm_ScO/
     1  5.59512884d+00, 5.52185551d+00, 5.39174989d+00, 5.23167239d+00, ScO
     2  5.05539351d+00, 4.87258556d+00, 4.68859609d+00, 4.50428272d+00, ScO
     3  4.31799343d+00, 4.12857387d+00, 3.93690933d+00, 3.74555188d+00, ScO
     4  3.55783229d+00, 3.37717072d+00, 3.20622489d+00, 3.04562772d+00, ScO
     5  2.89273853d+00, 2.74271001d+00, 2.58853619d+00, 2.43529990d+00, ScO
     6  2.27576309d+00, 2.11213131d+00, 1.93877112d+00, 1.77287249d+00, ScO
     7  1.60526750d+00, 1.42615239d+00, 1.24577988d+00, 1.07432549d+00, ScO
     8  8.93531427d-01, 7.07094457d-01, 5.62425349d-01, 4.44680081d-01, ScO
     9  3.21453908d-01, 1.97963171d-01, 1.47919688d-01, 3.01030520d-01/ ScO
      DATA Qm_TiO/
     1  6.12255387d+00, 6.04923119d+00, 5.91822944d+00, 5.75500971d+00, TiO
     2  5.57171780d+00, 5.37704785d+00, 5.17696947d+00, 4.97446135d+00, TiO
     3  4.77014048d+00, 4.56409989d+00, 4.35741537d+00, 4.15206213d+00, TiO
     4  3.94989554d+00, 3.75193648d+00, 3.55800141d+00, 3.36597992d+00, TiO
     5  3.17103843d+00, 2.96763270d+00, 2.75120004d+00, 2.53367971d+00, TiO
     6  2.31864474d+00, 2.11721102d+00, 1.92011325d+00, 1.75308784d+00, TiO
     7  1.59187259d+00, 1.40753010d+00, 1.22553266d+00, 1.05844963d+00, TiO
     8  8.75957565d-01, 6.88665980d-01, 5.47979132d-01, 4.35632362d-01, TiO
     9  3.15017194d-01, 1.94093563d-01, 1.46522767d-01, 3.01030273d-01/ TiO
      DATA Qm_TiS/
     1  6.64006381d+00, 6.57488741d+00, 6.45780102d+00, 6.31083904d+00, TiS
     2  6.14421785d+00, 5.96476833d+00, 5.77673049d+00, 5.58212154d+00, TiS
     3  5.38135386d+00, 5.17412689d+00, 4.96044249d+00, 4.74119206d+00, TiS
     4  4.51805610d+00, 4.29315095d+00, 4.06871254d+00, 3.84662116d+00, TiS
     5  3.62735429d+00, 3.40950939d+00, 3.18618271d+00, 2.96784821d+00, TiS
     6  2.74988957d+00, 2.54442228d+00, 2.34364743d+00, 2.17382354d+00, TiS
     7  2.01044814d+00, 1.82382239d+00, 1.63808695d+00, 1.46482758d+00, TiS
     8  1.27295183d+00, 1.07060910d+00, 9.04620630d-01, 7.47030657d-01, TiS
     9  5.53340726d-01, 3.40447830d-01, 2.03389876d-01, 3.04980694d-01/ TiS
      DATA Qm_VO/
     1  5.82740434d+00, 5.75910069d+00, 5.63699211d+00, 5.48485098d+00, VO
     2  5.31423823d+00, 5.13361768d+00, 4.94886865d+00, 4.76280364d+00, VO
     3  4.57557974d+00, 4.38654710d+00, 4.19618058d+00, 4.00663347d+00, VO
     4  3.82109999d+00, 3.64290223d+00, 3.47446568d+00, 3.31602511d+00, VO
     5  3.16454491d+00, 3.01511386d+00, 2.86113146d+00, 2.70788957d+00, VO
     6  2.54839022d+00, 2.38480601d+00, 2.21150972d+00, 2.04573113d+00, VO
     7  1.87830479d+00, 1.69943941d+00, 1.51946492d+00, 1.34864768d+00, VO
     8  1.16885212d+00, 9.84026192d-01, 8.42115717d-01, 7.29528544d-01, VO
     9  6.14624692d-01, 5.00340814d-01, 4.55537985d-01, 6.02060184d-01/ VO
      DATA Qm_CrO/
     1  6.21634269d+00, 6.15594791d+00, 6.04668304d+00, 5.90853675d+00, CrO
     2  5.75134869d+00, 5.58223349d+00, 5.40550441d+00, 5.22283136d+00, CrO
     3  5.03472034d+00, 4.84208884d+00, 4.64683317d+00, 4.45167065d+00, CrO
     4  4.25985926d+00, 4.07484447d+00, 3.89955349d+00, 3.73518630d+00, CrO
     5  3.57975630d+00, 3.42862556d+00, 3.27414087d+00, 3.12104020d+00, CrO
     6  2.96155691d+00, 2.79798396d+00, 2.62470321d+00, 2.45884520d+00, CrO
     7  2.29129831d+00, 2.11232740d+00, 1.93214822d+00, 1.76096610d+00, CrO
     8  1.58063675d+00, 1.39493821d+00, 1.25150540d+00, 1.13608595d+00, CrO
     9  1.01665113d+00, 8.97371394d-01, 8.49767303d-01, 1.00000034d+00/ CrO
      DATA Qm_FeO/
     1  6.23435069d+00, 6.17492990d+00, 6.06751110d+00, 5.93120650d+00, FeO
     2  5.77439868d+00, 5.60277955d+00, 5.42036195d+00, 5.22998289d+00, FeO
     3  5.03404357d+00, 4.83529110d+00, 4.63703580d+00, 4.44255685d+00, FeO
     4  4.25445057d+00, 4.07469724d+00, 3.90483237d+00, 3.74505393d+00, FeO
     5  3.59266687d+00, 3.44288578d+00, 3.28883420d+00, 3.13564463d+00, FeO
     6  2.97614808d+00, 2.81254240d+00, 2.63919879d+00, 2.47330962d+00, FeO
     7  2.30570604d+00, 2.12658399d+00, 1.94619542d+00, 1.77471172d+00, FeO
     8  1.59386833d+00, 1.40735052d+00, 1.26254163d+00, 1.14453336d+00, FeO
     9  1.02088105d+00, 8.96915027d-01, 8.46588759d-01, 1.00000055d+00/ FeO
      DATA Qm_YO/
     1  5.65923907d+00, 5.59476746d+00, 5.48039006d+00, 5.33919016d+00, YO
     2  5.18179284d+00, 5.01464532d+00, 4.84087113d+00, 4.66118213d+00, YO
     3  4.47551478d+00, 4.28451607d+00, 4.09003694d+00, 3.89485733d+00, YO
     4  3.70232056d+00, 3.51598996d+00, 3.33900224d+00, 3.17290362d+00, YO
     5  3.01611927d+00, 2.86423976d+00, 2.70938927d+00, 2.55618376d+00, YO
     6  2.39652875d+00, 2.23273801d+00, 2.05915728d+00, 1.89283368d+00, YO
     7  1.72459162d+00, 1.54458440d+00, 1.36278855d+00, 1.18905585d+00, YO
     8  1.00470312d+00, 8.12571611d-01, 6.58251917d-01, 5.22592083d-01, YO
     9  3.70261235d-01, 2.13120449d-01, 1.39142151d-01, 3.01048879d-01/ YO
      DATA Qm_ZrO/
     1  6.16061629d+00, 6.08908401d+00, 5.96037070d+00, 5.79792397d+00, ZrO
     2  5.61195024d+00, 5.40907889d+00, 5.19317519d+00, 4.96532856d+00, ZrO
     3  4.72448539d+00, 4.46880195d+00, 4.19709375d+00, 3.91039640d+00, ZrO
     4  3.61472398d+00, 3.32515206d+00, 3.06540854d+00, 2.85287096d+00, ZrO
     5  2.68046398d+00, 2.52624825d+00, 2.37186530d+00, 2.21836204d+00, ZrO
     6  2.05875128d+00, 1.89499297d+00, 1.72142664d+00, 1.55526708d+00, ZrO
     7  1.38725323d+00, 1.20745217d+00, 1.02603395d+00, 8.52957503d-01, ZrO
     8  6.69563348d-01, 4.78975396d-01, 3.27292876d-01, 1.96448146d-01, ZrO
     9  5.18864173d-02,-9.59872707d-02,-1.62445592d-01, 7.04188250d-06/ ZrO
      DATA Qm_LaO/
     1  5.96131674d+00, 5.88853337d+00, 5.75699985d+00, 5.59018407d+00, LaO
     2  5.39895853d+00, 5.19220739d+00, 4.97823171d+00, 4.76416353d+00, LaO
     3  4.55447665d+00, 4.35009272d+00, 4.14948864d+00, 3.95131278d+00, LaO
     4  3.75631887d+00, 3.56721958d+00, 3.38718826d+00, 3.21817557d+00, LaO
     5  3.05916993d+00, 2.90611699d+00, 2.75079292d+00, 2.59762599d+00, LaO
     6  2.43793502d+00, 2.27410401d+00, 2.10047546d+00, 1.93401843d+00, LaO
     7  1.76558097d+00, 1.58533314d+00, 1.40313905d+00, 1.22875430d+00, LaO
     8  1.04340393d+00, 8.49677045d-01, 6.92679072d-01, 5.52199415d-01, LaO
     9  3.92086036d-01, 2.25363778d-01, 1.42688202d-01, 3.01082140d-01/ LaO
C
C Equilibrium constants
C
      DATA Kp_H2p/
     1  8.77008105d+00, 8.67869030d+00, 8.50147737d+00, 8.25021225d+00, H2+
     2  7.91549511d+00, 7.47697477d+00, 6.90183299d+00, 6.14053333d+00, H2+
     3  5.12175022d+00, 3.74746302d+00, 1.88476758d+00,-6.51545961d-01, H2+
     4 -4.12680821d+00,-8.92548062d+00,-1.56093360d+01,-2.50074459d+01, H2+
     5 -3.83536754d+01,-5.74948691d+01,-8.52154540d+01,-1.25720114d+02, H2+
     6 -1.85451966d+02,-2.74295722d+02,-4.07564739d+02,-6.09076670d+02, H2+
     7 -9.16232567d+02,-1.38808924d+03,-2.11851356d+03,-3.25778937d+03, H2+
     8 -5.04816132d+03,-7.88242206d+03,-1.24016429d+04,-1.96592360d+04, H2+
     9 -3.13968079d+04,-5.05117370d+04,-8.18544981d+04,-1.33595888d+05/ H2+
      DATA Kp_H2/
     1  9.01858146d+00, 8.87784403d+00, 8.60111064d+00, 8.20109459d+00, H2
     2  7.65685320d+00, 6.92946534d+00, 5.96101222d+00, 4.67067585d+00, H2
     3  2.94757498d+00, 6.37837965d-01,-2.47576540d+00,-6.70343321d+00, H2
     4 -1.24922406d+01,-2.04927241d+01,-3.16590123d+01,-4.73988824d+01, H2
     5 -6.98019718d+01,-1.01990621d+02,-1.48667426d+02,-2.16944190d+02, H2
     6 -3.17702311d+02,-4.67646492d+02,-6.92621740d+02,-1.03285762d+03, H2
     7 -1.55150867d+03,-2.34832675d+03,-3.58195396d+03,-5.50630200d+03, H2
     8 -8.53054863d+03,-1.33182790d+04,-2.09525744d+04,-3.32129739d+04, H2
     9 -5.30415932d+04,-8.53331691d+04,-1.38282047d+05,-2.25691895d+05/ H2
      DATA Kp_H2m/
     1  9.91747880d+00, 9.86255426d+00, 9.75561706d+00, 9.60276114d+00, H2-
     2  9.39656020d+00, 9.12213747d+00, 8.75721080d+00, 8.27180717d+00, H2-
     3  7.62699673d+00, 6.77084157d+00, 5.63119103d+00, 4.10559692d+00, H2-
     4  2.04648792d+00,-7.61844876d-01,-4.63770458d+00,-1.00524214d+01, H2-
     5 -1.77064747d+01,-2.86466564d+01,-4.44522856d+01,-6.75018528d+01, H2-
     6 -1.01452812d+02,-1.51907491d+02,-2.27529234d+02,-3.41774558d+02, H2-
     7 -5.15794181d+02,-7.82994247d+02,-1.19651864d+03,-1.84140634d+03, H2-
     8 -2.85473330d+03,-4.45878323d+03,-7.01628281d+03,-1.11233486d+04, H2-
     9 -1.77655459d+04,-2.85824244d+04,-4.63186554d+04,-7.55977287d+04/ H2-
      DATA Kp_CH/
     1  8.83588004d+00, 8.73303773d+00, 8.52937624d+00, 8.23140063d+00, CH
     2  7.82033778d+00, 7.26467925d+00, 6.52002451d+00, 5.52481700d+00, CH
     3  4.19248931d+00, 2.40108315d+00,-2.05995618d-02,-3.31435530d+00, CH
     4 -7.82642185d+00,-1.40591028d+01,-2.27484776d+01,-3.49819098d+01, CH
     5 -5.23759279d+01,-7.73473234d+01,-1.13536967d+02,-1.66446555d+02, CH
     6 -2.44497488d+02,-3.60611349d+02,-5.34792844d+02,-7.98168009d+02, CH
     7 -1.19961787d+03,-1.81633374d+03,-2.77102718d+03,-4.26013743d+03, CH
     8 -6.60029312d+03,-1.03049679d+04,-1.62121879d+04,-2.56989046d+04, CH
     9 -4.10416357d+04,-6.60277352d+04,-1.06997623d+05,-1.74631963d+05/ CH
      DATA Kp_CHm/
     1  8.53485004d+00, 8.43200773d+00, 8.22834625d+00, 7.93037064d+00, CH-
     2  7.51930779d+00, 6.96364925d+00, 6.21899451d+00, 5.22378700d+00, CH-
     3  3.89145931d+00, 2.10005316d+00,-3.21629557d-01,-3.61538529d+00, CH-
     4 -8.12745185d+00,-1.43601328d+01,-2.30495076d+01,-3.52829398d+01, CH-
     5 -5.26769579d+01,-7.76483534d+01,-1.13837997d+02,-1.66747585d+02, CH-
     6 -2.44798518d+02,-3.60912379d+02,-5.35093874d+02,-7.98469039d+02, CH-
     7 -1.19991890d+03,-1.81663477d+03,-2.77132821d+03,-4.26043846d+03, CH-
     8 -6.60059415d+03,-1.03052689d+04,-1.62124889d+04,-2.56992056d+04, CH-
     9 -4.10419367d+04,-6.60280362d+04,-1.06997924d+05,-1.74632264d+05/ CH-
      DATA Kp_C2/
     1  8.42890535d+00, 8.24225926d+00, 7.87343649d+00, 7.33757462d+00, C2
     2  6.60677082d+00, 5.63122989d+00, 4.33675180d+00, 2.61668664d+00, C2
     3  3.19465062d-01,-2.76931325d+00,-6.95310351d+00,-1.26642942d+01, C2
     4 -2.05251337d+01,-3.14394536d+01,-4.67313819d+01,-6.83588649d+01, C2
     5 -9.92407874d+01,-1.43739264d+02,-2.08374935d+02,-3.03014802d+02, C2
     6 -4.42699791d+02,-6.50591153d+02,-9.62528733d+02,-1.43429964d+03, C2
     7 -2.15347368d+03,-3.25835503d+03,-4.96888394d+03,-7.63716492d+03, C2
     8 -1.18306947d+04,-1.84697196d+04,-2.90562815d+04,-4.60580838d+04, C2
     9 -7.35549836d+04,-1.18334714d+05,-1.91760776d+05,-3.12975396d+05/ C2
      DATA Kp_C2m/
     1  7.66671368d+00, 7.41075895d+00, 6.90363624d+00, 6.16349824d+00, C2-
     2  5.14884348d+00, 3.79008074d+00, 1.98888659d+00,-3.93756672d-01, C2-
     3 -3.55943345d+00,-7.80060272d+00,-1.35361228d+01,-2.13630448d+01, C2-
     4 -3.21375378d+01,-4.70987392d+01,-6.80556818d+01,-9.76715921d+01, C2-
     5 -1.39899272d+02,-2.00648359d+02,-2.88818108d+02,-4.17879886d+02, C2-
     6 -6.08418555d+02,-8.92049358d+02,-1.31770657d+03,-1.96159418d+03, C2-
     7 -2.94331266d+03,-4.45176324d+03,-6.78731793d+03,-1.04308103d+04, C2-
     8 -1.61571022d+04,-2.52227970d+04,-3.96790071d+04,-6.28954685d+04, C2-
     9 -1.00443366d+05,-1.61591610d+05,-2.61857571d+05,-4.27380821d+05/ C2-
      DATA Kp_CN/
     1  7.66305372d+00, 7.42687696d+00, 6.96281156d+00, 6.29322373d+00, CN
     2  5.38452728d+00, 4.17270504d+00, 2.56036149d+00, 4.08748229d-01, CN
     3 -2.47487131d+00,-6.35841738d+00,-1.16199379d+01,-1.88005207d+01, CN
     4 -2.86807369d+01,-4.23929495d+01,-6.15914300d+01,-8.87129203d+01, CN
     5 -1.27375324d+02,-1.82987442d+02,-2.63693771d+02,-3.81821428d+02, CN
     6 -5.56208784d+02,-8.15786980d+02,-1.20533645d+03,-1.79458852d+03, CN
     7 -2.69298881d+03,-4.07339880d+03,-6.21068525d+03,-9.54485518d+03, CN
     8 -1.47849874d+04,-2.30809957d+04,-3.63098378d+04,-5.75551356d+04, CN
     9 -9.19150476d+04,-1.47871512d+05,-2.39624375d+05,-3.91093796d+05/ CN
      DATA Kp_CNm/
     1  6.86275642d+00, 6.56099547d+00, 5.96346007d+00, 5.09279572d+00, CN-
     2  3.90053414d+00, 2.30043869d+00, 1.64710483d-01,-2.68776067d+00, CN-
     3 -6.51011482d+00,-1.16568481d+01,-1.86297431d+01,-2.81468695d+01, CN-
     4 -4.12435007d+01,-5.94217743d+01,-8.48778114d+01,-1.20848976d+02, CN-
     5 -1.72143760d+02,-2.45950639d+02,-3.53090304d+02,-5.09943453d+02, CN-
     6 -7.41533941d+02,-1.08629925d+03,-1.60373148d+03,-2.38648617d+03, CN-
     7 -3.57997436d+03,-5.41386556d+03,-8.25335849d+03,-1.26830561d+04, CN-
     8 -1.96450222d+04,-3.06670314d+04,-4.82428393d+04,-7.64693639d+04, CN-
     9 -1.22120044d+05,-1.96464047d+05,-3.18367457d+05,-5.19610889d+05/ CN-
      DATA Kp_NH/
     1  8.82957790d+00, 8.72116433d+00, 8.50745530d+00, 8.19697742d+00, NH
     2  7.77194252d+00, 7.20124873d+00, 6.44036535d+00, 5.42765208d+00, NH
     3  4.07749904d+00, 2.26993743d+00,-1.64678420d-01,-3.46806027d+00, NH
     4 -7.98722509d+00,-1.42252984d+01,-2.29194203d+01,-3.51584172d+01, NH
     5 -5.25602714d+01,-7.75433099d+01,-1.13750726d+02,-1.66688249d+02, NH
     6 -2.44784671d+02,-3.60975692d+02,-5.35294011d+02,-7.98917262d+02, NH
     7 -1.20081602d+03,-1.81832180d+03,-2.77432653d+03,-4.26556087d+03, NH
     8 -6.60908961d+03,-1.03191094d+04,-1.62348525d+04,-2.57352576d+04, NH
     9 -4.11001278d+04,-6.61222817d+04,-1.07151289d+05,-1.74883224d+05/ NH
      DATA Kp_N2/
     1  7.28668052d+00, 6.97249560d+00, 6.35866005d+00, 5.48086617d+00, N2
     2  4.30196395d+00, 2.74594067d+00, 6.94459925d-01,-2.02343325d+00, N2
     3 -5.64832739d+00,-1.05180844d+01,-1.71108508d+01,-2.61091900d+01, N2
     4 -3.84949329d+01,-5.56912466d+01,-7.97780153d+01,-1.13819980d+02, N2
     5 -1.62367058d+02,-2.32219898d+02,-3.33616938d+02,-4.82059010d+02, N2
     6 -7.01229378d+02,-1.02750796d+03,-1.51720857d+03,-2.25804226d+03, N2
     7 -3.38767197d+03,-5.12353157d+03,-7.81132986d+03,-1.20044587d+04, N2
     8 -1.85946424d+04,-2.90280531d+04,-4.56652643d+04,-7.23843960d+04, N2
     9 -1.15597152d+05,-1.85970884d+05,-3.01364140d+05,-4.91860297d+05/ N2
      DATA Kp_OH/
     1  8.73565022d+00, 8.59937397d+00, 8.32916057d+00, 7.93467537d+00, OH
     2  7.39409179d+00, 6.67024126d+00, 5.70928398d+00, 4.43481442d+00, OH
     3  2.73885688d+00, 4.69196401d-01,-2.58929189d+00,-6.74276262d+00, OH
     4 -1.24310368d+01,-2.02925846d+01,-3.12629492d+01,-4.67224812d+01, OH
     5 -6.87198771d+01,-1.00315024d+02,-1.46114872d+02,-2.13084436d+02, OH
     6 -3.11884551d+02,-4.58888075d+02,-6.79453997d+02,-1.01304965d+03, OH
     7 -1.52165340d+03,-2.30313845d+03,-3.51304511d+03,-5.40038437d+03, OH
     8 -8.36647661d+03,-1.30621347d+04,-2.05496020d+04,-3.25741975d+04, OH
     9 -5.20214581d+04,-8.36919723d+04,-1.35622472d+05,-2.21351123d+05/ OH
      DATA Kp_OHm/
     1  9.00672467d+00, 8.84835009d+00, 8.53792284d+00, 8.09172461d+00, OH-
     2  7.48964985d+00, 6.69337900d+00, 5.64501162d+00, 4.26147206d+00, OH-
     3  2.42548102d+00,-2.79918207d-02,-3.33203501d+00,-7.81834514d+00, OH-
     4 -1.39634788d+01,-2.24592303d+01,-3.43190301d+01,-5.10376538d+01, OH-
     5 -7.48332433d+01,-1.09018721d+02,-1.58581780d+02,-2.31064389d+02, OH-
     6 -3.38007981d+02,-4.97138743d+02,-7.35910113d+02,-1.09705293d+03, OH-
     7 -1.64766523d+03,-2.49370650d+03,-3.80357500d+03,-5.84686545d+03, OH-
     8 -9.05806607d+03,-1.41417787d+04,-2.22480491d+04,-3.52664425d+04, OH-
     9 -5.63209775d+04,-9.06090122d+04,-1.46831531d+05,-2.39645658d+05/ OH-
      DATA Kp_BO/
     1  7.59132147d+00, 7.34457388d+00, 6.85553199d+00, 6.14186738d+00, BO
     2  5.16436283d+00, 3.85556636d+00, 2.11650634d+00,-1.94430176d-01, BO
     3 -3.27860236d+00,-7.42172161d+00,-1.30299801d+01,-2.06838681d+01, BO
     4 -3.12184134d+01,-4.58433746d+01,-6.63250723d+01,-9.52638500d+01, BO
     5 -1.36518801d+02,-1.95858055d+02,-2.81966110d+02,-4.07987842d+02, BO
     6 -5.94015595d+02,-8.70916585d+02,-1.28647870d+03,-1.91511516d+03, BO
     7 -2.87362346d+03,-4.34647130d+03,-6.62693363d+03,-1.01844987d+04, BO
     8 -1.57757464d+04,-2.46276396d+04,-3.87429132d+04,-6.14118284d+04, BO
     9 -9.80741751d+04,-1.57780260d+05,-2.55681471d+05,-4.17300878d+05/ BO
      DATA Kp_CO/
     1  6.89095124d+00, 6.55234401d+00, 5.88324786d+00, 4.91333600d+00, CO
     2  3.59505700d+00, 1.84058597d+00,-4.82471355d-01,-3.56461569d+00, CO
     3 -7.67631270d+00,-1.32006681d+01,-2.06817815d+01,-3.08966952d+01, CO
     4 -4.49629765d+01,-6.45003244d+01,-9.18750962d+01,-1.30572648d+02, CO
     5 -1.85765679d+02,-2.65184266d+02,-3.80464084d+02,-5.49221738d+02, CO
     6 -7.98373254d+02,-1.16927090d+03,-1.72592430d+03,-2.56801889d+03, CO
     7 -3.85199891d+03,-5.82495945d+03,-8.87979663d+03,-1.36454556d+04, CO
     8 -2.11354524d+04,-3.29934374d+04,-5.19023162d+04,-8.22697509d+04, CO
     9 -1.31382934d+05,-2.11365781d+05,-3.42515359d+05,-5.59022821d+05/ CO
      DATA Kp_NOp/
     1  7.07693495d+00, 6.74598520d+00, 6.09242672d+00, 5.14531345d+00, NO+
     2  3.85771979d+00, 2.14301591d+00,-1.29234108d-01,-3.14634898d+00, NO+
     3 -7.17358469d+00,-1.25858544d+01,-1.99151838d+01,-2.99214733d+01, NO+
     4 -4.36979182d+01,-6.28291732d+01,-8.96309180d+01,-1.27514476d+02, NO+
     5 -1.81542480d+02,-2.59279072d+02,-3.72107782d+02,-5.37259472d+02, NO+
     6 -7.81055218d+02,-1.14392927d+03,-1.68848140d+03,-2.51219701d+03, NO+
     7 -3.76813594d+03,-5.69806068d+03,-8.68634962d+03,-1.33482776d+04, NO+
     8 -2.06752716d+04,-3.22751971d+04,-5.07725682d+04,-8.04791169d+04, NO+
     9 -1.28523437d+05,-2.06765594d+05,-3.35060919d+05,-5.46856448d+05/ NO+
      DATA Kp_NO/
     1  8.19987666d+00, 7.98933203d+00, 7.57718872d+00, 6.98651348d+00, NO
     2  6.19208295d+00, 5.14325209d+00, 3.76141900d+00, 1.93281102d+00, NO
     3 -5.03078671d-01,-3.77189083d+00,-8.19292715d+00,-1.42217922d+01, NO
     4 -2.25134201d+01,-3.40161749d+01,-5.01143393d+01,-7.28459673d+01, NO
     5 -1.05236694d+02,-1.51810495d+02,-2.19381621d+02,-3.18263223d+02, NO
     6 -4.64227592d+02,-6.81498301d+02,-1.00756894d+03,-1.50082457d+03, NO
     7 -2.25290519d+03,-3.40855290d+03,-5.19790767d+03,-7.98935787d+03, NO
     8 -1.23765276d+04,-1.93221358d+04,-3.03975961d+04,-4.81845778d+04, NO
     9 -7.69513556d+04,-1.23799142d+05,-2.00616254d+05,-3.27428988d+05/ NO
      DATA Kp_O2/
     1  9.35740706d+00, 9.20950522d+00, 8.91284889d+00, 8.47497479d+00, O2
     2  7.87082228d+00, 7.05867599d+00, 5.97658028d+00, 4.53506097d+00, O2
     3  2.60760444d+00, 1.65226041d-02,-3.48947146d+00,-8.26946054d+00, O2
     4 -1.48406026d+01,-2.39521505d+01,-3.66967621d+01,-5.46807336d+01, O2
     5 -8.02859398d+01,-1.17071862d+02,-1.70399280d+02,-2.48375584d+02, O2
     6 -3.63410399d+02,-5.34570650d+02,-7.91390519d+02,-1.17983324d+03, O2
     7 -1.77207232d+03,-2.68209700d+03,-4.09110286d+03,-6.28915509d+03, O2
     8 -9.74369333d+03,-1.52127601d+04,-2.39336775d+04,-3.79392823d+04, O2
     9 -6.05904539d+04,-9.74786931d+04,-1.57964851d+05,-2.57817625d+05/ O2
      DATA Kp_HF/
     1  8.49902959d+00, 8.32039582d+00, 7.96577127d+00, 7.44758098d+00, HF
     2  6.73711474d+00, 5.78520588d+00, 4.52021243d+00, 2.84048424d+00, HF
     3  6.02607816d-01,-2.39574208d+00,-6.44110098d+00,-1.19420292d+01, HF
     4 -1.94862861d+01,-2.99275380d+01,-4.45159567d+01,-6.50945687d+01, HF
     5 -9.43982473d+01,-1.36514174d+02,-1.97600787d+02,-2.86975533d+02, HF
     6 -4.18892456d+02,-6.15232478d+02,-9.09871438d+02,-1.35555649d+03, HF
     7 -2.03509241d+03,-3.07924834d+03,-4.69588664d+03,-7.21777195d+03, HF
     8 -1.11811780d+04,-1.74557791d+04,-2.74610766d+04,-4.35293076d+04, HF
     9 -6.95163591d+04,-1.11837225d+05,-1.81231395d+05,-2.95789932d+05/ HF
      DATA Kp_NaH/
     1  9.28269416d+00, 9.13911652d+00, 8.88293944d+00, 8.57774326d+00, NaH
     2  8.26031674d+00, 7.92113339d+00, 7.51071153d+00, 6.97308584d+00, NaH
     3  6.25020242d+00, 5.27215323d+00, 3.94602202d+00, 2.14015003d+00, NaH
     4 -3.35917727d-01,-3.75834847d+00,-8.52915012d+00,-1.52376810d+01, NaH
     5 -2.47563393d+01,-3.83898374d+01,-5.81099839d+01,-8.68934368d+01, NaH
     6 -1.29310892d+02,-1.92368510d+02,-2.86922326d+02,-4.29849662d+02, NaH
     7 -6.47683479d+02,-9.82326538d+02,-1.50037224d+03,-2.30842811d+03, NaH
     8 -3.57830713d+03,-5.58856849d+03,-8.79385893d+03,-1.39412745d+04, NaH
     9 -2.22660048d+04,-3.58229712d+04,-5.80522205d+04,-9.47486234d+04/ NaH
      DATA Kp_MgH/
     1  8.77216209d+00, 8.69584363d+00, 8.56871733d+00, 8.41672900d+00, MgH
     2  8.23596180d+00, 8.00909023d+00, 7.71134975d+00, 7.31157062d+00, MgH
     3  6.77202732d+00, 6.04496481d+00, 5.06438249d+00, 3.73555530d+00, MgH
     4  1.92212703d+00,-5.72458239d-01,-4.03325105d+00,-8.87793670d+00, MgH
     5 -1.57264664d+01,-2.55086209d+01,-3.96317957d+01,-6.02152438d+01, MgH
     6 -9.05200189d+01,-1.35539333d+02,-2.03014056d+02,-3.04967349d+02, MgH
     7 -4.60316577d+02,-6.98936059d+02,-1.06829318d+03,-1.64437416d+03, MgH
     8 -2.54964664d+03,-3.98264750d+03,-6.26739839d+03,-9.93642960d+03, MgH
     9 -1.58701825d+04,-2.55333281d+04,-4.13777587d+04,-6.75337649d+04/ MgH
      DATA Kp_MgO/
     1  8.23988246d+00, 8.09313842d+00, 7.82924281d+00, 7.48398601d+00, MgO
     2  7.04854679d+00, 6.49391213d+00, 5.77576501d+00, 4.83180697d+00, MgO
     3  3.57551302d+00, 1.88495743d+00,-4.16478828d-01,-3.58954241d+00, MgO
     4 -8.01729211d+00,-1.42437302d+01,-2.30258860d+01,-3.54550572d+01, MgO
     5 -5.31610286d+01,-7.85936515d+01,-1.15448283d+02,-1.69315956d+02, MgO
     6 -2.48762106d+02,-3.66943268d+02,-5.44239728d+02,-8.12359940d+02, MgO
     7 -1.22110988d+03,-1.84915031d+03,-2.82151085d+03,-4.33834487d+03, MgO
     8 -6.72221305d+03,-1.04961943d+04,-1.65140745d+04,-2.61786306d+04, MgO
     9 -4.18090625d+04,-6.72637635d+04,-1.09002156d+05,-1.77905374d+05/ MgO
      DATA Kp_AlH/
     1  8.80825270d+00, 8.71254780d+00, 8.53521954d+00, 8.29115018d+00, AlH
     2  7.96246016d+00, 7.51281344d+00, 6.89199432d+00, 6.03781322d+00, AlH
     3  4.87365349d+00, 3.29755022d+00, 1.16336604d+00,-1.74220597d+00, AlH
     4 -5.72861568d+00,-1.12440542d+01,-1.89414314d+01,-2.97806274d+01, AlH
     5 -4.51865663d+01,-6.72909467d+01,-9.93062440d+01,-1.46082172d+02, AlH
     6 -2.15040123d+02,-3.17569565d+02,-4.71334618d+02,-7.03808700d+02, AlH
     7 -1.05819066d+03,-1.60270076d+03,-2.44571478d+03,-3.76074481d+03, AlH
     8 -5.82742162d+03,-9.09914336d+03,-1.43159596d+04,-2.26938908d+04, AlH
     9 -3.62433852d+04,-5.83091148d+04,-9.44903877d+04,-1.54219437d+05/ AlH
      DATA Kp_AlO/
     1  8.34884675d+00, 8.18041596d+00, 7.86045687d+00, 7.41251565d+00, AlO
     2  6.81321704d+00, 6.01554501d+00, 4.95105831d+00, 3.52380438d+00, AlO
     3  1.59842421d+00,-1.01507134d+00,-4.58085328d+00,-9.46837422d+00, AlO
     4 -1.62056949d+01,-2.55616284d+01,-3.86635048d+01,-5.71707068d+01, AlO
     5 -8.35408310d+01,-1.21442193d+02,-1.76397612d+02,-2.56760923d+02, AlO
     6 -3.75308437d+02,-5.51667121d+02,-8.16255742d+02,-1.21641283d+03, AlO
     7 -1.82650737d+03,-2.76399055d+03,-4.21550734d+03,-6.47987314d+03, AlO
     8 -1.00386452d+04,-1.56727271d+04,-2.46567744d+04,-3.90849978d+04, AlO
     9 -6.24197384d+04,-1.00421252d+05,-1.62732968d+05,-2.65599605d+05/ AlO
      DATA Kp_AlF/
     1  7.86789935d+00, 7.65555520d+00, 7.24296900d+00, 6.65107872d+00, AlF
     2  5.84578424d+00, 4.76658601d+00, 3.32818632d+00, 1.41215993d+00, AlF
     3 -1.14844775d+00,-4.59067241d+00,-9.25225030d+00,-1.56167635d+01, AlF
     4 -2.43808552d+01,-3.65547173d+01,-5.36132002d+01,-7.77247617d+01, AlF
     5 -1.12101567d+02,-1.61539713d+02,-2.33261941d+02,-3.38200891d+02, AlF
     6 -4.93070987d+02,-7.23540153d+02,-1.06936789d+03,-1.59245685d+03, AlF
     7 -2.39002066d+03,-3.61559547d+03,-5.51321830d+03,-8.47356602d+03, AlF
     8 -1.31262057d+04,-2.04920952d+04,-3.22377377d+04,-5.11010810d+04, AlF
     9 -8.16087588d+04,-1.31291743d+05,-2.12757914d+05,-3.47245710d+05/ AlF
      DATA Kp_Al2/
     1  9.60358240d+00, 9.51547913d+00, 9.37163200d+00, 9.20491702d+00, Al2
     2  9.01278888d+00, 8.77600266d+00, 8.46662911d+00, 8.04822244d+00, Al2
     3  7.47274707d+00, 6.67934463d+00, 5.59136407d+00, 4.10373057d+00, Al2
     4  2.06186323d+00,-7.63041060d-01,-4.70711778d+00,-1.02648068d+01, Al2
     5 -1.81672442d+01,-2.95017283d+01,-4.58963039d+01,-6.97997798d+01, Al2
     6 -1.04944693d+02,-1.57071150d+02,-2.35128111d+02,-3.52997733d+02, Al2
     7 -5.32600577d+02,-8.08538399d+02,-1.23565729d+03,-1.90182855d+03, Al2
     8 -2.94873940d+03,-4.60604961d+03,-7.24860448d+03,-1.14923870d+04, Al2
     9 -1.83557954d+04,-2.95330222d+04,-4.78603106d+04,-7.81154158d+04/ Al2
      DATA Kp_SiH/
     1  8.81868077d+00, 8.71779549d+00, 8.52381193d+00, 8.24648060d+00, SiH
     2  7.86808518d+00, 7.35902203d+00, 6.67977656d+00, 5.77765707d+00, SiH
     3  4.57901090d+00, 2.97824169d+00, 8.23658783d-01,-2.10200965d+00, SiH
     4 -6.11008139d+00,-1.16492062d+01,-1.93724181d+01,-3.02407387d+01, SiH
     5 -4.56815419d+01,-6.78284228d+01,-9.98931181d+01,-1.46722168d+02, SiH
     6 -2.15735250d+02,-3.18322815d+02,-4.72140771d+02,-7.04645322d+02, SiH
     7 -1.05904061d+03,-1.60355998d+03,-2.44659003d+03,-3.76163951d+03, SiH
     8 -5.82832288d+03,-9.10004716d+03,-1.43168622d+04,-2.26947918d+04, SiH
     9 -3.62442883d+04,-5.83100203d+04,-9.44912938d+04,-1.54220339d+05/ SiH
      DATA Kp_SiHm/
     1  8.94881568d+00, 8.83727119d+00, 8.62431238d+00, 8.32485539d+00, SiH-
     2  7.92612962d+00, 7.40363154d+00, 6.72154754d+00, 5.82901242d+00, SiH-
     3  4.65373215d+00, 3.09265440d+00, 9.98662012d-01,-1.83827787d+00, SiH-
     4 -5.71857827d+00,-1.10748547d+01,-1.85369595d+01,-2.90323995d+01, SiH-
     5 -4.39390261d+01,-6.53157553d+01,-9.62612479d+01,-1.41450413d+02, SiH-
     6 -2.08041015d+02,-3.07020907d+02,-4.55422962d+02,-6.79732519d+02, SiH-
     7 -1.02163060d+03,-1.54694476d+03,-2.36023592d+03,-3.62888817d+03, SiH-
     8 -5.62264723d+03,-8.77891749d+03,-1.38116257d+04,-2.18938798d+04, SiH-
     9 -3.49651790d+04,-5.62521449d+04,-9.11564541d+04,-1.48777424d+05/ SiH-
      DATA Kp_SiC/
     1  8.64856138d+00, 8.49159706d+00, 8.18957475d+00, 7.76154326d+00, SiC
     2  7.18777813d+00, 6.43067348d+00, 5.43417837d+00, 4.11765356d+00, SiC
     3  2.36618370d+00, 1.69842334d-02,-3.16023519d+00,-7.49306310d+00, SiC
     4 -1.34520682d+01,-2.17194241d+01,-3.32910148d+01,-4.96302303d+01, SiC
     5 -7.29034614d+01,-1.06343868d+02,-1.54818771d+02,-2.25687731d+02, SiC
     6 -3.30206886d+02,-4.85660750d+02,-7.18820740d+02,-1.07133893d+03, SiC
     7 -1.60868175d+03,-2.43424184d+03,-3.71234504d+03,-5.70607576d+03, SiC
     8 -8.83945568d+03,-1.38000534d+04,-2.17101423d+04,-3.44135876d+04, SiC
     9 -5.49588363d+04,-8.84175375d+04,-1.43280291d+05,-2.33849829d+05/ SiC
      DATA Kp_SiN/
     1  8.83482244d+00, 8.69155423d+00, 8.41349361d+00, 8.01506516d+00, SiN
     2  7.47498317d+00, 6.75424466d+00, 5.79443539d+00, 4.51216098d+00, SiN
     3  2.79176380d+00, 4.73824749d-01,-2.66523158d+00,-6.94566224d+00, SiN
     4 -1.28309553d+01,-2.09941660d+01,-3.24170642d+01,-4.85417580d+01, SiN
     5 -7.15039369d+01,-1.04492817d+02,-1.52311019d+02,-2.22219323d+02, SiN
     6 -3.25323794d+02,-4.78680275d+02,-7.08712298d+02,-1.05653796d+03, SiN
     7 -1.58679502d+03,-2.40156987d+03,-3.66307768d+03,-5.63101082d+03, SiN
     8 -8.72387599d+03,-1.36203374d+04,-2.14281558d+04,-3.39673558d+04, SiN
     9 -5.42469650d+04,-8.72730452d+04,-1.41426386d+05,-2.30824741d+05/ SiN
      DATA Kp_SiO/
     1  7.85947717d+00, 7.62283230d+00, 7.14578721d+00, 6.43708774d+00, SiO
     2  5.45642499d+00, 4.14063373d+00, 2.39555973d+00, 8.24691045d-02, SiO
     3 -2.99845218d+00,-7.13195853d+00,-1.27240955d+01,-2.03562473d+01, SiO
     4 -3.08651003d+01,-4.54624316d+01,-6.59162137d+01,-9.48258406d+01, SiO
     5 -1.36043425d+02,-1.95322068d+02,-2.81321122d+02,-4.07144246d+02, SiO
     6 -5.92819360d+02,-8.69120755d+02,-1.28370885d+03,-1.91078434d+03, SiO
     7 -2.86689178d+03,-4.33610374d+03,-6.61099237d+03,-1.01599211d+04, SiO
     8 -1.57376295d+04,-2.45680846d+04,-3.86491704d+04,-6.12632377d+04, SiO
     9 -9.78369754d+04,-1.57398797d+05,-2.55063516d+05,-4.16292600d+05/ SiO
      DATA Kp_SiF/
     1  8.24896315d+00, 8.07692187d+00, 7.73925153d+00, 7.24947427d+00, SiF
     2  6.57978048d+00, 5.68362935d+00, 4.49455327d+00, 2.91831757d+00, SiF
     3  8.20652971d-01,-1.99016452d+00,-5.78837480d+00,-1.09671679d+01, SiF
     4 -1.80923888d+01,-2.79831332d+01,-4.18336885d+01,-6.13984504d+01, SiF
     5 -8.92745370d+01,-1.29338167d+02,-1.87426811d+02,-2.72371978d+02, SiF
     6 -3.97685162d+02,-5.84114066d+02,-8.63794236d+02,-1.28674408d+03, SiF
     7 -1.93156409d+03,-2.92239512d+03,-4.45652250d+03,-6.84977550d+03, SiF
     8 -1.06111129d+04,-1.65658967d+04,-2.60613483d+04,-4.13108851d+04, SiF
     9 -6.59739357d+04,-1.06138672d+05,-1.71997503d+05,-2.80719925d+05/ SiF
      DATA Kp_Si2/
     1  9.26870957d+00, 9.15601358d+00, 8.94296800d+00, 8.64282820d+00, Si2
     2  8.23807314d+00, 7.70016742d+00, 6.99131874d+00, 6.05956337d+00, Si2
     3  4.83024120d+00, 3.19467706d+00, 9.94994214d-01,-1.99702242d+00, Si2
     4 -6.11072042d+00,-1.18237243d+01,-1.98344524d+01,-3.11705855d+01, Si2
     5 -4.73487482d+01,-7.06151654d+01,-1.04325317d+02,-1.53551505d+02, Si2
     6 -2.26045101d+02,-3.33740255d+02,-4.95147797d+02,-7.39037359d+02, Si2
     7 -1.11075140d+03,-1.68189552d+03,-2.56613620d+03,-3.94547782d+03, Si2
     8 -6.11326813d+03,-9.54515169d+03,-1.50175001d+04,-2.38059299d+04, Si2
     9 -3.80194119d+04,-6.11665733d+04,-9.91212908d+04,-1.61778327d+05/ Si2
      DATA Kp_SH/
     1  8.76784894d+00, 8.65687618d+00, 8.43724426d+00, 8.11640597d+00, SH
     2  7.67551156d+00, 7.08352632d+00, 6.29694486d+00, 5.25481686d+00, SH
     3  3.87035183d+00, 2.01975220d+00,-4.72868776d-01,-3.85744456d+00, SH
     4 -8.49130711d+00,-1.48903994d+01,-2.38079088d+01,-3.63542896d+01, SH
     5 -5.41796490d+01,-7.97528030d+01,-1.16799540d+02,-1.70950413d+02, SH
     6 -2.50835543d+02,-3.69688815d+02,-5.48002763d+02,-8.17670721d+02, SH
     7 -1.22878323d+03,-1.86045319d+03,-2.83844047d+03,-4.36402992d+03, SH
     8 -6.76157555d+03,-1.05571172d+04,-1.66092353d+04,-2.63286588d+04, SH
     9 -4.20477497d+04,-6.76467681d+04,-1.09621675d+05,-1.78915147d+05/ SH
      DATA Kp_SHm/
     1  8.95299485d+00, 8.83408460d+00, 8.59928259d+00, 8.25724736d+00, SH-
     2  7.78830753d+00, 7.15946081d+00, 6.32413770d+00, 5.21712395d+00, SH-
     3  3.74576057d+00, 1.77802037d+00,-8.73676491d-01,-4.47600559d+00, SH-
     4 -9.41043515d+00,-1.62279498d+01,-2.57329765d+01,-3.91114354d+01, SH-
     5 -5.81253824d+01,-8.54108970d+01,-1.24945465d+02,-1.82741083d+02, SH-
     6 -2.68010034d+02,-3.94880914d+02,-5.85230886d+02,-8.73110854d+02, SH-
     7 -1.31199634d+03,-1.98634689d+03,-3.03042418d+03,-4.65912212d+03, SH-
     8 -7.21872271d+03,-1.12708291d+04,-1.77320734d+04,-2.81085536d+04, SH-
     9 -4.48903050d+04,-7.22199153d+04,-1.17032520d+05,-1.91010605d+05/ SH-
      DATA Kp_CS/
     1  8.24711670d+00, 8.04063845d+00, 7.62067128d+00, 6.99179572d+00, CS
     2  6.11706646d+00, 4.94061338d+00, 3.38023276d+00, 1.31399764d+00, CS
     3 -1.43577755d+00,-5.12375350d+00,-1.01126803d+01,-1.69203394d+01, CS
     4 -2.62905757d+01,-3.93002305d+01,-5.75200485d+01,-8.32591590d+01, CS
     5 -1.19940955d+02,-1.72683878d+02,-2.49204611d+02,-3.61183957d+02, CS
     6 -5.26485960d+02,-7.72532919d+02,-1.14176906d+03,-1.70028523d+03, CS
     7 -2.55181370d+03,-3.86019187d+03,-5.88593861d+03,-9.04609274d+03, CS
     8 -1.40127227d+04,-2.18757180d+04,-3.44140676d+04,-5.45504987d+04, CS
     9 -8.71171353d+04,-1.40153198d+05,-2.27117435d+05,-3.70681607d+05/ CS
      DATA Kp_NS/
     1  8.71311806d+00, 8.55617970d+00, 8.24729800d+00, 7.80241690d+00, NS
     2  7.20311088d+00, 6.41335546d+00, 5.37671823d+00, 4.01028418d+00, NS
     3  2.19582151d+00,-2.33864991d-01,-3.51583952d+00,-7.98818898d+00, NS
     4 -1.41364149d+01,-2.26629132d+01,-3.45915237d+01,-5.14266806d+01, NS
     5 -7.53991432d+01,-1.09847130d+02,-1.59809992d+02,-2.32912196d+02, NS
     6 -3.40816911d+02,-5.01412767d+02,-7.42394977d+02,-1.10688948d+03, NS
     7 -1.66260880d+03,-2.51649969d+03,-3.83858799d+03,-5.90104312d+03, NS
     8 -9.14245675d+03,-1.42741020d+04,-2.24569590d+04,-3.55984690d+04, NS
     9 -5.68521877d+04,-9.14646471d+04,-1.48219214d+05,-2.41911794d+05/ NS
      DATA Kp_SO/
     1  8.73293711d+00, 8.57112374d+00, 8.24986973d+00, 7.78138232d+00, SO
     2  7.14160719d+00, 6.28770792d+00, 5.15496617d+00, 3.64964391d+00, SO
     3  1.63876173d+00,-1.06492472d+00,-4.72583748d+00,-9.72065919d+00, SO
     4 -1.65911782d+01,-2.61226385d+01,-3.94603452d+01,-5.82864667d+01, SO
     5 -8.50925268d+01,-1.23602061d+02,-1.79431201d+02,-2.61077640d+02, SO
     6 -3.81554020d+02,-5.60836202d+02,-8.29857108d+02,-1.23676884d+03, SO
     7 -1.85717315d+03,-2.81047596d+03,-4.28649912d+03,-6.58910865d+03, SO
     8 -1.02079748d+04,-1.59372004d+04,-2.50729755d+04,-3.97448723d+04, SO
     9 -6.34736960d+04,-1.02116983d+05,-1.65481019d+05,-2.70084854d+05/ SO
      DATA Kp_MgS/
     1  8.80391364d+00, 8.65367472d+00, 8.38223921d+00, 8.02392447d+00, MgS
     2  7.56677292d+00, 6.97811287d+00, 6.21030965d+00, 5.19846069d+00, MgS
     3  3.85382258d+00, 2.05279622d+00,-3.79625326d-01,-3.69338628d+00, MgS
     4 -8.24785925d+00,-1.45635573d+01,-2.33998232d+01,-3.58723002d+01, MgS
     5 -5.36326209d+01,-7.91460677d+01,-1.16127078d+02,-1.70195760d+02, MgS
     6 -2.49960005d+02,-3.68633638d+02,-5.46677943d+02,-8.15936473d+02, MgS
     7 -1.22642276d+03,-1.85712690d+03,-2.83361400d+03,-4.35688796d+03, MgS
     8 -6.75087517d+03,-1.05408772d+04,-1.65843050d+04,-2.62898851d+04, MgS
     9 -4.19866635d+04,-6.75494429d+04,-1.09465123d+05,-1.78661125d+05/ MgS
      DATA Kp_AlS/
     1  9.13485311d+00, 8.99948162d+00, 8.74453958d+00, 8.39035561d+00, AlS
     2  7.91877184d+00, 7.29405557d+00, 6.46689264d+00, 5.37073564d+00, AlS
     3  3.91240041d+00, 1.95901054d+00,-6.79222565d-01,-4.27368941d+00, AlS
     4 -9.21453412d+00,-1.60666793d+01,-2.56539771d+01,-3.91863731d+01, AlS
     5 -5.84543242d+01,-8.61292433d+01,-1.26234779d+02,-1.84855619d+02, AlS
     6 -2.71306243d+02,-3.99884642d+02,-5.92754097d+02,-8.84397861d+02, AlS
     7 -1.32901202d+03,-2.01218581d+03,-3.06991152d+03,-4.71991726d+03, AlS
     8 -7.31310085d+03,-1.14184671d+04,-1.79647768d+04,-2.84779873d+04, AlS
     9 -4.54809591d+04,-7.31709326d+04,-1.18574623d+05,-1.93528792d+05/ AlS
      DATA Kp_SiS/
     1  8.51932571d+00, 8.32390073d+00, 7.93557206d+00, 7.36661804d+00, SiS
     2  6.58656416d+00, 5.54577146d+00, 4.17089993d+00, 2.35465850d+00, SiS
     3 -5.77213547d-02,-3.28778817d+00,-7.65223809d+00,-1.36044587d+01, SiS
     4 -2.17962974d+01,-3.31712166d+01,-4.91051569d+01,-7.16199063d+01, SiS
     5 -1.03708976d+02,-1.49840453d+02,-2.16742507d+02,-3.14596861d+02, SiS
     6 -4.58977728d+02,-6.73801973d+02,-9.96109137d+02,-1.48355901d+03, SiS
     7 -2.22674035d+03,-3.36872973d+03,-5.13692275d+03,-7.89534667d+03, SiS
     8 -1.22306208d+04,-1.90940589d+04,-3.00384778d+04,-4.76150510d+04, SiS
     9 -7.60416277d+04,-1.22335469d+05,-1.98244426d+05,-3.23558193d+05/ SiS
      DATA Kp_S2/
     1  9.16241910d+00, 9.03806004d+00, 8.78774028d+00, 8.41521233d+00, S2
     2  7.89658733d+00, 7.19605148d+00, 6.26352081d+00, 5.02643253d+00, S2
     3  3.37853935d+00, 1.16648379d+00,-1.82781383d+00,-5.91363754d+00, S2
     4 -1.15333931d+01,-1.93270983d+01,-3.02284506d+01,-4.56086762d+01, S2
     5 -6.74972408d+01,-9.89258739d+01,-1.44471826d+02,-2.11062548d+02, S2
     6 -3.09319109d+02,-4.55529867d+02,-6.74912481d+02,-1.00672222d+03, S2
     7 -1.51259776d+03,-2.28989173d+03,-3.49337070d+03,-5.37078015d+03, S2
     8 -8.32135921d+03,-1.29925555d+04,-2.04411635d+04,-3.24034631d+04, S2
     9 -5.17500432d+04,-8.32566919d+04,-1.34918694d+05,-2.20204400d+05/ S2
      DATA Kp_HCl/
     1  8.81949164d+00, 8.68341013d+00, 8.41313207d+00, 8.01779831d+00, HCl
     2  7.47519357d+00, 6.74769327d+00, 5.78078987d+00, 4.49727023d+00, HCl
     3  2.78815826d+00, 4.99661244d-01,-2.58587080d+00,-6.77808020d+00, HCl
     4 -1.25208146d+01,-2.04565422d+01,-3.15253840d+01,-4.71161904d+01, HCl
     5 -6.92952817d+01,-1.01153838d+02,-1.47347832d+02,-2.14915175d+02, HCl
     6 -3.14621964d+02,-4.62995020d+02,-6.85624283d+02,-1.02234645d+03, HCl
     7 -1.53571642d+03,-2.32453323d+03,-3.54585395d+03,-5.45106195d+03, HCl
     8 -8.44523971d+03,-1.31853666d+04,-2.07437520d+04,-3.28822400d+04, HCl
     9 -5.25136916d+04,-8.44841707d+04,-1.36906537d+05,-2.23447195d+05/ HCl
      DATA Kp_LiCl/
     1  8.18381416d+00, 7.99226373d+00, 7.63724533d+00, 7.16137812d+00, LiCl
     2  6.55380293d+00, 5.77158623d+00, 4.74713277d+00, 3.39014044d+00, LiCl
     3  1.57919183d+00,-8.54127036d-01,-4.14747224d+00,-8.63990519d+00, LiCl
     4 -1.48194181d+01,-2.33936374d+01,-3.53963804d+01,-5.23497060d+01, LiCl
     5 -7.65112895d+01,-1.11253174d+02,-1.61650492d+02,-2.35383630d+02, LiCl
     6 -3.44197765d+02,-5.06133626d+02,-7.49123318d+02,-1.11665204d+03, LiCl
     7 -1.67699863d+03,-2.53800342d+03,-3.87110631d+03,-5.95074430d+03, LiCl
     8 -9.21916561d+03,-1.43935688d+04,-2.26446076d+04,-3.58956176d+04, LiCl
     9 -5.73264341d+04,-9.22273127d+04,-1.49454824d+05,-2.43928189d+05/ LiCl
      DATA Kp_NaCl/
     1  8.57661088d+00, 8.35504692d+00, 7.94955622d+00, 7.43524426d+00, NaCl
     2  6.84060695d+00, 6.13666138d+00, 5.24400076d+00, 4.06609481d+00, NaCl
     3  2.48994969d+00, 3.67262048d-01,-2.50852558d+00,-6.43228630d+00, NaCl
     4 -1.18293780d+01,-1.93173922d+01,-2.97994564d+01,-4.46065111d+01, NaCl
     5 -6.57144578d+01,-9.60753460d+01,-1.40128770d+02,-2.04590495d+02, NaCl
     6 -2.99718164d+02,-4.41275124d+02,-6.53672491d+02,-9.74910607d+02, NaCl
     7 -1.46466495d+03,-2.21718653d+03,-3.38230510d+03,-5.19986674d+03, NaCl
     8 -8.05638451d+03,-1.25786631d+04,-1.97898015d+04,-3.13707244d+04, NaCl
     9 -5.01005033d+04,-8.06026657d+04,-1.30617591d+05,-2.13184231d+05/ NaCl
      DATA Kp_AlCl/
     1  8.34659229d+00, 8.18686070d+00, 7.87953485d+00, 7.44137405d+00, AlCl
     2  6.84488155d+00, 6.04218219d+00, 4.96831568d+00, 3.53530402d+00, AlCl
     3  1.61977468d+00,-9.53786713d-01,-4.43584750d+00,-9.18512272d+00, AlCl
     4 -1.57183202d+01,-2.47848777d+01,-3.74799860d+01,-5.54164834d+01, AlCl
     5 -8.09869292d+01,-1.17763698d+02,-1.71117710d+02,-2.49172502d+02, AlCl
     6 -3.64335552d+02,-5.35674505d+02,-7.92736398d+02,-1.18151039d+03, AlCl
     7 -1.77424651d+03,-2.68505013d+03,-4.09525566d+03,-6.29517245d+03, AlCl
     8 -9.75264825d+03,-1.52263586d+04,-2.39546718d+04,-3.79721798d+04, AlCl
     9 -6.06426674d+04,-9.75624539d+04,-1.58100524d+05,-2.58039277d+05/ AlCl
      DATA Kp_CaH/
     1  8.86559294d+00, 8.73698020d+00, 8.50853162d+00, 8.22991142d+00, CaH
     2  7.92131269d+00, 7.57711722d+00, 7.16909133d+00, 6.65455638d+00, CaH
     3  5.97939263d+00, 5.07715795d+00, 3.86118566d+00, 2.21062801d+00, CaH
     4 -4.75762350d-02,-3.16329919d+00,-7.49958321d+00,-1.35884049d+01, CaH
     5 -2.22178430d+01,-3.45677042d+01,-5.24219622d+01,-7.84715657d+01, CaH
     6 -1.16850664d+02,-1.73894277d+02,-2.59419823d+02,-3.88686205d+02, CaH
     7 -5.85687252d+02,-8.88314856d+02,-1.35678410d+03,-2.08749471d+03, CaH
     8 -3.23581631d+03,-5.05364326d+03,-7.95209310d+03,-1.26067221d+04, CaH
     9 -2.01344494d+04,-3.23934627d+04,-5.24944355d+04,-8.56773940d+04/ CaH
      DATA Kp_CaF/
     1  7.85157654d+00, 7.61063787d+00, 7.15785071d+00, 6.54979264d+00, CaF
     2  5.78970281d+00, 4.84514751d+00, 3.64897465d+00, 2.09906903d+00, CaF
     3  5.09036902d-02,-2.69327775d+00,-6.40608573d+00,-1.14728434d+01, CaF
     4 -1.84473083d+01,-2.81325027d+01,-4.17011538d+01,-6.08779348d+01, CaF
     5 -8.82178541d+01,-1.27535255d+02,-1.84573435d+02,-2.68027738d+02, CaF
     6 -3.91198992d+02,-5.74514687d+02,-8.49599604d+02,-1.26569265d+03, CaF
     7 -1.90009888d+03,-2.87491577d+03,-4.38425448d+03,-6.73884187d+03, CaF
     8 -1.04393991d+04,-1.62979579d+04,-2.56399710d+04,-4.06430792d+04, CaF
     9 -6.49075730d+04,-1.04423268d+05,-1.69217916d+05,-2.76183630d+05/ CaF
      DATA Kp_CaCl/
     1  8.20330207d+00, 8.00633025d+00, 7.64005375d+00, 7.15637515d+00, CaCl
     2  6.56314407d+00, 5.83719463d+00, 4.92688494d+00, 3.75432565d+00, CaCl
     3  2.21005451d+00, 1.45232545d-01,-2.64406342d+00,-6.44504957d+00, CaCl
     4 -1.16699034d+01,-1.89162117d+01,-2.90574730d+01,-4.33807705d+01, CaCl
     5 -6.37967884d+01,-9.31600001d+01,-1.35763117d+02,-1.98099207d+02, CaCl
     6 -2.90086455d+02,-4.26966436d+02,-6.32342672d+02,-9.42956925d+02, CaCl
     7 -1.41651006d+03,-2.14413416d+03,-3.27069947d+03,-5.02811333d+03, CaCl
     8 -7.79009720d+03,-1.21627100d+04,-1.91351867d+04,-3.03328188d+04, CaCl
     9 -4.84426994d+04,-7.79353292d+04,-1.26294897d+05,-2.06128792d+05/ CaCl
      DATA Kp_ScO/
     1  8.41256844d+00, 8.17292754d+00, 7.70582422d+00, 7.04015234d+00, ScO
     2  6.15194890d+00, 4.99299270d+00, 3.48754103d+00, 1.51940324d+00, ScO
     3 -1.08433683d+00,-4.57045294d+00,-9.28565904d+00,-1.57211604d+01, ScO
     4 -2.45818217d+01,-3.68888369d+01,-5.41327466d+01,-7.85046525d+01, ScO
     5 -1.13249571d+02,-1.63211051d+02,-2.35679308d+02,-3.41687940d+02, ScO
     6 -4.98112911d+02,-7.30894546d+02,-1.08021262d+03,-1.60861065d+03, ScO
     7 -2.41427873d+03,-3.65229832d+03,-5.56919576d+03,-8.55961864d+03, ScO
     8 -1.32595218d+04,-2.07002397d+04,-3.25652091d+04,-5.16201895d+04, ScO
     9 -8.24377991d+04,-1.32625527d+05,-2.14919353d+05,-3.50773506d+05/ ScO
      DATA Kp_TiO/
     1  8.32708851d+00, 8.09510734d+00, 7.64160031d+00, 6.99288535d+00, TiO
     2  6.12353618d+00, 4.98313825d+00, 3.49388361d+00, 1.54129278d+00, TiO
     3 -1.04049666d+00,-4.48985043d+00,-9.14796863d+00,-1.55025038d+01, TiO
     4 -2.42519822d+01,-3.64046663d+01,-5.34300227d+01,-7.74877392d+01, TiO
     5 -1.11776236d+02,-1.61070438d+02,-2.32563270d+02,-3.37146405d+02, TiO
     6 -4.91496646d+02,-7.21232565d+02,-1.06601245d+03,-1.58758023d+03, TiO
     7 -2.38284043d+03,-3.60484911d+03,-5.49696100d+03,-8.44872292d+03, TiO
     8 -1.30878535d+04,-2.04323586d+04,-3.21439099d+04,-5.09524997d+04, TiO
     9 -8.13716139d+04,-1.30910373d+05,-2.12140064d+05,-3.46237489d+05/ TiO
      DATA Kp_TiS/
     1  9.21674922d+00, 9.03554376d+00, 8.68494970d+00, 8.19029219d+00, TiS
     2  7.53768866d+00, 6.69602774d+00, 5.61532870d+00, 4.21913025d+00, TiS
     3  2.39347616d+00,-2.73343496d-02,-3.28105188d+00,-7.70723224d+00, TiS
     4 -1.37913459d+01,-2.22327502d+01,-3.40491559d+01,-5.07345493d+01, TiS
     5 -7.44993843d+01,-1.08641582d+02,-1.58130728d+02,-2.30492017d+02, TiS
     6 -3.37264379d+02,-4.96165781d+02,-7.34619771d+02,-1.09531324d+03, TiS
     7 -1.64524556d+03,-2.49023836d+03,-3.79855070d+03,-5.83951612d+03, TiS
     8 -9.04715114d+03,-1.41253182d+04,-2.22228969d+04,-3.52274418d+04, TiS
     9 -5.62596472d+04,-9.05114184d+04,-1.46674683d+05,-2.39391310d+05/ TiS
      DATA Kp_VO/
     1  9.00539863d+00, 8.78583664d+00, 8.35829777d+00, 7.74949427d+00, VO
     2  6.93625182d+00, 5.87007193d+00, 4.47452708d+00, 2.63758516d+00, VO
     3  2.00202143d-01,-3.06118355d+00,-7.46217884d+00,-1.34523409d+01, VO
     4 -2.16782826d+01,-3.30791180d+01,-4.90289011d+01,-7.15479003d+01, VO
     5 -1.03623476d+02,-1.49709880d+02,-2.16516035d+02,-3.14195122d+02, VO
     6 -4.58299430d+02,-6.72719006d+02,-9.94461364d+02,-1.48112359d+03, VO
     7 -2.22314331d+03,-3.36335519d+03,-5.12879734d+03,-7.88292929d+03, VO
     8 -1.22114582d+04,-1.90642174d+04,-2.99916080d+04,-4.75408381d+04, VO
     9 -7.59231931d+04,-1.22144992d+05,-1.97935764d+05,-3.23054338d+05/ VO
      DATA Kp_CrO/
     1  9.10622598d+00, 8.91783392d+00, 8.56093864d+00, 8.07217459d+00, CrO
     2  7.44613363d+00, 6.65524269d+00, 5.64871909d+00, 4.35099838d+00, CrO
     3  2.65597834d+00, 4.11598404d-01,-2.60277958d+00,-6.70314759d+00, CrO
     4 -1.23396379d+01,-2.01598707d+01,-3.11068044d+01,-4.65654951d+01, CrO
     5 -6.85862642d+01,-1.00229205d+02,-1.46103469d+02,-2.13183263d+02, CrO
     6 -3.12145801d+02,-4.59390762d+02,-6.80317689d+02,-1.01445637d+03, CrO
     7 -1.52388323d+03,-2.30664041d+03,-3.51857675d+03,-5.40917827d+03, CrO
     8 -8.38049436d+03,-1.30845209d+04,-2.05854800d+04,-3.26318598d+04, CrO
     9 -5.21144257d+04,-8.38425104d+04,-1.35867554d+05,-2.21752503d+05/ CrO
      DATA Kp_FeO/
     1  9.44933058d+00, 9.28596728d+00, 8.97412755d+00, 8.54113513d+00, FeO
     2  7.97649754d+00, 7.25052753d+00, 6.31345503d+00, 5.09136007d+00, FeO
     3  3.47906356d+00, 1.32777593d+00,-1.57377476d+00,-5.52574917d+00, FeO
     4 -1.09560633d+01,-1.84806874d+01,-2.89958096d+01,-4.38169654d+01, FeO
     5 -6.48929249d+01,-9.51384650d+01,-1.38955668d+02,-2.03000789d+02, FeO
     6 -2.97477631d+02,-4.38041145d+02,-6.48937579d+02,-9.67899929d+02, FeO
     7 -1.45418300d+03,-2.20137255d+03,-3.35823355d+03,-5.16291127d+03, FeO
     8 -7.99918091d+03,-1.24894025d+04,-1.96494219d+04,-3.11482521d+04, FeO
     9 -4.97452631d+04,-8.00311830d+04,-1.29691475d+05,-2.11672583d+05/ FeO
      DATA Kp_YO/
     1  8.20797845d+00, 7.96388599d+00, 7.48260592d+00, 6.78779527d+00, YO
     2  5.85172088d+00, 4.62511146d+00, 3.03296513d+00, 9.58235971d-01, YO
     3 -1.77841089d+00,-5.43686093d+00,-1.03825368d+01,-1.71313398d+01, YO
     4 -2.64218168d+01,-3.93226167d+01,-5.73920528d+01,-8.29191951d+01, YO
     5 -1.19294260d+02,-1.71583719d+02,-2.47429231d+02,-3.58398390d+02, YO
     6 -5.22191280d+02,-7.65981242d+02,-1.13184114d+03,-1.68528277d+03, YO
     7 -2.52914240d+03,-3.82584350d+03,-5.83361208d+03,-8.96580941d+03, YO
     8 -1.38885330d+04,-2.16820206d+04,-3.41095303d+04,-5.40679471d+04, YO
     9 -8.63466841d+04,-1.38913940d+05,-2.25109577d+05,-3.67405085d+05/ YO
      DATA Kp_ZrO/
     1  7.83108622d+00, 7.58427527d+00, 7.09718944d+00, 6.39139484d+00, ZrO
     2  5.43239598d+00, 4.15824591d+00, 2.47755673d+00, 2.60063884d-01, ZrO
     3 -2.67941196d+00,-6.60517909d+00,-1.18970484d+01,-1.91049544d+01, ZrO
     4 -2.90284281d+01,-4.28314420d+01,-6.22106418d+01,-8.96435247d+01, ZrO
     5 -1.28776687d+02,-1.85055210d+02,-2.66702561d+02,-3.86175250d+02, ZrO
     6 -5.62530282d+02,-8.25028053d+02,-1.21897211d+03,-1.81490851d+03, ZrO
     7 -2.72357189d+03,-4.11986130d+03,-6.28184034d+03,-9.65462474d+03, ZrO
     8 -1.49554776d+04,-2.33476179d+04,-3.67297580d+04,-5.82213117d+04, ZrO
     9 -9.29795978d+04,-1.49584911d+05,-2.42401852d+05,-3.95628124d+05/ ZrO
      DATA Kp_LaO/
     1  7.84593402d+00, 7.58571517d+00, 7.07215399d+00, 6.32781963d+00, LaO
     2  5.31574181d+00, 3.96906654d+00, 2.18833537d+00,-1.68494357d-01, LaO
     3 -3.30145610d+00,-7.49119010d+00,-1.31341879d+01,-2.07981275d+01, LaO
     4 -3.13059056d+01,-4.58594160d+01,-6.62195735d+01,-9.49772183d+01, LaO
     5 -1.35973093d+02,-1.94943323d+02,-2.80518687d+02,-4.05762523d+02, LaO
     6 -5.90641665d+02,-8.65834036d+02,-1.27883548d+03,-1.90360689d+03, LaO
     7 -2.85624332d+03,-4.32011019d+03,-6.58673111d+03,-1.01227691d+04, LaO
     8 -1.56802071d+04,-2.44785720d+04,-3.85084898d+04,-6.10403734d+04, LaO
     9 -9.74811900d+04,-1.56826588d+05,-2.54136536d+05,-4.14780107d+05/ LaO
C
      DATA TEMP0/-1./,FIRST/.TRUE./
C
C Compute 2nd derivatives for spline interpolation
C
      IF(FIRST) THEN
        DO 1 I=1,MSPEC
          CALL SPL_INIT(TH,Qm(1,I),Qm2(1,I),U,NTH)
          CALL SPL_INIT(TH,Kp(1,I),Kp2(1,I),U,NTH)
  1     CONTINUE
        FIRST=.FALSE.
      ENDIF
C
C Check if we have new temperature. If we do, find the bracketing indices
C
      IF(ABS(TEMP-TEMP0).GT.1.0) THEN
        THETA=5040.d0/TEMP
        TEMP0=TEMP
        KHI=NTH
        KLO=1
  2     CONTINUE
        I=(KLO+KHI)/2
        A=TH(I)
        IF(A.GT.THETA) THEN
          KHI=I
        ELSE IF(A.LE.THETA) THEN
          KLO=I
        END IF
        IF(KHI-KLO.GT.1) GO TO 2
      ENDIF
C
C Find species name
C
      DO 3 I=1,MSPEC
        ISPEC=I
        IF(SPLIST(I).EQ.SPNAME) GO TO 4
c        IF(SPLIST(I).EQ.SPNAME.AND.
c     *     ISPEC.NE.1.AND.ISPEC.NE.16) GO TO 4
   3  CONTINUE
C
C Species was not found
C
      Qm_spln=-1.D0
      Kp_spln=-1.D0
      RETURN
C
C Do the interpolation
C
   4  Qm_spln=SPL_INTERP(KLO,KHI,TH,Qm(1,ISPEC),Qm2(1,ISPEC),NTH,THETA)
      Qm_spln=10.d0**Qm_spln
      Kp_spln=SPL_INTERP(KLO,KHI,TH,Kp(1,ISPEC),Kp2(1,ISPEC),NTH,THETA)
      Kp_spln=10.d0**(Kp_spln+1.D0)
      RETURN
C
      END

      SUBROUTINE SPL_INIT(X,Y,Y2,U,N)
C
C  Computes second derivative approximations for cubic spline interpolation
C
      IMPLICIT NONE
      INTEGER N
      REAL*8 X(N),Y(N),Y2(N),U(N)
      INTEGER I
      REAL*8 SIG,P,YY1,YY2,YY3
C
C  Natural lower boundary condition
C
      Y2(1)=0.D0
      U(1)=0.D0
      DO 1 I=2,N-1
      SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
      P=SIG*Y2(I-1)+2.D0
      Y2(I)=(SIG-1.D0)/P
      YY1=Y(I-1)
      YY2=Y(I  )
      YY3=Y(I+1)
      U(I)=(6.D0*((YY3-YY2)/(X(I+1)-X(I))-(YY2-YY1)/
     /     (X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
  1   CONTINUE
C
C  Natural upper boundary condition
C
      Y2(N)=0.D0
      DO 2 I=N-1,1,-1
  2   Y2(I)=Y2(I)*Y2(I+1)+U(I)
C
      RETURN
      END

      REAL*8 FUNCTION SPL_INTERP(KLO,KHI,XA,YA,Y2A,N,X)
C
C  Performs cubic spline interpolation
C
      IMPLICIT NONE
      INTEGER KLO,KHI,N
      REAL*8 XA(N),YA(N),Y2A(N),X
      REAL*8 A,B,H,Y1,Y2
C
      H=XA(KHI)-XA(KLO)
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y1=YA(KLO)
      Y2=YA(KHI)
      SPL_INTERP=A*Y1+B*Y2+((A*A-1.D0)*A*Y2A(KLO)+
     +                      (B*B-1.D0)*B*Y2A(KHI))*(H*H)/6.D0
C
      RETURN
      END

      SUBROUTINE NEGION(NUM,TEMP,QPRD,IT,PART,POTION)
C
C                                     (3/2)              Eion
C   1     N(A)*N(e)        (2*Pi*m*kT)       2*U(A)    - ----
C  -- =   --------- = kT * -----------    *  ------ * e   kT
C  IT       N(A-)              h^3           U(A-)
C
C  QPRD is 2*U(A)
C                        (3/2)
C  Const = (2*Pi*m*k/h^2)
C
      DOUBLE PRECISION QPRD,IT,PART,Const,kBoleV,TkeV,Tlog
      REAL      POTION,TEMP
      PARAMETER (Const=0.3333984D0,kBoleV=8.6173175D-5)
      PARAMETER (NTABLE=6,NCOEF=11)
      REAL      TABLE(NCOEF,NTABLE)
      SAVE      TABLE
C
      DATA TABLE/
C     Atom  Partition function                           Electr. affin  Species
     *  1, 2.0,8*0.,                                        -0.754190,  H-
     *  6, 4.0,8*0.,                                        -1.262114,  C-
     *  8, 2.0,8*0.,                                        -1.461112,  O-
     *  9, 2.0,8*0.,                                        -3.401200,  F-
     * 14, 2.0,8*0.,                                        -1.389517,  Si-
     * 17, 6.0,8*0.,                                        -3.614400/  Cl-
C
      DO 1 I=1,NTABLE
        IF(NUM.EQ.INT(TABLE(1,I)+0.5)) THEN
          J=I
          GO TO 2
        ENDIF
   1  CONTINUE
      WRITE(*,*) 'NEGION: negative ion for element',NUM,' is unknown'
      STOP
C
   2  TkeV=kBoleV*TEMP
      IT=Const*QPRD*EXP(TABLE(NCOEF,J)/TkeV)*SQRT(TEMP)*TEMP*TEMP
      Tlog=LOG(TEMP)
      PART=TABLE(NCOEF-1,J)
      DO 3 I=NCOEF-2,2,-1
    3 PART=Tlog*PART+TABLE(I,J)
c      write(*,*) '->',QPRD,EXP(TABLE(NCOEF,J)/TkeV),PART
      IT=1.D0/IT
      POTION=-TABLE(NCOEF,J)
C
      RETURN
      END


      SUBROUTINE XSAHA(IEL,TT,XNELEC,XNATOM,MAXION,POTI,FRCT,
     *                 IONSIZ,MODE)
C
C     MODE=1 returns ionization fractions/partition functions
C     MODE=2 returns ionization fractions
C     MODE=3 returns partition functions
C     MODE=4 returns total number of electrons produced
C     MODE=5 returns in MAXION(!) the number of ionization stages
C            available in XSAHA
C
C     ALL OF THE ABOVE IS FOR ALL IONIZATION STAGES UP TO MAXION
C
C  Parameters:
C     IEL    - (input) element atomic number (Hydrogen: 1)
C     TT     - (input) temperature (Kelvins)
C     XNELEC - (input) electron number density (cm^-3)
C     XNATOM - (input) particle number density (excluding electrons) (cm^-3)
C     MAXION - (input/output) size of the output arrays
C     POTI   - (output array of MAXION) ionization potential (eV)
C     FRCT   - (output array of MAXION) results according to MODE
C     MODE   - (input) see above
C
      INTEGER ELESIZ
      PARAMETER (ELESIZ=100)
      DOUBLE PRECISION F(IONSIZ),FEXARG,FRCT(MAXION),CF
      REAL IP(IONSIZ),PART(IONSIZ),POTLO(IONSIZ),SCALE(4),
     *     POTI(MAXION),TT
      INTEGER LOCZ(ELESIZ+1)
      LOGICAL FIRST

      INTEGER SIZ_H ,SIZ_He,SIZ_Li,SIZ_Be,SIZ_B ,SIZ_C ,SIZ_N ,SIZ_O ,
     1        SIZ_F ,SIZ_Ne,SIZ_Na,SIZ_Mg,SIZ_Al,SIZ_Si,SIZ_P ,SIZ_S ,
     2        SIZ_Cl,SIZ_Ar,SIZ_K ,SIZ_Ca,SIZ_Sc,SIZ_Ti,SIZ_V ,SIZ_Cr,
     3        SIZ_Mn,SIZ_Fe,SIZ_Co,SIZ_Ni,SIZ_Cu,SIZ_Zn,SIZ_Ga,SIZ_Ge,
     4        SIZ_As,SIZ_Se,SIZ_Br,SIZ_Kr,SIZ_Rb,SIZ_Sr,SIZ_Y ,SIZ_Zr,
     5        SIZ_Nb,SIZ_Mo,SIZ_Tc,SIZ_Ru,SIZ_Rh,SIZ_Pd,SIZ_Ag,SIZ_Cd,
     6        SIZ_In,SIZ_Sn,SIZ_Sb,SIZ_Te,SIZ_I ,SIZ_Xe,SIZ_Cs,SIZ_Ba,
     7        SIZ_La,SIZ_Ce,SIZ_Pr,SIZ_Nd,SIZ_Pm,SIZ_Sm,SIZ_Eu,SIZ_Gd,
     8        SIZ_Tb,SIZ_Dy,SIZ_Ho,SIZ_Er,SIZ_Tm,SIZ_Yb,SIZ_Lu,SIZ_Hf,
     9        SIZ_Ta,SIZ_W ,SIZ_Re,SIZ_Os,SIZ_Ir,SIZ_Pt,SIZ_Au,SIZ_Hg,
     A        SIZ_Tl,SIZ_Pb,SIZ_Bi,SIZ_Po,SIZ_At,SIZ_Rn,SIZ_Fr,SIZ_Ra,
     B        SIZ_Ac,SIZ_Th,SIZ_Pa,SIZ_U ,SIZ_Np,SIZ_Pu,SIZ_Am,SIZ_Cm,
     C        SIZ_Bk,SIZ_Cf,SIZ_Es
      INTEGER OFF_H ,OFF_He,OFF_Li,OFF_Be,OFF_B ,OFF_C ,OFF_N ,OFF_O ,
     1        OFF_F ,OFF_Ne,OFF_Na,OFF_Mg,OFF_Al,OFF_Si,OFF_P ,OFF_S ,
     2        OFF_Cl,OFF_Ar,OFF_K ,OFF_Ca,OFF_Sc,OFF_Ti,OFF_V ,OFF_Cr,
     3        OFF_Mn,OFF_Fe,OFF_Co,OFF_Ni,OFF_Cu,OFF_Zn,OFF_Ga,OFF_Ge,
     4        OFF_As,OFF_Se,OFF_Br,OFF_Kr,OFF_Rb,OFF_Sr,OFF_Y ,OFF_Zr,
     5        OFF_Nb,OFF_Mo,OFF_Tc,OFF_Ru,OFF_Rh,OFF_Pd,OFF_Ag,OFF_Cd,
     6        OFF_In,OFF_Sn,OFF_Sb,OFF_Te,OFF_I ,OFF_Xe,OFF_Cs,OFF_Ba,
     7        OFF_La,OFF_Ce,OFF_Pr,OFF_Nd,OFF_Pm,OFF_Sm,OFF_Eu,OFF_Gd,
     8        OFF_Tb,OFF_Dy,OFF_Ho,OFF_Er,OFF_Tm,OFF_Yb,OFF_Lu,OFF_Hf,
     9        OFF_Ta,OFF_W ,OFF_Re,OFF_Os,OFF_Ir,OFF_Pt,OFF_Au,OFF_Hg,
     A        OFF_Tl,OFF_Pb,OFF_Bi,OFF_Po,OFF_At,OFF_Rn,OFF_Fr,OFF_Ra,
     B        OFF_Ac,OFF_Th,OFF_Pa,OFF_U ,OFF_Np,OFF_Pu,OFF_Am,OFF_Cm,
     C        OFF_Bk,OFF_Cf,OFF_Es
C
C In order to add data for another ionization stage to a perticular element
C one has to do two things: increase the value of SIZ_<elname> and add the
C data line(s) in the DATA NNN_<elname>
C
      PARAMETER (SIZ_H = 2, OFF_H = 1)
      INTEGER NNN_H (8*SIZ_H )
      PARAMETER (SIZ_He= 3, OFF_He=OFF_H +SIZ_H  )
      INTEGER NNN_He(8*SIZ_He)
      PARAMETER (SIZ_Li= 4, OFF_Li=OFF_He+SIZ_He)
      INTEGER NNN_Li(8*SIZ_Li)
      PARAMETER (SIZ_Be= 4, OFF_Be=OFF_Li+SIZ_Li)
      INTEGER NNN_Be(8*SIZ_Be)
      PARAMETER (SIZ_B = 4, OFF_B =OFF_Be+SIZ_Be)
      INTEGER NNN_B (8*SIZ_B )
      PARAMETER (SIZ_C = 6, OFF_C =OFF_B +SIZ_B )
      INTEGER NNN_C (8*SIZ_C )
      PARAMETER (SIZ_N = 6, OFF_N =OFF_C +SIZ_C )
      INTEGER NNN_N (8*SIZ_N )
      PARAMETER (SIZ_O = 6, OFF_O =OFF_N +SIZ_N )
      INTEGER NNN_O (8*SIZ_O )
      PARAMETER (SIZ_F = 6, OFF_F =OFF_O +SIZ_O )
      INTEGER NNN_F (8*SIZ_F )
      PARAMETER (SIZ_Ne= 6, OFF_Ne=OFF_F +SIZ_F )
      INTEGER NNN_Ne(8*SIZ_Ne)
      PARAMETER (SIZ_Na= 6, OFF_Na=OFF_Ne+SIZ_Ne)
      INTEGER NNN_Na(8*SIZ_Na)
      PARAMETER (SIZ_Mg= 6, OFF_Mg=OFF_Na+SIZ_Na)
      INTEGER NNN_Mg(8*SIZ_Mg)
      PARAMETER (SIZ_Al= 6, OFF_Al=OFF_Mg+SIZ_Mg)
      INTEGER NNN_Al(8*SIZ_Al)
      PARAMETER (SIZ_Si= 6, OFF_Si=OFF_Al+SIZ_Al)
      INTEGER NNN_Si(8*SIZ_Si)
      PARAMETER (SIZ_P = 6, OFF_P =OFF_Si+SIZ_Si)
      INTEGER NNN_P (8*SIZ_P )
      PARAMETER (SIZ_S = 6, OFF_S =OFF_P +SIZ_P )
      INTEGER NNN_S (8*SIZ_S )
      PARAMETER (SIZ_Cl= 5, OFF_Cl=OFF_S +SIZ_S )
      INTEGER NNN_Cl(8*SIZ_Cl)
      PARAMETER (SIZ_Ar= 5, OFF_Ar=OFF_Cl+SIZ_Cl)
      INTEGER NNN_Ar(8*SIZ_Ar)
      PARAMETER (SIZ_K = 5, OFF_K =OFF_Ar+SIZ_Ar)
      INTEGER NNN_K (8*SIZ_K )
      PARAMETER (SIZ_Ca= 5, OFF_Ca=OFF_K +SIZ_K )
      INTEGER NNN_Ca(8*SIZ_Ca)
      PARAMETER (SIZ_Sc= 5, OFF_Sc=OFF_Ca+SIZ_Ca)
      INTEGER NNN_Sc(8*SIZ_Sc)
      PARAMETER (SIZ_Ti= 5, OFF_Ti=OFF_Sc+SIZ_Sc)
      INTEGER NNN_Ti(8*SIZ_Ti)
      PARAMETER (SIZ_V = 5, OFF_V =OFF_Ti+SIZ_Ti)
      INTEGER NNN_V (8*SIZ_V )
      PARAMETER (SIZ_Cr= 5, OFF_Cr=OFF_V +SIZ_V )
      INTEGER NNN_Cr(8*SIZ_Cr)
      PARAMETER (SIZ_Mn= 5, OFF_Mn=OFF_Cr+SIZ_Cr)
      INTEGER NNN_Mn(8*SIZ_Mn)
      PARAMETER (SIZ_Fe= 5, OFF_Fe=OFF_Mn+SIZ_Mn)
      INTEGER NNN_Fe(8*SIZ_Fe)
      PARAMETER (SIZ_Co= 5, OFF_Co=OFF_Fe+SIZ_Fe)
      INTEGER NNN_Co(8*SIZ_Co)
      PARAMETER (SIZ_Ni= 5, OFF_Ni=OFF_Co+SIZ_Co)
      INTEGER NNN_Ni(8*SIZ_Ni)
      PARAMETER (SIZ_Cu= 3, OFF_Cu=OFF_Ni+SIZ_Ni)
      INTEGER NNN_Cu(8*SIZ_Cu)
      PARAMETER (SIZ_Zn= 3, OFF_Zn=OFF_Cu+SIZ_Cu)
      INTEGER NNN_Zn(8*SIZ_Zn)
      PARAMETER (SIZ_Ga= 3, OFF_Ga=OFF_Zn+SIZ_Zn)
      INTEGER NNN_Ga(8*SIZ_Ga)
      PARAMETER (SIZ_Ge= 3, OFF_Ge=OFF_Ga+SIZ_Ga)
      INTEGER NNN_Ge(8*SIZ_Ge)
      PARAMETER (SIZ_As= 3, OFF_As=OFF_Ge+SIZ_Ge)
      INTEGER NNN_As(8*SIZ_As)
      PARAMETER (SIZ_Se= 3, OFF_Se=OFF_As+SIZ_As)
      INTEGER NNN_Se(8*SIZ_Se)
      PARAMETER (SIZ_Br= 3, OFF_Br=OFF_Se+SIZ_Se)
      INTEGER NNN_Br(8*SIZ_Br)
      PARAMETER (SIZ_Kr= 3, OFF_Kr=OFF_Br+SIZ_Br)
      INTEGER NNN_Kr(8*SIZ_Kr)
      PARAMETER (SIZ_Rb= 3, OFF_Rb=OFF_Kr+SIZ_Kr)
      INTEGER NNN_Rb(8*SIZ_Rb)
      PARAMETER (SIZ_Sr= 3, OFF_Sr=OFF_Rb+SIZ_Rb)
      INTEGER NNN_Sr(8*SIZ_Sr)
      PARAMETER (SIZ_Y = 3, OFF_Y =OFF_Sr+SIZ_Sr)
      INTEGER NNN_Y (8*SIZ_Y )
      PARAMETER (SIZ_Zr= 3, OFF_Zr=OFF_Y +SIZ_Y )
      INTEGER NNN_Zr(8*SIZ_Zr)
      PARAMETER (SIZ_Nb= 3, OFF_Nb=OFF_Zr+SIZ_Zr)
      INTEGER NNN_Nb(8*SIZ_Nb)
      PARAMETER (SIZ_Mo= 3, OFF_Mo=OFF_Nb+SIZ_Nb)
      INTEGER NNN_Mo(8*SIZ_Mo)
      PARAMETER (SIZ_Tc= 3, OFF_Tc=OFF_Mo+SIZ_Mo)
      INTEGER NNN_Tc(8*SIZ_Tc)
      PARAMETER (SIZ_Ru= 3, OFF_Ru=OFF_Tc+SIZ_Tc)
      INTEGER NNN_Ru(8*SIZ_Ru)
      PARAMETER (SIZ_Rh= 3, OFF_Rh=OFF_Ru+SIZ_Ru)
      INTEGER NNN_Rh(8*SIZ_Rh)
      PARAMETER (SIZ_Pd= 3, OFF_Pd=OFF_Rh+SIZ_Rh)
      INTEGER NNN_Pd(8*SIZ_Pd)
      PARAMETER (SIZ_Ag= 3, OFF_Ag=OFF_Pd+SIZ_Pd)
      INTEGER NNN_Ag(8*SIZ_Ag)
      PARAMETER (SIZ_Cd= 3, OFF_Cd=OFF_Ag+SIZ_Ag)
      INTEGER NNN_Cd(8*SIZ_Cd)
      PARAMETER (SIZ_In= 3, OFF_In=OFF_Cd+SIZ_Cd)
      INTEGER NNN_In(8*SIZ_In)
      PARAMETER (SIZ_Sn= 3, OFF_Sn=OFF_In+SIZ_In)
      INTEGER NNN_Sn(8*SIZ_Sn)
      PARAMETER (SIZ_Sb= 3, OFF_Sb=OFF_Sn+SIZ_Sn)
      INTEGER NNN_Sb(8*SIZ_Sb)
      PARAMETER (SIZ_Te= 3, OFF_Te=OFF_Sb+SIZ_Sb)
      INTEGER NNN_Te(8*SIZ_Te)
      PARAMETER (SIZ_I = 3, OFF_I =OFF_Te+SIZ_Te)
      INTEGER NNN_I (8*SIZ_I )
      PARAMETER (SIZ_Xe= 3, OFF_Xe=OFF_I +SIZ_I )
      INTEGER NNN_Xe(8*SIZ_Xe)
      PARAMETER (SIZ_Cs= 3, OFF_Cs=OFF_Xe+SIZ_Xe)
      INTEGER NNN_Cs(8*SIZ_Cs)
      PARAMETER (SIZ_Ba= 3, OFF_Ba=OFF_Cs+SIZ_Cs)
      INTEGER NNN_Ba(8*SIZ_Ba)
      PARAMETER (SIZ_La= 3, OFF_La=OFF_Ba+SIZ_Ba)
      INTEGER NNN_La(8*SIZ_La)
      PARAMETER (SIZ_Ce= 3, OFF_Ce=OFF_La+SIZ_La)
      INTEGER NNN_Ce(8*SIZ_Ce)
      PARAMETER (SIZ_Pr= 4, OFF_Pr=OFF_Ce+SIZ_Ce)
      INTEGER NNN_Pr(8*SIZ_Pr)
      PARAMETER (SIZ_Nd= 4, OFF_Nd=OFF_Pr+SIZ_Pr)
      INTEGER NNN_Nd(8*SIZ_Nd)
      PARAMETER (SIZ_Pm= 3, OFF_Pm=OFF_Nd+SIZ_Nd)
      INTEGER NNN_Pm(8*SIZ_Pm)
      PARAMETER (SIZ_Sm= 3, OFF_Sm=OFF_Pm+SIZ_Pm)
      INTEGER NNN_Sm(8*SIZ_Sm)
      PARAMETER (SIZ_Eu= 4, OFF_Eu=OFF_Sm+SIZ_Sm)
      INTEGER NNN_Eu(8*SIZ_Eu)
      PARAMETER (SIZ_Gd= 3, OFF_Gd=OFF_Eu+SIZ_Eu)
      INTEGER NNN_Gd(8*SIZ_Gd)
      PARAMETER (SIZ_Tb= 3, OFF_Tb=OFF_Gd+SIZ_Gd)
      INTEGER NNN_Tb(8*SIZ_Tb)
      PARAMETER (SIZ_Dy= 3, OFF_Dy=OFF_Tb+SIZ_Tb)
      INTEGER NNN_Dy(8*SIZ_Dy)
      PARAMETER (SIZ_Ho= 3, OFF_Ho=OFF_Dy+SIZ_Dy)
      INTEGER NNN_Ho(8*SIZ_Ho)
      PARAMETER (SIZ_Er= 3, OFF_Er=OFF_Ho+SIZ_Ho)
      INTEGER NNN_Er(8*SIZ_Er)
      PARAMETER (SIZ_Tm= 3, OFF_Tm=OFF_Er+SIZ_Er)
      INTEGER NNN_Tm(8*SIZ_Tm)
      PARAMETER (SIZ_Yb= 3, OFF_Yb=OFF_Tm+SIZ_Tm)
      INTEGER NNN_Yb(8*SIZ_Yb)
      PARAMETER (SIZ_Lu= 3, OFF_Lu=OFF_Yb+SIZ_Yb)
      INTEGER NNN_Lu(8*SIZ_Lu)
      PARAMETER (SIZ_Hf= 3, OFF_Hf=OFF_Lu+SIZ_Lu)
      INTEGER NNN_Hf(8*SIZ_Hf)
      PARAMETER (SIZ_Ta= 3, OFF_Ta=OFF_Hf+SIZ_Hf)
      INTEGER NNN_Ta(8*SIZ_Ta)
      PARAMETER (SIZ_W = 3, OFF_W =OFF_Ta+SIZ_Ta)
      INTEGER NNN_W (8*SIZ_W )
      PARAMETER (SIZ_Re= 3, OFF_Re=OFF_W +SIZ_W )
      INTEGER NNN_Re(8*SIZ_Re)
      PARAMETER (SIZ_Os= 3, OFF_Os=OFF_Re+SIZ_Re)
      INTEGER NNN_Os(8*SIZ_Os)
      PARAMETER (SIZ_Ir= 3, OFF_Ir=OFF_Os+SIZ_Os)
      INTEGER NNN_Ir(8*SIZ_Ir)
      PARAMETER (SIZ_Pt= 3, OFF_Pt=OFF_Ir+SIZ_Ir)
      INTEGER NNN_Pt(8*SIZ_Pt)
      PARAMETER (SIZ_Au= 3, OFF_Au=OFF_Pt+SIZ_Pt)
      INTEGER NNN_Au(8*SIZ_Au)
      PARAMETER (SIZ_Hg= 3, OFF_Hg=OFF_Au+SIZ_Au)
      INTEGER NNN_Hg(8*SIZ_Hg)
      PARAMETER (SIZ_Tl= 3, OFF_Tl=OFF_Hg+SIZ_Hg)
      INTEGER NNN_Tl(8*SIZ_Tl)
      PARAMETER (SIZ_Pb= 3, OFF_Pb=OFF_Tl+SIZ_Tl)
      INTEGER NNN_Pb(8*SIZ_Pb)
      PARAMETER (SIZ_Bi= 3, OFF_Bi=OFF_Pb+SIZ_Pb)
      INTEGER NNN_Bi(8*SIZ_Bi)
      PARAMETER (SIZ_Po= 3, OFF_Po=OFF_Bi+SIZ_Bi)
      INTEGER NNN_Po(8*SIZ_Po)
      PARAMETER (SIZ_At= 3, OFF_At=OFF_Po+SIZ_Po)
      INTEGER NNN_At(8*SIZ_At)
      PARAMETER (SIZ_Rn= 3, OFF_Rn=OFF_At+SIZ_At)
      INTEGER NNN_Rn(8*SIZ_Rn)
      PARAMETER (SIZ_Fr= 3, OFF_Fr=OFF_Rn+SIZ_Rn)
      INTEGER NNN_Fr(8*SIZ_Fr)
      PARAMETER (SIZ_Ra= 3, OFF_Ra=OFF_Fr+SIZ_Fr)
      INTEGER NNN_Ra(8*SIZ_Ra)
      PARAMETER (SIZ_Ac= 3, OFF_Ac=OFF_Ra+SIZ_Ra)
      INTEGER NNN_Ac(8*SIZ_Ac)
      PARAMETER (SIZ_Th= 3, OFF_Th=OFF_Ac+SIZ_Ac)
      INTEGER NNN_Th(8*SIZ_Th)
      PARAMETER (SIZ_Pa= 3, OFF_Pa=OFF_Th+SIZ_Th)
      INTEGER NNN_Pa(8*SIZ_Pa)
      PARAMETER (SIZ_U = 3, OFF_U =OFF_Pa+SIZ_Pa)
      INTEGER NNN_U (8*SIZ_U )
      PARAMETER (SIZ_Np= 3, OFF_Np=OFF_U +SIZ_U )
      INTEGER NNN_Np(8*SIZ_Np)
      PARAMETER (SIZ_Pu= 3, OFF_Pu=OFF_Np+SIZ_Np)
      INTEGER NNN_Pu(8*SIZ_Pu)
      PARAMETER (SIZ_Am= 3, OFF_Am=OFF_Pu+SIZ_Pu)
      INTEGER NNN_Am(8*SIZ_Am)
      PARAMETER (SIZ_Cm= 3, OFF_Cm=OFF_Am+SIZ_Am)
      INTEGER NNN_Cm(8*SIZ_Cm)
      PARAMETER (SIZ_Bk= 3, OFF_Bk=OFF_Cm+SIZ_Cm)
      INTEGER NNN_Bk(8*SIZ_Bk)
      PARAMETER (SIZ_Cf= 3, OFF_Cf=OFF_Bk+SIZ_Bk)
      INTEGER NNN_Cf(8*SIZ_Cf)
      PARAMETER (SIZ_Es= 3, OFF_Es=OFF_Cf+SIZ_Cf)
      INTEGER NNN_Es(8*SIZ_Es)

      PARAMETER (NTABLE=OFF_Es+SIZ_Es-1)
      INTEGER NNNPFN(8,NTABLE)

      EQUIVALENCE (NNNPFN(1,OFF_H ),NNN_H (1))
      EQUIVALENCE (NNNPFN(1,OFF_He),NNN_He(1))
      EQUIVALENCE (NNNPFN(1,OFF_Li),NNN_Li(1))
      EQUIVALENCE (NNNPFN(1,OFF_Be),NNN_Be(1))
      EQUIVALENCE (NNNPFN(1,OFF_B ),NNN_B (1))
      EQUIVALENCE (NNNPFN(1,OFF_C ),NNN_C (1))
      EQUIVALENCE (NNNPFN(1,OFF_N ),NNN_N (1))
      EQUIVALENCE (NNNPFN(1,OFF_O ),NNN_O (1))
      EQUIVALENCE (NNNPFN(1,OFF_F ),NNN_F (1))
      EQUIVALENCE (NNNPFN(1,OFF_Ne),NNN_Ne(1))
      EQUIVALENCE (NNNPFN(1,OFF_Na),NNN_Na(1))
      EQUIVALENCE (NNNPFN(1,OFF_Mg),NNN_Mg(1))
      EQUIVALENCE (NNNPFN(1,OFF_Al),NNN_Al(1))
      EQUIVALENCE (NNNPFN(1,OFF_Si),NNN_Si(1))
      EQUIVALENCE (NNNPFN(1,OFF_P ),NNN_P (1))
      EQUIVALENCE (NNNPFN(1,OFF_S ),NNN_S (1))
      EQUIVALENCE (NNNPFN(1,OFF_Cl),NNN_Cl(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ar),NNN_Ar(1))
      EQUIVALENCE (NNNPFN(1,OFF_K ),NNN_K (1))
      EQUIVALENCE (NNNPFN(1,OFF_Ca),NNN_Ca(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sc),NNN_Sc(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ti),NNN_Ti(1))
      EQUIVALENCE (NNNPFN(1,OFF_V ),NNN_V (1))
      EQUIVALENCE (NNNPFN(1,OFF_Cr),NNN_Cr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Mn),NNN_Mn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Fe),NNN_Fe(1))
      EQUIVALENCE (NNNPFN(1,OFF_Co),NNN_Co(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ni),NNN_Ni(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cu),NNN_Cu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Zn),NNN_Zn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ga),NNN_Ga(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ge),NNN_Ge(1))
      EQUIVALENCE (NNNPFN(1,OFF_As),NNN_As(1))
      EQUIVALENCE (NNNPFN(1,OFF_Se),NNN_Se(1))
      EQUIVALENCE (NNNPFN(1,OFF_Br),NNN_Br(1))
      EQUIVALENCE (NNNPFN(1,OFF_Kr),NNN_Kr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Rb),NNN_Rb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sr),NNN_Sr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Y ),NNN_Y (1))
      EQUIVALENCE (NNNPFN(1,OFF_Zr),NNN_Zr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Nb),NNN_Nb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Mo),NNN_Mo(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tc),NNN_Tc(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ru),NNN_Ru(1))
      EQUIVALENCE (NNNPFN(1,OFF_Rh),NNN_Rh(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pd),NNN_Pd(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ag),NNN_Ag(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cd),NNN_Cd(1))
      EQUIVALENCE (NNNPFN(1,OFF_In),NNN_In(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sn),NNN_Sn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sb),NNN_Sb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Te),NNN_Te(1))
      EQUIVALENCE (NNNPFN(1,OFF_I ),NNN_I (1))
      EQUIVALENCE (NNNPFN(1,OFF_Xe),NNN_Xe(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cs),NNN_Cs(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ba),NNN_Ba(1))
      EQUIVALENCE (NNNPFN(1,OFF_La),NNN_La(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ce),NNN_Ce(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pr),NNN_Pr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Nd),NNN_Nd(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pm),NNN_Pm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Sm),NNN_Sm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Eu),NNN_Eu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Gd),NNN_Gd(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tb),NNN_Tb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Dy),NNN_Dy(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ho),NNN_Ho(1))
      EQUIVALENCE (NNNPFN(1,OFF_Er),NNN_Er(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tm),NNN_Tm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Yb),NNN_Yb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Lu),NNN_Lu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Hf),NNN_Hf(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ta),NNN_Ta(1))
      EQUIVALENCE (NNNPFN(1,OFF_W ),NNN_W (1))
      EQUIVALENCE (NNNPFN(1,OFF_Re),NNN_Re(1))
      EQUIVALENCE (NNNPFN(1,OFF_Os),NNN_Os(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ir),NNN_Ir(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pt),NNN_Pt(1))
      EQUIVALENCE (NNNPFN(1,OFF_Au),NNN_Au(1))
      EQUIVALENCE (NNNPFN(1,OFF_Hg),NNN_Hg(1))
      EQUIVALENCE (NNNPFN(1,OFF_Tl),NNN_Tl(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pb),NNN_Pb(1))
      EQUIVALENCE (NNNPFN(1,OFF_Bi),NNN_Bi(1))
      EQUIVALENCE (NNNPFN(1,OFF_Po),NNN_Po(1))
      EQUIVALENCE (NNNPFN(1,OFF_At),NNN_At(1))
      EQUIVALENCE (NNNPFN(1,OFF_Rn),NNN_Rn(1))
      EQUIVALENCE (NNNPFN(1,OFF_Fr),NNN_Fr(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ra),NNN_Ra(1))
      EQUIVALENCE (NNNPFN(1,OFF_Ac),NNN_Ac(1))
      EQUIVALENCE (NNNPFN(1,OFF_Th),NNN_Th(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pa),NNN_Pa(1))
      EQUIVALENCE (NNNPFN(1,OFF_U ),NNN_U (1))
      EQUIVALENCE (NNNPFN(1,OFF_Np),NNN_Np(1))
      EQUIVALENCE (NNNPFN(1,OFF_Pu),NNN_Pu(1))
      EQUIVALENCE (NNNPFN(1,OFF_Am),NNN_Am(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cm),NNN_Cm(1))
      EQUIVALENCE (NNNPFN(1,OFF_Bk),NNN_Bk(1))
      EQUIVALENCE (NNNPFN(1,OFF_Cf),NNN_Cf(1))
      EQUIVALENCE (NNNPFN(1,OFF_Es),NNN_Es(1))
      SAVE NNNPFN,LOCZ,SCALE,FIRST
C      ( 1)( 2)  ( 3)( 4)  ( 5)( 6)  ( 7)( 8)  ( 9)(10)   ( IP )   G Ion REF
      DATA NNN_H/
     1 200020001,200020011,201620881,231228281,378953411, 1359502, 1,00,D+F H  1
     2 100010001,100010001,100010001,100010001,100010001, 1359500, 1,01/G   H  2
      DATA NNN_He/
     1 100010001,100010011,102111241,145022061,363059451, 2458104, 2,00,D+F He 1
     2 200020001,200020071,208524971,382669341,128222452, 5440302, 2,01,D+F He 2
     3 100010001,100010001,100010001,100010001,100010001, 5440300, 2,02/G   He 3
      DATA NNN_Li/
     1 200020011,201220481,212922881,258731081,394251691,  538901, 3,00,D+F Li 1
     2 100010001,100010201,126225521, 67216512,351165562, 7561907, 3,01,D+F Li 2
     3 200020001,200020211,227936571, 69610342,137217102,12241800, 3,02,D+F Li 3
     4 100010001,100010001,100010001,100010001,100010001,12241800, 3,03/G   Li 4
      DATA NNN_Be/
     1 100010051,104311441,131615641,190623681,298037691,  931900, 4,00,AEL Be 1
     2 200120231,211422771,249627631,309034911,398545051, 1820600, 4,01,AEL Be 2
     3 100010001,100010201,126225521, 67216512,351165562,15385000, 4,02,AEL Be 3
     4 200020001,200020011,201220661,223426161,332644691,21765700, 4,03/AEL Be 4
      DATA NNN_B/
     1 600060001,600560281,608761991,637466191,693973361,  829500, 5,00,AEL B  1
     2 100310831,132016901,214226411,315736741,419147071, 2514900, 5,01,AEL B  2
     3 200721061,233526401,297533311,369040481,440747651, 3792000, 5,02,AEL B  3
     4 100010001,100010001,100010001,100010001,100010001,25929800, 5,03/G   B  4
      DATA NNN_C/
     1 893292271, 96110042,105311262,126315202,196126432, 1125508, 6,00,D+F C  1
     2 595060251,620865751,713280191, 95712292,167623542, 2437501, 6,01,D+F C  2
     3 105513201,180324851,341851341, 88416332,296550722, 4787101, 6,02,D+F C  3
     4 204922771,262630421,350941931,494556971,644872001, 6447600, 6,03,D+F C  4
     5 100010001,100010001,100010001,100010001,100010001,39207700, 6,04,G   C  5
     6 200020001,200020001,200020001,200020001,200020001,48998100, 6,05/G   C  6
      DATA NNN_N/
     1 403141851,457051681,594071181, 92913362,203331152, 1452915, 7,00,D+F N  1
     2 919899541,107211512,124914302,182526232,403762662, 2959202, 7,01,D+F N  2
     3 596862721,684177081, 88110342,128317062,239334312, 4742501, 7,02,D+F N  3
     4 112816481,240733751,462068491,116419932,283736822, 7744900, 7,03,D+F N  4
     5 210124681,293634211,391145791,539862151,703178471, 9786200, 7,04,D+F N  5
     6 100010001,100010001,100010001,100010001,100010001,55205700, 7,05/G   N  6
      DATA NNN_O/
     1 874789691,924795711, 99410492,115213492,169022242, 1361307, 8,00,D+F O  1
     2 424151091,622874781, 91312832,221842502, 79914013, 3510711, 8,01,D+F O  2
     3  95610702,118113032,149619922,329761642,101914173, 5488500, 8,02,D+F O  3
     4 603567171,775391141,106612482,143716252,181420032, 7739300, 8,03,D+F O  4
     5 124420321,306943181,606281181,101712232,142916342,11387300, 8,04,D+F O  5
     6 215026541,323137551,421546491,508255151,594863811,13807900, 8,05/AEL O  6
      DATA NNN_F/
     1 575958511,589859231,595860671,636470031,815199581, 1741802, 9,00,D+F F  1
     2 900296401,102610802,113912542,152921152,318348952, 3498003, 9,01,D+F F  2
     3 469162651,791295541,121419552,402686872,154822203, 6264500, 9,02,D+F F  3
     4  99511422,129214572,170523002,320140922,498458762, 8713900, 9,03,D+F F  4
     5 615472711, 87710602,127215002,172919582,218624152,11421300, 9,04,D+F F  5
     6 135324181,377252001,661580261, 94410852,122613672,15711700, 9,05/AEL F  6
      DATA NNN_Ne/
     1 100010001,100010051,105313051,210239461, 74013022, 2155808,10,00,D+F Ne 1
     2 580158751,591759741,642687101,159332652, 64111533, 4106907,10,01,D+F Ne 2
     3  93510272,110411662,127116062,257647882, 75110223, 6350000,10,02,D+F Ne 3
     4 529774371, 94611322,135816202,188221442,240626682, 9701900,10,03,D+F Ne 4
     5 103312152,140616092,181320182,222224262,263128352,12630000,10,04,AEL Ne 5
     6 629178711, 98311802,136715512,173619202,210422892,15790900,10,05/AEL Ne 6
      DATA NNN_Na/
     1 200020001,200320211,207322131,253031421,417657451,  513802,11,00,D+F Na 1
     2 100010001,100010161,119621261, 50711872,246445382, 4728901,11,01,D+F Na 2
     3 580158751,591860351, 71813142,321968812,106014333, 7165000,11,02,D+F Na 3
     4  96910772,116012242,130714232,153916552,177118872, 9888000,11,03,D+F Na 4
     5 601386081,108812932,148916832,187820722,226624612,13836900,11,04,AEL Na 5
     6 105712442,144616652,189221182,234425702,279630222,17209000,11,05/AEL Na 6
      DATA NNN_Mg/
     1 100010011,101410621,118414581,204831781,509479731,  764404,12,00,D+F Mg 1
     2 200120051,202921001,226926901,368457091, 92814872, 1503101,12,01,D+F Mg 2
     3 100010001,100110611,177455431,176546012, 99718753, 8011905,12,02,D+F Mg 3
     4 579758751,591459501,600560591,611461681,622362781,10928900,12,03,AEL Mg 4
     5 100611232,120612752,134214102,147815462,161416822,14122900,12,04,AEL Mg 5
     6 674896701,121814462,167018942,211723412,256527892,18648900,12,05/AEL Mg 6
      DATA NNN_Al/
     1 558857701,583558761,593260591,635969541,796790971,  598400,13,00,D+F Al 1
     2 100310211,110313021,172828201, 55311252,215637942, 1882203,13,01,D+F Al 2
     3 200320201,208622331,250530971,410251081,611571211, 2844000,13,02,D+F Al 3
     4 100010001,100210881,207436531,523168101,838999681,11996000,13,03,D+F Al 4
     5 577758651,591259631,604461351,622563161,640764981,15377000,13,04,AEL Al 5
     6 103511582,124713242,140014772,155316292,170517812,19042000,13,05/AEL Al 6
      DATA NNN_Si/
     1 825189211, 95210052,106211532,134317202,237934082,  814913,14,00,D+F Si 1
     2 563057761,588160311,631768671,791097651,127817282, 1634000,14,01,D+F Si 2
     3 101110771,126716471,232438081, 71914052,262045302, 3346001,14,02,D+F Si 3
     4 200720521,217224081,284439171,551370951, 86810262, 4513000,14,03,D+F Si 4
     5 100010001,100210881,207436531,523168101,838999681,16672900,14,04,FAK Si 5
     6 575458521,591459851,610063201,672674071,843698661,20510900,14,05/AEL Si 6
      DATA NNN_P/
     1 402643441,496757481,658274401,833492941,103511532, 1048300,15,00,AEL P  1
     2 874497931,106011282,119812802,138415142,164717802, 1972000,15,01,AEL P  2
     3 564058061,604164611,709579551, 90410172,112912422, 3015500,15,02,AEL P  3
     4 100811411,149720221,280936121,441552181,602168241, 5135400,15,03,AEL P  4
     5 200420781,227025361,281430911,336936471,392542021, 6500700,15,04,AEL P  5
     6 100010001,100010001,100010001,100010001,100010001,22041300,15,05/G   P  6
      DATA NNN_S/
     1 822887891,930697831,102610932,121614492,185124742, 1035708,16,00,D+F S  1
     2 443056011,694982961, 96911522,144218572,227326892, 2339900,16,01,D+F S  2
     3  91610392,113512242,136416942,233429882,364242962, 3500000,16,02,D+F S  3
     4 560058861,633871081, 82410062,123314602,168619132, 4728900,16,03,D+F S  4
     5 104512901,177025421,375163021,122420462,286036742, 7250000,16,04,D+F S  5
     6 202321571,241428261,358355061, 78310152,124814802, 8802800,16,05/D+F S  6
      DATA NNN_Cl/
     1 538155931,571657911,598067191, 89013782,227737172, 1300916,17,00,D+F Cl 1
     2 873396771,104411072,118513532,175525872,406763932, 2379903,17,01,D+F Cl 2
     3 506569571, 87610522,134421682,439092662,182132573, 3990006,17,02,D+F Cl 3
     4  95110872,120013232,154921252,345149322,641378942, 5350000,17,03,D+F Cl 4
     5 558960371,677779341, 95311692,138816082,182720472, 6780000,17,04/D+F Cl 5
      DATA NNN_Ar/
     1 100010001,100010051,106913911,240147261, 90716112, 1575411,18,00,D+F Ar 1
     2 550256831,578158781,636585461,151530162, 58010303, 2762007,18,01,D+F Ar 2
     3  92110362,112412002,133216772,254443722, 76512833, 4090003,18,02,D+F Ar 3
     4 582082081,103112292,149920212,309750502,720793642, 5978900,18,03,D+F Ar 4
     5  97111072,123213982,172625622,463976582,106413633, 7500000,18,04/D+F Ar 5
      DATA NNN_K/
     1 200020011,200720361,211923291,280137141,525575741,  433803,19,00,D+F K  1
     2 100010001,100110341,135929551, 79119282,405274892, 3180905,19,01,D+F K  2
     3 554657081,581260301, 73012702,285363872,129023363, 4600005,19,02,D+F K  3
     4  96010862,118413212,180836632, 90321023,416863253, 6090000,19,03,D+F K  4
     5 657793361,119515082,195826322,352944302,533162332, 8259900,19,04/D+F K  5
      DATA NNN_Ca/
     1 100110061,104311741,145919971,294345051, 69010322,  611003,20,00,D+F Ca 1
     2 205822781,279234761,427553061,688994901,136319772, 1186701,20,01,D+F Ca 2
     3 100010001,100510821,168744821,130232522, 69012813, 5121003,20,02,D+F Ca 3
     4 555157161,585662471, 82816862, 42510013,168423663, 6700000,20,03,D+F Ca 4
     5  99411262,123814062,182930402,484766392, 84310223, 8438900,20,04/D+F Ca 5
      DATA NNN_Sc/
     1 924696691,105212282,151219062,240530032,368944512,  653900,21,00,AEL Sc 1
     2 190424662,297634542,391743752,482952832,573761912, 1280000,21,01,AEL Sc 2
     3 976799291,101110322,105810882,111911502,118112122, 2475000,21,02,AEL Sc 3
     4 100010001,100510821,168744821,130232522, 69012813, 7390000,21,03,FAK Sc 4
     5 555157161,585662471, 82816862, 42510013,168423663, 9200000,21,04/FAK Sc 5
      DATA NNN_Ti/
     1 181021172,260333222,430155582,710089242,110213293,  681900,22,00,D+F Ti 1
     2 474659872,721284672, 98211413,134515623,177919963, 1356900,22,01,D+F Ti 2
     3 228327012,308134272,381143862,534563472,734983512, 2747000,22,02,D+F Ti 3
     4 971498311, 99210032,102610572,108711172,114711782, 4324000,22,03,D+F Ti 4
     5 100010001,100510821,168744821,130232522, 69012813, 9980000,22,04/FAK Ti 5
      DATA NNN_V/
     1 272835172,425851532,632278322, 97212013,146817723,  674000,23,00,AEL V  1
     2 373954132,743597002,121414713,173920143,229225713, 1464900,23,01,AEL V  2
     3 323142642,519660272,679975352,824789522, 96610363, 2930900,23,02,AEL V  3
     4 248329302,324234952,373439752,421744582,469949412, 4800000,23,03,AEL V  4
     5 970698231,990699881,100710152,102410322,104010482, 6500000,23,04/AEL V  5
      DATA NNN_Cr/
     1 717277611, 92911652,152620872,295141952,550468122,  676400,24,00,D+F Cr 1
     2  71611552,205635512,558281952,115315823,205625293, 1649000,24,01,D+F Cr 2
     3 280639822,538369722, 87610823,129115003,170919183, 3095000,24,02,D+F Cr 3
     4 377150952,616070292,791788382, 97610683,116012523, 5000000,24,03,D+F Cr 4
     5 264730962,341436462,394042872,463549832,533056782, 7300000,24,04/D+F Cr 5
      DATA NNN_Mn/
     1 600060321,629270891, 86911302,151020222,267534752,  743100,25,00,AEL Mn 1
     2 739594821,139921212,309342852,567372412, 97112553, 1563600,25,01,AEL Mn 2
     3  98417472,265535782,454754842,641973532,828792212, 3369000,25,02,AEL Mn 3
     4 328847052,586668342,771785912, 94710343,112112093, 5300000,25,03,AEL Mn 4
     5 422055132,636770792,779285062,921999322,106411363, 7600000,25,04/AEL Mn 5
      DATA NNN_Fe/
     1 197023222,274433302,416753952,723799822,139419053,  787038,26,00,D+F Fe 1
     2 409453722,686687452,110213823,174322233,286437043, 1617902,26,01,D+F Fe 2
     3 262136422,501167232, 87911303,138916483,190721673, 3064300,26,02,D+F Fe 3
     4  98723522,420363072, 87011423,145117913,215925463, 5700000,26,03,AEL Fe 4
     5 388854482,666275742,846693572,102511143,120312923, 7900000,26,04/D+F Fe 5
      DATA NNN_Co/
     1 199427202,335740022,474957182,708090462,118315403,  786000,27,00,D+F Co 1
     2 279739202,490858232,684582472,104713233,159818733, 1704900,27,01,D+F Co 2
     3 279836622,461857562,720693022,124915873,192522633, 3349000,27,02,D+F Co 3
     4 262136422,501167232, 87911303,138916483,190821673, 5300000,27,03,FAK Co 4
     5  98723522,420363072, 87011423,145117913,215925463, 8300000,27,04/FAK Co 5
      DATA NNN_Ni/
     1 227027622,306233052,356839222,446052912,652382292,  763314,28,00,D+F Ni 1
     2 108416342,222428472,353944332,577378932,110314303, 1814900,28,01,D+F Ni 2
     3 198724282,293236452,468362702, 86511123,136016073, 3516000,28,02,D+F Ni 3
     4 279836622,461857562,720693022,124915873,192522633, 5600000,28,03,FAK Ni 4
     5 262136422,501167232, 87911303,138916483,190721673, 7900000,28,04/FAK Ni 5
      DATA NNN_Cu/
     1 201620781,231026761,314737361,450555381,692386911,  772301,29,00,D+F Cu 1
     2 109415761,247938311, 58910042,190937022, 68311693, 2028903,29,01,D+F Cu 2
     3 897195961,107212972,165021182,260230862,356940532, 3682900,29,02/D+F Cu 3
      DATA NNN_Zn/
     1 100010001,100410231,108712611,167124841,388460411,  939102,30,00,D+F Zn 1
     2 200020021,201620761,223726341,351352061, 80812472, 1796001,30,01,D+F Zn 2
     3 100610471,122617301,300566361,149924112,332342352, 3970000,30,02/D+F Zn 3
      DATA NNN_Ga/
     1 403245601,493151431,529654331,559358091,611065171,  600000,31,00,AEL Ga 1
     2  99710051,104511541,135016501,208226431,321837921, 2050900,31,01,AEL Ga 2
     3 199820071,204521391,229124761,266028451,302932131, 3070000,31,02/AEL Ga 3
      DATA NNN_Ge/
     1 502665261,755183501,901496201,102410942,117912812,  787900,32,00,AEL Ge 1
     2 422848161,512153401,557458941,636270361,794489061, 1593000,32,01,AEL Ge 2
     3 100010261,114613921,175221251,249828711,324436181, 3421000,32,02/AEL Ge 3
      DATA NNN_As/
     1 403143241,491856701,649173781,840396751,113013392,  981000,33,00,AEL As 1
     2 593676641,884697521,105911572,129515012,180322212, 1858700,33,01,AEL As 2
     3 484470541, 91510972,125614082,157017612,199722912, 2829900,33,02/AEL As 3
      DATA NNN_Se/
     1 630172361,799686381,919797221,102810942,117712832,  975000,34,00,AEL Se 1
     2 438055511,691582151, 94510732,121413672,152016732, 2150000,34,01,AEL Se 2
     3 651982921, 94610382,113212492,139515462,169718482, 3200000,34,02/AEL Se 3
      DATA NNN_Br/
     1 437347431,498951671,538559501, 74710812,169126672, 1183910,35,00,D+F Br 1
     2 705183611, 93510092,111614162,222932532,427652992, 2160000,35,01,D+F Br 2
     3 510869921, 87410312,123116552,236530712,377744832, 3590000,35,02/D+F Br 3
      DATA NNN_Kr/
     1 100010001,100010051,105012781,198535971, 65911422, 1399507,36,00,D+F Kr 1
     2 461049811,522254261,609088131,168935052, 68612253, 2455908,36,01,D+F Kr 2
     3 759990901,101911142,129017782,302856642, 99414333, 3690000,36,02/D+F Kr 3
      DATA NNN_Rb/
     1 200020011,200720361,211523021,269434141,459163351,  417502,37,00,D+F Rb 1
     2 100010001,100110321,129524961, 61014202,291753192, 2750004,37,01,D+F Rb 2
     3 473650891,533156051, 66810932,232950852, 99915303, 4000000,37,02/D+F Rb 3
      DATA NNN_Sr/
     1 100110041,104111741,146019721,281941411,607785251,  569202,38,00,D+F Sr 1
     2 202621931,255331271,384347931,624085761,122417632, 1102600,38,01,D+F Sr 2
     3 100010001,100110321,129524961, 61014202,291753192, 4300000,38,02/FAK Sr 3
      DATA NNN_Y/
     1 791587851,100012192,155119942,254031782,389946932,  637900,39,00,AEL Y  1
     2 118217102,220827002,319036792,416646512,513256072, 1223000,39,01,AEL Y  2
     3  92510012,104710862,112311612,120212472,132814282, 2050000,39,02/AEL Y  3
      DATA NNN_Zr/
     1 141320802,291439702,531170262, 92712273,162521053,  684000,40,00,D+F Zr 1
     2 354454352,724689652,107212643,148517093,193321573, 1312900,40,01,D+F Zr 2
     3 209727032,324537052,415446282,510255752,604965222, 2298000,40,02/D+F Zr 3
      DATA NNN_Nb/
     1 256636022,465759302,749693962,116514243,171520333,  687900,41,00,AEL Nb 1
     2 335157222, 84511463,147718363,221826083,299933893, 1431900,41,01,AEL Nb 2
     3 223725352,280830972,340937362,406844002,473150632, 2503900,41,02/AEL Nb 3
      DATA NNN_Mo/
     1 703972941, 82610822,154822682,327244912,571469372,  709900,42,00,D+F Mo 1
     2  75714552,274347322,718897632,123414913,174920063, 1614900,42,01,D+F Mo 2
     3 267645462,669890262,115514323,173620673,242528083, 2714900,42,02/AEL Mo 3
      DATA NNN_Tc/
     1  90613732,184823562,291735332,419949102,565764332,  728000,43,00,AEL Tc 1
     2 131318312,227126932,311735452,397644072,483852692, 1525900,43,01,AEL Tc 2
c     3 204721673,234725733,284031463,348738613,426546943, 3000000,43,02/AEL Tc 3
C corrected following ATLAS12
     3 600460071,607964351,731488341,110013702,168420312, 3000000,43,02/AEL Tc 3
      DATA NNN_Ru/
     1 176824122,318941082,515263202,761790472,106112303,  736400,44,00,AEL Ru 1
     2 221934642,501968372, 88911173,136316243,189221613, 1675900,44,01,AEL Ru 2
     3 210622722,241025422,267928262,297731272,327834282, 2846000,44,02/AEL Ru 3
      DATA NNN_Rh/
     1 148520202,255230902,364942462,489656082,638872352,  746000,45,00,AEL Rh 1
     2 153421292,288137912,484660322,720187062,101011483, 1807000,45,01,AEL Rh 2
     3 254537212,492362292,770592182,107312243,137615273, 3104900,45,02/AEL Rh 3
      DATA NNN_Pd/
     1 115919651,320746011,607576761, 95011642,141817172,  832900,46,00,AEL Pd 1
     2 755087211,105913442,173122222,282034722,412247732, 1941900,46,01,AEL Pd 2
     3 180223462,289735212,414247632,538460052,662672472, 3292000,46,02/AEL Pd 3
      DATA NNN_Ag/
     1 200020001,200220141,206422141,257633021,455164681,  757403,47,00,D+F Ag 1
     2 100810581,125817401,260641031, 66210072,135316982, 2148000,47,01,D+F Ag 2
     3 795887491, 97711762,156620252,248329422,340038582, 3481900,47,02/D+F Ag 3
      DATA NNN_Cd/
     1 100010001,100410241,109212891,176827421,444268771,  899003,48,00,D+F Cd 1
     2 200020021,201720921,233329881,451475371,127520782, 1690301,48,01,D+F Cd 2
     3 100310281,114815371,246138311,519265531,791492761, 3747000,48,02/D+F Cd 3
      DATA NNN_In/
     1 252431921,368440461,433746521,512259221,723389021,  578400,49,00,D+F In 1
     2 100110071,104611651,146118581,225426511,304734431, 1886000,49,01,D+F In 2
     3 200120111,205021611,243628031,317035371,390442701, 2802900,49,02/D+F In 3
      DATA NNN_Sn/
     1 232637101,488058571,669074381,816189091, 97210632,  734200,50,00,AEL Sn 1
     2 286335941,408144471,479351961,571862901,686274341, 1462700,50,01,AEL Sn 2
     3 100010251,114013811,175321601,256829751,338337901, 3049000,50,02/AEL Sn 3
      DATA NNN_Sb/
     1 404043481,494656811,646772781,813490751,101411372,  863900,51,00,AEL Sb 1
     2 303147981,618472951,827392621,103711702,131214532, 1650000,51,01,AEL Sb 2
     3 313037601,429347901,536260591,689477591,862494881, 2529900,51,02/AEL Sb 3
      DATA NNN_Te/
     1 526258801,657372351,784284071,897095741,102711082,  900900,52,00,AEL Te 1
     2 440855541,686481251, 93810792,125414792,176321132, 1860000,52,01,AEL Te 2
     3 349054751,699883081, 96611302,134216202,197724212, 2800000,52,02/AEL Te 3
      DATA NNN_I/
     1 405342041,438645621,475751071,587974491,102214572, 1045404,53,00,D+F I  1
     2 568567471,773485861, 94510362,112712182,130914002, 1909000,53,01,D+F I  2
     3 514269581, 86910562,130716652,215327742,351843662, 3200000,53,02/AEL I  3
      DATA NNN_Xe/
     1 100010001,100010091,109515351,291060661,119621482, 1212716,54,00,D+F Xe 1
     2 414844131,465649111,538464651, 87112232,158019362, 2120000,54,01,D+F Xe 2
     3 615475101,867797531,112213462,157618062,203622662, 3209900,54,02/D+F Xe 3
      DATA NNN_Cs/
     1 200020001,201020501,215623871,283536181,462756261,  389300,55,00,D+F Cs 1
     2 100010001,100310371,119016501,269146361, 77912412, 2510000,55,01,D+F Cs 2
     3 424445601,481750061,516953311,549356551,581759791, 3500000,55,02/D+F Cs 3
      DATA NNN_Ba/
     1 101210791,135119351,282340571,574580391,111015062,  521002,56,00,D+F Ba 1
     2 262638611,504160621,698579371, 91010692,129115952, 1000000,56,01,D+F Ba 2
     3 100010001,100310351,118416321,264945521, 76512182, 3700000,56,02/FAK Ba 3
      DATA NNN_La/
     1  71111992,172323592,312540402,510763182,765791012,  557700,57,00,AEL La 1
     2 204529582,383647882,582469262,807992692,104911723, 1106000,57,01,AEL La 2
     3  94712552,148416582,179819212,203621522,227424042, 1917700,57,02/AEL La 3
      DATA NNN_Ce/
     1 295959132,103515693,215527593,335939413,449650223,  553870,58,00,AEL Ce 1
     2  80118633,304342383,541765723,769387773, 98210814, 1085000,58,01,MZH Ce 2
     3 506183092,108612923,146416133,174418603,196520603, 2020000,58,02/CCB Ce 3
      DATA NNN_Pr/
C     1 460693672,158523823,327242303,519661563,709379783,  546400,59,00,FAK Pr 1
C     2 455480232,114014653,178521013,240927073,299232633, 1055000,59,01,AEL Pr 2
C     3  46410533,183826893,354443773,518459633,674375243, 2162400,59,02/AEL Pr 3
C Pr I corrected following MOOG (Sneden)
     1 146526632,508289352,142720943,287237333,465456163,  547300,59,00,Sne Pr 1
C Pr II corrected following ISAN's theoretical calculations of Pr II atomic levels
     2  53014313,298951373, 77710794,141117664,213725224, 1055000,59,01,ISA Pr 2
C Pr III calculated using new energy levels from Wyart et al. 2006, Phys. Scr. (in press)
     3 421093902,165924663,331041793,507660143,700980743, 2162400,59,02,ISA Pr 3
C Pr IV calculated using NIST energy levels
     4 373649462,593368882,785988552, 98810923,119813043, 3900000,59,03/AEL Pr 4
      DATA NNN_Nd/
C     1 139623042,364860002, 96114603,209828633,373446973,  552500,60,00,AEL Nd 1
C     2 460493692,158523823,327142303,519661563,709279783, 1073000,60,01,AEL Nd 2
C     3 455480232,114014653,178521013,240927073,299232633, 2218000,60,02/BOR Nd 3
C Nd I corrected following MOOG (Sneden)
     1 145623072,410172132,120218793,276138313,505263693,  552500,60,00,Sne Nd 1
C Nd II corrected following ISAN's theoretical calculations of Nd II atomic levels (Mashonkina et al. 2005, A&A, 441, 309
     2  47511303,223037433,559777223,100512564,151817894, 1073000,60,01,ISA Nd 2
C Nd III corrected following ISAN's theoretical calculations of Nd III atomic levels (Ryabchikova et al. 2006, A&A, 456, 329
     3 432699302,204835193,525971403, 90710984,128314614, 2218000,60,02,ISA Nd 3
C Nd IY calculated following using energy levels from Wyart et al. 2006, J. Phys. B39, L77
     4 104717683,241529543,339937663,407343323,455447453, 4042000,60,03/Wyt Nd 4
      DATA NNN_Pm/
     1 131720482,280535692,441254492,676583972,103412583,  555400,61,00,AEL Pm 1
     2 139623042,364860002, 96114603,209828633,373446973, 1090000,61,01,FAK Pm 2
     3 460493682,158523823,327142303,519661563,709279783, 2230000,61,02/FAK Pm 3
      DATA NNN_Sm/
     1  92915672,222431062,444763802, 89612173,159520253,  564370,62,00,AEL Sm 1
     2 315059662, 97114563,204627093,342541693,490556383, 1106900,62,01,AEL Sm 2
     3 269037812,520270372, 91111273,133915483,172719093, 2340000,62,02/AEL Sm 3
      DATA NNN_Eu/
     1 800080571,851699301,127617362,240433032,444958442,  567045,63,00,AEL Eu 1
     2 125416052,211828182,375549622,644381732,101112213, 1124100,63,01,AEL Eu 2
C     3 800080571,851699301,127617362,240433032,444958442, 2492000,63,02,FAK Eu 3
C Eu III calculated using new energy levels from Wyart et al. 2006, Phys. Scr. (in press)
     3  82514782, 47913863,315459503, 98114674,204226924, 2492000,63,02,ISA Eu 3
C Eu IV g-factor fot the ground level
     4 490049002,490049002,490049002,490049002,490049002, 4000000,63,03/TAR Eu 4
      DATA NNN_Gd/
C     1 240432982,427555202,708489962,112613853,167319843,  615000,64,00,AEL Gd 1
C Pr I corrected following MOOG (Sneden)
     1 244232982,441460242,826089963,112613853,167319843,  615000,64,00,Sne Gd 1
     2 534793262,139219123,247730843,371043333,495055893, 1209000,64,01,AEL Gd 2
     3 364145232,514756362,604864112,673870372,732276072, 2063000,64,02/AEL Gd 3
      DATA NNN_Tb/
C     1 480767202, 89011393,144118243,230028753,354142883,  586390,65,00,AEL Tb 1
C     2 480767192, 89011393,144118243,230028753,354142883, 1151900,65,01,FAK Tb 2
C     3 480767202, 89011393,144118243,230028753,354142883, 2191000,65,02/FAK Tb 3
C Tb I corrected following MOOG (Sneden)
     1 546880382,113515623,209227313,347543173,524362333,  586390,65,00,Sne Tb 1
C Tb II corrected following MOOG (Sneden)
     2  56510823,163922043,279234353,417550623,615575303, 1151900,65,01,Sne Tb 2
C Tb III calculated using new energy levels from Wyart et al. 2006, Phys. Scr. (in press)
     3  53713323,276551143, 85012894,181224014,304037114, 2191000,65,02/ISA Tb 3
      DATA NNN_Dy/
C     1 343147532,645887152,115314793,183322063,257729373,  593890,66,00,FAK Dy 1
C     2 343147532,645887142,115314793,183322063,257729373, 1167000,66,01,AEL Dy 2
C     3 343147532,645887142,115314793,183322063,257729373, 2280000,66,02/FAK Dy 3
C Dy I, II corrected following MOOG (Sneden)
     1 175219662,262038952,604693903,142320733,288338103,  593890,66,00,Sne Dy 1
     2 347359162,108619003,300742453,533359923,606555733, 1167000,66,01,Sne Dy 2
c  Dy III partition function coefficients (Irwin 1981)
     3 370353332, 75410543,144819483,256633063,416951503, 2280000,66,02/Sne Dy 3
      DATA NNN_Ho/
     1 222635002,542276772,100312353,145716713,187020703,  602160,67,00,FAK Ho 1
C     2 222635002,542276772,100312353,145716713,187020703, 1180000,67,01,FAK Ho 2
C Ho II corrected following Bord & Cowley, 2002, Sol. Phys., 211, 3
     2 321455092,112322203,401966563,102014674,200226144, 1180000,67,01,Bor Ho 2
     3 222635002,542276772,100312353,145716713,187020703, 2284000,67,02/AEL Ho 3
      DATA NNN_Er/
C     1 133715382,209130152,429859382, 79410293,129815983,  610780,68,00,AEL Er 1
C     2 265934782,497877532,120517733,245032063,400448073, 1193000,68,01,AEL Er 2
C     3 265934782,497877532,120517733,245032063,400448073, 2274000,68,02/FAK Er 3
C Er I corrected following MOOG (Sneden)
     1 131715322,213632462,504577482,115416533,226829683,  610780,68,00,Sne Er 1
C Er II corrected following MOOG (Sneden)
     2 282946962, 81713443,201827463,339638403,399938623, 1193000,68,01,Sne Er 2
c  Er III partition function coefficients (Irwin 1981)
     3 801281851, 91511592,166126662,472591362,190642503, 2274000,68,02/FAK Er 3
      DATA NNN_Tm/
     1 800381111, 87510702,147621462,310343462,585475982,  618436,69,00,AEL Tm 1
     2 156718872,279244452,678196342,128316243,197823443, 1205000,69,01,AEL Tm 2
     3  93517192,364666132,103414613,192624193,293334613, 2368000,69,02/AEL Tm 3
      DATA NNN_Yb/
C     1 100010011,101310651,118613951,169120661,250629971,  625394,70,00,AEL Yb 1
C Yb I corrected following MOOG (Sneden)
     1 104410001,100011021,142920191,299545391, 68910342,  625394,70,00,Sne Yb 1
     2 200120901,270345231, 81714042,223533112,461959862, 1218400,70,01,AEL Yb 2
     3 100312561,250851931, 91914182,198626022,323638692, 2505000,70,02/AEL Yb 3
      DATA NNN_Lu/
     1 514664441,759086851, 99211442,133315612,182721252,  542589,71,00,AEL Lu 1
     2 125924831,438667801, 98714112,199727872,380850742, 1389900,71,01,AEL Lu 2
C     2 112718911,335853801,742987841,895879721,626944081, 1389900,71,01,Sne Lu 2
     3 323948621,661297271,158626482,426865032, 93712843, 2095960,71,02/AEL Lu 3
      DATA NNN_Hf/
     1 659294081,128016962,222528952,372047062,585171462,  700000,72,00,AEL Hf 1
     2  99117882,274638812,520867322, 84410313,123314453, 1489900,72,01,AEL Hf 2
     3 187427702,343739872,448049452,539358282,625266642, 2329900,72,02/AEL Hf 3
      DATA NNN_Ta/
     1  65210892,171325762,373552252,705192012,116414343,  787900,73,00,AEL Ta 1
     2 192837842,600784802,111113823,165419233,218524383, 1620000,73,01,AEL Ta 2
     3  99117872,274638812,520867312, 84410313,123314453, 2400000,73,02/FAK Ta 3
      DATA NNN_W/
     1 398981651,130019172,273438022,516168382, 88411163,  797900,74,00,AEL W  1
     2 131429482,523279952,111414623,183422233,262130233, 1770000,74,01,AEL W  2
     3 192837842,600784792,111113823,165419233,218524383, 2500000,74,02/FAK W  3
      DATA NNN_Re/
     1 600963001, 75910412,150121572,301940972,539168952,  787000,75,00,AEL Re 1
     2  73710852,190731262,464964142, 83810503,127315053, 1660000,75,01,AEL Re 2
     3 131429482,523279952,111414623,183422233,262130233, 2600000,75,02/FAK Re 3
      DATA NNN_Os/
     1 110815502,216829732,398752322,672484682,104612673,  850000,76,00,AEL Os 1
     2 168225972,362046562,566766422,757484612, 93010103, 1700000,76,01,AEL Os 2
     3  73710852,190731262,464964142, 83810503,127315053, 2700000,76,02/FAK Os 3
      DATA NNN_Ir/
     1 129117892,239430882,388748292,596173252, 89510843,  910000,77,00,AEL Ir 1
     2 110815502,216829732,398752322,672484682,104612673, 2000000,77,01,FAK Ir 2
     3 168225972,362046562,566766422,757484612, 93010103, 2800000,77,02/FAK Ir 3
      DATA NNN_Pt/
     1 158918512,207523002,254328242,316335762,407246582,  900000,78,00,AEL Pt 1
     2  98115462,224930742,401150612,623475412, 89910583, 1855900,78,01,AEL Pt 2
     3 110815502,216829732,398752322,672484682,104612673, 2900000,78,02/FAK Pt 3
      DATA NNN_Au/
     1 203222611,265731251,364042301,494958601,702084731,  922000,79,00,AEL Au 1
     2 120521331,357753801, 75310062,130516572,206925452, 2050000,79,01,AEL Au 2
     3 651780821,108814772,195925252,316338622,460853882, 3000000,79,02/AEL Au 3
      DATA NNN_Hg/
     1 100010001,100110111,105211851,152122101,341552811, 1043002,80,00,D+F Hg 1
     2 200320211,210023021,268834231,480472341,111416912, 1875000,80,01,D+F Hg 2
     3 104012871,186129471,458664151, 82410072,119013732, 3420000,80,02/D+F Hg 3
      DATA NNN_Tl/
     1 200420711,222424271,265429161,325637371,442853911,  610500,81,00,AEL Tl 1
     2 100010021,101910801,121414641,189525811,358949721, 2041900,81,01,AEL Tl 2
     3 200020311,216624611,296337451,489064791, 85711212, 2979900,81,02/AEL Tl 3
      DATA NNN_Pb/
     1 103411711,147819101,244331781,434862751, 93113762,  741404,82,00,D+F Pb 1
     2 204122231,248227841,311535621,429153941,651976431, 1502800,82,01,D+F Pb 2
     3 100210131,106812201,154522671,381665951, 95512512, 3192900,82,02/D+F Pb 3
      DATA NNN_Bi/
     1 400140351,416944121,474851591,564362181,690477231,  728700,83,00,AEL Bi 1
     2 106814451,204427341,350744811,586879131,108314772, 1667900,83,01,AEL Bi 2
     3 205523051,264830231,345439921,469156001,675281671, 2555900,83,02/AEL Bi 3
      DATA NNN_Po/
     1 500950661,518153561,559058941,628968071,748483501,  843000,84,00,AEL Po 1
     2 443756241,696282451, 95411012,128615262,182922012, 1900000,84,01,FAK Po 2
     3 336953201,682481011, 93810882,127915272,184622442, 2700000,84,02/FAK Po 3
      DATA NNN_At/
     1 402841621,431544771,463148311,520059491,734896851,  930000,85,00,FAK At 1
     2 576168741,788387631, 96910642,116012552,135014462, 2000000,85,01,FAK At 2
     3 490265341,812797201,116614322,179622692,285035302, 2900000,85,02/FAK At 3
      DATA NNN_Rn/
     1 100010001,100010031,102311051,133018071,264539391, 1074500,86,00,AEL Rn 1
     2 402841621,431544771,463148311,520059491,734996851, 2000000,86,01,FAK Rn 2
     3 576168741,788387631, 96910642,116012552,135014462, 3000000,86,02/FAK Rn 3
      DATA NNN_Fr/
     1 200020011,201220591,218124481,296538611,488859141,  400000,87,00,FAK Fr 1
     2 100010001,100010031,102311051,133018071,264539401, 2200000,87,01,FAK Fr 2
     3 421645151,477449611,511852711,542455761,572958821, 3300000,87,02/FAK Fr 3
      DATA NNN_Ra/
     1 100010041,105212131,153220271,270435641,460258111,  527600,88,00,AEL Ra 1
     2 201221791,258131471,381645781,546365131,777592781, 1014400,88,01,AEL Ra 2
     3 100010001,100010031,102311051,133018071,264539391, 3400000,88,02/FAK Ra 3
      DATA NNN_Ac/
     1 510064491, 82710872,142718412,232328712,348341572,  690000,89,00,AEL Ac 1
     2 228951571, 88513232,183324132,305537492,448152402, 1210000,89,01,AEL Ac 2
     3 723989131,103511752,130814352,155416652,177018682, 2000000,89,02/AEL Ac 3
      DATA NNN_Th/
C     1 620099241,162725772,391457072, 80110833,141818023,  600000,90,00,AEL Th 1
C     2 620099241,162725772,391457072, 80110833,141818023, 1200000,90,01,FAK Th 2
C     3 620099251,162725772,391457072, 80110833,141818023, 2000000,90,02/FAK Th 3
C Th I, II corrected following MOOG (Sneden)
     1  63810522,177929162,457168312, 97513353,175722323,  630670,90,00,Sne Th 1
     2 167341852, 79712803,185025013,323740733,503061383, 1190000,90,01,Sne Th 2
c  Th III partition function coefficients (Irwin 1981)
     3 583391701,143022052,322043402,533459552,604355902, 1830000,90,02/Sne Th 3
      DATA NNN_Pa/
     1 347877992,129318323,240730533,380546863,570368573,  600000,91,00,AEL Pa 1
     2 347877992,129318323,240730533,380546863,570368573, 1200000,91,01,FAK Pa 2
     3 347777992,129318323,240730533,380546863,570368573, 2000000,91,02/FAK Pa 3
      DATA NNN_U/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,92,00,AEL U  1
C  U II corrected following MOOG (Sneden)
     2  51311613,230239873,615986563,112513714,158317444, 1060000,92,01,Sne U  2
c  U III partition function coefficients (Irwin 1981)
     3 211130612,456267402, 94912483,151817063,177417123, 2000000,92,02/Sne U  3
C     2 209530092,450866762, 96613623,186524763,318839893, 1200000,92,01,FAK U  2
C     3 209530092,450866762, 96613623,186524763,318839893, 2000000,92,02/FAK U  3
      DATA NNN_Np/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,93,00,FAK Np 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,93,01,FAK Np 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,93,02/FAK Np 3
      DATA NNN_Pu/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,94,00,FAK Pu 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,94,01,FAK Pu 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,94,02/FAK Pu 3
      DATA NNN_Am/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,95,00,FAK Am 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,95,01,FAK Am 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,95,02/FAK Am 3
      DATA NNN_Cm/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,96,00,FAK Cm 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,96,01,FAK Cm 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,96,02/FAK Cm 3
      DATA NNN_Bk/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,97,00,FAK Bk 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,97,01,FAK Bk 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,97,02/FAK Bk 3
      DATA NNN_Cf/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,98,00,FAK Cf 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,98,01,FAK Cf 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,98,02/FAK Cf 3
      DATA NNN_Es/
     1 209530092,450866762, 96613623,186524763,318839893,  600000,99,00,FAK Es 1
     2 209530092,450866762, 96613623,186524763,318839893, 1200000,99,01,FAK Es 2
     3 209530092,450866762, 96613623,186524763,318839893, 2000000,99,02/FAK Es 3
      DATA SCALE/0.001,0.01,0.1,1.0/,FIRST/.TRUE./
C
C  First time XSAHA is called find the starting locations for each element
C
      IF(FIRST) THEN
        FIRST=.FALSE.
        IZ=0
        DO 1 N=1,NTABLE
        IF(NNNPFN(7,N).NE.IZ.AND.IZ.LE.ELESIZ) THEN
          IZ=NNNPFN(7,N)
          LOCZ(IZ)=N
        ENDIF
   1    CONTINUE
        LOCZ(IZ+1)=NTABLE+1
      ENDIF
C
C  Find starting row in the partition table and the number of ionization
C  stages available for a given element IEL
C
      N=LOCZ(IEL)
      NIONS=LOCZ(IEL+1)-N
C
C  For MODE=5 return the number of ionizations available for IEL
C
      IF(MODE.EQ.5) THEN
        MAXION=NIONS
        RETURN
      ENDIF
C
C  Compute T and kT in eV
C
      TTKEV=8.6171E-5*TT
      TV=TTKEV
      TTK=1.38065E-16*TT
C
C  Lowering of the ionization potential in Volts for unit Zeff
C
      CHARGE=2.*XNELEC
      EXCESS=XNELEC-XNATOM
C
C  Special allowance for doubly ionized Helium
C
      IF(EXCESS.GT.0.) CHARGE=CHARGE-EXCESS+4.*(2.*EXCESS)
      DEBYE=SQRT(TTK/(2.8965E-18*CHARGE))
      POTLOW=MIN(1.,1.44E-7/DEBYE)
C
C  Solve the Saha equation
C
      NION2=NIONS
      N=N-1
      DO 2 IONN=1,NION2
      Z=IONN
      POTLO(IONN)=POTLOW*Z
C      write(*,*) IP(IONN)-POTLO(IONN)
      N=N+1
      NNN100=NNNPFN(6,N)/100
      IP(IONN)=FLOAT(NNN100)/1000.
      G=NNNPFN(6,N)-NNN100*100
      IF(N.EQ.1) THEN
        PART(1)=2.
c        IF(TT.LT.9000.) GO TO 2
        PART(1)=PART(1)+8.*EXP(-10.196/TV)+18.*EXP(-12.084/TV)+32.*
     *          EXP(-12.745/TV)+50.*EXP(-13.051/TV)+72.*EXP(-13.217/TV)
        D1=13.595/6.5/6.5/TV
        D2=POTLO(1)/TV
      ELSE
        T2000=IP(IONN)*2000./11.
        IT=MAX(1,MIN(9,INT(TT/T2000-.5)))
        DT=TT/T2000-FLOAT(IT)-.5
        PMIN=1.
        I=(IT+1)/2
        K1=NNNPFN(I,N)/100000
        K2=NNNPFN(I,N)-K1*100000
        K3=K2/10
        KSCALE=K2-K3*10
        IF(MOD(IT,2).EQ.0) THEN
          P1=K3*SCALE(KSCALE)
          K1=NNNPFN(I+1,N)/100000
          KSCALE=MOD(NNNPFN(I+1,N),10)
          P2=K1*SCALE(KSCALE)
        ELSE
          P1=K1*SCALE(KSCALE)
          P2=K3*SCALE(KSCALE)
          IF(DT.LT.0.AND.KSCALE.LE.1) KP1=P1
          IF(DT.LT.0.AND.KSCALE.LE.1.AND.KP1.EQ.INT(P2+.5)) PMIN=KP1
        END IF
        PART(IONN)=MAX(PMIN,P1+(P2-P1)*DT)
        IF(G.EQ.0.0.OR.POTLO(IONN).LT.0.1.OR.TT.LT.T2000*4.0) GO TO 2
        IF(TT.GT.(T2000*11.)) TV=(T2000*11.)*8.6171E-5
        D1=.1/TV
      END IF
      D2=POTLO(IONN)/TV
      PART(IONN)=PART(IONN)+G*EXP(-IP(IONN)/TV)*
     *           (SQRT(13.595*Z*Z/TV/D2)**3*
     *           (1./3.+(1.-(.5+(1./18.+D2/120.)*D2)*D2)*D2)-
     -           SQRT(13.595*Z*Z/TV/D1)**3*
     *           (1./3.+(1.-(.5+(1./18.+D1/120.)*D1)*D1)*D1))
      TV=TTKEV
   2  CONTINUE
C
      N=N-NION2
      CF=2.*2.4148D15*TT*SQRT(TT)/XNELEC
      DO 3 IONN=2,NION2
      N=N+1
C
C  IF is to avoid annoying floating point underflows
C
      FEXARG=(IP(IONN-1)-POTLO(IONN-1))/TV
      IF(FEXARG.GT.80.) THEN
        F(IONN)=0.
      ELSE
        F(IONN)=CF*PART(IONN)/PART(IONN-1)*EXP(-FEXARG)
      END IF
   3  CONTINUE
      F(1)=1.
      DO 4 IONN=NION2,2,-1
   4  F(1)=1.+F(IONN)*F(1)
      F(1)=1./F(1)
      DO 5 IONN=2,NION2
   5  F(IONN)=F(IONN-1)*F(IONN)
C
C  Clear ionization stages for which infromation is not available
C
      DO 6 IONN=1,MAXION
      IF(MODE.EQ.3) THEN
        FRCT(IONN)=0.
      ELSE
        FRCT(IONN)=1.
      END IF
   6  CONTINUE
C
C  Formulate the answer according to MODE
C
      NIONS=MIN(MAXION,NION2)
      IF(MODE.EQ.1) THEN
        FRCT(1)=F(1)/PART(1)
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 7 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
   7      FRCT(IONN)=F(IONN)/PART(IONN)
        END IF
      ELSE IF(MODE.EQ.2) THEN
        FRCT(1)=F(1)
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 8 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
   8      FRCT(IONN)=F(IONN)
        END IF
      ELSE IF(MODE.EQ.3) THEN
        FRCT(1)=PART(1)
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 9 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
   9      FRCT(IONN)=PART(IONN)
        END IF
      ELSE IF(MODE.EQ.4) THEN
        FRCT(1)=0
        POTI(1)=IP(1)
        IF(NIONS.GT.1) THEN
          DO 10 IONN=2,NIONS
          POTI(IONN)=IP(IONN)
  10      FRCT(1)=FRCT(1)+F(IONN)*(IONN-1)
        END IF
      END IF
C
      RETURN
      END
