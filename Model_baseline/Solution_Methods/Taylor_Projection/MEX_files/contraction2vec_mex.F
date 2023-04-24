C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [NC] = contraction2vec_mex(l,m1,n1,n2,IA,JA,NA,...
C            IB1,JB1,NB1,IB2,JB2,NB2,...
C            IC_1,JC_1,IC_2,JC_2,...
C            maxload,dosum)
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxGetM, mxGetN
      mwPointer mxCreateDoubleMatrix
      integer*4 complexflag
      real*8 mxGetScalar

C     Pointers to input/output mxArrays:
      mwPointer l,IA,JA
      mwPointer m1,n1,n2
      mwPointer IB1,JB1,IB2,JB2
      mwPointer IC_1,IC_2
      mwPointer JC_1,JC_2
      mwPointer maxload,dosum
      mwPointer NA,NB1,NB2

      mwPointer NC

C     Array information:
      mwSize IArows,JArows
      mwSize IB1rows, JB1rows
      mwSize IB2rows, JB2rows
      mwSize ICrows_1,ICrows_2
      mwSize JCrows_1,JCrows_2

      mwSize s,sC,sA,sB1,sB2
      real*8 do_sum

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 19) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '19 inputs required.')
      endif
      if(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 1 outputs required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(5))
      JArows = mxGetM(prhs(6))

      IB1rows = mxGetM(prhs(8))
      JB1rows = mxGetM(prhs(9))
      IB2rows = mxGetM(prhs(11))
      JB2rows = mxGetM(prhs(12))

      ICrows_1 = mxGetM(prhs(14))
      JCrows_1 = mxGetM(prhs(15))

      ICrows_2 = mxGetM(prhs(16))
      JCrows_2 = mxGetM(prhs(17))

C     Create Fortran array from the input arguments.
      l = mxGetPr(prhs(1))
      m1 = mxGetPr(prhs(2))
      n1 = mxGetPr(prhs(3))
      n2 = mxGetPr(prhs(4))
      IA = mxGetPr(prhs(5))
      JA = mxGetPr(prhs(6))
      NA = mxGetPr(prhs(7))
      IB1 = mxGetPr(prhs(8))
      JB1 = mxGetPr(prhs(9))
      NB1 = mxGetPr(prhs(10))
      IB2 = mxGetPr(prhs(11))
      JB2 = mxGetPr(prhs(12))
      NB2 = mxGetPr(prhs(13))
      IC_1 = mxGetPr(prhs(14))
      JC_1 = mxGetPr(prhs(15))
      IC_2 = mxGetPr(prhs(16))
      JC_2 = mxGetPr(prhs(17))
      maxload = mxGetPr(prhs(18))

C     Get number of states
      sA = mxGetM(prhs(7))
      sB1 = mxGetM(prhs(10))
      sB2 = mxGetM(prhs(13))

C check compatibility of states
      s=max(sA,sB1,sB2)
      if ((sA.ne.s .and. sA.ne.1) .or. (sB1.ne.s .and. sB1.ne.1)
     +   .or. (sB1.ne.sB2)) then
         call mexErrMsgIdAndTxt ('MATLAB:states',
     +   'incompatible states.')
      endif

C     sum across rows of C if dosum==1
      do_sum=mxGetScalar(prhs(19))
 
      if (do_sum==1) then
       sC=1
      else
       sC=s
      endif

C     Create matrix for NC.
      complexflag = 0
      plhs(1) = mxCreateDoubleMatrix(sC,JCrows_1, complexflag)
      NC = mxGetPr(plhs(1))

C     Call the computational routine.
      call compute(%VAL(l),%VAL(m1),
     + %VAL(n1),%VAL(n2),
     + %VAL(IA),%VAL(JA),%VAL(NA),
     + %VAL(IB1),%VAL(JB1),%VAL(NB1),
     + %VAL(IB2),%VAL(JB2),%VAL(NB2),
     + %VAL(IC_1),%VAL(JC_1),%VAL(IC_2),%VAL(JC_2),
     + %VAL(maxload),
     + %VAL(NC),            
     + IArows,JArows,IB1rows,JB1rows,IB2rows,JB2rows,
     + ICrows_1,JCrows_1,ICrows_2,JCrows_2,
     + sC,sA,sB1,sB2)
 
      return
      end

C     Computational subroutine
      subroutine compute(l,m1,n1,n2,IA,JA,NA,
     + IB1,JB1,NB1,IB2,JB2,NB2,
     + IC_1,JC_1,IC_2,JC_2,
     + maxload,
     + NC,
     + IArows,JArows,
     + IB1rows,JB1rows,IB2rows,JB2rows,
     + ICrows_1,JCrows_1,ICrows_2,JCrows_2,
     + sC,sA,sB1,sB2)

C     Declarations
      implicit none

      mwSize  l,m1,n1,n2,sC,sA,sB1,sB2
      mwSize  IArows,JArows
      mwSize  IB1rows,JB1rows,IB2rows,JB2rows
      mwSize  ICrows_1,ICrows_2
      mwSize  JCrows_1,JCrows_2

      mwSize  s,sB,sC2,rows
      mwSize  maxload,n_vecload,loadi,load_

      real*8, target ::  NA(sA,JArows)
      real*8, target ::  NB1(sB1,JB1rows)
      real*8, target ::  NB2(sB2,JB2rows)

      real*8  NC(sC,JCrows_1)

      real*8, pointer :: pNA(:,:),pNB(:,:)
      mwSize, pointer :: pIA(:),pIB(:),pIC(:)
      mwSize, pointer :: pJA(:,:),pJC(:,:)
      mwSize, pointer :: pJB(:)

       mwSize, allocatable :: vecload(:)

      real*8, allocatable :: x(:,:,:) 
      logical, allocatable :: w(:,:)

C one temporary crs var is needed, 2
      real*8, target, allocatable :: tempNC2(:,:)

      mwSize, target ::  IA(IArows)
      mwSize, target ::  JA(JArows,2)

      mwSize, target ::  IB1(IB1rows)
      mwSize, target ::  JB1(JB1rows)
      mwSize, target ::  IB2(IB2rows)
      mwSize, target ::  JB2(JB2rows)

      mwSize, target ::  IC_1(ICrows_1)
      mwSize, target ::  JC_1(JCrows_1,2)

      mwSize, target ::  IC_2(ICrows_2)
      mwSize, target ::  JC_2(JCrows_2,2)

      mwSize i,jp,j,kp,k1,k2,vp,v1,v2
      mwSize is,js,nz

      mwSize jpstart, jpend
      mwSize kpstart, kpend
      mwSize vpstart, vpend
 
      integer AllocateStatus

C vecload
      s=max(sA,sB1,sB2)

      n_vecload=(s-1)/maxload+2

      allocate (vecload(n_vecload),STAT = AllocateStatus)
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:OutofMemory',
     +                           'Out of memory1.')
      endif

      vecload(1)=1
      vecload(n_vecload)=s+1
      do i=2,n_vecload-1
        vecload(i)=vecload(i-1)+maxload
      end do

      sB=sB1 ! sB1 and sB2 must all be equal

C first contraction
      sC2=max(sA,sB2)
      rows=IC_2(ICrows_2)-1
      allocate (w(m1,n2),
     +   x(min(sC2,maxload),m1,n2),
     +   tempNC2(sC2,rows),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory2.')
      endif

      w=.false.
      nz=1

      pIA=>IA
      pJA=>JA

      pIB=>IB2        
      pJB=>JB2

      pIC=>IC_2
      pJC=>JC_2

      do loadi=1,n_vecload-1
        load_=vecload(loadi+1)-vecload(loadi)
       
      if (sA>1) then
        pNA=>NA(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNA=>NA
      endif

      if (sB2>1) then
        pNB=>NB2(vecload(loadi):vecload(loadi+1)-1,:)       
      else 
        pNB=>NB2
      endif

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        j=pJA(jp,2)
        k1=pJA(jp,1)
      do kp=pIB(j),pIB(j+1)-1
        k2=pJB(kp)
C initialize x to zero
        if (w(k1,k2)==.false.) then
         w(k1,k2)=.true.
         x(:,k1,k2)=0
        endif

C Choose among three options:
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(1,jp)*pNB(is,kp)
       end do
      endif

      end do
      end do

C     assign the result in tempNC2
      do vp=pIC(i),pIC(i+1)-1
         v1=pJC(vp,1)
         v2=pJC(vp,2)
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         tempNC2(js,vp)=x(is,v1,v2)
         x(is,v1,v2)=0.0
      end do
      end do

      end do

      end do

      deallocate (w,x)

C second contraction
      rows=IC_1(ICrows_1)-1
      allocate (w(n1,n2),
     +   x(min(sC,maxload),n1,n2),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory4.')
      endif

      w=.false.
      nz=1

      pIA=>IC_2
      pJA=>JC_2

      pIB=>IB1        
      pJB=>JB1

      pIC=>IC_1
      pJC=>JC_1

      sA=sC2

      do loadi=1,n_vecload-1
        load_=vecload(loadi+1)-vecload(loadi)
       
      if (sC2>1) then
        pNA=>tempNC2(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNA=>tempNC2
      endif

      if (sB1>1) then
        pNB=>NB1(vecload(loadi):vecload(loadi+1)-1,:)       
      else 
        pNB=>NB1
      endif

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        j=pJA(jp,1)
        k2=pJA(jp,2)
      do kp=pIB(j),pIB(j+1)-1
        k1=pJB(kp)
C initialize x to zero
        if (w(k1,k2)==.false.) then
         w(k1,k2)=.true.
         x(:,k1,k2)=0
        endif

      if (sC>1) then
C Choose among three options and do not sum across rows
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(1,jp)*pNB(is,kp)
       end do
      endif
      else
C Choose among three options and sum across rows.
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2)=x(1,k1,k2)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2)=x(1,k1,k2)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2)=x(1,k1,k2)+pNA(1,jp)*pNB(is,kp)
       end do
      endif
      endif

      end do
      end do

C     assign the result in NC
      do vp=pIC(i),pIC(i+1)-1
         v1=pJC(vp,1)
         v2=pJC(vp,2)
      if (sC>1) then
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         NC(js,vp)=x(is,v1,v2)
         x(is,v1,v2)=0.0
      end do
      else
         NC(1,vp)=NC(1,vp)+x(1,v1,v2)
         x(1,v1,v2)=0.0
      endif
      end do

      end do

      end do

      deallocate (w,x)
      deallocate (tempNC2)

      return
      end

