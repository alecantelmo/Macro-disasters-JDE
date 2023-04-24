C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [NC] = contraction1IC_mex(l,n1,IA,JA,NA,IB1,JB1,NB1,IC_1,JC_1,maxload,dosum)
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetM, mxGetN
      real*8 mxGetScalar
      integer*4 complexflag

C     Pointers to input/output mxArrays:
      mwPointer l,n1,IA,JA,NA,IB1,JB1,NB1
      mwPointer IC_1,JC_1
      mwPointer maxload, dosum

      mwPointer NC

C     Array information:
      mwSize IArows,JArows
      mwSize IB1rows,IB1cols,JB1rows
      mwSize ICrows_1,ICcols_1,JCrows_1

      mwSize sC,sA,sB
      real*8 do_sum

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 12) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '12 inputs required.')
      endif
      if(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'One output required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(3))
      JArows = mxGetM(prhs(4))

      IB1rows = mxGetM(prhs(6))
      IB1cols = mxGetN(prhs(6))
      JB1rows = mxGetM(prhs(7))

      ICrows_1 = mxGetM(prhs(9))
      ICcols_1 = mxGetN(prhs(9))
      JCrows_1 = mxGetM(prhs(10))

C     Create Fortran array from the input arguments.
      l = mxGetPr(prhs(1))
      n1 = mxGetPr(prhs(2))
      IA = mxGetPr(prhs(3))
      JA = mxGetPr(prhs(4))
      NA = mxGetPr(prhs(5))
      IB1 = mxGetPr(prhs(6))
      JB1 = mxGetPr(prhs(7))
      NB1 = mxGetPr(prhs(8))
      IC_1 = mxGetPr(prhs(9))
      JC_1 = mxGetPr(prhs(10))
      maxload = mxGetPr(prhs(11))

C     Get number of states
      sA = mxGetM(prhs(5))
      sB = mxGetM(prhs(8))

C     sum across rows of C if dosum==1
      do_sum=mxGetScalar(prhs(12))
 
      if (do_sum==1) then
       sC=1
      else
       sC=max(sA,sB)
      endif
C     Create matrix for the return argument.
      complexflag=0
      plhs(1) = mxCreateDoubleMatrix(sC,JCrows_1, complexflag)
      NC = mxGetPr(plhs(1))

C     Call the computational routine.
      call compute(%VAL(l), %VAL(n1),
     +             %VAL(IA), %VAL(JA),%VAL(NA),
     +             %VAL(IB1), %VAL(JB1),%VAL(NB1),
     +             %VAL(IC_1), %VAL(JC_1),
     +             %VAL(maxload),
     +             IArows,JArows,
     +             IB1rows,IB1cols,JB1rows,
     +             ICrows_1,ICcols_1,JCrows_1,
     +             %VAL(NC),sC,sA,sB)

      return
      end

C     Computational subroutine
      subroutine compute(l,n1,
     +              IA,JA,NA,
     +              IB1,JB1,NB1,
     +              IC_1,JC_1,
     +              maxload,
     +              IArows,JArows,
     +              IB1rows,IB1cols,JB1rows,
     +              ICrows_1,ICcols_1,JCrows_1,
     +              NC,sC,sA,sB)

C     Declarations
      implicit none

      mwSize  l,n1,s,sC,sA,sB,maxload,n_vecload,loadi,load_
      mwSize  IArows,JArows
      mwSize  IB1rows,IB1cols,JB1rows
      mwSize  ICrows_1,ICcols_1,JCrows_1

      real*8, target ::  NA(sA,JArows)
      real*8, target ::  NB1(sB,JB1rows)
      real*8  NC(sC,JCrows_1)
 
      real*8, pointer :: pNA(:,:)
      mwSize, pointer :: pIA(:)
      mwSize, pointer :: pJA(:)
      mwSize, pointer :: pIC(:,:)
      mwSize, pointer :: pJC(:)
      real*8, pointer :: pNB(:,:)
      mwSize, pointer :: pIB(:,:)
 
      mwSize, allocatable :: vecload(:)

      real*8, allocatable :: x1(:,:)
     
      mwSize, target ::  IA(IArows)
      mwSize, target ::  JA(JArows)

      mwSize, target ::  IB1(IB1rows,IB1cols)
      mwSize  JB1(JB1rows)

      mwSize, target ::  IC_1(ICrows_1,ICcols_1)
      mwSize, target ::  JC_1(JCrows_1)

      mwSize i,jp,j,kp,k1,vp,v1,is,js,i2


      mwSize jpstart, jpend
      mwSize kpstart, kpend
      mwSize vpstart, vpend

      integer AllocateStatus

      s=max(sA,sB)

      if ((sA>sB .and. sB.ne.1) .or. (sA<sB .and. sA.ne.1)) then
         call mexErrMsgIdAndTxt ('MATLAB:states',
     +   'A.vals and B.vals must have the same number of rows.')
      endif

      if (sC>1) then
       allocate (x1(min(sC,maxload),n1),
     +   STAT = AllocateStatus)
      else
       allocate (x1(1,n1),
     +   STAT = AllocateStatus)
      endif

      IF (AllocateStatus .ne. 0) then
          call mexErrMsgIdAndTxt ('MATLAB:compress:OutofMemory',
     +                           'Out of memory.')
      endif

      x1=0.0

C     NC=0.0, mxCreateDoubleMatrix already initializes to zero

      n_vecload=(s-1)/maxload+2

      allocate (vecload(n_vecload),STAT = AllocateStatus)
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:OutofMemory',
     +                           'Out of memory.')
      endif

      vecload(1)=1
      vecload(n_vecload)=s+1
      do i=2,n_vecload-1
        vecload(i)=vecload(i-1)+maxload
      end do

      do loadi=1,n_vecload-1
 
        load_=vecload(loadi+1)-vecload(loadi)
      if (sA>1) then
        pNA=>NA(vecload(loadi):vecload(loadi+1)-1,:)
      else 
        pNA=>NA
      endif
        pIA=>IA
        pJA=>JA

        pIC=>IC_1
        pJC=>JC_1
        
      if (sB>1) then
        pNB=>NB1(vecload(loadi):vecload(loadi+1)-1,:)       
      else 
        pNB=>NB1
      endif
        pIB=>IB1        
        
      do i=1,l
      do i2=1,IB1cols-1  
      do jp=pIA(i),pIA(i+1)-1
        j=pJA(jp)
      do kp=pIB(j,i2),pIB(j,i2+1)-1
        k1=JB1(kp)
      if (sC>1) then
C Choose among three options and do not sum across rows: 1. A and B are both vectorized. 2. Only A is vectorized. 3. only B is vectorized.
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
         x1(is,k1)=x1(is,k1)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
         x1(is,k1)=x1(is,k1)+pNA(is,jp)*pNB(1,kp)
       end do
      else 
      !DEC$ simd
       do is=1,load_
         x1(is,k1)=x1(is,k1)+pNA(1,jp)*pNB(is,kp)
       end do
      endif
      else
C Choose among three options and sum across rows.
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
         x1(1,k1)=x1(1,k1)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
         x1(1,k1)=x1(1,k1)+pNA(is,jp)*pNB(1,kp)
       end do
      else 
      !DEC$ simd
       do is=1,load_
         x1(1,k1)=x1(1,k1)+pNA(1,jp)*pNB(is,kp)
       end do
      endif
      endif
      end do
      end do

C     assign the result in NC
      do vp=pIC(i,i2),pIC(i,i2+1)-1
         v1=pJC(vp)
      if (sC>1) then
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         NC(js,vp)=x1(is,v1)
         x1(is,v1)=0.0
      end do
      else
         NC(1,vp)=NC(1,vp)+x1(1,v1)
         x1(1,v1)=0.0
      endif

      end do
      end do
      end do
      end do
      deallocate (vecload)
      deallocate (x1)

      return
      end

