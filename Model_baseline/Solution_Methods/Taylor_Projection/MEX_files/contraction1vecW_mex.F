C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C     Gateway subroutine
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)

C     Declarations
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
      mwPointer pr_in1,pr_in2
      mwPointer pr_in3,pr_in4,pr_in5,pr_in6
      mwPointer pr_in7,pr_in8,pr_in9,pr_in10
      mwPointer pr_in11, dosum

      mwPointer pr_out

C     Array information:
      mwSize IArows,JArows
      mwSize IB1rows, JB1rows
      mwSize ICrows_1,JCrows_1

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
      JB1rows = mxGetM(prhs(7))

      ICrows_1 = mxGetM(prhs(9))
      JCrows_1 = mxGetM(prhs(10))



C     Create Fortran array from the input arguments.
      pr_in1 = mxGetPr(prhs(1))
      pr_in2 = mxGetPr(prhs(2))
      pr_in3 = mxGetPr(prhs(3))
      pr_in4 = mxGetPr(prhs(4))
      pr_in5 = mxGetPr(prhs(5))
      pr_in6 = mxGetPr(prhs(6))
      pr_in7 = mxGetPr(prhs(7))
      pr_in8 = mxGetPr(prhs(8))
      pr_in9 = mxGetPr(prhs(9))
      pr_in10 = mxGetPr(prhs(10))
      pr_in11 = mxGetPr(prhs(11))

C     Get number of states
      sA = mxGetM(prhs(5))
      sB = 1 ! B1 is a W matrix

C     sum across rows of C if dosum==1
      do_sum=mxGetScalar(prhs(12))
 
      if (do_sum==1) then
       sC=1
      else
       sC=max(sA,sB)
      endif
C     Create matrix for the return argument.
      complexflag = 0
      plhs(1) = mxCreateDoubleMatrix(sC,JCrows_1, complexflag)
      pr_out = mxGetPr(plhs(1))

C     Call the computational routine.
      call compute(%VAL(pr_in1), %VAL(pr_in2),
     +             %VAL(pr_in3), %VAL(pr_in4),
     +             %VAL(pr_in5), %VAL(pr_in6),
     +             %VAL(pr_in7), %VAL(pr_in8),
     +             %VAL(pr_in9), %VAL(pr_in10),
     +             %VAL(pr_in11),
     +             IArows,JArows,
     +             IB1rows,JB1rows,
     +             ICrows_1,JCrows_1,
     +             %VAL(pr_out),sC,sA,sB)

      return
      end

C     Computational subroutine
      subroutine compute(l,n1,
     +              IA,JA,NA,
     +              IB1,JB1,NB1,
     +              IC_1,JC1_1,
     +              maxload,
     +              IArows,JArows,
     +              IB1rows,JB1rows,
     +              ICrows_1,JCrows_1,
     +              NC,sC,sA,sB)

C     Declarations
      implicit none

      mwSize  l,n1,s,sC,sA,sB,maxload,n_vecload,loadi,load_
      mwSize  IArows,JArows
      mwSize  IB1rows,JB1rows
      mwSize  ICrows_1,JCrows_1

      real*8, target ::  NA(sA,JArows)
      real*8, target ::  NB1(sB,JB1rows)
      real*8  NC(sC,JCrows_1)
 
      real*8, pointer :: pNA(:,:)
      mwSize, pointer :: pIA(:)
      mwSize, pointer :: pJA(:)
      mwSize, pointer :: IC(:)
      mwSize, pointer :: JC1(:)
      real*8, pointer :: NB(:,:)
      mwSize, pointer :: IB(:)
 
      mwSize, allocatable :: vecload(:)

      real*8, allocatable :: x1(:,:)
     
      mwSize, target ::  IA(IArows)
      mwSize, target ::  JA(JArows)

      mwSize, target ::  IB1(IB1rows)
      mwSize  JB1(JB1rows)

      mwSize, target ::  IC_1(ICrows_1)
      mwSize, target ::  JC1_1(JCrows_1)

      mwSize i,jp,j,kp,k1,vp,v1,is,js,nz


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
       allocate (x1(maxload,n1),
     +   STAT = AllocateStatus)
      else
       allocate (x1(1,n1),
     +   STAT = AllocateStatus)
      endif

      IF (AllocateStatus .ne. 0) then
          call mexErrMsgIdAndTxt ('MATLAB:compress:OutofMemory',
     +                           'Out of memory.')
      endif

C     x1=0.0

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

        IC=>IC_1
        JC1=>JC1_1
        
        IB=>IB1        
      nz=1  
      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        j=pJA(jp)
      do kp=IB(j),IB(j+1)-1
C       k1=JB1(kp)
      if (sC>1) then
C Choose among three options and do not sum across rows: 1. A and B are both vectorized. 2. Only A is vectorized. 3. only B is vectorized.
      if (sA.ge.sB) then
      !DEC$ simd
       do is=1,load_
         NC(is,nz)=pNA(is,jp)
       end do
      else 
      !DEC$ simd
       do is=1,load_
         NC(is,nz)=pNA(1,jp)
       end do
      endif
      else
C Choose among three options and sum across rows.
      if (sA.ge.sB) then
      !DEC$ simd
       do is=1,load_
         NC(1,nz)=NC(1,nz)+pNA(is,jp)
       end do
      else 
      !DEC$ simd
       do is=1,load_
         NC(1,nz)=NC(1,nz)+pNA(1,jp)
       end do
      endif
      endif
      nz=nz+1
      end do
      end do
#if 0
C     assign the result in NC
      do vp=IC(i),IC(i+1)-1
         v1=JC1(vp)
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
#endif

      end do
      end do
      deallocate (vecload)
      deallocate (x1)

      return
      end

