C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [IC_1] = contraction1IC_mex(l,n1,IA,JA,IB1,JB1)
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxGetM, mxGetN
      mwPointer mxCreateNumericArray
      integer*4 mxGetClassID

C     Arguments for mxCreateNumericArray
      integer*4 classid
      integer*4 complexflag

C     Pointers to input/output mxArrays:
      mwPointer l,n1,IA,JA,IB1,JB1
      mwPointer IC_1

C     Array information:
      mwSize IArows,JArows
      mwSize IB1rows, JB1rows, ndim

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 6) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '6 inputs required.')
      endif
      if(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           '1 output required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(3))
      JArows = mxGetM(prhs(4))

      IB1rows = mxGetM(prhs(5))
      JB1rows = mxGetM(prhs(6))

C     Create Fortran array from the input arguments.
      l = mxGetPr(prhs(1))
      n1 = mxGetPr(prhs(2))
      IA = mxGetPr(prhs(3))
      JA = mxGetPr(prhs(4))
      IB1 = mxGetPr(prhs(5))
      JB1 = mxGetPr(prhs(6))

C     Create matrix for the return arguments.

      classid = mxGetClassID(prhs(3))
      complexflag = 0
      ndim=1
      plhs(1) = mxCreateNumericArray(ndim, IArows, 
     +          classid, 
     +          complexflag)

      IC_1 = mxGetPr(plhs(1))


C     Call the computational routine.
      call compute(%VAL(l),%VAL(n1),
     + %VAL(IA),%VAL(JA),%VAL(IB1),%VAL(JB1),
     + %VAL(IC_1),            
     + IArows,JArows,IB1rows,JB1rows)

      return
      end

C     Computational subroutine
      subroutine compute(l,n1,IA,JA,IB1,JB1,
     + IC_1,         
     + IArows,JArows,IB1rows,JB1rows)

C     Declarations
      implicit none

      mwSize  l,n1
      mwSize  IArows,JArows
      mwSize  IB1rows,JB1rows

      mwSize, pointer :: pIA(:)
      mwSize, pointer :: pJA(:)
      mwSize, pointer :: pIB(:)
 
      mwSize, allocatable :: w(:)
     
      mwSize, target :: IA(IArows)
      mwSize, target :: JA(JArows)

      mwSize, target :: IB1(IB1rows)
      mwSize  JB1(JB1rows)

      mwSize, target :: IC_1(IArows)

      mwSize i,jp,j,kp,k1,vp,v1,nz

      integer AllocateStatus

      allocate (w(n1),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      w=0
      nz=1

      pIA=>IA
      pJA=>JA

      pIB=>IB1        
       
      do i=1,l
        IC_1(i)=nz
      do jp=pIA(i),pIA(i+1)-1
        j=pJA(jp)
      do kp=pIB(j),pIB(j+1)-1
        k1=JB1(kp)
        if (w(k1)<i) then
         w(k1)=i
         nz=nz+1
        endif
      end do
      end do
      end do

      IC_1(IArows)=nz
      
      deallocate (w)

      return
      end

