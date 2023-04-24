C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [IC_1] = contraction1IC2D_mex(l,n1,IA,JA,IB1,JB1)
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
      mwSize IB1rows,IB1cols,JB1rows,ndim,dims(2)

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
      IB1cols = mxGetN(prhs(5))
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
      ndim=2
      dims(1)=IArows-1
      dims(2)=IB1cols
      plhs(1) = mxCreateNumericArray(ndim, dims, 
     +          classid, 
     +          complexflag)

      IC_1 = mxGetPr(plhs(1))


C     Call the computational routine.
      call compute(%VAL(l),%VAL(n1),
     + %VAL(IA),%VAL(JA),%VAL(IB1),%VAL(JB1),
     + %VAL(IC_1),            
     + IArows,JArows,IB1rows,IB1cols,JB1rows)

      return
      end

C     Computational subroutine
      subroutine compute(l,n1,IA,JA,IB1,JB1,
     + IC_1,         
     + IArows,JArows,IB1rows,IB1cols,JB1rows)

C     Declarations
      implicit none

      mwSize  l,n1
      mwSize  IArows,JArows
      mwSize  IB1rows,IB1cols,JB1rows

      mwSize, pointer :: pIA(:)
      mwSize, pointer :: pJA(:)
      mwSize, pointer :: pIB(:,:)
 
      mwSize, allocatable :: w(:)
     
      mwSize, target :: IA(IArows)
      mwSize, target :: JA(JArows)

      mwSize, target :: IB1(IB1rows,IB1cols)
      mwSize  JB1(JB1rows)

      mwSize, target :: IC_1(IArows-1,IB1cols)

      mwSize i,jp,j,kp,k1,vp,v1,nz,i2

      integer AllocateStatus

      allocate (w(n1),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      nz=1

      pIA=>IA
      pJA=>JA

      pIB=>IB1        
       
      do i=1,l
        w=0
      do i2=1,IB1cols-1  
        IC_1(i,i2)=nz
      do jp=pIA(i),pIA(i+1)-1
        j=pJA(jp)
      do kp=pIB(j,i2),pIB(j,i2+1)-1
        k1=JB1(kp)
        if (w(k1)<i2) then
         w(k1)=i2
         nz=nz+1
        endif
      end do
      end do
      end do
        IC_1(i,IB1cols)=nz
      end do

      
      deallocate (w)

      return
      end

