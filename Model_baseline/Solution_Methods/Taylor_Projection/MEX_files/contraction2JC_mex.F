C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [JC_1,JC_2] = contraction2JC_mex(l,m1,n1,n2,IA,JA,IB1,JB1,IB2,JB2,JCrows_1,JCrows_2,IC_2)
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
      real*8 mxGetScalar

C     Arguments for mxCreateNumericArray
      integer*4 classid
      integer*4 complexflag

C     Pointers to input/output mxArrays:
      mwPointer l,m1,n1,n2,IA,JA
      mwPointer IB1,JB1,IB2,JB2
      mwPointer IC_2

      mwPointer JCrows_1,JCrows_2

      mwPointer JC_1,JC_2

C     Array information:
      mwSize IArows,JArows
      mwSize IB1rows, JB1rows
      mwSize IB2rows, JB2rows

      mwSize rows1, rows2
      mwSize JC_1dims(2)
      mwSize JC_2dims(2)

      mwSize ICrows_2,ndim

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 13) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '13 inputs required.')
      endif
      if(nlhs .gt. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 2 outputs required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(5))
      JArows = mxGetM(prhs(6))

      IB1rows = mxGetM(prhs(7))
      JB1rows = mxGetM(prhs(8))
      IB2rows = mxGetM(prhs(9))
      JB2rows = mxGetM(prhs(10))

      ICrows_2 = mxGetM(prhs(13))

C     Create Fortran array from the input arguments.
      l = mxGetPr(prhs(1))
      m1 = mxGetPr(prhs(2))
      n1 = mxGetPr(prhs(3))
      n2 = mxGetPr(prhs(4))
      IA = mxGetPr(prhs(5))
      JA = mxGetPr(prhs(6))
      IB1 = mxGetPr(prhs(7))
      JB1 = mxGetPr(prhs(8))
      IB2 = mxGetPr(prhs(9))
      JB2 = mxGetPr(prhs(10))
      JCrows_1 = mxGetPr(prhs(11))
      JCrows_2 = mxGetPr(prhs(12))
      IC_2 = mxGetPr(prhs(13))

C     Create matrix for the return arguments.
      rows1=mxGetScalar(prhs(11))
      rows2=mxGetScalar(prhs(12))

      classid = mxGetClassID(prhs(6))

      complexflag = 0
      ndim=2
      JC_1dims(1)=rows1
      JC_1dims(2)=2
      plhs(1) = mxCreateNumericArray(ndim, JC_1dims, 
     +          classid, 
     +          complexflag)

      JC_1 = mxGetPr(plhs(1))

      JC_2dims(1)=rows2
      JC_2dims(2)=2
      plhs(2) = mxCreateNumericArray(ndim, JC_2dims, 
     +          classid, 
     +          complexflag)

      JC_2 = mxGetPr(plhs(2))

C     Call the computational routine.
      call compute(%VAL(l),%VAL(m1),%VAL(n1),%VAL(n2),
     + %VAL(IA),%VAL(JA),
     + %VAL(IB1),%VAL(JB1),%VAL(IB2),%VAL(JB2),
     + %VAL(JCrows_1),%VAL(JCrows_2),            
     + %VAL(IC_2),
     + %VAL(JC_1),%VAL(JC_2),
     + IArows,JArows,IB1rows,JB1rows,IB2rows,JB2rows,
     + ICrows_2)

      return
      end

C     Computational subroutine
      subroutine compute(l,m1,n1,n2,IA,JA,
     + IB1,JB1,IB2,JB2,
     + JCrows_1,JCrows_2,
     + IC_2,
     + JC_1,JC_2,
     + IArows,JArows,
     + IB1rows,JB1rows,IB2rows,JB2rows,
     + ICrows_2)

C     Declarations
      implicit none

      mwSize  l,m1,n1,n2
      mwSize  IArows,JArows
      mwSize  IB1rows,JB1rows,IB2rows,JB2rows
      mwSize  JCrows_1,JCrows_2
      mwSize  ICrows_2

      mwSize, pointer :: pIA(:)
      mwSize, pointer :: pJA(:,:)
      mwSize, pointer :: pIB(:)
 
      mwSize, allocatable :: w(:,:)
     
      mwSize, target :: IA(IArows)
      mwSize, target :: JA(JArows,2)

      mwSize, target ::  IB1(IB1rows)
      mwSize  JB1(JB1rows)
      mwSize, target ::  IB2(IB2rows)
      mwSize  JB2(JB2rows)


      mwSize, target :: JC_1(JCrows_1,2)
      mwSize, target :: JC_2(JCrows_2,2)

      mwSize, target :: IC_2(ICrows_2)

      mwSize i,jp,j,kp,k1,k2,vp,v1,v2,nz

      integer AllocateStatus

C first contraction

      allocate (w(m1,n2),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      w=0
      nz=1

      pIA=>IA
      pJA=>JA

      pIB=>IB2        
       
      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        j=pJA(jp,2)
        k1=pJA(jp,1)
      do kp=pIB(j),pIB(j+1)-1
        k2=JB2(kp)
        if (w(k1,k2)<i) then
         w(k1,k2)=i
         JC_2(nz,1)=k1
         JC_2(nz,2)=k2
         nz=nz+1
        endif
      end do
      end do
      end do

      deallocate (w)


C second contraction
      allocate (w(n1,n2),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      w=0
      nz=1

      pIA=>IC_2
      pJA=>JC_2

      pIB=>IB1

      do i=1,l

      do jp=pIA(i),pIA(i+1)-1
        j=pJA(jp,1)
        k2=pJA(jp,2)
      do kp=pIB(j),pIB(j+1)-1
        k1=JB1(kp)
        if (w(k1,k2)<i) then
         w(k1,k2)=i
         JC_1(nz,1)=k1
         JC_1(nz,2)=k2
         nz=nz+1
        endif
      end do
      end do

      end do

      deallocate (w)

      return
      end

