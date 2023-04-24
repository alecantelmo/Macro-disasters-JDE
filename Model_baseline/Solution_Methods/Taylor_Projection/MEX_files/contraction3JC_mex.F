C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [JC_1,JC_2,JC_3] = contraction3JC_mex(l,m1,m2,n1,n2,n3,IA,JA,...
C            IB1,JB1,IB2,JB2,IB3,JB3,...
C            JCrows_1,JCrows_2,JCrows_3,...
C            IC_2,IC_3)
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
      mwPointer l,IA,JA
      mwPointer m1,m2,n1,n2,n3
      mwPointer IB1,JB1,IB2,JB2,IB3,JB3
      mwPointer IC_2,IC_3

      mwPointer JCrows_1,JCrows_2,JCrows_3

      mwPointer JC_1,JC_2,JC_3

C     Array information:
      mwSize IArows,JArows
      mwSize IB1rows, JB1rows
      mwSize IB2rows, JB2rows
      mwSize IB3rows, JB3rows

      mwSize rows1, rows2, rows3
      mwSize JC_1dims(2)
      mwSize JC_2dims(2)
      mwSize JC_3dims(2)

      mwSize ICrows_2,ICrows_3,ndim

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 19) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '19 inputs required.')
      endif
      if(nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 3 outputs required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(7))
      JArows = mxGetM(prhs(8))

      IB1rows = mxGetM(prhs(9))
      JB1rows = mxGetM(prhs(10))
      IB2rows = mxGetM(prhs(11))
      JB2rows = mxGetM(prhs(12))
      IB3rows = mxGetM(prhs(13))
      JB3rows = mxGetM(prhs(14))

      ICrows_2 = mxGetM(prhs(18))
      ICrows_3 = mxGetM(prhs(19))

C     Create Fortran array from the input arguments.
      l = mxGetPr(prhs(1))
      m1 = mxGetPr(prhs(2))
      m2 = mxGetPr(prhs(3))
      n1 = mxGetPr(prhs(4))
      n2 = mxGetPr(prhs(5))
      n3 = mxGetPr(prhs(6))
      IA = mxGetPr(prhs(7))
      JA = mxGetPr(prhs(8))
      IB1 = mxGetPr(prhs(9))
      JB1 = mxGetPr(prhs(10))
      IB2 = mxGetPr(prhs(11))
      JB2 = mxGetPr(prhs(12))
      IB3 = mxGetPr(prhs(13))
      JB3 = mxGetPr(prhs(14))
      JCrows_1 = mxGetPr(prhs(15))
      JCrows_2 = mxGetPr(prhs(16))
      JCrows_3 = mxGetPr(prhs(17))
      IC_2 = mxGetPr(prhs(18))
      IC_3 = mxGetPr(prhs(19))

C     Create matrix for the return arguments.
      rows1=mxGetScalar(prhs(15))
      rows2=mxGetScalar(prhs(16))
      rows3=mxGetScalar(prhs(17))

      classid = mxGetClassID(prhs(8))

      complexflag = 0
      JC_1dims(1)=rows1
      JC_1dims(2)=3
      ndim=2

      plhs(1) = mxCreateNumericArray(ndim, JC_1dims, 
     +          classid, 
     +          complexflag)

      JC_1 = mxGetPr(plhs(1))

      JC_2dims(1)=rows2
      JC_2dims(2)=3
      plhs(2) = mxCreateNumericArray(ndim, JC_2dims, 
     +          classid, 
     +          complexflag)

      JC_2 = mxGetPr(plhs(2))

      JC_3dims(1)=rows3
      JC_3dims(2)=3
      plhs(3) = mxCreateNumericArray(ndim, JC_3dims, 
     +          classid, 
     +          complexflag)

      JC_3 = mxGetPr(plhs(3))
 
C     Call the computational routine.
      call compute(%VAL(l),%VAL(m1),%VAL(m2),
     + %VAL(n1),%VAL(n2),%VAL(n3),
     + %VAL(IA),%VAL(JA),
     + %VAL(IB1),%VAL(JB1),%VAL(IB2),%VAL(JB2),
     + %VAL(IB3),%VAL(JB3),
     + %VAL(JCrows_1),%VAL(JCrows_2),%VAL(JCrows_3),            
     + %VAL(IC_2),%VAL(IC_3),
     + %VAL(JC_1),%VAL(JC_2),%VAL(JC_3),
     + IArows,JArows,IB1rows,JB1rows,IB2rows,JB2rows,
     + IB3rows,JB3rows,
     + ICrows_2,ICrows_3)
 
      return
      end

C     Computational subroutine
      subroutine compute(l,m1,m2,n1,n2,n3,IA,JA,
     + IB1,JB1,IB2,JB2,IB3,JB3,
     + JCrows_1,JCrows_2,JCrows_3,
     + IC_2,IC_3,
     + JC_1,JC_2,JC_3,
     + IArows,JArows,
     + IB1rows,JB1rows,IB2rows,JB2rows,
     + IB3rows,JB3rows,
     + ICrows_2,ICrows_3)

C     Declarations
      implicit none

      mwSize  l,m1,m2,n1,n2,n3
      mwSize  IArows,JArows
      mwSize  IB1rows,JB1rows,IB2rows,JB2rows
      mwSize  IB3rows,JB3rows
      mwSize  JCrows_1,JCrows_2,JCrows_3
      mwSize  ICrows_2,ICrows_3

      mwSize, pointer :: pIA(:)
      mwSize, pointer :: pJA(:,:)
      mwSize, pointer :: pIB(:)
 
      mwSize, allocatable :: w(:,:,:)
     
      mwSize, target :: IA(IArows)
      mwSize, target :: JA(JArows,3)

      mwSize, target ::  IB1(IB1rows)
      mwSize  JB1(JB1rows)
      mwSize, target ::  IB2(IB2rows)
      mwSize  JB2(JB2rows)
      mwSize, target ::  IB3(IB3rows)
      mwSize  JB3(JB3rows)

      mwSize, target :: JC_1(JCrows_1,3)
      mwSize, target :: JC_2(JCrows_2,3)
      mwSize, target :: JC_3(JCrows_3,3)

      mwSize, target :: IC_2(ICrows_2)
      mwSize, target :: IC_3(ICrows_3)

      mwSize i,jp,j,kp,k1,k2,k3,vp,v1,v2,v3,nz

      integer AllocateStatus

C first contraction
      allocate (w(m1,m2,n3),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      w=0
      nz=1

      pIA=>IA
      pJA=>JA

      pIB=>IB3        

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        k1=pJA(jp,1)
        k2=pJA(jp,2)
        j=pJA(jp,3)
      do kp=pIB(j),pIB(j+1)-1
        k3=JB3(kp)
        if (w(k1,k2,k3)<i) then
         w(k1,k2,k3)=i
         JC_3(nz,1)=k1
         JC_3(nz,2)=k2
         JC_3(nz,3)=k3
         nz=nz+1
        endif
      end do
      end do
      end do

      deallocate (w)

C second contraction
      allocate (w(m1,n2,n3),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      w=0
      nz=1

      pIA=>IC_3
      pJA=>JC_3

      pIB=>IB2       
       
      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        k1=pJA(jp,1)
        j=pJA(jp,2)
        k3=pJA(jp,3)
      do kp=pIB(j),pIB(j+1)-1
        k2=JB2(kp)
        if (w(k1,k2,k3)<i) then
         w(k1,k2,k3)=i
         JC_2(nz,1)=k1
         JC_2(nz,2)=k2
         JC_2(nz,3)=k3
         nz=nz+1
        endif
      end do
      end do
      end do
      deallocate (w)

C third contraction
      allocate (w(n1,n2,n3),
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
        k3=pJA(jp,3)
      do kp=pIB(j),pIB(j+1)-1
        k1=JB1(kp)
        if (w(k1,k2,k3)<i) then
         w(k1,k2,k3)=i
         JC_1(nz,1)=k1
         JC_1(nz,2)=k2
         JC_1(nz,3)=k3
         nz=nz+1
        endif
      end do
      end do
      end do

      deallocate (w)

      return
      end

