C � Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [JC_1,JC_2,JC_3] = AkronBU3JC_mex(l,m,n,IA,JA,...
C            IB,JB,...
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
      mwPointer l,m,n,IA,JA
      mwPointer IB,JB
      mwPointer IC_2,IC_3

      mwPointer JCrows_1,JCrows_2,JCrows_3

      mwPointer JC_1,JC_2,JC_3

C     Array information:
      mwSize IArows,JArows
      mwSize IBrows, JBrows


      mwSize rows1, rows2, rows3
      mwSize JC_1dims(2)
      mwSize JC_2dims(2)
      mwSize JC_3dims(2)

      mwSize ICrows_2,ICrows_3,ndim

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 12) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '12 inputs required.')
      endif
      if(nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 3 outputs required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(4))
      JArows = mxGetM(prhs(5))

      IBrows = mxGetM(prhs(6))
      JBrows = mxGetM(prhs(7))

      ICrows_2 = mxGetM(prhs(11))
      ICrows_3 = mxGetM(prhs(12))

C     Create Fortran array from the input arguments.
      l = mxGetPr(prhs(1))
      m = mxGetPr(prhs(2))
      n = mxGetPr(prhs(3))
      IA = mxGetPr(prhs(4))
      JA = mxGetPr(prhs(5))
      IB = mxGetPr(prhs(6))
      JB = mxGetPr(prhs(7))

      JCrows_1 = mxGetPr(prhs(8))
      JCrows_2 = mxGetPr(prhs(9))
      JCrows_3 = mxGetPr(prhs(10))
      IC_2 = mxGetPr(prhs(11))
      IC_3 = mxGetPr(prhs(12))

C     Create matrix for the return arguments.
      rows1=mxGetScalar(prhs(8))
      rows2=mxGetScalar(prhs(9))
      rows3=mxGetScalar(prhs(10))

      classid = mxGetClassID(prhs(5))

      complexflag = 0
      JC_1dims(1)=rows1
      JC_1dims(2)=4
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
      call compute(%VAL(l),%VAL(m),%VAL(n),
     + %VAL(IA),%VAL(JA),
     + %VAL(IB),%VAL(JB),
     + %VAL(JCrows_1),%VAL(JCrows_2),%VAL(JCrows_3),            
     + %VAL(IC_2),%VAL(IC_3),
     + %VAL(JC_1),%VAL(JC_2),%VAL(JC_3),
     + IArows,JArows,IBrows,JBrows,
     + ICrows_2,ICrows_3)
 
      return
      end

C     Computational subroutine
      subroutine compute(l,m,n,IA,JA,
     + IB,JB,
     + JCrows_1,JCrows_2,JCrows_3,
     + IC_2,IC_3,
     + JC_1,JC_2,JC_3,
     + IArows,JArows,
     + IBrows,JBrows,
     + ICrows_2,ICrows_3)

C     Declarations
      implicit none

      mwSize  l,m,n
      mwSize  IArows,JArows
      mwSize  IBrows,JBrows

      mwSize  JCrows_1,JCrows_2,JCrows_3
      mwSize  ICrows_2,ICrows_3

      mwSize, pointer :: pIA(:)
      mwSize, pointer :: pJA(:,:)
      mwSize, pointer :: pIB(:)
 
      mwSize, allocatable :: w(:,:,:)
     
      mwSize, target :: IA(IArows)
      mwSize, target :: JA(JArows,3)

      mwSize, target ::  IB(IBrows)
      mwSize  JB(JBrows)

      mwSize, target :: JC_1(JCrows_1,4)
      mwSize, target :: JC_2(JCrows_2,3)
      mwSize, target :: JC_3(JCrows_3,3)

      mwSize, target :: IC_2(ICrows_2)
      mwSize, target :: IC_3(ICrows_3)

      mwSize i,jp,j,kp,k1,k2,k3,vp,v1,v2,v3,nz

      mwSize mati,k,kpstart,loc,j2
      parameter (k=3)
      real*8 mat(k,k), fact(k), prod1

      integer AllocateStatus

C prepare matrices for factorials and products

      mat(1,1:k)=(/ (mati, mati = n+1, n+k) /)
      mat(2,1:k)=(/ (mati, mati = n, n+k-1) /)
      do j=2,k-1
        mat(j+1,1:k-j+1)=(/ (mati, mati = n+k-j+1, n+1, -1) /)
      end do
      fact(1)=1
      do j=2,k
       fact(j)=fact(j-1)*j
      end do
      prod1=product(mat(2,:))

C first contraction
      allocate (w(m,m,n),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      w=0
      nz=1

      pIA=>IA
      pJA=>JA

      pIB=>IB

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        k1=pJA(jp,1)
        k2=pJA(jp,2)
        j=pJA(jp,3)
      do kp=pIB(j),pIB(j+1)-1
        k3=JB(kp)
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
      allocate (w(m,n,n),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      w=0
      nz=1

      pIA=>IC_3
      pJA=>JC_3

      pIB=>IB       
       
      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        k1=pJA(jp,1)
        j=pJA(jp,2)
        k3=pJA(jp,3)
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (JB(pIB(j+1)-1).ge.k3)) then
        kpstart=pIB(j)
      do while (JB(kpstart)<k3)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k2=JB(kp)
        if (w(k1,k2,k3)<i) then
         w(k1,k2,k3)=i
         JC_2(nz,1)=k1
         JC_2(nz,2)=k2
         JC_2(nz,3)=k3
         nz=nz+1
        endif
      end do
      endif
      end do
      end do
      deallocate (w)

C third contraction
      allocate (w(n,n,n),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      w=0
      nz=1

      pIA=>IC_2
      pJA=>JC_2

      pIB=>IB

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        j=pJA(jp,1)
        k2=pJA(jp,2)
        k3=pJA(jp,3)
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (JB(pIB(j+1)-1).ge.k2)) then
        kpstart=pIB(j)
      do while (JB(kpstart)<k2)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k1=JB(kp)
        if (w(k1,k2,k3)<i) then
         w(k1,k2,k3)=i
C determine the location of (k1,k2) in a compressed matrix and assign in JC_1
         loc=(prod1-product(mat(1,:)-k3))/fact(k)+1
         loc=loc+(product(mat(3,1:k-1)-k3)-product(mat(3,1:k-1)-k2))/
     +       fact(k-1)
         loc=loc+k1-k2
         JC_1(nz,1)=loc
         JC_1(nz,2)=k1
         JC_1(nz,3)=k2
         JC_1(nz,4)=k3
         nz=nz+1
        endif
      end do
      endif
      end do
      end do

      deallocate (w)

      return
      end

