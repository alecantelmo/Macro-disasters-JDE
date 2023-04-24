C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [IC_1,IC_2,IC_3] = AkronBU3IC_mex(l,m,n,IA,JA,IB,JB)
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
      mwPointer l,m,n,IA,JA
      mwPointer IB,JB
      mwPointer IC_1,IC_2,IC_3

C     Array information:
      mwSize IArows,JArows
      mwSize IBrows, JBrows,ndim


C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 7) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '14 inputs required.')
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

C     Create Fortran array from the input arguments.
      l = mxGetPr(prhs(1))
      m = mxGetPr(prhs(2))
      n = mxGetPr(prhs(3))
      IA = mxGetPr(prhs(4))
      JA = mxGetPr(prhs(5))
      IB = mxGetPr(prhs(6))
      JB = mxGetPr(prhs(7))

C     Create matrix for the return arguments.
      classid = mxGetClassID(prhs(4))

      complexflag = 0
      ndim=1
      plhs(1) = mxCreateNumericArray(ndim, IArows, 
     +          classid, 
     +          complexflag)

      IC_1 = mxGetPr(plhs(1))

      plhs(2) = mxCreateNumericArray(ndim, IArows, 
     +          classid, 
     +          complexflag)

      IC_2 = mxGetPr(plhs(2))

      plhs(3) = mxCreateNumericArray(ndim, IArows, 
     +          classid, 
     +          complexflag)

      IC_3 = mxGetPr(plhs(3))

C     Call the computational routine.
      call compute(%VAL(l),%VAL(m),%VAL(n),
     + %VAL(IA),%VAL(JA),
     + %VAL(IB),%VAL(JB),
     + %VAL(IC_1),%VAL(IC_2),%VAL(IC_3),            
     + IArows,JArows,IBrows,JBrows)

      return
      end

C     Computational subroutine
      subroutine compute(l,m,n,IA,JA,
     + IB,JB,
     + IC_1,IC_2,IC_3,
     + IArows,JArows,
     + IBrows,JBrows)

C     Declarations
      implicit none

      mwSize  l,m,n
      mwSize  IArows,JArows
      mwSize  IBrows,JBrows

      mwSize, pointer :: pIA(:)
      mwSize, pointer :: pJA(:,:)
      mwSize, pointer :: pIB(:)
 
      mwSize, allocatable :: w(:,:,:)
     
      mwSize, target :: IA(IArows)
      mwSize, target :: JA(JArows,3)

      mwSize, target ::  IB(IBrows)
      mwSize  JB(JBrows)

      mwSize, target :: IC_1(IArows)
      mwSize, target :: IC_2(IArows)
      mwSize, target :: IC_3(IArows)

      mwSize, target,allocatable :: JC_2(:,:)
      mwSize, target,allocatable :: JC_3(:,:)

      mwSize i,jp,j,kp,k1,k2,k3,vp,v1,v2,v3,nz

      mwSize kpstart

      integer AllocateStatus

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
        IC_3(i)=nz
      do jp=pIA(i),pIA(i+1)-1
        k1=pJA(jp,1)
        k2=pJA(jp,2)
        j=pJA(jp,3)
      do kp=pIB(j),pIB(j+1)-1
        k3=JB(kp)
        if (w(k1,k2,k3)<i) then
         w(k1,k2,k3)=i
         nz=nz+1
        endif
      end do
      end do
      end do

      IC_3(IArows)=nz
      
      w=0
      nz=1

C calculate temporary JC_3
      allocate (JC_3(IC_3(IArows)-1,3),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

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
        IC_2(i)=nz
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
         nz=nz+1
        endif
      end do
      endif
      end do
      end do

      IC_2(IArows)=nz
      
      w=0
      nz=1

C calculate temporary JC_2
      allocate (JC_2(IC_2(IArows)-1,3),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

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
        IC_1(i)=nz
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
         nz=nz+1
        endif
      end do
      endif
      end do
      end do

      IC_1(IArows)=nz        

      
      deallocate (JC_2)
      deallocate (JC_3)
      deallocate (w)

      return
      end

