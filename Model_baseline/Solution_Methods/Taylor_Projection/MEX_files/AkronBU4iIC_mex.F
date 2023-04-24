C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [IC_1,IC_2,IC_3,IC_4] = AkronBU4iIC_mex(l,m,n,IA,JA,IB,JB,convertind,maxind)
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
      mwPointer IC_1,IC_2,IC_3,IC_4
      mwPointer convertind,maxind

C     Array information:
      mwSize IArows,JArows
      mwSize IBrows, JBrows,ndim
C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 9) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '9 inputs required.')
      endif
      if(nlhs .gt. 4) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 4 outputs required.')
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
      convertind = mxGetPr(prhs(8))
      maxind = mxGetPr(prhs(9))

C     Create matrix for the return arguments.
      classid = mxGetClassID(prhs(4))
      ndim=1
      complexflag = 0
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

      plhs(4) = mxCreateNumericArray(ndim, IArows, 
     +          classid, 
     +          complexflag)

      IC_4 = mxGetPr(plhs(4))

C     Call the computational routine.
      call compute(%VAL(l),%VAL(m),%VAL(n),
     + %VAL(IA),%VAL(JA),
     + %VAL(IB),%VAL(JB),
     + %VAL(IC_1),%VAL(IC_2),%VAL(IC_3),%VAL(IC_4),
     + %VAL(convertind),%VAL(maxind),
     + IArows,JArows,IBrows,JBrows)

      return
      end

C     Computational subroutine
      subroutine compute(l,m,n,IA,JA,
     + IB,JB,
     + IC_1,IC_2,IC_3,IC_4,
     + convertind,maxind,
     + IArows,JArows,
     + IBrows,JBrows)

C     Declarations
      implicit none

      mwSize  l,m,n,smallm
      mwSize  IArows,JArows
      mwSize  IBrows,JBrows

      mwSize, pointer :: pIA(:)
      mwSize, pointer :: pJA(:,:)
      mwSize, pointer :: pIB(:)
 
      mwSize, allocatable :: w(:,:,:,:)
     
      mwSize, target :: IA(IArows)
      mwSize, target :: JA(JArows,2)

      mwSize, target ::  IB(IBrows)
      mwSize  JB(JBrows)

      mwSize, target :: IC_1(IArows)
      mwSize, target :: IC_2(IArows)
      mwSize, target :: IC_3(IArows)
      mwSize, target :: IC_4(IArows)

      mwSize, target,allocatable :: JC_2(:,:)
      mwSize, target,allocatable :: JC_3(:,:)
      mwSize, target,allocatable :: JC_4(:,:)

      mwSize i,jp,j,kp,k1,k2,k3,k4,vp,v1,v2,v3,v4,nz

      mwSize kpstart,maxind
      mwSize convertind(l,m)

      integer AllocateStatus

      smallm=maxind ! here is the memory saving

C first contraction
      allocate (w(smallm,smallm,smallm,n),
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
        IC_4(i)=nz
      do jp=pIA(i),pIA(i+1)-1
        k1=convertind(i,pJA(jp,1))
        k2=convertind(i,pJA(jp,2)) 
        k3=convertind(i,pJA(jp,3)) 
        j=pJA(jp,4)
      do kp=pIB(j),pIB(j+1)-1
        k4=JB(kp)
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
         nz=nz+1
        endif
      end do
      end do
      end do

      IC_4(IArows)=nz
      
      w=0
      nz=1

C calculate temporary JC_4
      allocate (JC_4(IC_4(IArows)-1,4),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        k1=convertind(i,pJA(jp,1))
        k2=convertind(i,pJA(jp,2)) 
        k3=convertind(i,pJA(jp,3))
        j=pJA(jp,4)
      do kp=pIB(j),pIB(j+1)-1
        k4=JB(kp)
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
         JC_4(nz,1)=pJA(jp,1)
         JC_4(nz,2)=pJA(jp,2)
         JC_4(nz,3)=pJA(jp,3)
         JC_4(nz,4)=k4
         nz=nz+1
        endif
      end do
      end do
      end do

      deallocate (w)

C second contraction
      allocate (w(smallm,smallm,n,n),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      w=0
      nz=1

      pIA=>IC_4
      pJA=>JC_4

      pIB=>IB        
       
      do i=1,l
        IC_3(i)=nz
      do jp=pIA(i),pIA(i+1)-1
        k1=convertind(i,pJA(jp,1))
        k2=convertind(i,pJA(jp,2))
        j=pJA(jp,3)
        k4=pJA(jp,4)
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (JB(pIB(j+1)-1).ge.k4)) then
        kpstart=pIB(j)
      do while (JB(kpstart)<k4)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k3=JB(kp)
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
         nz=nz+1
        endif
      end do
      endif
      end do
      end do

      IC_3(IArows)=nz
      
      w=0
      nz=1

C calculate temporary JC_3
      allocate (JC_3(IC_3(IArows)-1,4),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        k1=convertind(i,pJA(jp,1))
        k2=convertind(i,pJA(jp,2))
        j=pJA(jp,3)
        k4=pJA(jp,4)
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (JB(pIB(j+1)-1).ge.k4)) then
        kpstart=pIB(j)
      do while (JB(kpstart)<k4)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k3=JB(kp)
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
         JC_3(nz,1)=pJA(jp,1)
         JC_3(nz,2)=pJA(jp,2)
         JC_3(nz,3)=k3
         JC_3(nz,4)=k4
         nz=nz+1
        endif
      end do
      endif
      end do
      end do

      deallocate (w)

C third contraction
      allocate (w(smallm,n,n,n),
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
        k1=convertind(i,pJA(jp,1))
        j=pJA(jp,2)
        k3=pJA(jp,3)
        k4=pJA(jp,4)
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (JB(pIB(j+1)-1).ge.k3)) then
        kpstart=pIB(j)
      do while (JB(kpstart)<k3)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k2=JB(kp)
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
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
      allocate (JC_2(IC_2(IArows)-1,4),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        k1=convertind(i,pJA(jp,1))
        j=pJA(jp,2)
        k3=pJA(jp,3)
        k4=pJA(jp,4)
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (JB(pIB(j+1)-1).ge.k3)) then
        kpstart=pIB(j)
      do while (JB(kpstart)<k3)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k2=JB(kp)
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
         JC_2(nz,1)=pJA(jp,1)
         JC_2(nz,2)=k2
         JC_2(nz,3)=k3
         JC_2(nz,4)=k4
         nz=nz+1
        endif
      end do
      endif
      end do
      end do
      deallocate (w)

C fourth contraction
      allocate (w(n,n,n,n),
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
        k4=pJA(jp,4)
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (JB(pIB(j+1)-1).ge.k2)) then
        kpstart=pIB(j)
      do while (JB(kpstart)<k2)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k1=JB(kp)
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
         nz=nz+1
        endif
      end do
      endif
      end do
      end do

      IC_1(IArows)=nz        

      
      deallocate (JC_2)
      deallocate (JC_3)
      deallocate (JC_4)
      deallocate (w)

      return
      end

