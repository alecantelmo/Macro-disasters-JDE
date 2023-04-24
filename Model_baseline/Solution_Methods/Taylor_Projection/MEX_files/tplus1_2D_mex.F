C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [NC,IC,JC] = tplus1_2D_mex(l1,l2,m1,...
C            IB1,JB1,NB1,IB2,JB2,NB2,
C            maxload)
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxGetM, mxGetN
      mwPointer mxCreateDoubleMatrix
      mwPointer mxCreateNumericArray
      integer*4 mxGetClassID

C     Arguments for mxCreateNumericArray
      integer*4 classid
      integer*4 complexflag

C     Pointers to input/output mxArrays:
      mwPointer l1,l2,m1
      mwPointer IB1,JB1,IB2,JB2
      mwPointer IC
      mwPointer JC
      mwPointer maxload
      mwPointer NB1,NB2

      mwPointer NC

C     Array information:
      mwSize IB1rows, IB1cols, JB1rows
      mwSize IB2rows, JB2rows
      mwSize ICrows,ICcols
      mwSize JCrows

      mwSize s,sB1,sB2,ndim
      mwSize dim(2)

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 10) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '10 inputs required.')
      endif
      if(nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 3 outputs required.')
      endif

C     Get the size of the inputs array.

      IB1rows = mxGetM(prhs(4))
      IB1cols = mxGetN(prhs(4))
      JB1rows = mxGetM(prhs(5))
      IB2rows = mxGetM(prhs(7))
      JB2rows = mxGetM(prhs(8))
      
      ICrows=IB1rows
      ICcols=IB1cols

      JCrows=JB1rows+JB2rows ! temporary size, should be adjusted after calling the mex function.

C     Create Fortran array from the input arguments.
      l1 = mxGetPr(prhs(1))
      l2 = mxGetPr(prhs(2))
      m1 = mxGetPr(prhs(3))
      IB1 = mxGetPr(prhs(4))
      JB1 = mxGetPr(prhs(5))
      NB1 = mxGetPr(prhs(6))
      IB2 = mxGetPr(prhs(7))
      JB2 = mxGetPr(prhs(8))
      NB2 = mxGetPr(prhs(9))
      maxload = mxGetPr(prhs(10))

C     Get number of states
      sB1 = mxGetM(prhs(6))
      sB2 = mxGetM(prhs(9))

      s=sB1

C check compatibility of states
      if (sB1.ne.sB2) then
         call mexErrMsgIdAndTxt ('MATLAB:states',
     +   'incompatible states.')
      endif

C     Create matrix for NC.
      complexflag = 0

      plhs(1) = mxCreateDoubleMatrix(s,JCrows,complexflag)
      NC = mxGetPr(plhs(1))

C     Create matrix for IC and JC
      classid = mxGetClassID(prhs(4))
      ndim=2
      dim(1)=ICrows
      dim(2)=ICcols
      plhs(2) = mxCreateNumericArray(ndim, dim, 
     +          classid, 
     +          complexflag)

      IC = mxGetPr(plhs(2))

      classid = mxGetClassID(prhs(5))
      ndim=2
      dim(1)=JCrows
      dim(2)=1
      plhs(3) = mxCreateNumericArray(ndim, dim, 
     +          classid, 
     +          complexflag)

      JC = mxGetPr(plhs(3))

C     Call the computational routine.
      call compute(%VAL(l1),%VAL(l2),%VAL(m1),
     + %VAL(IB1),%VAL(JB1),%VAL(NB1),
     + %VAL(IB2),%VAL(JB2),%VAL(NB2),
     + %VAL(maxload),
     + %VAL(NC),%VAL(IC),%VAL(JC),            
     + IB1rows,JB1rows,IB2rows,JB2rows,
     + ICrows,JCrows,
     + s)

      return
      end

C     Computational subroutine
      subroutine compute(l1,l2,m1,
     + IB1,JB1,NB1,IB2,JB2,NB2,
     + maxload,
     + NC,IC,JC,
     + IB1rows,JB1rows,IB2rows,JB2rows,
     + ICrows,JCrows,
     + s)

C     Declarations
      implicit none

      mwSize  l1,l2,m1,s
      mwSize  IB1rows,JB1rows,IB2rows,JB2rows
      mwSize  ICrows,ICcols
      mwSize  JCrows

      mwSize  maxload,n_vecload,loadi,load_
      logical doind
      real*8, target ::  NB1(s,JB1rows)
      real*8, target ::  NB2(s,JB2rows)

      real*8, pointer :: pNB1(:,:)
      real*8, pointer :: pNB2(:,:)

      real*8  NC(s,JCrows)

      mwSize, allocatable :: vecload(:)

      real*8, allocatable :: x(:,:) 
      mwSize, allocatable :: w(:)

      mwSize, target ::  IB1(l1,l2+1)
      mwSize, target ::  JB1(JB1rows,1)
      mwSize, target ::  IB2(l1,l2+1)
      mwSize, target ::  JB2(JB2rows,1)

      mwSize, target ::  IC(l1,l2+1)
      mwSize, target ::  JC(JCrows,1)

      mwSize i,i1,i2,jp,j,kp,k1,k2,vp,v1,v2
      mwSize is,js,nz

      mwSize jpstart, jpend
      mwSize kpstart, kpend
      mwSize vpstart, vpend
 
      integer AllocateStatus

C vecload
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

C allocate
      allocate (w(m1),
     +   x(min(s,maxload),m1),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory.')
      endif

      nz=1

C run first time and calculate indices
      loadi=1
      load_=vecload(loadi+1)-vecload(loadi)
       
      if (s>1) then
        pNB1=>NB1(vecload(loadi):vecload(loadi+1)-1,:)
        pNB2=>NB2(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNB1=>NB1
        pNB2=>NB2
      endif

      do i1=1,l1
       w=0
      do i2=1,l2
        IC(i1,i2)=nz
C run over rows of B1
      do jp=IB1(i1,i2),IB1(i1,i2+1)-1
        j=JB1(jp,1)
        !DEC$ simd
        do is=1,load_
         x(is,j)=pNB1(is,jp)
        end do
        JC(nz,1)=j
        w(j)=i2
        nz=nz+1
       end do
C run over rows of B2
      do jp=IB2(i1,i2),IB2(i1,i2+1)-1
        j=JB2(jp,1)
         if (w(j)<i2) then
         !DEC$ simd
         do is=1,load_
          x(is,j)=pNB2(is,jp)
         end do
         JC(nz,1)=j
         w(j)=i2
         nz=nz+1
         else
         !DEC$ simd
          do is=1,load_
           x(is,j)=x(is,j)+pNB2(is,jp)
          end do
         endif
       end do
C assign result in NC
      do vp=IC(i1,i2),nz-1
       v1=JC(vp,1)
       !DEC$ simd
       do is=1,load_
         js=vecload(loadi)+is-1
         NC(js,vp)=x(is,v1)
         x(is,v1)=0.0
       end do
      end do
      end do
      IC(i1,l2+1)=nz
       
      end do

C if needed run again, but don't calcualte indices
      if (n_vecload>2) then
        w=0
      endif

      do loadi=2,n_vecload-1
        load_=vecload(loadi+1)-vecload(loadi)
       
      if (s>1) then
        pNB1=>NB1(vecload(loadi):vecload(loadi+1)-1,:)
        pNB2=>NB2(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNB1=>NB1
        pNB2=>NB2
      endif

      do i1=1,l1
      do i2=1,l2
C run over rows of B1
      do jp=IB1(i1,i2),IB1(i1,i2+1)-1
        j=JB1(jp,1)
        !DEC$ simd
        do is=1,load_
         x(is,j)=pNB1(is,jp)
        end do
       end do
C run over rows of B2
      do jp=IB2(i1,i2),IB2(i1,i2+1)-1
       j=JB2(jp,1)
       !DEC$ simd
       do is=1,load_
        x(is,j)=x(is,j)+pNB2(is,jp)
       end do
      end do
C assign result in NC
      do vp=IC(i1,i2),IC(i1,i2+1)-1
       v1=JC(vp,1)
       !DEC$ simd
       do is=1,load_
         js=vecload(loadi)+is-1
         NC(js,vp)=x(is,v1)
         x(is,v1)=0.0
       end do
      end do
      end do

      end do
      end do
   
      deallocate (w,x)

      return
      end

