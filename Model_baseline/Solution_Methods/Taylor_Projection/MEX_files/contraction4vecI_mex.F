C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [NC,JC_1] = contraction4vecI_mex(l,m1,m2,m3,n1,n2,n3,n4,IA,JA,NA,...
C            IB1,JB1,NB1,IB2,JB2,NB2,IB3,JB3,NB3,IB4,JB4,NB4,...
C            IC_1,IC_2,IC_3,IC_4,JCrows_1,...
C            maxload,dosum) 
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
      real*8 mxGetScalar

C     Arguments for mxCreateNumericArray
      integer*4 classid
      integer*4 complexflag

C     Pointers to input/output mxArrays:
      mwPointer l,IA,JA
      mwPointer m1,m2,m3,n1,n2,n3,n4
      mwPointer IB1,JB1,IB2,JB2,IB3,JB3,IB4,JB4
      mwPointer IC_1,IC_2,IC_3,IC_4
      mwPointer maxload,dosum
      mwPointer NA,NB1,NB2,NB3,NB4

      mwPointer NC,JC_1

C     Array information:
      mwSize IArows,JArows
      mwSize IB1rows, JB1rows
      mwSize IB2rows, JB2rows
      mwSize IB3rows, JB3rows
      mwSize IB4rows, JB4rows
      mwSize ICrows_1,ICrows_2,ICrows_3,ICrows_4
      mwSize JCrows_1

      mwSize s,sC,sA,sB1,sB2,sB3,sB4
      real*8 do_sum

      mwSize rows
      mwSize JCdims(2),ndim

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 30) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '30 inputs required.')
      endif
      if(nlhs .gt. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 2 outputs required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(9))
      JArows = mxGetM(prhs(10))

      IB1rows = mxGetM(prhs(12))
      JB1rows = mxGetM(prhs(13))
      IB2rows = mxGetM(prhs(15))
      JB2rows = mxGetM(prhs(16))
      IB3rows = mxGetM(prhs(18))
      JB3rows = mxGetM(prhs(19))
      IB4rows = mxGetM(prhs(21))
      JB4rows = mxGetM(prhs(22))

      ICrows_1 = mxGetM(prhs(24))
      ICrows_2 = mxGetM(prhs(25))
      ICrows_3 = mxGetM(prhs(26))
      ICrows_4 = mxGetM(prhs(27))

      rows = mxGetScalar(prhs(28))

C     Create Fortran array from the input arguments.
      l = mxGetPr(prhs(1))
      m1 = mxGetPr(prhs(2))
      m2 = mxGetPr(prhs(3))
      m3 = mxGetPr(prhs(4))
      n1 = mxGetPr(prhs(5))
      n2 = mxGetPr(prhs(6))
      n3 = mxGetPr(prhs(7))
      n4 = mxGetPr(prhs(8))
      IA = mxGetPr(prhs(9))
      JA = mxGetPr(prhs(10))
      NA = mxGetPr(prhs(11))
      IB1 = mxGetPr(prhs(12))
      JB1 = mxGetPr(prhs(13))
      NB1 = mxGetPr(prhs(14))
      IB2 = mxGetPr(prhs(15))
      JB2 = mxGetPr(prhs(16))
      NB2 = mxGetPr(prhs(17))
      IB3 = mxGetPr(prhs(18))
      JB3 = mxGetPr(prhs(19))
      NB3 = mxGetPr(prhs(20))
      IB4 = mxGetPr(prhs(21))
      JB4 = mxGetPr(prhs(22))
      NB4 = mxGetPr(prhs(23))
      IC_1 = mxGetPr(prhs(24))
      IC_2 = mxGetPr(prhs(25))
      IC_3 = mxGetPr(prhs(26))
      IC_4 = mxGetPr(prhs(27))
      JCrows_1 = mxGetPr(prhs(28))
      maxload = mxGetPr(prhs(29))

C     Get number of states
      sA = mxGetM(prhs(11))
      sB1 = mxGetM(prhs(14))
      sB2 = mxGetM(prhs(17))
      sB3 = mxGetM(prhs(20))
      sB4 = mxGetM(prhs(23))

C check compatibility of states
      s=max(sA,sB1,sB2,sB3,sB4)
      if ((sA.ne.s .and. sA.ne.1) .or. (sB1.ne.s .and. sB1.ne.1)
     +   .or. (sB1.ne.sB2) .or. (sB1.ne.sB3) .or. (sB1.ne.sB4)) then
         call mexErrMsgIdAndTxt ('MATLAB:states',
     +   'incompatible states.')
      endif

C     sum across rows of C if dosum==1
      do_sum=mxGetScalar(prhs(30))
 
      if (do_sum==1) then
       sC=1
      else
       sC=s
      endif

C     Create matrix for NC.

      plhs(1) = mxCreateDoubleMatrix(sC,rows, 0)
      NC = mxGetPr(plhs(1))

C     Create matrix for JC_1.

      classid = mxGetClassID(prhs(10))

      complexflag = 0
      JCdims(1)=rows
      JCdims(2)=4
      ndim=2
      plhs(2) = mxCreateNumericArray(ndim, JCdims, 
     +          classid, 
     +          complexflag)

      JC_1 = mxGetPr(plhs(2))
  
C     Call the computational routine.
      call compute(%VAL(l),%VAL(m1),%VAL(m2),%VAL(m3),
     + %VAL(n1),%VAL(n2),%VAL(n3),%VAL(n4),
     + %VAL(IA),%VAL(JA),%VAL(NA),
     + %VAL(IB1),%VAL(JB1),%VAL(NB1),
     + %VAL(IB2),%VAL(JB2),%VAL(NB2),
     + %VAL(IB3),%VAL(JB3),%VAL(NB3),
     + %VAL(IB4),%VAL(JB4),%VAL(NB4),
     + %VAL(IC_1),%VAL(IC_2),%VAL(IC_3),%VAL(IC_4),
     + %VAL(JCrows_1),%VAL(maxload),
     + %VAL(NC),%VAL(JC_1),            
     + IArows,JArows,IB1rows,JB1rows,IB2rows,JB2rows,
     + IB3rows,JB3rows,IB4rows,JB4rows,
     + ICrows_1,ICrows_2,ICrows_3,ICrows_4,
     + sC,sA,sB1,sB2,sB3,sB4)
 
      return
      end

C     Computational subroutine
      subroutine compute(l,m1,m2,m3,n1,n2,n3,n4,IA,JA,NA,
     + IB1,JB1,NB1,IB2,JB2,NB2,IB3,JB3,NB3,IB4,JB4,NB4,
     + IC_1,IC_2,IC_3,IC_4,
     + JCrows_1,maxload,
     + NC,JC_1,
     + IArows,JArows,
     + IB1rows,JB1rows,IB2rows,JB2rows,
     + IB3rows,JB3rows,IB4rows,JB4rows,
     + ICrows_1,ICrows_2,ICrows_3,ICrows_4,
     + sC,sA,sB1,sB2,sB3,sB4)

C     Declarations
      implicit none

      mwSize  l,m1,m2,m3,n1,n2,n3,n4,sC,sA,sB1,sB2,sB3,sB4
      mwSize  IArows,JArows
      mwSize  IB1rows,JB1rows,IB2rows,JB2rows
      mwSize  IB3rows,JB3rows,IB4rows,JB4rows
      mwSize  ICrows_1,ICrows_2,ICrows_3,ICrows_4
      mwSize  JCrows_1,JCrows_2,JCrows_3,JCrows_4

      mwSize  s,sB,sC4,sC3,sC2,rows
      mwSize  maxload,n_vecload,loadi,load_

      real*8, target ::  NA(sA,JArows)
      real*8, target ::  NB1(sB1,JB1rows)
      real*8, target ::  NB2(sB2,JB2rows)
      real*8, target ::  NB3(sB3,JB3rows)
      real*8, target ::  NB4(sB4,JB4rows)

      real*8  NC(sC,JCrows_1)

      real*8, pointer :: pNA(:,:),pNB(:,:)
      mwSize, pointer :: pIA(:),pIB(:),pIC(:)
      mwSize, pointer :: pJA(:,:),pJC(:,:)
      mwSize, pointer :: pJB(:)

       mwSize, allocatable :: vecload(:)

      real*8, allocatable :: x(:,:,:,:,:) 
      mwSize, allocatable :: w(:,:,:,:)

C three temporary crs vars are needed
      real*8, target, allocatable :: tempNC4(:,:)
      mwSize, target, allocatable :: tempJC4(:,:)
      real*8, target, allocatable :: tempNC3(:,:)
      mwSize, target, allocatable :: tempJC3(:,:)
      real*8, target, allocatable :: tempNC2(:,:)
      mwSize, target, allocatable :: tempJC2(:,:)

      mwSize, target ::  IA(IArows)
      mwSize, target ::  JA(JArows,4)

      mwSize, target ::  IB1(IB1rows)
      mwSize, target ::  JB1(JB1rows)
      mwSize, target ::  IB2(IB2rows)
      mwSize, target ::  JB2(JB2rows)
      mwSize, target ::  IB3(IB3rows)
      mwSize, target ::  JB3(JB3rows)
      mwSize, target ::  IB4(IB4rows)
      mwSize, target ::  JB4(JB4rows)

      mwSize, target ::  IC_1(ICrows_1)
      mwSize, target ::  JC_1(JCrows_1,4)

      mwSize, target ::  IC_2(ICrows_2)
      mwSize, target ::  IC_3(ICrows_3)
      mwSize, target ::  IC_4(ICrows_4)

      mwSize i,jp,j,kp,k1,k2,k3,k4,vp,v1,v2,v3,v4
      mwSize is,js,nz

      mwSize jpstart, jpend
      mwSize kpstart, kpend
      mwSize vpstart, vpend
 
      integer AllocateStatus

C vecload
      s=max(sA,sB1,sB2,sB3,sB4)

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

      sB=sB1 ! sB1, sB2, sB3 and sB4 must all be equal

C first contraction
      sC4=max(sA,sB4)
      rows=IC_4(ICrows_4)-1
      allocate (w(m1,m2,m3,n4),
     +   x(min(sC4,maxload),m1,m2,m3,n4),
     +   tempNC4(sC4,rows),
     +   tempJC4(rows,4),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory2.')
      endif

      w=0
      nz=1

      pIA=>IA
      pJA=>JA

      pIB=>IB4        
      pJB=>JB4

      pIC=>IC_4
      pJC=>tempJC4

      do loadi=1,n_vecload-1
        load_=vecload(loadi+1)-vecload(loadi)
       
      if (sA>1) then
        pNA=>NA(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNA=>NA
      endif

      if (sB4>1) then
        pNB=>NB4(vecload(loadi):vecload(loadi+1)-1,:)       
      else 
        pNB=>NB4
      endif

      do i=1,l

      do jp=pIA(i),pIA(i+1)-1
        k1=pJA(jp,1)
        k2=pJA(jp,2)
        k3=pJA(jp,3)
        j=pJA(jp,4)
      do kp=pIB(j),pIB(j+1)-1
        k4=JB4(kp)

C collect information for pJC and initialize x to zero
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
         pJC(nz,1)=k1
         pJC(nz,2)=k2
         pJC(nz,3)=k3
         pJC(nz,4)=k4
         nz=nz+1
         x(:,k1,k2,k3,k4)=0
        endif

C Choose among three options:
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(1,jp)*pNB(is,kp)
       end do
      endif

      end do
      end do

C     assign the result in tempNC4
      do vp=pIC(i),pIC(i+1)-1
         v1=pJC(vp,1)
         v2=pJC(vp,2)
         v3=pJC(vp,3)
         v4=pJC(vp,4)
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         tempNC4(js,vp)=x(is,v1,v2,v3,v4)
         x(is,v1,v2,v3,v4)=0.0
      end do
      end do

      end do

      end do

      deallocate (w,x)

C second contraction
      sC3=max(sC4,sB3)
      rows=IC_3(ICrows_3)-1
      allocate (w(m1,m2,n3,n4),
     +   x(min(sC3,maxload),m1,m2,n3,n4),
     +   tempNC3(sC3,rows),
     +   tempJC3(rows,4),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory2.')
      endif

      w=0
      nz=1

      pIA=>IC_4
      pJA=>tempJC4

      pIB=>IB3        
      pJB=>JB3

      pIC=>IC_3
      pJC=>tempJC3

      sA=sC4

      do loadi=1,n_vecload-1
        load_=vecload(loadi+1)-vecload(loadi)
       
      if (sC4>1) then
        pNA=>tempNC4(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNA=>tempNC4
      endif

      if (sB3>1) then
        pNB=>NB3(vecload(loadi):vecload(loadi+1)-1,:)       
      else 
        pNB=>NB3
      endif

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        k1=pJA(jp,1)
        k2=pJA(jp,2)
        j=pJA(jp,3)
        k4=pJA(jp,4)
      do kp=pIB(j),pIB(j+1)-1
        k3=JB3(kp)
C collect information for pJC and initialize x to zero
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
         pJC(nz,1)=k1
         pJC(nz,2)=k2
         pJC(nz,3)=k3
         pJC(nz,4)=k4
         nz=nz+1
         x(:,k1,k2,k3,k4)=0
        endif

C Choose among three options:
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(1,jp)*pNB(is,kp)
       end do
      endif

      end do
      end do

C     assign the result in tempNC3
      do vp=pIC(i),pIC(i+1)-1
         v1=pJC(vp,1)
         v2=pJC(vp,2)
         v3=pJC(vp,3)
         v4=pJC(vp,4)
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         tempNC3(js,vp)=x(is,v1,v2,v3,v4)
         x(is,v1,v2,v3,v4)=0.0
      end do
      end do

      end do

      end do

      deallocate (w,x)
      deallocate (tempNC4,tempJC4)

C third contraction
      sC2=max(sC3,sB2)
      rows=IC_2(ICrows_2)-1
      allocate (w(m1,n2,n3,n4),
     +   x(min(sC2,maxload),m1,n2,n3,n4),
     +   tempNC2(sC2,rows),
     +   tempJC2(rows,4),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory3.')
      endif

      w=0
      nz=1

      pIA=>IC_3
      pJA=>tempJC3

      pIB=>IB2        
      pJB=>JB2

      pIC=>IC_2
      pJC=>tempJC2

      do loadi=1,n_vecload-1
        load_=vecload(loadi+1)-vecload(loadi)
       
      if (sC3>1) then
        pNA=>tempNC3(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNA=>tempNC3
      endif

      if (sB2>1) then
        pNB=>NB2(vecload(loadi):vecload(loadi+1)-1,:)       
      else 
        pNB=>NB2
      endif

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        k1=pJA(jp,1)
        j=pJA(jp,2)
        k3=pJA(jp,3)
        k4=pJA(jp,4)
      do kp=pIB(j),pIB(j+1)-1
        k2=JB2(kp)
C collect information for pJC and initializes x to zero
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
         pJC(nz,1)=k1
         pJC(nz,2)=k2
         pJC(nz,3)=k3
         pJC(nz,4)=k4
         nz=nz+1
         x(:,k1,k2,k3,k4)=0
        endif

C Choose among three options:
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(1,jp)*pNB(is,kp)
       end do
      endif

      end do
      end do

C     assign the result in tempNC2
      do vp=pIC(i),pIC(i+1)-1
         v1=pJC(vp,1)
         v2=pJC(vp,2)
         v3=pJC(vp,3)
         v4=pJC(vp,4)
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         tempNC2(js,vp)=x(is,v1,v2,v3,v4)
         x(is,v1,v2,v3,v4)=0.0
      end do
      end do

      end do

      end do

      deallocate (w,x)
      deallocate (tempNC3,tempJC3)

C fourth contraction
      rows=IC_1(ICrows_1)-1
      allocate (w(n1,n2,n3,n4),
     +   x(min(sC,maxload),n1,n2,n3,n4),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory4.')
      endif

      w=0
      nz=1

      pIA=>IC_2
      pJA=>tempJC2

      pIB=>IB1        
      pJB=>JB1

      pIC=>IC_1
      pJC=>JC_1

      do loadi=1,n_vecload-1
        load_=vecload(loadi+1)-vecload(loadi)
       
      if (sC2>1) then
        pNA=>tempNC2(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNA=>tempNC2
      endif

      if (sB1>1) then
        pNB=>NB1(vecload(loadi):vecload(loadi+1)-1,:)       
      else 
        pNB=>NB1
      endif

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        j=pJA(jp,1)
        k2=pJA(jp,2)
        k3=pJA(jp,3)
        k4=pJA(jp,4)
      do kp=pIB(j),pIB(j+1)-1
        k1=JB1(kp)
C collect information for pJC and initializes x to zero
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
         pJC(nz,1)=k1
         pJC(nz,2)=k2
         pJC(nz,3)=k3
         pJC(nz,4)=k4
         nz=nz+1
         x(:,k1,k2,k3,k4)=0
        endif

      if (sC>1) then
C Choose among three options and do not sum across rows
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3,k4)=x(is,k1,k2,k3,k4)+pNA(1,jp)*pNB(is,kp)
       end do
      endif
      else
C Choose among three options and sum across rows.
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2,k3,k4)=x(1,k1,k2,k3,k4)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2,k3,k4)=x(1,k1,k2,k3,k4)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2,k3,k4)=x(1,k1,k2,k3,k4)+pNA(1,jp)*pNB(is,kp)
       end do
      endif
      endif

      end do
      end do

C     assign the result in NC
      do vp=pIC(i),pIC(i+1)-1
         v1=pJC(vp,1)
         v2=pJC(vp,2)
         v3=pJC(vp,3)
         v4=pJC(vp,4)
      if (sC>1) then
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         NC(js,vp)=x(is,v1,v2,v3,v4)
         x(is,v1,v2,v3,v4)=0.0
      end do
      else
         NC(1,vp)=NC(1,vp)+x(1,v1,v2,v3,v4)
         x(1,v1,v2,v3,v4)=0.0
      endif
      end do

      end do

      end do

      deallocate (w,x)
      deallocate (tempNC2,tempJC2)
       
      return
      end

