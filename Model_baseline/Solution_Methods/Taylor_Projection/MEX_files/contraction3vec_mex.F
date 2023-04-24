C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [NC] = contraction3vec_mex(l,m1,m2,n1,n2,n3,IA,JA,NA,...
C            IB1,JB1,NB1,IB2,JB2,NB2,IB3,JB3,NB3,...
C            IC_1,JC_1,IC_2,JC_2,IC_3,JC_3,...
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
      real*8 mxGetScalar
      integer*4 complexflag

C     Pointers to input/output mxArrays:
      mwPointer l,IA,JA
      mwPointer m1,m2,n1,n2,n3
      mwPointer IB1,JB1,IB2,JB2,IB3,JB3
      mwPointer IC_1,IC_2,IC_3
      mwPointer JC_1,JC_2,JC_3
      mwPointer maxload,dosum
      mwPointer NA,NB1,NB2,NB3

      mwPointer NC

C     Array information:
      mwSize IArows,JArows
      mwSize IB1rows, JB1rows
      mwSize IB2rows, JB2rows
      mwSize IB3rows, JB3rows
      mwSize ICrows_1,ICrows_2,ICrows_3
      mwSize JCrows_1,JCrows_2,JCrows_3

      mwSize s,sC,sA,sB1,sB2,sB3
      real*8 do_sum

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 26) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '26 inputs required.')
      endif
      if(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 1 output required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(7))
      JArows = mxGetM(prhs(8))

      IB1rows = mxGetM(prhs(10))
      JB1rows = mxGetM(prhs(11))
      IB2rows = mxGetM(prhs(13))
      JB2rows = mxGetM(prhs(14))
      IB3rows = mxGetM(prhs(16))
      JB3rows = mxGetM(prhs(17))

      ICrows_1 = mxGetM(prhs(19))
      JCrows_1 = mxGetM(prhs(20))

      ICrows_2 = mxGetM(prhs(21))
      JCrows_2 = mxGetM(prhs(22))

      ICrows_3 = mxGetM(prhs(23))
      JCrows_3 = mxGetM(prhs(24))

C     Create Fortran array from the input arguments.
      l = mxGetPr(prhs(1))
      m1 = mxGetPr(prhs(2))
      m2 = mxGetPr(prhs(3))
      n1 = mxGetPr(prhs(4))
      n2 = mxGetPr(prhs(5))
      n3 = mxGetPr(prhs(6))
      IA = mxGetPr(prhs(7))
      JA = mxGetPr(prhs(8))
      NA = mxGetPr(prhs(9))
      IB1 = mxGetPr(prhs(10))
      JB1 = mxGetPr(prhs(11))
      NB1 = mxGetPr(prhs(12))
      IB2 = mxGetPr(prhs(13))
      JB2 = mxGetPr(prhs(14))
      NB2 = mxGetPr(prhs(15))
      IB3 = mxGetPr(prhs(16))
      JB3 = mxGetPr(prhs(17))
      NB3 = mxGetPr(prhs(18))
      IC_1 = mxGetPr(prhs(19))
      JC_1 = mxGetPr(prhs(20))
      IC_2 = mxGetPr(prhs(21))
      JC_2 = mxGetPr(prhs(22))
      IC_3 = mxGetPr(prhs(23))
      JC_3 = mxGetPr(prhs(24))
      maxload = mxGetPr(prhs(25))

C     Get number of states
      sA = mxGetM(prhs(9))
      sB1 = mxGetM(prhs(12))
      sB2 = mxGetM(prhs(15))
      sB3 = mxGetM(prhs(18))

C check compatibility of states
      s=max(sA,sB1,sB2,sB3)
      if ((sA.ne.s .and. sA.ne.1) .or. (sB1.ne.s .and. sB1.ne.1)
     +   .or. (sB1.ne.sB2) .or. (sB1.ne.sB3)) then
         call mexErrMsgIdAndTxt ('MATLAB:states',
     +   'incompatible states.')
      endif

C     sum across rows of C if dosum==1
      do_sum=mxGetScalar(prhs(26))
 
      if (do_sum==1) then
       sC=1
      else
       sC=s
      endif

C     Create matrix for NC.
      complexflag=0
      plhs(1) = mxCreateDoubleMatrix(sC,JCrows_1, complexflag)
      NC = mxGetPr(plhs(1))

C     Call the computational routine.
      call compute(%VAL(l),%VAL(m1),%VAL(m2),
     + %VAL(n1),%VAL(n2),%VAL(n3),
     + %VAL(IA),%VAL(JA),%VAL(NA),
     + %VAL(IB1),%VAL(JB1),%VAL(NB1),
     + %VAL(IB2),%VAL(JB2),%VAL(NB2),
     + %VAL(IB3),%VAL(JB3),%VAL(NB3),
     + %VAL(IC_1),%VAL(JC_1),%VAL(IC_2),%VAL(JC_2),
     + %VAL(IC_3),%VAL(JC_3),%VAL(maxload),
     + %VAL(NC),          
     + IArows,JArows,IB1rows,JB1rows,IB2rows,JB2rows,
     + IB3rows,JB3rows,
     + ICrows_1,JCrows_1,ICrows_2,JCrows_2,
     + ICrows_3,JCrows_3,
     + sC,sA,sB1,sB2,sB3)
 
      return
      end

C     Computational subroutine
      subroutine compute(l,m1,m2,n1,n2,n3,IA,JA,NA,
     + IB1,JB1,NB1,IB2,JB2,NB2,IB3,JB3,NB3,
     + IC_1,JC_1,IC_2,JC_2,IC_3,JC_3,
     + maxload,
     + NC,
     + IArows,JArows,
     + IB1rows,JB1rows,IB2rows,JB2rows,
     + IB3rows,JB3rows,
     + ICrows_1,JCrows_1,ICrows_2,JCrows_2,
     + ICrows_3,JCrows_3,
     + sC,sA,sB1,sB2,sB3)


C     Declarations
      implicit none

      mwSize  l,m1,m2,n1,n2,n3,sC,sA,sB1,sB2,sB3
      mwSize  IArows,JArows
      mwSize  IB1rows,JB1rows,IB2rows,JB2rows
      mwSize  IB3rows,JB3rows
      mwSize  ICrows_1,ICrows_2,ICrows_3
      mwSize  JCrows_1,JCrows_2,JCrows_3

      mwSize  s,sB,sC3,sC2,rows
      mwSize  maxload,n_vecload,loadi,load_

      real*8, target ::  NA(sA,JArows)
      real*8, target ::  NB1(sB1,JB1rows)
      real*8, target ::  NB2(sB2,JB2rows)
      real*8, target ::  NB3(sB3,JB3rows)

      real*8  NC(sC,JCrows_1)

      real*8, pointer :: pNA(:,:),pNB(:,:)
      mwSize, pointer :: pIA(:),pIB(:),pIC(:)
      mwSize, pointer :: pJA(:,:),pJC(:,:)
      mwSize, pointer :: pJB(:)

       mwSize, allocatable :: vecload(:)

      real*8, allocatable :: x(:,:,:,:) 
      logical, allocatable :: w(:,:,:)

C two temporary crs vars are needed, 2 and 3
      real*8, target, allocatable :: tempNC3(:,:)
      real*8, target, allocatable :: tempNC2(:,:)

      mwSize, target ::  IA(IArows)
      mwSize, target ::  JA(JArows,3)

      mwSize, target ::  IB1(IB1rows)
      mwSize, target ::  JB1(JB1rows)
      mwSize, target ::  IB2(IB2rows)
      mwSize, target ::  JB2(JB2rows)
      mwSize, target ::  IB3(IB3rows)
      mwSize, target ::  JB3(JB3rows)

      mwSize, target ::  IC_1(ICrows_1)
      mwSize, target ::  JC_1(JCrows_1,3)

      mwSize, target ::  IC_2(ICrows_2)
      mwSize, target ::  JC_2(JCrows_2,3)

      mwSize, target ::  IC_3(ICrows_3)
      mwSize, target ::  JC_3(JCrows_3,3)

      mwSize i,jp,j,kp,k1,k2,k3,vp,v1,v2,v3
      mwSize is,js,nz

      mwSize jpstart, jpend
      mwSize kpstart, kpend
      mwSize vpstart, vpend
 
      integer AllocateStatus

C vecload
      s=max(sA,sB1,sB2,sB3)

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

      sB=sB1 ! sB1, sB2 and sB3 must all be equal

C first contraction
      sC3=max(sA,sB3)
      rows=IC_3(ICrows_3)-1
      allocate (w(m1,m2,n3),
     +   x(min(sC3,maxload),m1,m2,n3),
     +   tempNC3(sC3,rows),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory2.')
      endif

      w=.false.
      nz=1

      pIA=>IA
      pJA=>JA

      pIB=>IB3        
      pJB=>JB3

      pIC=>IC_3
      pJC=>JC_3

      do loadi=1,n_vecload-1
        load_=vecload(loadi+1)-vecload(loadi)
       
      if (sA>1) then
        pNA=>NA(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNA=>NA
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
      do kp=pIB(j),pIB(j+1)-1
        k3=JB3(kp)
C initialize x to zero
      if (w(k1,k2,k3)==.false.) then
        w(k1,k2,k3)=.true.
        x(:,k1,k2,k3)=0
      endif 

C Choose among three options:
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3)=x(is,k1,k2,k3)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3)=x(is,k1,k2,k3)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3)=x(is,k1,k2,k3)+pNA(1,jp)*pNB(is,kp)
       end do
      endif

      end do
      end do

C     assign the result in tempNC3
      do vp=pIC(i),pIC(i+1)-1
         v1=pJC(vp,1)
         v2=pJC(vp,2)
         v3=pJC(vp,3)
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         tempNC3(js,vp)=x(is,v1,v2,v3)
         x(is,v1,v2,v3)=0.0
      end do
      end do

      end do

      end do

      deallocate (w,x)

C second contraction
      sC2=max(sC3,sB2)
      rows=IC_2(ICrows_2)-1
      allocate (w(m1,n2,n3),
     +   x(min(sC2,maxload),m1,n2,n3),
     +   tempNC2(sC2,rows),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory3.')
      endif

      w=.false.
      nz=1

      pIA=>IC_3
      pJA=>JC_3

      pIB=>IB2        
      pJB=>JB2

      pIC=>IC_2
      pJC=>JC_2

      sA=sC3

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
      do kp=pIB(j),pIB(j+1)-1
        k2=JB2(kp)
C initialize x to zero
      if (w(k1,k2,k3)==.false.) then
        w(k1,k2,k3)=.true.
        x(:,k1,k2,k3)=0
      endif 

C Choose among three options:
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3)=x(is,k1,k2,k3)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3)=x(is,k1,k2,k3)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3)=x(is,k1,k2,k3)+pNA(1,jp)*pNB(is,kp)
       end do
      endif

      end do
      end do

C     assign the result in tempNC2
      do vp=pIC(i),pIC(i+1)-1
         v1=pJC(vp,1)
         v2=pJC(vp,2)
         v3=pJC(vp,3)
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         tempNC2(js,vp)=x(is,v1,v2,v3)
         x(is,v1,v2,v3)=0.0
      end do
      end do

      end do

      end do

      deallocate (w,x)
      deallocate (tempNC3)

C third contraction
      rows=IC_1(ICrows_1)-1
      allocate (w(n1,n2,n3),
     +   x(min(sC,maxload),n1,n2,n3),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory4.')
      endif

      w=.false.
      nz=1

      pIA=>IC_2
      pJA=>JC_2

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
      do kp=pIB(j),pIB(j+1)-1
        k1=JB1(kp)
C initialize x to zero
      if (w(k1,k2,k3)==.false.) then
        w(k1,k2,k3)=.true.
        x(:,k1,k2,k3)=0
      endif 

      if (sC>1) then
C Choose among three options and do not sum across rows
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3)=x(is,k1,k2,k3)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3)=x(is,k1,k2,k3)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2,k3)=x(is,k1,k2,k3)+pNA(1,jp)*pNB(is,kp)
       end do
      endif
      else
C Choose among three options and sum across rows.
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2,k3)=x(1,k1,k2,k3)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2,k3)=x(1,k1,k2,k3)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2,k3)=x(1,k1,k2,k3)+pNA(1,jp)*pNB(is,kp)
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
      if (sC>1) then
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         NC(js,vp)=x(is,v1,v2,v3)
         x(is,v1,v2,v3)=0.0
      end do
      else
         NC(1,vp)=NC(1,vp)+x(1,v1,v2,v3)
         x(1,v1,v2,v3)=0.0
      endif
      end do

      end do

      end do

      deallocate (w,x)
      deallocate (tempNC2)

      return
      end

