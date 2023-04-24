C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [NC] = AkronBU3vec_mex(l,m,n,IA,JA,NA,...
C            IB,JB,NB,...
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
      real*8 mxGetScalar
      integer*4 complexflag

C     Pointers to input/output mxArrays:
      mwPointer l,IA,JA
      mwPointer m,n
      mwPointer IB,JB
      mwPointer IC_1,IC_2,IC_3
      mwPointer JC_1,JC_2,JC_3
      mwPointer maxload,dosum
      mwPointer NA,NB

      mwPointer NC

C     Array information:
      mwSize IArows,JArows
      mwSize IBrows, JBrows
      mwSize ICrows_1,ICrows_2,ICrows_3
      mwSize JCrows_1,JCrows_2,JCrows_3

      mwSize s,sC,sA,sB
      real*8 do_sum

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 17) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '17 inputs required.')
      endif
      if(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 1 output required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(4))
      JArows = mxGetM(prhs(5))

      IBrows = mxGetM(prhs(7))
      JBrows = mxGetM(prhs(8))

      ICrows_1 = mxGetM(prhs(10))
      JCrows_1 = mxGetM(prhs(11))

      ICrows_2 = mxGetM(prhs(12))
      JCrows_2 = mxGetM(prhs(13))

      ICrows_3 = mxGetM(prhs(14))
      JCrows_3 = mxGetM(prhs(15))

C     Create Fortran array from the input arguments.
      l = mxGetPr(prhs(1))
      m = mxGetPr(prhs(2))
      n = mxGetPr(prhs(3))
      IA = mxGetPr(prhs(4))
      JA = mxGetPr(prhs(5))
      NA = mxGetPr(prhs(6))
      IB = mxGetPr(prhs(7))
      JB = mxGetPr(prhs(8))
      NB = mxGetPr(prhs(9))
      IC_1 = mxGetPr(prhs(10))
      JC_1 = mxGetPr(prhs(11))
      IC_2 = mxGetPr(prhs(12))
      JC_2 = mxGetPr(prhs(13))
      IC_3 = mxGetPr(prhs(14))
      JC_3 = mxGetPr(prhs(15))
      maxload = mxGetPr(prhs(16))

C     Get number of states
      sA = mxGetM(prhs(6))
      sB = mxGetM(prhs(9))

C check compatibility of states
      s=max(sA,sB)
      if ((sA.ne.s .and. sA.ne.1) .or. (sB.ne.s .and. sB.ne.1))
     +   then
         call mexErrMsgIdAndTxt ('MATLAB:states',
     +   'incompatible states.')
      endif

C     sum across rows of C if dosum==1
      do_sum=mxGetScalar(prhs(17))
 
      if (do_sum==1) then
       sC=1
      else
       sC=s
      endif

C     Create matrix for NC.
      complexflag=0
      plhs(1) = mxCreateDoubleMatrix(sC,JCrows_1,complexflag)
      NC = mxGetPr(plhs(1))

C     Call the computational routine.
      call compute(%VAL(l),%VAL(m),
     + %VAL(n),
     + %VAL(IA),%VAL(JA),%VAL(NA),
     + %VAL(IB),%VAL(JB),%VAL(NB),
     + %VAL(IC_1),%VAL(JC_1),%VAL(IC_2),%VAL(JC_2),
     + %VAL(IC_3),%VAL(JC_3),%VAL(maxload),
     + %VAL(NC),          
     + IArows,JArows,IBrows,JBrows,
     + ICrows_1,JCrows_1,ICrows_2,JCrows_2,
     + ICrows_3,JCrows_3,
     + sC,sA,sB)
 
      return
      end

C     Computational subroutine
      subroutine compute(l,m,n,IA,JA,NA,
     + IB,JB,NB,
     + IC_1,JC_1,IC_2,JC_2,IC_3,JC_3,
     + maxload,
     + NC,
     + IArows,JArows,
     + IBrows,JBrows,
     + ICrows_1,JCrows_1,ICrows_2,JCrows_2,
     + ICrows_3,JCrows_3,
     + sC,sA,sB)


C     Declarations
      implicit none

      mwSize  l,m,n,sC,sA,sB
      mwSize  IArows,JArows
      mwSize  IBrows,JBrows
      mwSize  ICrows_1,ICrows_2,ICrows_3
      mwSize  JCrows_1,JCrows_2,JCrows_3

      mwSize  s,sC3,sC2,rows
      mwSize  maxload,n_vecload,loadi,load_

      real*8, target ::  NA(sA,JArows)
      real*8, target ::  NB(sB,JBrows)

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

      mwSize, target ::  IB(IBrows)
      mwSize, target ::  JB(JBrows)

      mwSize, target ::  IC_1(ICrows_1)
      mwSize, target ::  JC_1(JCrows_1,4)

      mwSize, target ::  IC_2(ICrows_2)
      mwSize, target ::  JC_2(JCrows_2,3)

      mwSize, target ::  IC_3(ICrows_3)
      mwSize, target ::  JC_3(JCrows_3,3)

      mwSize i,jp,j,kp,k1,k2,k3,vp,v1,v2,v3
      mwSize is,js,nz

      mwSize mati,k,kpstart,loc,j2
      real*8 mat(2,2), fact(2), prod1

 
      integer AllocateStatus

C vecload
      s=max(sA,sB)

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
      sC3=max(sA,sB)
      rows=IC_3(ICrows_3)-1
      allocate (w(m,m,n),
     +   x(min(sC3,maxload),m,m,n),
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

      pIB=>IB       
      pJB=>JB

      pIC=>IC_3
      pJC=>JC_3

      do loadi=1,n_vecload-1
        load_=vecload(loadi+1)-vecload(loadi)
       
      if (sA>1) then
        pNA=>NA(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNA=>NA
      endif

      if (sB>1) then
        pNB=>NB(vecload(loadi):vecload(loadi+1)-1,:)       
      else 
        pNB=>NB
      endif

      do i=1,l
      do jp=pIA(i),pIA(i+1)-1
        k1=pJA(jp,1)
        k2=pJA(jp,2)
        j=pJA(jp,3)
      do kp=pIB(j),pIB(j+1)-1
        k3=pJB(kp)
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
      sC2=max(sC3,sB)
      rows=IC_2(ICrows_2)-1
      allocate (w(m,n,n),
     +   x(min(sC2,maxload),m,n,n),
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

      pIB=>IB        
      pJB=>JB

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

      if (sB>1) then
        pNB=>NB(vecload(loadi):vecload(loadi+1)-1,:)       
      else 
        pNB=>NB
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
        k2=pJB(kp)
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
      endif
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
      allocate (w(n,n,n),
     +   x(min(sC,maxload),n,n,n),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory4.')
      endif

      w=.false.
      nz=1

      pIA=>IC_2
      pJA=>JC_2

      pIB=>IB        
      pJB=>JB

      pIC=>IC_1
      pJC=>JC_1

      do loadi=1,n_vecload-1
        load_=vecload(loadi+1)-vecload(loadi)
       
      if (sC2>1) then
        pNA=>tempNC2(vecload(loadi):vecload(loadi+1)-1,:)
      else
        pNA=>tempNC2
      endif

      if (sB>1) then
        pNB=>NB(vecload(loadi):vecload(loadi+1)-1,:)       
      else 
        pNB=>NB
      endif

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
        k1=pJB(kp)
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
      endif
      end do

C     assign the result in NC
      do vp=pIC(i),pIC(i+1)-1
         v1=pJC(vp,2)
         v2=pJC(vp,3)
         v3=pJC(vp,4)
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

