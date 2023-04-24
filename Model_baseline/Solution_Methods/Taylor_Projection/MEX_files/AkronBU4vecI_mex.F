C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [NC,JC_1] = AkronBU4vecI_mex(l,m,n,IA,JA,NA,...
C            IB,JB,NB,...
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
      mwPointer m,n
      mwPointer IB,JB
      mwPointer IC_1,IC_2,IC_3,IC_4
      mwPointer maxload,dosum
      mwPointer NA,NB

      mwPointer NC,JC_1

C     Array information:
      mwSize IArows,JArows
      mwSize IBrows, JBrows
      mwSize ICrows_1,ICrows_2,ICrows_3,ICrows_4
      mwSize JCrows_1

      mwSize s,sC,sA,sB
      real*8 do_sum

      mwSize rows
      mwSize JCdims(2),ndim

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 16) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '16 inputs required.')
      endif
      if(nlhs .gt. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 2 outputs required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(4))
      JArows = mxGetM(prhs(5))

      IBrows = mxGetM(prhs(7))
      JBrows = mxGetM(prhs(8))

      ICrows_1 = mxGetM(prhs(10))
      ICrows_2 = mxGetM(prhs(11))
      ICrows_3 = mxGetM(prhs(12))
      ICrows_4 = mxGetM(prhs(13))

      rows = mxGetScalar(prhs(14))

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
      IC_2 = mxGetPr(prhs(11))
      IC_3 = mxGetPr(prhs(12))
      IC_4 = mxGetPr(prhs(13))
      JCrows_1 = mxGetPr(prhs(14))
      maxload = mxGetPr(prhs(15))

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
      do_sum=mxGetScalar(prhs(16))
 
      if (do_sum==1) then
       sC=1
      else
       sC=s
      endif

C     Create matrix for NC.
      complexflag = 0

      plhs(1) = mxCreateDoubleMatrix(sC,rows,complexflag)
      NC = mxGetPr(plhs(1))

C     Create matrix for JC_1.

      classid = mxGetClassID(prhs(5))

      JCdims(1)=rows
      JCdims(2)=1
      ndim=2
      plhs(2) = mxCreateNumericArray(ndim, JCdims, 
     +          classid, 
     +          complexflag)

      JC_1 = mxGetPr(plhs(2))
  
C     Call the computational routine.
      call compute(%VAL(l),%VAL(m),
     + %VAL(n),
     + %VAL(IA),%VAL(JA),%VAL(NA),
     + %VAL(IB),%VAL(JB),%VAL(NB),
     + %VAL(IC_1),%VAL(IC_2),%VAL(IC_3),%VAL(IC_4),
     + %VAL(JCrows_1),%VAL(maxload),
     + %VAL(NC),%VAL(JC_1),            
     + IArows,JArows,IBrows,JBrows,
     + ICrows_1,ICrows_2,ICrows_3,ICrows_4,
     + sC,sA,sB)
 
      return
      end

C     Computational subroutine
      subroutine compute(l,m,n,IA,JA,NA,
     + IB,JB,NB,
     + IC_1,IC_2,IC_3,IC_4,
     + JCrows_1,maxload,
     + NC,JC_1,
     + IArows,JArows,
     + IBrows,JBrows,
     + ICrows_1,ICrows_2,ICrows_3,ICrows_4,
     + sC,sA,sB)

C     Declarations
      implicit none

      mwSize  l,m,n,sC,sA,sB
      mwSize  IArows,JArows
      mwSize  IBrows,JBrows
      mwSize  ICrows_1,ICrows_2,ICrows_3,ICrows_4
      mwSize  JCrows_1,JCrows_2,JCrows_3,JCrows_4

      mwSize  s,sC4,sC3,sC2,rows
      mwSize  maxload,n_vecload,loadi,load_

      real*8, target ::  NA(sA,JArows)
      real*8, target ::  NB(sB,JBrows)

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

      mwSize, target, allocatable :: tempJC1(:,:)

      mwSize, target ::  IA(IArows)
      mwSize, target ::  JA(JArows,4)

      mwSize, target ::  IB(IBrows)
      mwSize, target ::  JB(JBrows)

      mwSize, target ::  IC_1(ICrows_1)
      mwSize, target ::  JC_1(JCrows_1,1)

      mwSize, target ::  IC_2(ICrows_2)
      mwSize, target ::  IC_3(ICrows_3)
      mwSize, target ::  IC_4(ICrows_4)

      mwSize i,jp,j,kp,k1,k2,k3,k4,vp,v1,v2,v3,v4
      mwSize is,js,nz

      mwSize mati,k,kpstart,loc,j2
      parameter (k=4)
      real*8 mat(k,k), fact(k), prod1
 
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
      sC4=max(sA,sB)
      rows=IC_4(ICrows_4)-1
      allocate (w(m,m,m,n),
     +   x(min(sC4,maxload),m,m,m,n),
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

      pIB=>IB  
      pJB=>JB

      pIC=>IC_4
      pJC=>tempJC4

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
        k3=pJA(jp,3)
        j=pJA(jp,4)
      do kp=pIB(j),pIB(j+1)-1
        k4=JB(kp)

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
      sC3=max(sC4,sB)
      rows=IC_3(ICrows_3)-1
      allocate (w(m,m,n,n),
     +   x(min(sC3,maxload),m,m,n,n),
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

      pIB=>IB        
      pJB=>JB

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
        k4=pJA(jp,4)
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (JB(pIB(j+1)-1).ge.k4)) then
        kpstart=pIB(j)
      do while (JB(kpstart)<k4)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k3=JB(kp)
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
      endif
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
      sC2=max(sC3,sB)
      rows=IC_2(ICrows_2)-1
      allocate (w(m,n,n,n),
     +   x(min(sC2,maxload),m,n,n,n),
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

      pIB=>IB
      pJB=>JB

      pIC=>IC_2
      pJC=>tempJC2

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
        k4=pJA(jp,4)
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (JB(pIB(j+1)-1).ge.k3)) then
        kpstart=pIB(j)
      do while (JB(kpstart)<k3)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k2=JB(kp)
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
      endif
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
      allocate (w(n,n,n,n),
     +   x(min(sC,maxload),n,n,n,n),
     +   tempJC1(JCrows_1,4),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:OutOfMemory',
     +                           'Out of memory4.')
      endif

      w=0
      nz=1

      pIA=>IC_2
      pJA=>tempJC2

      pIB=>IB
      pJB=>JB

      pIC=>IC_1
      pJC=>tempJC1

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
        k4=pJA(jp,4)
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (JB(pIB(j+1)-1).ge.k2)) then
        kpstart=pIB(j)
      do while (JB(kpstart)<k2)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k1=JB(kp)
C collect information for pJC and initializes x to zero
        if (w(k1,k2,k3,k4)<i) then
         w(k1,k2,k3,k4)=i
C determine the location of (k1,k2,k3,k4) in a compressed matrix and assign in JC_1
         loc=(prod1-product(mat(1,:)-k4))/fact(k)+1
         loc=loc+(product(mat(3,1:k-1)-k4)-product(mat(3,1:k-1)-k3))/
     +       fact(k-1)
         loc=loc+(product(mat(4,1:k-2)-k3)-product(mat(4,1:k-2)-k2))/
     +       fact(k-2)
         loc=loc+k1-k2
         pJC(nz,1)=k1
         pJC(nz,2)=k2
         pJC(nz,3)=k3
         pJC(nz,4)=k4
         JC_1(nz,1)=loc
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
      endif
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
      deallocate (tempNC2,tempJC2,tempJC1)
       
      return
      end

