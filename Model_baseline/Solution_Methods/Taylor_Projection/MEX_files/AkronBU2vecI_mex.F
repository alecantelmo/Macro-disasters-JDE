C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [NC,JC_1] = AkronBU2vecI_mex(l,m,n,IA,JA,NA,...
C            IB,JB,NB,...
C            IC_1,IC_2,JCrows_1,...
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
      mwPointer IC_1,IC_2
      mwPointer maxload,dosum
      mwPointer NA,NB

      mwPointer NC,JC_1

C     Array information:
      mwSize IArows,JArows
      mwSize IBrows, JBrows
      mwSize ICrows_1,ICrows_2
      mwSize JCrows_1

      mwSize s,sC,sA,sB
      real*8 do_sum

      mwSize rows
      mwSize JCdims(2),ndim

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 14) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '14 inputs required.')
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

      rows = mxGetScalar(prhs(12))

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
      JCrows_1 = mxGetPr(prhs(12))
      maxload = mxGetPr(prhs(13))

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
      do_sum=mxGetScalar(prhs(14))
 
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
     + %VAL(IC_1),%VAL(IC_2),
     + %VAL(JCrows_1),%VAL(maxload),
     + %VAL(NC),%VAL(JC_1),            
     + IArows,JArows,IBrows,JBrows,
     + ICrows_1,ICrows_2,
     + sC,sA,sB)
 
      return
      end

C     Computational subroutine
      subroutine compute(l,m,n,IA,JA,NA,
     + IB,JB,NB,
     + IC_1,IC_2,
     + JCrows_1,maxload,
     + NC,JC_1,
     + IArows,JArows,
     + IBrows,JBrows,
     + ICrows_1,ICrows_2,
     + sC,sA,sB)

C     Declarations
      implicit none

      mwSize  l,m,n,sC,sA,sB
      mwSize  IArows,JArows
      mwSize  IBrows,JBrows
      mwSize  ICrows_1,ICrows_2
      mwSize  JCrows_1,JCrows_2

      mwSize  s,sC2,rows
      mwSize  maxload,n_vecload,loadi,load_

      real*8, target ::  NA(sA,JArows)
      real*8, target ::  NB(sB,JBrows)

      real*8  NC(sC,JCrows_1)

      real*8, pointer :: pNA(:,:),pNB(:,:)
      mwSize, pointer :: pIA(:),pIB(:),pIC(:)
      mwSize, pointer :: pJA(:,:),pJC(:,:)
      mwSize, pointer :: pJB(:)

       mwSize, allocatable :: vecload(:)

      real*8, allocatable :: x(:,:,:) 
      mwSize, allocatable :: w(:,:)

C two temporary crs var is needed, 2
      real*8, target, allocatable :: tempNC2(:,:)
      mwSize, target, allocatable :: tempJC2(:,:)

      mwSize, target, allocatable :: tempJC1(:,:)

      mwSize, target ::  IA(IArows)
      mwSize, target ::  JA(JArows,2)

      mwSize, target ::  IB(IBrows)
      mwSize, target ::  JB(JBrows)

      mwSize, target ::  IC_1(ICrows_1)
      mwSize, target ::  JC_1(JCrows_1,1)

      mwSize, target ::  IC_2(ICrows_2)

      mwSize i,jp,j,kp,k1,k2,vp,v1,v2,v
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
      k=2
      mat(1,1:k)=(/ (mati, mati = n+1, n+k) /)
      mat(2,1:k)=(/ (mati, mati = n, n+k-1) /)

      fact(1)=1
      do j=2,k
       fact(j)=fact(j-1)*j
      end do
      prod1=product(mat(2,:))

C first contraction
      sC2=max(sA,sB)
      rows=IC_2(ICrows_2)-1
      allocate (w(m,n),
     +   x(min(sC2,maxload),m,n),
     +   tempNC2(sC2,rows),
     +   tempJC2(rows,2),
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

      pIC=>IC_2
      pJC=>tempJC2

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
        j=pJA(jp,2)
        k1=pJA(jp,1)
      do kp=pIB(j),pIB(j+1)-1
        k2=pJB(kp)
C collect information for pJC and initialize x to zero
        if (w(k1,k2)<i) then
         w(k1,k2)=i
         pJC(nz,1)=k1
         pJC(nz,2)=k2
         nz=nz+1
         x(:,k1,k2)=0
        endif

C Choose among three options:
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(1,jp)*pNB(is,kp)
       end do
      endif

      end do
      end do

C     assign the result in tempNC2
      do vp=pIC(i),pIC(i+1)-1
         v1=pJC(vp,1)
         v2=pJC(vp,2)
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         tempNC2(js,vp)=x(is,v1,v2)
         x(is,v1,v2)=0.0
      end do
      end do

      end do

      end do

      deallocate (w,x)

C second contraction
      rows=IC_1(ICrows_1)-1
      allocate (w(n,n),
     +   x(min(sC,maxload),n,n),
     +   tempJC1(JCrows_1,2),
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

      sA=sC2

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
      if ((pIB(j).le.(pIB(j+1)-1)) .and. (pJB(pIB(j+1)-1).ge.k2)) then
        kpstart=pIB(j)
      do while (pJB(kpstart)<k2)
        kpstart=kpstart+1
      end do
      do kp=kpstart,pIB(j+1)-1
        k1=pJB(kp)
C collect information for pJC and initialize x to zero
        if (w(k1,k2)<i) then
         w(k1,k2)=i
C determine the location of (k1,k2) in a compressed matrix and assign in pJC
         loc=(prod1-product(mat(1,:)-k2))/fact(k)+1
         loc=loc+k1-k2
         pJC(nz,1)=k1
         pJC(nz,2)=k2
         JC_1(nz,1)=loc
         nz=nz+1
         x(:,k1,k2)=0
        endif

      if (sC>1) then
C Choose among three options and do not sum across rows
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(is,k1,k2)=x(is,k1,k2)+pNA(1,jp)*pNB(is,kp)
       end do
      endif
      else
C Choose among three options and sum across rows.
      if (sA==sB) then
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2)=x(1,k1,k2)+pNA(is,jp)*pNB(is,kp)
       end do
      elseif (sA>sB) then
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2)=x(1,k1,k2)+pNA(is,jp)*pNB(1,kp)
       end do
      else
      !DEC$ simd
       do is=1,load_
        x(1,k1,k2)=x(1,k1,k2)+pNA(1,jp)*pNB(is,kp)
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
      if (sC>1) then
      !DEC$ simd
      do is=1,load_
         js=vecload(loadi)+is-1
         NC(js,vp)=x(is,v1,v2)
         x(is,v1,v2)=0.0
      end do
      else
         NC(1,vp)=NC(1,vp)+x(1,v1,v2)
         x(1,v1,v2)=0.0
      endif
      end do

      end do

      end do

      deallocate (w,x)
      deallocate (tempNC2,tempJC2,tempJC1)
      return
      end

