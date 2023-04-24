C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [IC,JC,NC] = takerows_mex(IA,JA,NA,l,takel)
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
      mwPointer l,takel,IA,JA,NA
      mwPointer IC,JC,NC

C     Array information:
      mwSize IArows,JArows,JAcols,NArows,NAcols

      mwSize sA,dims(2),ndims,n_takel


C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 5) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '5 inputs required.')
      endif
      if(nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 3 outputs required.')
      endif

C     Get the size of the inputs array.

      IArows = mxGetM(prhs(1))
      JArows = mxGetM(prhs(2))
      NArows = mxGetM(prhs(3))

      JAcols = mxGetN(prhs(2))
      NAcols = mxGetN(prhs(3))

      n_takel=mxGetM(prhs(5))*mxGetN(prhs(5))

C     Create Fortran array from the input arguments.
      IA = mxGetPr(prhs(1))
      JA = mxGetPr(prhs(2))
      NA = mxGetPr(prhs(3))
      l = mxGetPr(prhs(4))
      takel = mxGetPr(prhs(5))

C     Create return arguments.

      classid = mxGetClassID(prhs(1))
      complexflag = 0
      ndims=1
      plhs(1) = mxCreateNumericArray(ndims, IArows, 
     +          classid, 
     +          complexflag)

      IC = mxGetPr(plhs(1))

      classid = mxGetClassID(prhs(2))
      ndims=2
      dims(1)=JArows
      dims(2)=JAcols
      plhs(2) = mxCreateNumericArray(ndims, dims, 
     +          classid, 
     +          complexflag)

      JC = mxGetPr(plhs(2))

      plhs(3) = mxCreateDoubleMatrix(NArows,NAcols,complexflag)
      NC = mxGetPr(plhs(3))
 
C     Call the computational routine.
      call compute(%VAL(IA),%VAL(JA),%VAL(NA),
     + %VAL(l),%VAL(takel),
     + %VAL(IC),%VAL(JC),%VAL(NC),
     + IArows,JArows,JAcols,NArows,NAcols,n_takel)
 
      return
      end

C     Computational subroutine
      subroutine compute(IA,JA,NA,l,takel,
     + IC,JC,NC,
     + IArows,JArows,JAcols,NArows,NAcols,n_takel)

C     Declarations
      implicit none

      mwSize  IArows,JArows,JAcols,NArows,NAcols,n_takel
      mwSize  IA(IArows),IC(IArows)
      mwSize  JA(JArows,JAcols),JC(JArows,JAcols)
      real*8  NA(NArows,NAcols),NC(NArows,NAcols)

      mwSize  l,takel(n_takel),i,j,k,nz,startc,endc
 
C create IC,JC,NC
      nz=0
      k=0
      IC(1)=1
      do i=1,n_takel
       k=k+1
       j=takel(i)
       IC(k+1)=IC(k)+IA(j+1)-IA(j)
       JC(IC(k):IC(k+1)-1,:)=JA(IA(j):IA(j+1)-1,:)
       NC(:,IC(k):IC(k+1)-1)=NA(:,IA(j):IA(j+1)-1)

      enddo
 
      return
      end

