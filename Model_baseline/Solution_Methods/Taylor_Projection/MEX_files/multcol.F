C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [A] = multcol(mat,vec)
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
      implicit none

C     mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

C     Function declarations:
      mwPointer mxGetPr
      mwPointer mxGetM, mxGetN
      mwPointer mxCreateDoubleMatrix
      integer*4 complexflag

C     Pointers to input/output mxArrays:
      mwPointer mat,vec,A

C     Array information:
      mwSize matM,matN
      mwSize vecM,vecN


C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction2:nInput',
     +                           '2 inputs required.')
      endif
      if(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:contraction1:nOutput',
     +                           'Max 1 output required.')
      endif

C     Get the size of the inputs array.

      matM = mxGetM(prhs(1))
      matN = mxGetN(prhs(1))

      vecM = mxGetM(prhs(2))
      vecN = mxGetN(prhs(2))

C     Create Fortran array from the input arguments.
      mat = mxGetPr(prhs(1))
      vec = mxGetPr(prhs(2))

C     Create matrix for the return arguments.
      complexflag=0

      plhs(1) = mxCreateDoubleMatrix(matM,matN,complexflag)
      A = mxGetPr(plhs(1))

 
C     Call the computational routine.
      call compute(%VAL(mat),%VAL(vec),%VAL(A),
     + matM,matN,vecM,vecN)
 
      return
      end

C     Computational subroutine
      subroutine compute(mat,vec,A,
     + matM,matN,vecM,vecN)

C     Declarations
      implicit none

      mwSize  matM,matN,vecM,vecN,i,j
      real*8 mat(matM,matN),A(matM,matN),vec(vecM,vecN)

      if ((vecM>1) .and. (vecN>1)) then
         call mexErrMsgIdAndTxt ('MATLAB:multrow',
     +   'second argument is not a vector.')
      elseif ((vecM>1) .and. (matM.ne.vecM)) then
         call mexErrMsgIdAndTxt ('MATLAB:multrow',
     +   'incompatible dimensions.')
      elseif ((vecN>1) .and. (matM.ne.vecN)) then
         call mexErrMsgIdAndTxt ('MATLAB:multrow',
     +   'incompatible dimensions.')
      endif

      if (vecN==1) then
        do j=1,matN
         !DEC$ simd
         do i=1,matM
          A(i,j)=mat(i,j)*vec(i,1)
         end do
        end do
      else
        do j=1,matN
         !DEC$ simd
         do i=1,matM
          A(i,j)=mat(i,j)*vec(1,i)
         end do
        end do
      endif

      return
      end

