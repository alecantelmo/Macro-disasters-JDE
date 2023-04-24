C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"

C [ newcols ] = ptr2col_mex( ptr,n_vals )
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
      implicit none

      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxCreateNumericArray
      integer*4 mxGetClassID
      mwPointer mxGetM, mxGetN
      real*8 mxGetScalar

      mwPointer ptr
      mwPointer newcols

      mwSize ptrM

C     Arguments for mxCreateNumericArray
      integer*4 classid
      integer*4 complexflag
      mwSize ndim
      mwSize dims(1)
      mwSize l1,l2
      real*8 n_vals

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:nInput',
     +                           '2 inputs required.')
      endif
      if(nlhs .gt. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:nOutput',
     +                           '1 output required.')
      endif


C     Create Fortran array from the input arguments.

      ptr = mxGetPr(prhs(1))

      l1 = mxGetM(prhs(1))
      l2 = mxGetN(prhs(1))-1

C     Create matrix for the return arguments.

      n_vals=mxGetScalar(prhs(2))

      complexflag = 0

      classid = mxGetClassID(prhs(1))
      ndim = 1
      dims(1) = n_vals
      plhs(1) = mxCreateNumericArray(ndim, dims, 
     +          classid, 
     +          complexflag)

      newcols = mxGetPr(plhs(1))

C     Call the computational routine.
      call compute(%VAL(ptr),l1,l2,%VAL(newcols),n_vals)

      return
      end


      subroutine compute(ptr,l1,l2,newcols,n_vals)
      implicit none

      real*8 n_vals
      mwSize l1,l2
      mwSize ptr(l1,l2+1),newcols(n_vals)

      mwSize i1,i2,i
      do i1=1,l1
       do i2=1,l2
        do i=ptr(i1,i2),ptr(i1,i2+1)-1
         newcols(i)=i2
        enddo
       enddo
      enddo

      return
      end

