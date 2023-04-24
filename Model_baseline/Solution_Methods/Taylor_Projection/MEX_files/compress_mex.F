C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"

C [ newvals,ptr,newcols ] = compress( vals,row,cols,rowdim,n_vals )
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
      implicit none

      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxCreateNumericArray
      integer*4 mxClassIDFromClassName
      mwPointer mxGetM, mxGetN
      real*8 mxGetScalar
      integer*4 mexPrintf

      mwPointer vals,row,cols
      mwPointer rowdim, n_vals

      mwPointer newvals,ptr,newcols
      mwSize n_ptr

C     Arguments for mxCreateNumericArray
      integer*4 classid
      integer*4 complexflag
      mwSize ndim
      mwSize dims(2)
      
      mwSize valsM, valsN, colsM, colsN

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 5) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:nInput',
     +                           '5 inputs required.')
      endif
      if(nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:nOutput',
     +                           '3 outputs required.')
      endif


C     Create Fortran array from the input arguments.
      vals = mxGetPr(prhs(1))
      row = mxGetPr(prhs(2))
      cols = mxGetPr(prhs(3))
      rowdim = mxGetPr(prhs(4))
      n_vals = mxGetPr(prhs(5))

      valsM = mxGetM(prhs(1))
      valsN = mxGetN(prhs(1))

      colsM = mxGetM(prhs(3))
      colsN = mxGetN(prhs(3))

C     Create matrix for the return arguments.
      plhs(1) = mxCreateDoubleMatrix(valsM, valsN, 0)

      n_ptr=mxGetScalar(prhs(4))+1

      if (sizeof (n_ptr) == 4 ) then
         classid = mxClassIDFromClassName('int32')
      else if (sizeof (n_ptr) == 8 ) then
         classid = mxClassIDFromClassName('int64')
      else
          call mexErrMsgIdAndTxt ('MATLAB:compress:OutofMemory',
     +                           'failed')
      end if


      complexflag = 0
      plhs(2) = mxCreateNumericArray(1, n_ptr, 
     +          classid, 
     +          complexflag)

      ndim = 2
      dims(1) = colsM
      dims(2) = colsN
      plhs(3) = mxCreateNumericArray(ndim, dims, 
     +          classid, 
     +          complexflag)


      newvals = mxGetPr(plhs(1))
      ptr = mxGetPr(plhs(2))
      newcols = mxGetPr(plhs(3))

C     Call the computational routine.
      call compute(%VAL(vals),%VAL(row),%VAL(cols),
     +             %VAL(n_vals),%VAL(rowdim),
     +             valsM,valsN,colsM,colsN,
     +             %VAL(newvals),%VAL(ptr),%VAL(newcols))



      return
      end


      subroutine compute(vals,row,cols, n_vals, rowdim, 
     +      valsM,valsN,colsM,colsN,
     +      newvals,ptr,newcols)
      implicit none

      mwSize i, n, size, valsM, valsN, colsM, colsN
      mwSize n_vals, rowdim

      real*8 vals(valsM,valsN), newvals(valsM,valsN)
      mwSize cols(colsM,colsN), newcols(colsM,colsN)
      mwSize row(n_vals), ptr(rowdim+1)



      mwSize, allocatable :: countrow(:), z(:), loc(:), ind(:)

      integer AllocateStatus

      allocate (countrow(rowdim),z(rowdim+1),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:OutofMemory',
     +                           'Out of memory.')
      endif

      call countdim(countrow, row, n_vals, rowdim)


      z(1)=1
      z(2:rowdim+1)=countrow

      deallocate(countrow)

      call cumsum(ptr,z,rowdim+1)

      deallocate(z)
      allocate (loc(rowdim), ind(n_vals), 
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:OutofMemory',
     +                           'Out of memory.')
      endif
       
      loc=ptr(1:rowdim)
      
      do i=1,n_vals
        ind(i)=loc(row(i))
        loc(row(i))=loc(row(i))+1
      end do
      
      do i=1,n_vals
        newvals(1:valsM,ind(i))=vals(1:valsM,i)
        newcols(ind(i),1:colsN)=cols(i,1:colsN)
      end do

      deallocate(loc,ind,STAT = AllocateStatus)
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:OutofMemory',
     +                           'cannot deallocate.')
      endif

      return
      end

      subroutine countdim(countd, d, n_d, dim)
      implicit none

      mwSize i, n_d, dim 
      mwSize countd(dim), d(n_d)

      countd=0      
      do i=1,n_d
        countd(d(i))=countd(d(i))+1
      end do

      return
      end

      subroutine cumsum(cumx, x, n_x)
      implicit none

      mwSize i, n_x
      mwSize cumx(n_x), x(n_x)
      
      cumx(1)=x(1)
      do i=2,n_x
        cumx(i)=cumx(i-1)+x(i)
      end do

      return
      end
