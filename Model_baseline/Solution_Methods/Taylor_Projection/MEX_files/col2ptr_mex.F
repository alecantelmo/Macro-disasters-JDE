C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"

C [ newptr,newcols,newvals ] = col2ptr_mex( ptr,cols,vals,l,coldim,colptr,ptrcoldim )
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
      integer*4 mexPrintf

      mwPointer ptr,cols,vals,l,coldim,colptr
      mwPointer newcols,newvals

      mwSize ptrM,colsM,colsN,valsM,valsN,coldimM,coldimN
      mwSize newptrM,newptrN
      mwSize ptrcoldim

C     Arguments for mxCreateNumericArray
      integer*4 classid
      integer*4 complexflag
      mwSize ndim
      mwSize dims(2)
      mwSize newptr


C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 7) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:nInput',
     +                           '7 inputs required.')
      endif
      if(nlhs .gt. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:nOutput',
     +                           '3 outputs required.')
      endif


C     Create Fortran array from the input arguments.

      ptr = mxGetPr(prhs(1))
      cols = mxGetPr(prhs(2))
      vals = mxGetPr(prhs(3))
      l = mxGetPr(prhs(4))
      coldim = mxGetPr(prhs(5))
      colptr = mxGetPr(prhs(6))

      ptrM = mxGetM(prhs(1))

      colsM = mxGetM(prhs(2))
      colsN = mxGetN(prhs(2))

      valsM = mxGetM(prhs(3))
      valsN = mxGetN(prhs(3))

      coldimM = mxGetM(prhs(5))
      coldimN = mxGetN(prhs(5))
      coldimN = coldimM*coldimN

C     Create matrix for the return arguments.
C     2 dimensional ptr

      ptrcoldim=mxGetScalar(prhs(7)) ! second dimension of newptr


      newptrM=ptrM-1
      newptrN=ptrcoldim

      classid = mxGetClassID(prhs(1))
      ndim=2
      dims(1)=newptrM
      dims(2)=newptrN
      complexflag = 0
      
      plhs(1) = mxCreateNumericArray(ndim, dims, 
     +          classid, 
     +          complexflag)

      newptr = mxGetPr(plhs(1))

      classid = mxGetClassID(prhs(2))
      ndim = 2
      dims(1) = colsM
      dims(2) = colsN-1
      plhs(2) = mxCreateNumericArray(ndim, dims, 
     +          classid, 
     +          complexflag)

      newcols = mxGetPr(plhs(2))

      plhs(3) = mxCreateDoubleMatrix(valsM,valsN,complexflag)
      
      newvals = mxGetPr(plhs(3))

C     Call the computational routine.
      call compute(%VAL(ptr),%VAL(cols),%VAL(vals),
     +      %VAL(l),%VAL(coldim),%VAL(colptr),
     +      ptrM,colsM,colsN,valsM,valsN,coldimN,
     +      %VAL(newptr),%VAL(newcols),%VAL(newvals),
     +      newptrM,newptrN)

      return
      end


      subroutine compute(ptr,cols,vals,l,coldim,colptr, 
     +      ptrM,colsM,colsN,valsM,valsN,coldimN,
     +      newptr,newcols,newvals,
     +      newptrM,newptrN)
      implicit none

      mwSize ptrM,colsM,colsN,valsM,valsN,coldimN,newptrM,newptrN
      mwSize ptr(ptrM),coldim(coldimN)

      mwSize l,colptr,colptrdim,newcolsind(colsN-1)
      mwSize, target :: cols(colsM,colsN)
      mwSize, pointer :: pcol(:,:)

      mwSize newptr(newptrM,newptrN),newcols(colsM,colsN-1)
      real*8 vals(valsM,valsN),newvals(valsM,valsN)

      mwSize i,n_vals,j,temploc1,temploc2,n_valsi

      mwSize, allocatable :: countcol(:),z(:),loc(:),ind(:),tempptr(:)

      integer AllocateStatus

      colptrdim=coldim(colptr)
      n_vals=ptr(ptrM)-1
      newcolsind(1:colptr-1)=(/(i,i=1,colptr-1)/)
      newcolsind(colptr:colsN-1)=(/(i,i=colptr+1,colsN)/)

      allocate (countcol(colptrdim),z(colptrdim+1),
     +   tempptr(colptrdim+1),loc(colptrdim),ind(n_vals),
     +   STAT = AllocateStatus)      
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:OutofMemory',
     +                           'Out of memory.')
      endif
      do i=1,l
       n_valsi=ptr(i+1)-ptr(i)
C count cols
       countcol=0      
       do j=ptr(i),ptr(i+1)-1
        countcol(cols(j,colptr))=countcol(cols(j,colptr))+1
       end do

C create ptr for the column dimension
       z(1)=1
       z(2:colptrdim+1)=countcol
       call cumsum(tempptr,z,colptrdim+1)

C assign in newptr
       newptr(i,:)=tempptr+ptr(i)-1      

C create sort index
       loc=tempptr(1:colptrdim)
       do j=1,n_valsi
        temploc1=ptr(i)+j-1
        ind(j)=loc(cols(temploc1,colptr))
        loc(cols(temploc1,colptr))=loc(cols(temploc1,colptr))+1
       end do

C sort      
       do j=1,n_valsi
        temploc1=ptr(i)+ind(j)-1
        temploc2=ptr(i)+j-1
        newvals(:,temploc1)=vals(:,temploc2)
        newcols(temploc1,1:colsN-1)=cols(temploc2,newcolsind)
       end do
      end do

      deallocate(countcol,z,
     +   tempptr,loc,ind,
     +   STAT = AllocateStatus)
      IF (AllocateStatus .ne. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:OutofMemory',
     +                           'cannot deallocate.')
      endif

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
