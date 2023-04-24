C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"

C [ newcols,newvals ] = sortcol_mex( ptr,cols,vals,l,coldim,coli,colidim )
C the function sorts columns by coli
      subroutine mexfunction(nlhs, plhs, nrhs, prhs)
      implicit none

      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      mwPointer mxCreateNumericArray
      integer*4 mxGetClassID

      mwPointer mxGetM, mxGetN

      mwPointer ptr,cols,vals,l,coldim,coli
      mwPointer newcols,newvals

      mwSize ptrM,colsM,colsN,valsM,valsN,coldimM,coldimN
      mwSize colidim

C     Arguments for mxCreateNumericArray
      integer*4 classid
      integer*4 complexflag
      mwSize ndim
      mwSize dims(2)

C-----------------------------------------------------------------------
C     Check for proper number of arguments. 
      if(nrhs .ne. 7) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:nInput',
     +                           '7 inputs required.')
      endif
      if(nlhs .gt. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:compress:nOutput',
     +                           '2 outputs required.')
      endif


C     Create Fortran array from the input arguments.

      ptr = mxGetPr(prhs(1))
      cols = mxGetPr(prhs(2))
      vals = mxGetPr(prhs(3))
      l = mxGetPr(prhs(4))
      coldim = mxGetPr(prhs(5))
      coli = mxGetPr(prhs(6))

      ptrM = mxGetM(prhs(1))

      colsM = mxGetM(prhs(2))
      colsN = mxGetN(prhs(2))

      valsM = mxGetM(prhs(3))
      valsN = mxGetN(prhs(3))

      coldimM = mxGetM(prhs(5))
      coldimN = mxGetN(prhs(5))
      coldimN = coldimM*coldimN

C     Create matrix for the return arguments.

      classid = mxGetClassID(prhs(2))
      complexflag=0

      ndim = 2
      dims(1) = colsM
      dims(2) = colsN
      plhs(1) = mxCreateNumericArray(ndim, dims, 
     +          classid, 
     +          complexflag)

      newcols = mxGetPr(plhs(1))

      plhs(2) = mxCreateDoubleMatrix(valsM, valsN, complexflag)
      
      newvals = mxGetPr(plhs(2))

C     Call the computational routine.
      call compute(%VAL(ptr),%VAL(cols),%VAL(vals),
     +      %VAL(l),%VAL(coldim),%VAL(coli),
     +      ptrM,colsM,colsN,valsM,valsN,coldimN,
     +      %VAL(newcols),%VAL(newvals))

      return
      end

      subroutine compute(ptr,cols,vals,l,coldim,coli, 
     +      ptrM,colsM,colsN,valsM,valsN,coldimN,
     +      newcols,newvals)
      implicit none

      mwSize ptrM,colsM,colsN,valsM,valsN,coldimN
      mwSize ptr(ptrM),coldim(coldimN)

      mwSize l,coli,colidim
      mwSize cols(colsM,colsN)

      mwSize newcols(colsM,colsN)
      real*8 vals(valsM,valsN),newvals(valsM,valsN)

      mwSize i,n_vals,j,temploc1,temploc2,n_valsi

      mwSize, allocatable :: countcol(:),z(:),loc(:),ind(:),tempptr(:)

      integer AllocateStatus

      colidim=coldim(coli)
      n_vals=ptr(ptrM)-1

      allocate (countcol(colidim),z(colidim+1),
     +   tempptr(colidim+1),loc(colidim),ind(n_vals),
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
        countcol(cols(j,coli))=countcol(cols(j,coli))+1
       end do

C create ptr for the column dimension
       z(1)=1
       z(2:colidim+1)=countcol
       call cumsum(tempptr,z,colidim+1)

C create sort index
       loc=tempptr(1:colidim)
       do j=1,n_valsi
        temploc1=ptr(i)+j-1
        ind(j)=loc(cols(temploc1,coli))
        loc(cols(temploc1,coli))=loc(cols(temploc1,coli))+1
       end do

C sort      
       do j=1,n_valsi
        temploc1=ptr(i)+ind(j)-1
        temploc2=ptr(i)+j-1
        newvals(:,temploc1)=vals(:,temploc2)
        newcols(temploc1,:)=cols(temploc2,:)
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
