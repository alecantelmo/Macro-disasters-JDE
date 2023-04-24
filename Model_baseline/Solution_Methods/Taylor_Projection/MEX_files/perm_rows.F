C © Copyright, Oren Levintal, June 13, 2016.

#include "fintrf.h"
C [IC,JC,NC] = perm_rows(IA,JA,NA,l,perml)
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
      mwPointer l,perml,IA,JA,NA
      mwPointer IC,JC,NC

C     Array information:
      mwSize IArows,JArows,JAcols,NArows,NAcols

      mwSize sA,dims(2),ndims


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

C     Create Fortran array from the input arguments.
      IA = mxGetPr(prhs(1))
      JA = mxGetPr(prhs(2))
      NA = mxGetPr(prhs(3))
      l = mxGetPr(prhs(4))
      perml = mxGetPr(prhs(5))

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
      if (NArows==1) then
      call compute1(%VAL(IA),%VAL(JA),%VAL(NA),
     + %VAL(l),%VAL(perml),
     + %VAL(IC),%VAL(JC),%VAL(NC),
     + IArows,JArows,JAcols,NArows,NAcols)
      else
      call compute2(%VAL(IA),%VAL(JA),%VAL(NA),
     + %VAL(l),%VAL(perml),
     + %VAL(IC),%VAL(JC),%VAL(NC),
     + IArows,JArows,JAcols,NArows,NAcols)
      endif

      return
      end

C     Computational subroutine
      subroutine compute1(IA,JA,NA,l,perml,
     + IC,JC,NC,
     + IArows,JArows,JAcols,NArows,NAcols)

C     Declarations
      implicit none

      mwSize  IArows,JArows,JAcols,NArows,NAcols,l
      mwSize  IA(IArows),IC(IArows)
      mwSize  JA(JArows,JAcols),JC(JArows,JAcols)
      real*8  NA(NAcols),NC(NAcols)

      mwSize  perml(l),i,j,k,startj,endj,si
 
C create IC,JC,NC
      IC(1)=1
      do i=1,l
        startj=IA(perml(i))
        endj=IA(perml(i)+1)-1
        IC(i+1)=IC(i)+endj-startj+1
      !DEC$ simd
        do j=startj,endj
         NC(IC(i)+j-startj)=NA(j)
        enddo
        do k=1,JAcols
      !DEC$ simd
        do j=startj,endj
         JC(IC(i)+j-startj,k)=JA(j,k)
        enddo
        enddo
      enddo
      return
      end


      subroutine compute2(IA,JA,NA,l,perml,
     + IC,JC,NC,
     + IArows,JArows,JAcols,NArows,NAcols)

C     Declarations
      implicit none

      mwSize  IArows,JArows,JAcols,NArows,NAcols,l
      mwSize  IA(IArows),IC(IArows)
      mwSize  JA(JArows,JAcols),JC(JArows,JAcols)
      real*8  NA(NArows,NAcols),NC(NArows,NAcols)

      mwSize  perml(l),i,j,k,startj,endj,si
 
C create IC,JC,NC
      IC(1)=1
      do i=1,l
        startj=IA(perml(i))
        endj=IA(perml(i)+1)-1
        IC(i+1)=IC(i)+endj-startj+1
        do j=startj,endj
      !DEC$ simd
         do si=1,NArows
          NC(si,IC(i)+j-startj)=NA(si,j)
         enddo
        enddo
        do k=1,JAcols
      !DEC$ simd
        do j=startj,endj
         JC(IC(i)+j-startj,k)=JA(j,k)
        enddo
        enddo
      enddo
      return
      end

