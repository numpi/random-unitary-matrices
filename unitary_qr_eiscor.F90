#include <fintrf.h>

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
  implicit none

  ! There is a small technical difficulty here: since eiscor is usually
  ! compiled with the default integer type of the system, which will be
  ! integer*4 most of the time, we need to pass the integer variables in
  ! this format. This is not the same as MATLAB integers. 
  integer*4 :: eiscor_n, info
  
  integer :: nlhs, nrhs, n
  mwPointer :: plhs(*), prhs(*)
  mwPointer :: mxGetN, mxGetPr, mxCreateDoubleMatrix
  real*8, allocatable :: Q(:, :), D(:, :)
  integer*4, allocatable :: its(:)
  character(LEN=256) :: buf

  ! We expect the first input to be a 3 x (N - 1) double
  ! matrix containing the cosines and sines of the rotation
  ! generating the unitary matrix Q.
  n = mxGetN(prhs(2))

  allocate(Q(3, n-1))
  allocate(D(2, n))
  allocate(its(n-1))

  ! Copy the data in the matrices
  call mxCopyPtrToReal8(mxGetPr(prhs(1)), Q, 3 * (n - 1))
  call mxCopyPtrToReal8(mxGetPr(prhs(2)), D, 2 * n)

  ! Call the function in eiscor
  eiscor_n = n
  call z_unifact_qr(.false., .false., eiscor_n, Q, D, &
       eiscor_n, 0, its, info)

  if (info .ne. 0) then
     call mexPrintf("Warning: non-zero exit code from z_unifact_qr\n")
     write(buf, *) " :: UNITARY_QR_EISCOR // INFO = ", info, "\n"
     call mexPrintf(buf)
  end if

  ! Copy the data back in the output
  plhs(1) = mxCreateDoubleMatrix(2, n, 0)
  call mxCopyReal8ToPtr(D, mxGetPr(plhs(1)), 2 * n)

  ! Cleanup
  deallocate(Q, D, its)

end subroutine mexFunction

