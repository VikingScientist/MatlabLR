find_program(MATLAB_MEX NAMES mex
  HINTS "/usr/lib/matlab2011a/bin"
  "/usr/lib/matlab2011b/bin"
  "/usr/lib/matlab2012a/bin"
  "/usr/lib/matlab2012b/bin"
  "/usr/lib/matlab2013a/bin"
  "/usr/lib/matlab2013b/bin"
  "/usr/lib/matlab2014a/bin"
  "/usr/lib/matlab2014b/bin"
  "/usr/lib/matlab2015a/bin"
  "/usr/lib/matlab2015b/bin")

find_program(MATLAB_CC NAMES gcc-4.4 gcc)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Matlab DEFAULT_MSG MATLAB_MEX MATLAB_CC)
