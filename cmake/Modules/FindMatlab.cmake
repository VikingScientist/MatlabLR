find_program(MATLAB_MEX NAMES mex
  HINTS "/usr/lib/matlab2021a/bin"
        "/usr/lib/matlab2021b/bin"
        "/usr/lib/matlab2020a/bin"
        "/usr/lib/matlab2020b/bin"
        "/usr/lib/matlab2019a/bin"
        "/usr/lib/matlab2019b/bin"
        "/usr/lib/matlab2018a/bin"
        "/usr/lib/matlab2018b/bin"
        "/usr/lib/matlab2017a/bin"
        "/usr/lib/matlab2017b/bin"
        "/usr/lib/matlab2016a/bin"
        "/usr/lib/matlab2016b/bin"
        "/usr/lib/matlab2015a/bin"
        "/usr/lib/matlab2015b/bin"
        "/usr/lib/matlab2014a/bin"
        "/usr/lib/matlab2014b/bin"
        "/usr/lib/matlab2013a/bin"
        "/usr/lib/matlab2013b/bin"
        "/usr/lib/matlab2012a/bin"
        "/usr/lib/matlab2012b/bin")

find_program(MATLAB_CC NAMES gcc-4.4 gcc)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Matlab DEFAULT_MSG MATLAB_MEX MATLAB_CC)
