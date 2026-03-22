# CMake generated Testfile for 
# Source directory: /Users/youenn/Documents/0-Work/SimQJ_project/SimQJ_project/tests
# Build directory: /Users/youenn/Documents/0-Work/SimQJ_project/SimQJ_project/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
include("/Users/youenn/Documents/0-Work/SimQJ_project/SimQJ_project/tests/unit_tests-b12d07c_include.cmake")
add_test(AllTests "/Users/youenn/Documents/0-Work/SimQJ_project/SimQJ_project/tests/unit_tests")
set_tests_properties(AllTests PROPERTIES  _BACKTRACE_TRIPLES "/Users/youenn/Documents/0-Work/SimQJ_project/SimQJ_project/tests/CMakeLists.txt;17;add_test;/Users/youenn/Documents/0-Work/SimQJ_project/SimQJ_project/tests/CMakeLists.txt;0;")
add_test(MathTagGroup "/Users/youenn/Documents/0-Work/SimQJ_project/SimQJ_project/tests/unit_tests" "[math]")
set_tests_properties(MathTagGroup PROPERTIES  LABELS "MathTagGroup" _BACKTRACE_TRIPLES "/Users/youenn/Documents/0-Work/SimQJ_project/SimQJ_project/tests/CMakeLists.txt;26;add_test;/Users/youenn/Documents/0-Work/SimQJ_project/SimQJ_project/tests/CMakeLists.txt;0;")
