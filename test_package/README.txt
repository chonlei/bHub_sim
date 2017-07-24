README (last updated: 24/07/17)
===============================
These tests should serve as testing purpose for project bHub_sim and should be tested on ARCUS(-B).

Please compile folder test1 using NEURON, this can be done by:
$ cd /dir/to/test_package
$ cd test1 && nrnivmodl

Please run the test script runTests.sh by either submitting it as a job or change it into an appropriate job script. If either does not apply, please directly run the script locally by:
$ cd /dir/to/test_package
$ bash runTests.sh