#! /bin/bash
source CoreFlows.sh
mkdir CoreFlowsTests
cd CoreFlowsTests
$CoreFlows_INSTALL/bin/Executable/CoreFlowsExampleExe
a=$?
if [[ $a == 1 ]]; then
	testsCPP=OK
else
	testsCPP=KO
fi
if [ "$CoreFlows_PYTHON" == 'ON' ]; then
	mkdir CoreFlowsPythonTests
	cd CoreFlowsPythonTests
	$CoreFlows_INSTALL/share/examples/Python/main_tests.py
	b=$?
	if [[ $b == 1 ]]; then
		testsPYTHON=OK
	else
		testsPYTHON=KO
	fi
	cd ..
else
	cd ..
fi

