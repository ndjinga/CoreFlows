source CoreFlows.sh
rm -rf CoreFlows_build $CoreFlows_INSTALL
mkdir -p CoreFlows_build
mkdir -p $CoreFlows_INSTALL

cd CoreFlows_build
cmake $CoreFlows_ROOT/CoreFlows_src -DCMAKE_INSTALL_PREFIX=$CoreFlows_INSTALL -DCMAKE_BUILD_TYPE=Debug -G"Eclipse CDT4 - Unix Makefiles" -D_ECLIPSE_VERSION=4.3 -DCOREFLOWS_WITH_DOCUMENTATION=$CoreFlows_DOC -DCOREFLOWS_WITH_PYTHON=$CoreFlows_PYTHON -DCOREFLOWS_WITH_GUI=$CoreFlows_GUI -DCOREFLOWS_WITH_PACKAGE=OFF

make 

if [ $CoreFlows_DOC=ON ] 
then
 make doc 
fi

make install -j4 

cd ..
chmod -R 755 $CoreFlows_INSTALL/bin/CoreFlows_Python/* $CoreFlows_INSTALL/share/examples/Python/*
