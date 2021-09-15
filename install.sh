set -e

if [ $# -ne 1 ]; then
  echo "Usage: $0 INSTALL_PREFIX"
  exit 1
fi

INSTALL_PREFIX=$1

rm -rf build
mkdir -p build
pushd build/
cmake .. \
	-DCMAKE_INSTALL_PREFIX=${INSTALL_PREFIX} \
	-DCMAKE_BUILD_TYPE=Debug \
	-DCMAKE_EXPORT_COMPILE_COMMANDS=true

cmake --build . -- VERBOSE=1
cmake --build . --target install
popd

echo ""
echo "*****************************************************************************"
echo "  DONE: MBH built and installed to: "
echo "      ${INSTALL_PREFIX}"
echo "  To use: "
echo "      export PATH=${INSTALL_PREFIX}/bin:\${PATH}"
echo "*****************************************************************************"


