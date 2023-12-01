scriptDir=$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")
cd "$scriptDir" || exit 1
cd cpp
mkdir build
cd build
cmake ..
make
cd "$scriptDir" || exit 1
pip install -e .
# rm cpp/build/main
