module purge
module load cmake/3.28.3 gcc/latest-v13 openmpi/4.1.6  python/3.13.5  cgal/6.0.1  boost/1.88

build='../RELEASE_202604'
#install='../INSTALL_202604'
install='/home/user/omierka/nobackup/MyApplications/GenDieMesher'

#mkdir -p ${build}
#cmake -S . -B ${build} -DCMAKE_BUILD_TYPE=Release
cmake --build ${build} -j
cmake --install ${build} --prefix ${install}

