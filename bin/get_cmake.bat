set CMAKE_VERSION=cmake-3.5.2-win32-x86
wget_win\wget.exe https://cmake.org/files/v3.5/%CMAKE_VERSION%.zip
7z x -y %CMAKE_VERSION%.zip
del %CMAKE_VERSION%.zip
move %CMAKE_VERSION% cmake_win
