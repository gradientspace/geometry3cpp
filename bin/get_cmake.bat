set CMAKE_VERSION=cmake-3.12.2-win32-x86
%~dp0\wget_win\wget.exe https://cmake.org/files/v3.12/%CMAKE_VERSION%.zip
%~dp0\7z x -y %CMAKE_VERSION%.zip
del %CMAKE_VERSION%.zip
move %CMAKE_VERSION% %~dp0\cmake_win
