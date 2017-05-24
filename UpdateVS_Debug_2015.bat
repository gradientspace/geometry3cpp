pushd .

set SOURCE_DIR=%cd%
set BUILD_DIR=build\Win64_Debug

mkdir %BUILD_DIR%
cd %BUILD_DIR%

%SOURCE_DIR%\bin\cmake_win\bin\cmake.exe -G "Visual Studio 14 2015 Win64" %SOURCE_DIR%

call "%VS140COMNTOOLS%\vsvars32.bat"

if not exist geometry3.sln echo CMAKE FAILED

popd
