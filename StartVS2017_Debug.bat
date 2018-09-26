pushd .

set SOURCE_DIR=%cd%
set BUILD_DIR=build\Win64_Debug

mkdir %BUILD_DIR%
cd %BUILD_DIR%

%SOURCE_DIR%\bin\cmake_win\bin\cmake.exe -G "Visual Studio 15 2017 Win64" %SOURCE_DIR%


set INSTALLPATH=

if exist "%programfiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" (
  for /F "tokens=* USEBACKQ" %%F in (`"%programfiles(x86)%\Microsoft Visual Studio\Installer\vswhere.exe" -version 15.0 -property installationPath`) do set INSTALLPATH=%%F
)

echo INSTALLPATH is "%INSTALLPATH%"


if exist geometry3.sln start "%INSTALLPATH%\Common7\IDE\devenv.exe" geometry3.sln
if not exist geometry3.sln echo CMAKE FAILED

popd
