rem helper.bat - used to copy package warp.lib and warp.dll to Python installation
echo off
echo
echo "****************IMPORTANT Read First ********************************"
echo "First Make sure your <package> package doesn't exist already in Python27"
echo "This installs warp.pyd and warp.lib. Do not name package that of an existing one."
echo "This will overwrite existing warp.pyd and warp.lib in C:\_Python27\DLLs and c:\_Python27\libs"
echo "**********************************************************************"
echo
echo "Hit space to continue."
echo on
pause
copy warp.dll warp.pyd
Copy warp.pyd C:\_Python27\DLLs\warp.pyd
copy warp.lib C:\_Python27\libs\warp.lib
