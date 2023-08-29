echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v231\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v231\fluent\ntbin\win64\tell.exe" Cameo 57275 CLEANUP_EXITING
if /i "%LOCALHOST%"=="Cameo" (%KILL_CMD% 16836) 
if /i "%LOCALHOST%"=="Cameo" (%KILL_CMD% 15080) 
if /i "%LOCALHOST%"=="Cameo" (%KILL_CMD% 14448) 
if /i "%LOCALHOST%"=="Cameo" (%KILL_CMD% 13256)
del "C:\Users\harsh\Documents\Course 7th sem\internship\CAD\cleanup-fluent-Cameo-14448.bat"
