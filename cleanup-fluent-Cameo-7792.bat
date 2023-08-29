echo off
set LOCALHOST=%COMPUTERNAME%
set KILL_CMD="C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v231\fluent/ntbin/win64/winkill.exe"

"C:\PROGRA~1\ANSYSI~1\ANSYSS~1\v231\fluent\ntbin\win64\tell.exe" Cameo 58338 CLEANUP_EXITING
if /i "%LOCALHOST%"=="Cameo" (%KILL_CMD% 16168) 
if /i "%LOCALHOST%"=="Cameo" (%KILL_CMD% 10160) 
if /i "%LOCALHOST%"=="Cameo" (%KILL_CMD% 7792) 
if /i "%LOCALHOST%"=="Cameo" (%KILL_CMD% 13740)
del "C:\Users\harsh\Documents\Course 7th sem\internship\CAD\cleanup-fluent-Cameo-7792.bat"
