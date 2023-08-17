@echo off
set loopcount=844
:loop
echo Iter %loopcount%
"C:\Program Files\CloudCompare\CloudCompare.exe" -SILENT -AUTO_SAVE OFF -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\TroubleshootingV\tiles\BCI22Tiles\BCI22_%LOOPCOUNT%.laz" -APPLY_TRANS "D:\TroubleshootingV\tiles\BCI22Tiles_dec\BCI22mat_%LOOPCOUNT%.txt" -C_EXPORT_FMT LAS -SAVE_CLOUDS FILE "D:\TroubleshootingV\tiles\BC22Tiles_alignedto21Full\BCI22af_%LOOPCOUNT%.las" -CLEAR
set /a loopcount=loopcount-1
if %loopcount%==0 goto exitloop
goto loop
:exitloop
pause