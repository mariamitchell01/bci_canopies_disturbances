@echo off
set loopcount=844
:loop
echo Iter %loopcount%
"C:\Program Files\CloudCompare\CloudCompare.exe" -SILENT -AUTO_SAVE OFF -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\TroubleshootingV\tiles\BCI22Tiles_dec\BCI22d_%LOOPCOUNT%.laz" -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\TroubleshootingV\tiles\BCI21Tiles_alignedto20Trim\BCI21at_%LOOPCOUNT%.laz" -ICP -ITER 800 -OVERLAP 80 -C_EXPORT_FMT LAS -SAVE_CLOUDS FILE "D:\TroubleshootingV\tiles\BCI22Tiles_alignedto21\BCI22a_%LOOPCOUNT%.las D:\TroubleshootingV\tiles\BCI21Tiles_ref\BCI21r_%LOOPCOUNT%.las" -CLEAR
set /a loopcount=loopcount-1
if %loopcount%==0 goto exitloop
goto loop
:exitloop
pause