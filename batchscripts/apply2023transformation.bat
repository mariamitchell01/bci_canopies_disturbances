@echo off
set loopcount=844
:loop
echo Iter %loopcount%
"C:\Program Files\CloudCompare\CloudCompare.exe" -SILENT -AUTO_SAVE OFF -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\TroubleshootingV\tiles\BCI23FTiles\BCIf23_%LOOPCOUNT%.laz" -APPLY_TRANS "D:\TroubleshootingV\tiles\BCI23FTiles_dec\BCI23fmat_%LOOPCOUNT%.txt" -C_EXPORT_FMT LAS -SAVE_CLOUDS FILE "D:\TroubleshootingV\tiles\BCI23FTiles_alignedto22Full\BCI23faf_%LOOPCOUNT%.las" -CLEAR
set /a loopcount=loopcount-1
if %loopcount%==0 goto exitloop
goto loop
:exitloop
pause