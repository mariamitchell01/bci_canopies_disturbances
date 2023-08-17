@echo off
set loopcount=844
:loop
echo Iter %loopcount%
"C:\Program Files\CloudCompare\CloudCompare.exe" -SILENT -AUTO_SAVE OFF -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\TroubleshootingV\tiles\BCI23FTiles_dec\BCI23fd_%LOOPCOUNT%.laz" -O -GLOBAL_SHIFT -620000.00 -1000000.00 0 "D:\TroubleshootingV\tiles\BCI22Tiles_alignedto21Full\BCI22af_%LOOPCOUNT%.las" -ICP -ITER 800 -OVERLAP 80 -C_EXPORT_FMT LAS -SAVE_CLOUDS FILE "D:\TroubleshootingV\tiles\BCIF23Tiles_alignedto22\BCI23fa_%LOOPCOUNT%.las D:\TroubleshootingV\tiles\BCI22FTiles_ref\BCI22fr_%LOOPCOUNT%.las" -CLEAR
set /a loopcount=loopcount-1
if %loopcount%==0 goto exitloop
goto loop
:exitloop
pause