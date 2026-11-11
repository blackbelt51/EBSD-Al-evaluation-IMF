@echo off
rem ------------------------------------------------------------------------
rem  Neu erstelltes Programm (%LOCALAPPDATA%\EBSD-Build\dist\EBSD-Auswertung)
rem  in den Programmordner uebernehmen.
rem  Die bisherige Version 1.2 (Programm und Quellcode-Ordner) wird nach "alte Versionen\V1.2" verschoben.
rem  Vorher EBSD-Auswertung beenden!
rem ------------------------------------------------------------------------
setlocal
cd /d "%~dp0"
set "QUELLE=%LOCALAPPDATA%\EBSD-Build\dist\EBSD-Auswertung"
set "ZIEL=%~dp0..\.."
if not exist "%QUELLE%\EBSD-Auswertung.exe" (
    echo Zuerst build_windows.bat ausfuehren - %QUELLE% fehlt.
    goto :fehler
)
tasklist /fi "imagename eq EBSD-Auswertung.exe" | find /i "EBSD-Auswertung.exe" >nul && (
    echo EBSD-Auswertung laeuft noch - bitte beenden und erneut starten.
    goto :fehler
)
if exist "%ZIEL%\EBSD-Auswertung.exe" (
    if not exist "%ZIEL%\alte Versionen\V1.2" mkdir "%ZIEL%\alte Versionen\V1.2"
    move /y "%ZIEL%\EBSD-Auswertung.exe" "%ZIEL%\alte Versionen\V1.2\" >nul || goto :fehler
    if exist "%ZIEL%\_internal" move /y "%ZIEL%\_internal" "%ZIEL%\alte Versionen\V1.2\_internal" >nul || goto :fehler
    if exist "%ZIEL%\LIESMICH.txt" move /y "%ZIEL%\LIESMICH.txt" "%ZIEL%\alte Versionen\V1.2\" >nul
    if exist "%ZIEL%\README_EN.txt" del /q "%ZIEL%\README_EN.txt"
)
if exist "%ZIEL%\EBSD-Auswertung_1.2_Quellcode" if not exist "%ZIEL%\alte Versionen\V1.2\EBSD-Auswertung_1.2_Quellcode" (
    move /y "%ZIEL%\EBSD-Auswertung_1.2_Quellcode" "%ZIEL%\alte Versionen\V1.2\" >nul || echo Hinweis: Quellcode 1.2 konnte nicht verschoben werden.
)
robocopy "%QUELLE%" "%ZIEL%" /E /NFL /NDL /NJH /NJS >nul
if errorlevel 8 goto :fehler
echo.
echo Fertig: EBSD-Auswertung 1.3 liegt in %ZIEL%
echo Die Version 1.2 liegt in "%ZIEL%\alte Versionen\V1.2".
pause
exit /b 0

:fehler
echo.
echo FEHLER - siehe Meldung oben.
pause
exit /b 1
