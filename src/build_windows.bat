@echo off
rem ------------------------------------------------------------------------
rem  EBSD-Auswertung: Windows-Programmordner neu erzeugen
rem  Voraussetzung: Python 3.12 oder 3.13 (64 bit) von python.org (Python-Launcher "py")
rem
rem  Python-Umgebung, Zwischendateien und Ergebnis liegen in einem kurzen Pfad
rem  (%LOCALAPPDATA%\EBSD-Build), weil PySide6 sehr lange Dateipfade enthaelt und
rem  Windows ohne "Long Path Support" nur 260 Zeichen pro Pfad erlaubt.
rem  Ergebnis: %LOCALAPPDATA%\EBSD-Build\dist\EBSD-Auswertung\EBSD-Auswertung.exe
rem  Danach uebernehmen.bat ausfuehren.
rem ------------------------------------------------------------------------
setlocal
cd /d "%~dp0"
set "BUILD=%LOCALAPPDATA%\EBSD-Build"
if not exist "%BUILD%" mkdir "%BUILD%" || goto :fehler

rem Reste eines frueheren Versuchs im Quellcode-Ordner entfernen
if exist ".venv" (
    echo Entferne alte .venv im Quellcode-Ordner ...
    rmdir /s /q ".venv"
)

if not exist "%BUILD%\venv\Scripts\python.exe" (
    py -3.12 -m venv "%BUILD%\venv" 2>nul || py -3.13 -m venv "%BUILD%\venv" || goto :fehler
)
call "%BUILD%\venv\Scripts\activate.bat" || goto :fehler
python -m pip install --upgrade pip || goto :fehler
python -m pip install -r requirements.txt || goto :fehler

python -m PyInstaller --noconfirm --clean --workpath "%BUILD%\work" --distpath "%BUILD%\dist" EBSD-Auswertung.spec || goto :fehler
copy /y LIESMICH.txt "%BUILD%\dist\EBSD-Auswertung\" >nul
copy /y README_EN.txt "%BUILD%\dist\EBSD-Auswertung\" >nul
if exist EBSD-Analysis_1.3_Documentation.pdf copy /y EBSD-Analysis_1.3_Documentation.pdf "%BUILD%\dist\EBSD-Auswertung\" >nul
if exist ..\docs\EBSD-Analysis_1.3_Documentation.pdf copy /y ..\docs\EBSD-Analysis_1.3_Documentation.pdf "%BUILD%\dist\EBSD-Auswertung\" >nul

echo.
echo Fertig: %BUILD%\dist\EBSD-Auswertung\EBSD-Auswertung.exe
echo Jetzt uebernehmen.bat ausfuehren.
pause
exit /b 0

:fehler
echo.
echo FEHLER beim Erstellen - siehe Ausgabe oben.
pause
exit /b 1
