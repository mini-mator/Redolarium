@echo off
setlocal EnableDelayedExpansion
title Redolarium Ground Truth Validation Suite

set PYTHON_PATH="C:\Users\antog\AppData\Local\Programs\Python\Python312\python.exe"

:MENU
cls
echo ===============================================================================
echo                REDOLARIUM - MIBiG GROUND TRUTH VALIDATION SUITE                
echo ===============================================================================
echo.
echo [1] Harvest Golden Genomes (Downloads MIBiG 3.1 ^& 100 NCBI Genomes)
echo [2] Run Validation Pipeline (Executes Redolarium with Live Progress)
echo [3] Generate Truth Analytics (Computes Recall ^& Compiles FAIR Bundle)
echo [4] Exit
echo.
set /p "choice=Select an option (1-4): "

if "%choice%"=="1" goto HARVEST
if "%choice%"=="2" goto VALIDATE
if "%choice%"=="3" goto ANALYZE
if "%choice%"=="4" goto EOF

echo Invalid choice. Please press any key to try again.
pause >nul
goto MENU

:HARVEST
cls
echo ===============================================================================
echo                     PHASE 1: HARVESTING GOLDEN GENOMES
echo ===============================================================================
%PYTHON_PATH% "%~dp0harvest_golden_genomes.py"
echo.
pause
goto MENU

:VALIDATE
cls
echo ===============================================================================
echo                     PHASE 2: RUNNING VALIDATION SUITE
echo ===============================================================================
%PYTHON_PATH% "%~dp0run_validation_suite.py"
echo.
pause
goto MENU

:ANALYZE
cls
echo ===============================================================================
echo                     PHASE 3: GENERATING TRUTH ANALYTICS
echo ===============================================================================
%PYTHON_PATH% "%~dp0generate_truth_analytics.py"
echo.
pause
goto MENU

:EOF
exit /b 0
