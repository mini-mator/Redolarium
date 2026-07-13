@echo off
cd /d "%~dp0"

echo =========================================
echo Redolarium Benchmark Suite Runner
echo =========================================
echo 1) Run Phase A: Scrape 100 Genomes
echo 2) Run Phase B: Execute Benchmark Suite
echo 3) Run Phase C: Generate Analytics ^& Validation Graphs
echo =========================================
set /p choice="Enter 1, 2, or 3: "

set PYTHON_PATH="C:\Users\antog\AppData\Local\Programs\Python\Python312\python.exe"

if "%choice%"=="1" (
    echo.
    echo Starting the Scraper...
    %PYTHON_PATH% scrape_benchmark_genomes.py
    pause
) else if "%choice%"=="2" (
    echo.
    echo Starting the Benchmark Dashboard...
    %PYTHON_PATH% run_benchmark_suite.py
    pause
) else if "%choice%"=="3" (
    echo.
    echo Starting Post-Run Analytics...
    %PYTHON_PATH% analyze_benchmark_results.py
    pause
) else (
    echo Invalid choice.
    pause
)
