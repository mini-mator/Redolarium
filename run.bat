@echo off
cd /d "%~dp0"
echo Starting Redolarium Streamlit App...
streamlit run src/front_end/app.py
pause
