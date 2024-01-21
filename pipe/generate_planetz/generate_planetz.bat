@echo off
for /l %%i in (16,2,220) do (
echo %%i
echo %%i | python generate_planetz.py
)