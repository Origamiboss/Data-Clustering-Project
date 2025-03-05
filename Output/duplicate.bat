@echo off
setlocal enabledelayedexpansion

set "input_file=iris_bezdek.txt"
set "output_file=iris_bezdek_mod.txt"
set /a count=0

> "%output_file%" (
    for /f "tokens=*" %%a in (%input_file%) do (
        set /a count+=1
        echo %%a
        set /a mod=count%%5
        if !mod! == 0 (
            for /L %%i in (1,1,10) do echo %%a
        )
    )
)

echo Processing complete. Duplicated every 5th line 10 times.
