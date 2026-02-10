@echo off
set "R_HOME=D:\conda\envs\paper_fig_env\Lib\R"
set "PATH=D:\conda\envs\paper_fig_env\Lib\R\bin\x64;D:\conda\envs\paper_fig_env\Library\bin;D:\conda\envs\paper_fig_env\Library\mingw-w64\bin;D:\conda\envs\paper_fig_env\Scripts;%PATH%"
"D:\conda\envs\paper_fig_env\Lib\R\bin\x64\Rscript.exe" --encoding=UTF-8 figure5/fig5_threshold_sensitivity_table.R > fig5_run.log 2>&1
if %errorlevel% neq 0 (
    echo Error occurred with code %errorlevel% >> fig5_run.log
) else (
    echo Success >> fig5_run.log
)
