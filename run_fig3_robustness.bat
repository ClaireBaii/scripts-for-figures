@echo off
set "R_HOME=D:\conda\envs\paper_fig_env\Lib\R"
set "PATH=D:\conda\envs\paper_fig_env\Lib\R\bin\x64;D:\conda\envs\paper_fig_env\Library\bin;D:\conda\envs\paper_fig_env\Library\mingw-w64\bin;D:\conda\envs\paper_fig_env\Scripts;%PATH%"
"D:\conda\envs\paper_fig_env\Lib\R\bin\x64\Rscript.exe" --encoding=UTF-8 figure3/figure3_rarefaction_robustness.R > fig3_robustness.log 2>&1
if %errorlevel% neq 0 (
    echo Error occurred with code %errorlevel% >> fig3_robustness.log
) else (
    echo Success >> fig3_robustness.log
)
