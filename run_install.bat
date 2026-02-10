@echo off
set "R_HOME=D:\conda\envs\paper_fig_env\Lib\R"
set "PATH=D:\conda\envs\paper_fig_env\Lib\R\bin\x64;D:\conda\envs\paper_fig_env\Library\bin;D:\conda\envs\paper_fig_env\Library\mingw-w64\bin;D:\conda\envs\paper_fig_env\Scripts;%PATH%"
"D:\conda\envs\paper_fig_env\Lib\R\bin\x64\Rscript.exe" install_rlang.R > install_run.log 2>&1
