@ECHO OFF

:: Assign env variables necessary
:: %HOMEDRIVE% - C:
:: %HOMEPATH% - \Users\username
:: %ProgramFiles% - C:\Program Files\
:: This assumes that STK is installed in the default location
:: You will need to update the CONDA_ENVIRONMENT variable to the location of the stk_base env.

SET CONDA_ENVIRONMENT=C:\Users\asligar\Anaconda3\envs\farfields
SET CONDA_BASE=C:\Users\asligar\Anaconda3\envs\farfields
SET PYTHONHOME=%CONDA_ENVIRONMENT%
SET PYTHONPATH=%CONDA_ENVIRONMENT%\python.exe
SET PATH=%CONDA_BASE%;%CONDA_BASE%\Scripts;%CONDA_BASE%\Library\bin;%CONDA_BASE%\Library\mingw-w64\bin;%PATH%

start "" "%ProgramFiles%\AGI\STK 12\bin\AgUiApplication.exe" /pers "STK"