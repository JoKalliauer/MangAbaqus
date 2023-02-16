@ECHO OFF
::REM:: chmod u+rx *.sh
::REM:: chmod u+r *
::REM:: make clean
if exist *.o del *.o
if exist output del output
git config --global core.autocrlf input
git config --global core.eol lf
git config core.eol lf
git config pull.rebase false
git pull
REM git status
git add .
git reset -- *.o *.csv *.txt output .goutputstream-D4N0N1 *.eps Output/* Figures/* *.emf *.svg *.mat *.png *.asv
git status
set /p Answer="Do you wish to proceed? (y/n)"

if %Answer%==y (
 git commit -m "Update"
 git push
 ) else (
 echo please reply with "y" if you like to proceed
)
::wait

