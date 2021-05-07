#!/bin/sh

#find . -path ./.git -prune -false -o -name '*' -size +1M
find . -type d \( -path ./.git -o -path ./Output \) -prune -false -o -name '*' -size +70k #https://stackoverflow.com/a/4210072/6747994

git config --global core.eol lf
git config core.eol lf
git pull
#git status
git add .
git reset -- AbaqusRuns/* AnalysisResults/*  Output/* *.mat *.m~ *.mtx
git status


    read -p "Do you wish to proceed?" yn
    case $yn in
        [Yy]* ) git commit -m "Update";git push;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac


#git status

#git status

