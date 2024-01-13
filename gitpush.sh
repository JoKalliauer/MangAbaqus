#!/bin/sh

#find . -path ./.git -prune -false -o -name '*' -size +1M
find . -type d \( -path ./.git -o -path ./Output \) -prune -false -o -name '*' -size +200k #https://stackoverflow.com/a/4210072/6747994

if [ -f *.owncloud ]; then
    echo "Please download all files from owncloud before using git"
    read -p "Do you wish to proceed?" yn
    case $yn in
        [Yy]* ) echo "not shure if thats a good choise";;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac
fi


git config --global core.eol lf
git config core.eol lf
git pull
#git status
git add .
git reset -- AbaqusRuns/* AnalysisResults/*  Output/* *.mat *.m~ *.mtx *.png *.eps *.asv *.svg
git status


    read -p "Do you wish to proceed?" yn
    case $yn in
        [Yy]* ) git commit -m "Update";git push;;
        [Nn]* ) exit;;
        * ) echo "Please answer yes or no.";;
    esac


#git status

#git status

