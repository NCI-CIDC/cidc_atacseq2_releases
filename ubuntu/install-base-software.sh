#!/usr/bin/env bash
##############################################################################################################################
# pipeline_template
# 
# This program is free software that contains third party software subject to various licenses, 
# namely, the GNU General Public License version 3 (or later), the GNU Affero General Public License 
# version 3 (or later), and the LaTeX Project Public License v.1.3(c). A list of the software contained 
# in this program, including the applicable licenses, can be accessed here: 
#
# 
# You can redistribute and/or modify this program, including its components, only under the terms of 
# the applicable license(s).  
#
# This program is distributed in the hope that it will be useful, but "as is," WITHOUT ANY WARRANTY; 
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#
# Program:  install-software.sh
# Version:  v0.1
# Author:   Kevin Conway, Leigh Villarroel, Sami Cherikh
# Purpose:  install software onto Ubuntu machine
# Input:    N/A
# Output:   N/A
##############################################################################################################################

##Add Swap for extra memory if needed
sudo /bin/dd if=/dev/zero of=/var/swap.1 bs=1M count=1024
sudo /sbin/mkswap /var/swap.1
sudo chmod 600 /var/swap.1
sudo /sbin/swapon /var/swap.1



##Initial setup
sudo apt-get update



## create/cd into download/pre-installation dir - making on home/<myuser> currently due to permissions
mydir=~/software_install
if [ ! -d "$mydir" ]
then
	echo creating directory "$mydir"
	mkdir $mydir
	if [ ! -d "$mydir" ]
	then
		echo cannot create directory "$mydir"
		echo exiting with error code 1 ...
		exit 1
	fi
fi
echo cd-ing into directory "$mydir"
cd "$mydir"



## GCC and other Configure/Make Components
mygcc=/usr/bin/gcc
if [ ! -f "$mygcc" ]
then
	sudo apt-get -y install build-essential
fi

sudo apt-get -y install aptitude
sudo aptitude -y install libreadline-dev

echo gcc and other necessary compiler tools INSTALLED



## JAVA
#myjava=/usr/bin/java
#if [ ! -f "$myjava" ]
#then
#	sudo apt -y install default-jre
#fi
##sudo apt-get -y install openjdk-6-jdk
#echo java is INSTALLED



### PYTHON 2
#mypython27=/usr/bin/python2.7
#mypython=/usr/bin/python
#if [ ! -f "$mypython" ] && [ ! -f "$mypython27" ]
#then
#	sudo apt-get -y install python2.7 python-dev libpython2.7-dev
#	sudo cp /usr/bin/python2.7 /usr/bin/python2
#elif [ ! -f "$mypython" ] && [ -f "$mypython27" ]
#then
#	sudo apt-get -y install python-dev libpython2.7-dev
#	sudo cp /usr/bin/python2.7 /usr/bin/python2
#fi
#echo "Installing required python modules."
#sudo apt-get -y install python-dev libpython2.7-dev #Gives the Python.h header files for rseqc install
#sudo apt-get -y install libbz2-1.0 libbz2-dev libbz2-ocaml libbz2-ocaml-dev #for samtools/htslib-1.5
##Make python bin python2
#
#sudo mv /usr/bin/python /usr/bin/python2
#
#echo python INSTALLED

### PYTHON3
##Python3.10.12
#mypython3=/usr/bin/python3
#if [ ! -f "$mypython3" ]
#then
#	sudo apt-get update
#	sudo apt-get -y install python3.10-dev
##	sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.6 1
##	sudo update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.7 2
##	#Point python3.7 to python3
##	sudo update-alternatives --config python3
#fi
#
###RSeQC Needs Python3 so fluip the current python with the latest version
#sudo cp -a /usr/bin/python3 /usr/bin/python
#
##Needs another update before pip can be installed
#sudo apt-add-repository universe
#sudo apt-get update



## PIP
mypip=/usr/bin/pip
if [ ! -f "$mypip" ]
then
	sudo apt-get -y install python-pip
	sudo pip install --upgrade pip
fi
echo pip INSTALLED

## PIP3
mypip=/usr/bin/pip3
if [ ! -f "$mypip" ]
then
	sudo apt-get -y install python3-pip
fi
echo pip3 INSTALLED



## UNZIP
#myfile=/usr/bin/unzip
#if [ ! -f "$myfile" ]
#then
#	echo installing unzip
#	sudo apt-get -y install unzip
#    if [ ! -f "$myfile" ]
#    then
#		echo could not install unzip
#		echo exiting with error code 1 ...
#		exit 1
#	fi
#fi
#echo unzip INSTALLED



## GIT
#mygit=/usr/bin/git
#if [ ! -f "$mygit" ]
#then
#	sudo apt-get -y install git
#fi
#echo git INSTALLED



## Anaconda latest
sudo mkdir anaconda_install
cd anaconda_install
sudo wget https://repo.anaconda.com/archive/Anaconda3-2023.07-2-Linux-x86_64.sh
sudo bash Anaconda3-2023.07-2-Linux-x86_64.sh
echo conda INSTALLED
cd ..



## Snakemake latest
sudo pip3 install snakemake==7.31.1
echo snakemake INSTALLED


## R latest
#sudo apt install r-base-core
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
sudo apt update
sudo apt install r-base
echo R INSTALLED



# install openxlsx in R
myfile=/usr/local/lib/R/library/openxlsx/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package openxlsx is already installed, so do nothing
else
	echo installing package openxlsx in R
	echo install.packages\(\"openxlsx\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_openxlsx.R
	cat install_openxlsx.R
	sudo R CMD BATCH install_openxlsx.R
	cat install_openxlsx.Rout
	rm --interactive=never install_openxlsx.R install_openxlsx.Rout
fi
echo R package openxlsx INSTALLED



# install stringr in R
myfile=/usr/local/lib/R/library/stringr/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package stringr is already installed, so do nothing
else
	echo installing package stringr in R
	echo install.packages\(\"stringr\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_stringr.R
	cat install_stringr.R
	sudo R CMD BATCH install_stringr.R
	cat install_stringr.Rout
	rm --interactive=never install_stringr.R install_stringr.Rout
fi
echo R package stringr INSTALLED



# install yaml in R
myfile=/usr/local/lib/R/library/yaml/DESCRIPTION
if [ -f "$myfile" ]
then
	echo R package yaml is already installed, so do nothing
else
	echo installing package yaml in R
	echo install.packages\(\"yaml\"\, repos \= \"http\://cran.us.r\-project.org\"\) > install_yaml.R
	cat install_yaml.R
	sudo R CMD BATCH install_yaml.R
	cat install_yaml.Rout
	rm --interactive=never install_yaml.R install_yaml.Rout
fi
echo R package yaml INSTALLED



## pdfrop and other utilities
#sudo apt-get install texlive-extra-utils
#sudo apt-get install poppler-utils



## install LATEX FONTS
## Check for /usr/share/fonts/type1/texlive-fonts-recommended/ and /usr/share/doc/texlive-fonts-extra/
#myfonts1=/usr/share/fonts/type1/texlive-fonts-recommended/
#myfonts2=/usr/share/doc/texlive-fonts-extra/
#if [ ! -d "$myfonts1" ] || [ ! -d "$myfonts2" ]
#then
#	sudo apt-get -y install texlive-fonts-recommended && sudo apt-get -y install texlive-fonts-extra
#fi
#echo Latex fonts are INSTALLED



## CLEANUP
#Apt get
sudo apt-get clean
sudo apt-get autoclean
sudo apt autoremove



## Remove swap
sudo /sbin/swapoff /var/swap.1
sudo rm -f /var/swap.1

## may need to remove old distribute package..
#sudo rm -fr /usr/local/lib/python2.7/dist-packages/distribute*

## Delete software download/pre-installation dir to gain back space on server
rm -rf ~/software_install



echo Script Completed.
##end of script 
exit 0
