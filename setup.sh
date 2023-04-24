#!/bin/bash

: '
Bash script to install the following packages:
    - US-align
    - LGA
    - PHENIX
    - rna-tools
    - lDDT

If packages are already installed on the system, script will request directory to create the following symbolic links:
    - bins/us-align -> <path to US-align package>
    - bins/lga -> <path to LGA package>
    - bins/phenix -> <path to PHENIX package>
    - bins/lddt -> <path to lDDT package>

This script will create the following directories:
    - downloads
    - bins
    - bins/us-align
    - bins/lga
    - bins/phenix
    - bins/lddt

For LGA/GDT, script will apply the following patches:
    - patch_runlga.patch
    - unknown_residue_correspondence.patch
'

project_dir="$(pwd)"
mkdir -p downloads bins

mkdir -p bins/us-align
mkdir -p bins/lga

cd downloads

# Install US-align by the Zhang group
read -p "Do you have US-align installed? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
    read -p "Enter the path to the US-align package: " us_align_path
    ln -s $us_align_path $project_dir/bins/us-align
else
    # Install US-align by the Zhang group
    wget https://zhanggroup.org/US-align/bin/module/USalign.cpp
    g++ -O3 -ffast-math -o USalign USalign.cpp

    mv USalign ../bins/us-align
    rm USalign.cpp
fi

# Install LGA
read -p "Do you have LGA installed? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
    read -p "Enter the path to the LGA package: " lga_path
    ln -s $lga_path $project_dir/bins/lga
else
    lga_pkg="LGA_package_src"
    
    echo "$lga_pkg.tar exists!"
    tar -xvf "$lga_pkg.tar"
    
    cd $lga_pkg

    cd SRC
    make -f Makefile.linux_x86_64_msse2
    chmod +x lga
    mv lga ..
    cd ..

    patch -p1 runlga.mol_mol.pl < ../../patches/patch_runlga.patch
    sed -i".backup" "s|SED_BIN_PATH|\"$project_dir/bins/lga/LGA_package_src/\"|g" runlga.mol_mol.pl
    echo "Applied patch_runlga.patch"

    patch -p1 run_GDT_for_structures_with_unknown_residue_residue_correspondences.sh < ../../patches/unknown_residue_correspondence.patch
    sed -i".backup" "s|SED_BIN_PATH|\"$project_dir/bins/lga/LGA_package_src\"|g" run_GDT_for_structures_with_unknown_residue_residue_correspondences.sh
    echo "Applied unknown_residue.correspondence.patch"

    cd ..
    mv $lga_pkg ../bins/lga
fi


# Install PHENIX
read -p "Do you have PHENIX installed? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "Provide the path to to the /build/bin directory in the PHENIX package (e.g. ../phenix-1.18.2-3874/build/bin/)"
    read -p "Enter the path to the PHENIX package: " phenix_path
    ln -s $phenix_path $project_dir/bins/phenix
else
    phenix_pkg="phenix-installer-*.tar.gz"

    if [ -e $phenix_pkg]; then
        tar -xf $phenix_pkg
        ./install --prefix="$project_dir/bins/phenix"
    else
        printf "${phenix_pkg} does not exist.\nPlease download the phenix package and place in downloads/\n" >&2
    fi
fi

# Install rna-tools
if pip3 list | grep "rna-tools"; then
    echo "rna-tools is already installed"
else
    echo "rna-tools is not installed"
    read -p "Do you want to install rna-tools using pip3? (y/n) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]
    then
        pip3 install rna-tools
    fi
fi

# Install lDDT
# TODO: Currently using singularity container, but will use OpenStructure package later
read -p "Do you OpenStructure/lDDT installed? (y/n) " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]
then
    echo "Provide the path to to the lddt directory in the LGA package (e.g. ../lddt)"
    read -p "Enter the path to the lDDT package: " lddt_path
    ln -s $lddt_path $project_dir/bins/lddt
else
    echo "lDDT currently must be installed. Future revisions will enable installation of lDDT."
fi

echo "Finished running setup.sh"