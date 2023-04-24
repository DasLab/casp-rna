#!/bin/bash

phenix_location=PHENIX_LOCATION

pdb=$1

${phenix_location}/phenix.clashscore $pdb nuclear=True keep_hydrogens=True > "${pdb}.out.tmp"
score=$(awk '/^clashscore = /{print $NF}' "${pdb}.out.tmp")
rm "${pdb}.out.tmp"

echo " $score" > "clash_parallel/${pdb}.out"
