#Script to coarse-grain and run energy minimization

module load cmake/3.23.2 gcc/12.1.0 cuda/12.0.0

#Loop through IDPs
for i in A2 K25 A1 CoRNID ColNT FhuA Hst52 K19 PNt Sic1 aSyn Hst5 ACTR OPN htau40 FUS
do

if [[ "$i" == "A2" ]]
then
	salt=0.005
fi

if [[ "$i" == "A1" ]]
then
        salt=0.05
fi

if [[ "$i" == "htau40" ]]
then
	salt=0.1
fi

if [[ "$i" == "FhuA" || "$i" == "Hst52" || "$i" == "K19" || "$i" == "K25" || "$i" == "PNt" || "$i" == "Hst5" || "$i" == "FUS" || "$i" == "OPN" ]]
then
        salt=0.15
fi

if [[ "$i" == "CoRNID" || "$i" == "Sic1" || "$i" == "ACTR" || "$i" == "aSyn" ]]
then
        salt=0.2
fi

if [[ "$i" == "ColNT" ]]
then
        salt=0.4
fi

echo "$i has ${salt}M salt"

#Starting structure
pdb=${i}.pdb

#Make directory, cp .pdb and insane.py file there and cd there
dir=${i}
cp $insane $dir
cp $pdb $dir
cd $dir
gmx editconf -o min_AA.gro -f $pdb 

#Martinize
martinize2 -f $pdb -o PRO_topol.top -x PRO_CG.pdb -ff martini3001 -ff-dir /fsa/home/ww_zhuy/documents/martini_v300/martini_v3.0.0_proteins/force_fields/

#Put protein in box
gmx editconf -f PRO_CG.pdb -o PRO_CG.gro -bt dodecahedron -d 4.0 <<EOF
1
EOF

#Solvate using insane.py
source activate python27
python insane.py -f PRO_CG.gro -o PRO_SOL_IONS.gro -pbc keep -salt $salt -sol W -center -p PRO_topol_SOL_IONS.top
conda activate py3

#The next few blocks modify the toplogy file and molecule_0.itp file:
#Remove #include include martini.itp and substitute ion names in topology file
perl -pi -e's/#include "martini.itp"//g' PRO_topol_SOL_IONS.top
perl -pi -e's/NA\+/NA/g' PRO_topol_SOL_IONS.top
perl -pi -e's/CL-/CL/g' PRO_topol_SOL_IONS.top

#Add "#include .itp" lines to PRO_topol_SOL_IONS.top
cat <<EOF > others.top
#include "$ffdir/martini_v3.0.0.itp"
#include "PRO.itp"
#include "$ffdir/martini_v3.0.0_ions_v1.itp"
#include "$ffdir/martini_v3.0.0_solvents_v1.itp"
EOF
cat others.top PRO_topol_SOL_IONS.top >a
mv a PRO_topol_SOL_IONS.top

#Run energy minimization
gmx grompp -f ../minimization.mdp -p PRO_topol_SOL_IONS.top -c PRO_SOL_IONS.gro -o min.tpr -pp all_PRO.top -maxwarn 3 -r PRO_SOL_IONS.gro
gmx mdrun -deffnm min -v -ntomp 1 -ntmpi 1 &

#所有蛋白质和水珠之间的 LennardJones 势均按系数 λ 重新调整。 eg. λ =1.06
mkdir lambda_1.12
cd lambda_1.12
python ../PW_rescaling_martini3.py -i ../all_PRO.top -o all_PRO_lambda.top -l 1.12

#relax/ equilibration
#need relax_Martini_298K.mdp/relax_Martini_296K.mdp
module load cuda/12.0.0 openmpi/4.1.2-gcc11.2.0 gcc/12.1.0
mv all_PRO_lambda_1.12.top all_PRO_lambda.top
gmx grompp -f ../relax_Martini_298K.mdp -p all_PRO_lambda.top -c ../min.gro -o relax.tpr -maxwarn 2 -v 
#上传作业mdrun.sh，命令如下
gmx mdrun -s relax.tpr -deffnm relax -ntomp 40 -maxh 47.9 -v
#bsub<mdrun.sh
#gmx mdrun -s relax.tpr -dlb yes -nb gpu -deffnm relax -maxh 47.9 -v -ntmpi 6

#production_run
#need md_Martini_298K.mdp/md_Martini_296K.mdp,改一下模拟数据步长
gmx grompp -f ../md_Martini_298K.mdp -p all_PRO_lambda.top -c relax.gro -t relax.cpt -o prodrun.tpr -maxwarn 2 -v
#上传作业md_mdrun.sh，命令如下
gmx mdrun -deffnm prodrun -s prodrun.tpr -cpi prodrun.cpt -maxh 239.9 -v
#bsub<md_mdrun.sh
#gmx mdrun -s prodrun.tpr -dlb yes -nb gpu -deffnm prodrun -cpi prodrun.cpt -maxh 239.9 -v

#数据分析
gmx trjconv -s prodrun.tpr -f prodrun.xtc -o prodrun_nopbc.xtc -pbc mol -center <<EOF
1 
1 
EOF 
#计算回旋半径
gmx gyrate -s prodrun.tpr -f prodrun.xtc -o aSyn_1.00_gyrate.xvg
gmx gyrate -s prodrun.tpr -f prodrun.xtc -o aSyn_1.10_gyrate.xvg
gmx gyrate -s prodrun.tpr -f prodrun.xtc -o aSyn_1.12_gyrate.xvg

gmx gyrate -s prodrun.tpr -f prodrun.xtc -o K25_1.00_right.xvg
gmx gyrate -s prodrun.tpr -f prodrun.xtc -o K25_1.10_right.xvg
gmx gyrate -s prodrun.tpr -f prodrun.xtc -o K25_1.12_right.xvg


