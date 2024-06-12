
pmemd.cuda -O -i min1.in -p pro_solv.prmtop -c pro_solv.inpcrd -ref pro_solv.inpcrd -o min1.out -r min1.rst -AllowSmallBox
pmemd.cuda -O -i min2.in -p pro_solv.prmtop -c min1.rst -ref min1.rst -o min2.out -r min2.rst -AllowSmallBox
pmemd.cuda -O -i min3.in -p pro_solv.prmtop -c min2.rst -ref min2.rst -o min3.out -r min3.rst -AllowSmallBox
pmemd.cuda -O -i heat.in -p pro_solv.prmtop -c min3.rst -ref min3.rst -o heat.out -r heat.rst -AllowSmallBox
pmemd.cuda -O -i nvt.in -p pro_solv.prmtop -c heat.rst -ref heat.rst -o nvt.out -r nvt.rst -AllowSmallBox
pmemd.cuda -O -i npt.in -p pro_solv.prmtop -c nvt.rst -ref nvt.rst -o npt.out -r npt.rst -AllowSmallBox
pmemd.cuda -O -i prod.in -p pro_solv.prmtop -c npt.rst -ref npt.rst -o prod.out -r prod.rst -x prod.nc -AllowSmallBox

cpptraj -i cluster.in
