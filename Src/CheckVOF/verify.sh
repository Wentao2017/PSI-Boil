#!/bin/tcsh
### Files requred: *.txt *.gth *.stl, main*.cpp *.gnu *.lay

rm -rf Result_tmp
mkdir Result_tmp
mkdir Result_tmp/SrcVOF
cp ../Equation/Centered/VOF/*.cpp Result_tmp/SrcVOF
cp ../Equation/Centered/VOF/*.h   Result_tmp/SrcVOF
mkdir Result_tmp/SrcPhaseChange
cp ../Equation/Centered/PhaseChange/*.cpp Result_tmp/SrcPhaseChange
cp ../Equation/Centered/PhaseChange/*.h Result_tmp/SrcPhaseChange
cd ..

echo '###############'
echo '### Zalesak ###'
echo '###############'
mkdir CheckVOF/Result_tmp/Zalesak
cp CheckVOF/main-zalesak.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Zalesak
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpirun -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
# tecplot
cp ../../zalesak.gth .
gather.exe < zalesak.gth >& /dev/null
# cp ../../make_png.mcr .
# cp ../../zalesak.lay .
# foreach i (*.plt)
#   tec360 zalesak.lay -b make_png.mcr $i >& /dev/null
#   mv tmp.png $i:r.png
# end
cd ../../..

echo '##W##################'
echo '### Circle vortex ###'
echo '#####################'
mkdir CheckVOF/Result_tmp/CircleVortex
cp CheckVOF/main-circle-vortex.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/CircleVortex
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpirun -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
# tecplot
cp ../../circleVortex.gth .
gather.exe < circleVortex.gth >& /dev/null
cp ../../make_png.mcr .
cp ../../circleVortex.lay .
foreach i (*.plt)
  tec360 circleVortex.lay -b make_png.mcr $i >& /dev/null
  mv tmp.png $i:r.png
end
cd ../../..

echo '###################'
echo '### Corner flow ###'
echo '###################'
mkdir CheckVOF/Result_tmp/CornerFlow
cp CheckVOF/main-corner-flow.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/CornerFlow
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpirun -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
# tecplot
cp ../../cornerFlow.gth .
gather.exe < cornerFlow.gth >& /dev/null
cp ../../make_png.mcr .
cp ../../cornerFlow.lay .
foreach i (*.plt)
  tec360 cornerFlow.lay -b make_png.mcr $i >& /dev/null
  mv tmp.png $i:r.png
end
cd ../../..

echo '################################'
echo '### sliding sphare near wall ###'
echo '################################'
mkdir CheckVOF/Result_tmp/Sliding_Sphare_nearWall
cp CheckVOF/main-slide-sphare-nearWall.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Sliding_Sphare_nearWall
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpirun -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
# tecplot 1
cp ../../slidingSphare.gth .
gather.exe < slidingSphare.gth >& /dev/null
cp ../../make_png.mcr .
cp ../../slidingSphareWall-cont.lay .
tec360 slidingSphareWall-cont.lay -b make_png.mcr >& /dev/null
mv tmp.png slidingSphareWall-cont.png
# tecplot 2
cp ../../slidingSphareWall-iso.lay .
tec360 slidingSphareWall-iso.lay -b make_png.mcr >& /dev/null
mv tmp.png slidingSphareWall-iso.png
cd ../../..

echo '##############################'
echo '### sliding sphare near IB ###'
echo '##############################'
mkdir CheckVOF/Result_tmp/Sliding_Sphare_nearIB
cp CheckVOF/main-slide-sphare-nearIB.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Sliding_Sphare_nearIB
cp ../../../Boil .
cp ../../../main.cpp .
cp ../../floor.stl .
set sec0 = `date +%s`
mpirun -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
# tecplot
cp ../../slidingSphare.gth .
gather.exe < slidingSphare.gth >& /dev/null
cp ../../make_png.mcr .
cp ../../slidingSphareIB-cont.lay .
tec360 slidingSphareIB-cont.lay -b make_png.mcr >& /dev/null
mv tmp.png slidingSphareIB-cont.png
cp ../../slidingSphareIB-iso.lay .
tec360 slidingSphareIB-iso.lay -b make_png.mcr >& /dev/null
mv tmp.png slidingSphareIB-iso.png
cd ../../..

echo '#########################'
echo '### 1D Stefan problem ###'
echo '#########################'
mkdir CheckVOF/Result_tmp/1D-stefan
cp CheckVOF/main-phaseChange-stefan-JCP.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/1D-stefan
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
cat log.txt |grep x-min > front.out
foreach i (*.dat)
  preplot $i >& /dev/null
end
rm *.dat
cp ../../front.gnu .
gnuplot front.gnu >& /dev/null
cp ../../make_png.mcr .
cp ../../1D-stefan.lay .
tec360 1D-stefan.lay -b make_png.mcr >& /dev/null
mv tmp.png 1D-stefan.png
cd ../../..

echo '##########################'
echo '### 1D sucking problem ###'
echo '##########################'
mkdir CheckVOF/Result_tmp/1D-sucking
cp CheckVOF/main-phaseChange-sucking-JCP.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/1D-sucking
cp ../../../Boil .
cp ../../../main.cpp .
cp ../../input.txt .
set sec0 = `date +%s`
Boil < input.txt > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
cat log.txt |grep x-min > front.out
foreach i (*.dat)
  preplot $i >& /dev/null
end
rm *.dat
cp ../../front.gnu .
gnuplot front.gnu >& /dev/null
cp ../../make_png.mcr .
cp ../../1D-sucking.lay .
tec360 1D-sucking.lay -b make_png.mcr >& /dev/null
mv tmp.png 1D-sucking.png
cd ../../..

echo '###############'
echo '### Enright ###'
echo '###############'
mkdir CheckVOF/Result_tmp/Enright
cp CheckVOF/main-enright.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Enright
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpirun -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
# tecplot
cp ../../enright.gth .
gather.exe < enright.gth >& /dev/null
cp ../../make_png.mcr .
cp ../../enright.lay .
foreach i (*.plt)
  tec360 enright.lay -b make_png.mcr $i >& /dev/null
  mv tmp.png $i:r.png
end
cd ../../..

echo '##################'
echo '### TENSION 2D ###'
echo '##################'
mkdir CheckVOF/Result_tmp/Tension2D
cp CheckVOF/main-surfaceTension-2D.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Tension2D
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpirun -np 8 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
cat log.txt |grep velocity > vel.out
cp ../../vel.gnu .
gnuplot vel.gnu >& /dev/null
# tecplot
cp ../../tension2D.gth .
gather.exe < tension2D.gth >& /dev/null
cp ../../make_png.mcr .
cp ../../tension2D.lay .
tec360 tension2D.lay -b make_png.mcr >& /dev/null
mv tmp.png tension2D.png
cd ../../..

echo '##################'
echo '### TENSION 3D ###'
echo '##################'
mkdir CheckVOF/Result_tmp/Tension3D
cp CheckVOF/main-surfaceTension-3D.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/Tension3D
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpirun -np 16 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
cat log.txt |grep velocity > vel.out
cp ../../vel.gnu .
gnuplot vel.gnu >& /dev/null
cat log.txt |grep kappa_min_max > kappa.out
cp ../../kappa.gnu .
gnuplot kappa.gnu >& /dev/null
# tecplot
cp ../../tension3D.gth .
gather.exe < tension3D.gth >& /dev/null
cp ../../tension3D2.gth .
gather.exe < tension3D2.gth >& /dev/null
cp ../../make_png.mcr .
cp ../../tension3D.lay .
tec360 tension3D.lay -b make_png.mcr >& /dev/null
mv tmp.png tension3D.png
cd ../../..

echo '###########################'
echo '### RISING BUBBLE IJNMF ###'
echo '###########################'
mkdir CheckVOF/Result_tmp/RisingBubbleIJNMF
cp CheckVOF/main-risingBubble-IJNMF.cpp main.cpp
make >& /dev/null
cd CheckVOF/Result_tmp/RisingBubbleIJNMF
cp ../../../Boil .
cp ../../../main.cpp .
set sec0 = `date +%s`
mpirun -np 16 Boil > log.txt
set sec1 = `date +%s`
@ difft = $sec1 - $sec0
echo $difft sec.
# gnuplot
cat log.txt |grep totalvol > vol.out
cp ../../vol.gnu .
gnuplot vol.gnu >& /dev/null
cat log.txt |grep front > front.out
cp ../../front-z.gnu .
gnuplot front-z.gnu >& /dev/null
# tecplot
cp ../../risingBubble-uvw.gth .
gather.exe < risingBubble-uvw.gth >& /dev/null
cp ../../risingBubble-xyz.gth .
gather.exe < risingBubble-xyz.gth >& /dev/null
cp ../../make_png.mcr .
cp ../../risingBubble.lay .
tec360 risingBubble.lay -b make_png.mcr >& /dev/null
mv tmp.png risingBubble.png
cp ../../risingBubble-kappa.lay .
tec360 risingBubble-kappa.lay -b make_png.mcr >& /dev/null
mv tmp.png risingBubble-kappa.png
cp ../../risingBubble-side.lay .
tec360 risingBubble-side.lay -b make_png.mcr >& /dev/null
mv tmp.png risingBubble-side.png
cd ../../..


exit
