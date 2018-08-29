#!/bin/sh
OUTPUT=$1
binFile=$2
#This script downlds the ParSy(SC18) matrix set and extract it to output folder and then generates .
wget https://sparse.tamu.edu/MM/AMD/G3_circuit.tar.gz -O ${OUTPUT}G3_circuit.tar.gz
wget https://sparse.tamu.edu/MM/McRae/ecology2.tar.gz -O ${OUTPUT}ecology2.tar.gz
wget https://sparse.tamu.edu/MM/Schmid/thermal2.tar.gz -O ${OUTPUT}thermal2.tar.gz
wget https://sparse.tamu.edu/MM/GHS_psdef/apache2.tar.gz -O ${OUTPUT}apache2.tar.gz
wget https://sparse.tamu.edu/MM/Janna/StocF-1465.tar.gz -O ${OUTPUT}StocF-1465.tar.gz
wget https://sparse.tamu.edu/MM/Janna/Hook_1498.tar.gz -O ${OUTPUT}Hook_1498.tar.gz 
wget https://sparse.tamu.edu/MM/CEMW/tmt_sym.tar.gz -O ${OUTPUT}tmt_sym.tar.gz
wget https://sparse.tamu.edu/MM/Janna/PFlow_742.tar.gz -O ${OUTPUT}PFlow_742.tar.gz
wget https://sparse.tamu.edu/MM/Janna/Flan_1565.tar.gz -O ${OUTPUT}Flan_1565.tar.gz
wget https://sparse.tamu.edu/MM/GHS_psdef/audikw_1.tar.gz -O ${OUTPUT}audikw_1.tar.gz
wget https://sparse.tamu.edu/MM/Oberwolfach/bone010.tar.gz -O ${OUTPUT}bone010.tar.gz
wget https://sparse.tamu.edu/MM/Botonakis/thermomech_dM.tar.gz -O ${OUTPUT}thermomech_dM.tar.gz
wget https://sparse.tamu.edu/MM/Janna/Emilia_923.tar.gz -O ${OUTPUT}Emilia_923.tar.gz
wget https://sparse.tamu.edu/MM/Janna/Fault_639.tar.gz -O ${OUTPUT}Fault_639.tar.gz
wget https://sparse.tamu.edu/MM/GHS_psdef/bmwcra_1.tar.gz -O ${OUTPUT}bmwcra_1.tar.gz
wget https://sparse.tamu.edu/MM/ND/nd24k.tar.gz -O ${OUTPUT}nd24k.tar.gz
wget https://sparse.tamu.edu/MM/ND/nd12k.tar.gz -O ${OUTPUT}nd12k.tar.gz



for f in ${OUTPUT}*.gz; do tar xvf $f -C ${OUTPUT}; done;
#copying all matrices into root so it is easier to explore later.
cp ${OUTPUT}G3_circuit/G3_circuit.mtx ${OUTPUT}G3_circuit.mtx;
cp ${OUTPUT}ecology2/ecology2.mtx ${OUTPUT}ecology2.mtx;
cp ${OUTPUT}thermal2/thermal2.mtx ${OUTPUT}thermal2.mtx;
cp ${OUTPUT}apache2/apache2.mtx ${OUTPUT}apache2.mtx;
cp ${OUTPUT}StocF-1465/StocF-1465.mtx ${OUTPUT}StocF-1465.mtx;
cp ${OUTPUT}Hook_1498/Hook_1498.mtx ${OUTPUT}Hook_1498.mtx;
cp ${OUTPUT}tmt_sym/tmt_sym.mtx ${OUTPUT}tmt_sym.mtx;
cp ${OUTPUT}PFlow_742/PFlow_742.mtx ${OUTPUT}PFlow_742.mtx;
cp ${OUTPUT}Flan_1565/Flan_1565.mtx ${OUTPUT}Flan_1565.mtx;
cp ${OUTPUT}audikw_1/audikw_1.mtx ${OUTPUT}audikw_1.mtx;
cp ${OUTPUT}bone010/bone010.mtx ${OUTPUT}bone010.mtx;
cp ${OUTPUT}thermomech_dM/thermomech_dM.mtx ${OUTPUT}thermomech_dM.mtx;
cp ${OUTPUT}Emilia_923/Emilia_923.mtx ${OUTPUT}Emilia_923.mtx;
cp ${OUTPUT}Fault_639/Fault_639.mtx ${OUTPUT}Fault_639.mtx;
cp ${OUTPUT}bmwcra_1/bmwcra_1.mtx ${OUTPUT}bmwcra_1.mtx;
cp ${OUTPUT}nd24k/nd24k.mtx ${OUTPUT}nd24k.mtx;
cp ${OUTPUT}nd12k/nd12k.mtx ${OUTPUT}nd12k.mtx;

