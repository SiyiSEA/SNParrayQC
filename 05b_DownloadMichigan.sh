#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=05bdownloadMichigan.o
#SBATCH --error=05bdownloadMichigan.e
#SBATCH --job-name=downloadMichigan

echo "checking the arguments for config file----------------------------------------------------------------------------"
datapeth=$1

if [ -z "$1" ]
then
        echo "No argument supplied"
        echo "Please input the paht of the data folder as the first argument"
		exit 1 # fail
fi

echo "running the PostQCSanger at $datapeth"
source ${datapeth}/config

mv 05bdownloadMichigan.o ${JOBSDIR}/05bdownloadMichigan.o
mv 05bdownloadMichigan.e ${JOBSDIR}/05bdownloadMichigan.e

cd ${IMPUTEDIR}/ImputationOutputMichiganHRC || exit

# data for HRC
wget https://imputationserver.sph.umich.edu/share/results/8425007043dd7f5683ed7ac2214e256f83b0cddede3c794ec1994339f271b96e/chr_1.zip
wget https://imputationserver.sph.umich.edu/share/results/0bc56332f2cf31a8ecf17cdb74fbb77db74eed785b4f8e2a5d13d1155372a081/chr_10.zip
wget https://imputationserver.sph.umich.edu/share/results/bd3c74c1f63248a61be85ceb96934fe3e18e8bd876a46cb87f1a46e649c5502e/chr_11.zip
wget https://imputationserver.sph.umich.edu/share/results/da4dc39313c47fd3021584435c2e7705d1c182b459ae81d03030ed6feb78dbf3/chr_12.zip
wget https://imputationserver.sph.umich.edu/share/results/8e324abde7dd729dbde05a33c86b142f15198c2fa6d60548fc94c2a605d11c08/chr_13.zip
wget https://imputationserver.sph.umich.edu/share/results/7ff8bb5c12b6fd789746e58ef4bd6f00b5083710a77039bbf3ef7d69db76920d/chr_14.zip
wget https://imputationserver.sph.umich.edu/share/results/8d4ef845dbb7cc82c6af2c3220ccd7cf111c56608036b146076a7f913fdee8f0/chr_15.zip
wget https://imputationserver.sph.umich.edu/share/results/c7d65b8e917eecab9f636ff3ad47271ea75ffd67ee0207f5921a462f786444e6/chr_16.zip
wget https://imputationserver.sph.umich.edu/share/results/3dc6a52f1b323b9c8245b646a035f7d9773dc40a64e88c9c7834063973cff05d/chr_17.zip
wget https://imputationserver.sph.umich.edu/share/results/f08db44a3478ae2cd8b31812195915c56a65ae68bba0a8bdd40815243be0c294/chr_18.zip
wget https://imputationserver.sph.umich.edu/share/results/1c2bfe90a6e16ae1d12a8ac099856f30625bfbd9423dbe2f7cc2afd970fec5ac/chr_19.zip
wget https://imputationserver.sph.umich.edu/share/results/10557b158c60122dd718254ad910ec653ab169047206841a57b2f7b76084c98c/chr_2.zip
wget https://imputationserver.sph.umich.edu/share/results/78abea29a36d517aa7cbe0cfcb80626bb68faa65c83953ca06cf6be75f89ed2f/chr_20.zip
wget https://imputationserver.sph.umich.edu/share/results/20c9f22a106c3f07bc3afca470e421adc323de66b405b28e1ce38322252a6a53/chr_21.zip
wget https://imputationserver.sph.umich.edu/share/results/96c83550ba2684aaf56c624edeedd9c00d35f648da88a0346e0e63becb65f5c3/chr_22.zip
wget https://imputationserver.sph.umich.edu/share/results/5338db79ead45d05d6a694671200194ed0929d5882901ade0fed30f07ba348da/chr_3.zip
wget https://imputationserver.sph.umich.edu/share/results/c175595b7c95ce404e9b94e8e33ce5e55c77c1b6b911d69a12e36c135361ef92/chr_4.zip
wget https://imputationserver.sph.umich.edu/share/results/5ccfbf924d54aca4256df4c2b807d8350bae8177cb4d814d94a89f3018c00556/chr_5.zip
wget https://imputationserver.sph.umich.edu/share/results/71b70d7db2ff20d0ab8dad2c25c9ced974a440d067699f7ff9dd0299eb4e7988/chr_6.zip
wget https://imputationserver.sph.umich.edu/share/results/c79dcd064a8b897588aa1ccad6e4bc789b4b6c845465e0d0650c07a74fde51be/chr_7.zip
wget https://imputationserver.sph.umich.edu/share/results/1b0300443b379e01794a2b5a0431906bd6e5d79c0d4a22b86b7dbbda6c4df38f/chr_8.zip
wget https://imputationserver.sph.umich.edu/share/results/405bb4fc8af01ceb44a482e6fd1ce168eef229b419127896ad509aaf14095e4d/chr_9.zip
wget https://imputationserver.sph.umich.edu/share/results/288cff741c913cc4abe5e24b5f7bfe0b699a9f59948c3bf83739ff776df84690/qc_report.txt
wget https://imputationserver.sph.umich.edu/share/results/9e12efa8c820c1ec67df13a5799aea7f406cd9b576b03bc1845c6c654d9857aa/quality-control.html
wget https://imputationserver.sph.umich.edu/share/results/9f950ad82def2fd8ac861d595a522b95733d513a4d774b379412e2414e3fd792/statistics/chunks-excluded.txt

# passcode: ieZM%4]C7VxUht

