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
wget https://imputationserver.sph.umich.edu/share/results/99a390ebdb128fd6d833aa7a97ecd743e9ca36b79edfc98e5d442d278d0cad32/chr_1.zip
wget https://imputationserver.sph.umich.edu/share/results/4d332391af74cc5b26ba3514903ea44b34a97135c5d004d2d8a6e331734e56c2/chr_10.zip
wget https://imputationserver.sph.umich.edu/share/results/212ffb04b24129c0af84a5076150faed35427885ad1d1e8c0a9c261a16b89fbf/chr_11.zip
wget https://imputationserver.sph.umich.edu/share/results/d0af5b04d80cecad4ed8339a0f20f8c8abd0419137a290ad264bb4b56838d061/chr_12.zip
wget https://imputationserver.sph.umich.edu/share/results/c32e69f3f28b956b522daf8f1d504ce48215f9e5dd252c426911a17bc312a738/chr_13.zip
wget https://imputationserver.sph.umich.edu/share/results/a4a387b538de3f6c51ff8eadc4e541496d744d06a7156ef1f6c6f14abb030a32/chr_14.zip
wget https://imputationserver.sph.umich.edu/share/results/465736e9a93d41d8316682d5adef3e97dc599ee56aeb0539553d93fa2aad23eb/chr_15.zip
wget https://imputationserver.sph.umich.edu/share/results/7800aee7e236c1e08fe73354a5d1cde37c860c8188d7f250d60fcab14c34c196/chr_16.zip
wget https://imputationserver.sph.umich.edu/share/results/d29253e3a5a5d53f8fc0a93c4c66204fdf93e082065a994ded9c32c69dacf82f/chr_17.zip
wget https://imputationserver.sph.umich.edu/share/results/18f8cd3d12494d5fb33ca9885a971a237df1af95359f3bb5e0a9eb0d7679cf07/chr_18.zip
wget https://imputationserver.sph.umich.edu/share/results/258e87871611f1ba56105e3863124efd549f064b0a6116b49d98ab05463ef637/chr_19.zip
wget https://imputationserver.sph.umich.edu/share/results/9802e0960082268d91705f69ae92dfbf7060d074df1edd27068dd563181d8898/chr_2.zip
wget https://imputationserver.sph.umich.edu/share/results/78852b730be61991c7413e57f0f7743e08244eb3431a1119a64ee67919b07df5/chr_20.zip
wget https://imputationserver.sph.umich.edu/share/results/fa27f540a4eae714c3fc7d835f9ae781764182888fae2a0ce321b65fca470e1e/chr_21.zip
wget https://imputationserver.sph.umich.edu/share/results/5d11e366eec60d58c71967fc62ad90fa20e0d1bf700e3ea0bfd5729727420746/chr_22.zip
wget https://imputationserver.sph.umich.edu/share/results/1192defcef1134dba4b39392c43c69cb6c590c5de8ab2a03025775f8814e4479/chr_3.zip
wget https://imputationserver.sph.umich.edu/share/results/90c5cf629492d4b46ce85099a47fb81ba66bd3f2b7046a640805f9b76d2a9aad/chr_4.zip
wget https://imputationserver.sph.umich.edu/share/results/8f866b7b60bbc3d1c4a425fa3e099ccf60e1332fac4fce3ccd976f7f670ee4e8/chr_5.zip
wget https://imputationserver.sph.umich.edu/share/results/57a23f7a9ebd27b36e2085c5f630d08271bb26158193df5578f7df7def3baead/chr_6.zip
wget https://imputationserver.sph.umich.edu/share/results/3b81c54fae6323eaf4326b5badfce684f187c089b3ab6a5265d8e626a489bc40/chr_7.zip
wget https://imputationserver.sph.umich.edu/share/results/24129eeb85e4b1e1ad49766bb41a9a96cf4e2a76905296735673a80c06658759/chr_8.zip
wget https://imputationserver.sph.umich.edu/share/results/11350a6c393ee13030f20e9566d588378a1a268ef8de767f7d021169626dda5b/chr_9.zip
wget https://imputationserver.sph.umich.edu/share/results/db529a3389c2037301f26647b29908dfc292ba1ee527e165d626ef42de809a71/results.md5

wget https://imputationserver.sph.umich.edu/share/results/e76c69ddce5fcc17d7ca40faf2fa92f064fc8c6aa35d7d1262c2a6239f86caec/chunks-excluded.txt
wget https://imputationserver.sph.umich.edu/share/results/47e2333dfb905ff31693327a39c576ff6c8cd3301dd8e5a07e8c709fdbbb5359/snps-excluded.txt
wget https://imputationserver.sph.umich.edu/share/results/7d8431d8eca62d99c158887b58db999e14f4c1c392deccb06ac0ef29cd2a66e2/qcreport.html

wget https://imputationserver.sph.umich.edu/share/results/68ba95db8cbc68f9e5c5335c96ffa0e19b20d476b1431793751640c0f4abf15e/chr_1.log
wget https://imputationserver.sph.umich.edu/share/results/85682148e88eea0787cf53a14daa70de5e0260a2c6fa855396124484ee7046cb/chr_10.log
wget https://imputationserver.sph.umich.edu/share/results/e9dc3f49d64da0b67d411d251fc0119d769327417c0eab96d4b33216ebc30ab6/chr_11.log
wget https://imputationserver.sph.umich.edu/share/results/6831b6630523957e6164fce8b65e82012c73734712704dc5513a84db9d5e3008/chr_12.log
wget https://imputationserver.sph.umich.edu/share/results/228dc72eaf588459ed497435f99327303c521b6b58bf1546a3c6d7e29e823607/chr_13.log
wget https://imputationserver.sph.umich.edu/share/results/d5ad16eae68df5ab02a40fa58c28ff23ba91a994e08559a66a19403f5fdf6a1a/chr_14.log
wget https://imputationserver.sph.umich.edu/share/results/666043f81ade184ba73c611c83e41a053f3018aa4cb840137633899f09efc715/chr_15.log
wget https://imputationserver.sph.umich.edu/share/results/c9add6266cd64e33503abee3a1685e90455db7ca71180a6e32e9496b274daa90/chr_16.log
wget https://imputationserver.sph.umich.edu/share/results/2449639a0b68032c8dadd4c8dd9d558e855624f8a8bc63b16e0ee1e59f2dab33/chr_17.log
wget https://imputationserver.sph.umich.edu/share/results/13a724fb282179e059b2344d88f7cff8e16a0ab9bde615083eae6b5f32b8dd50/chr_18.log
wget https://imputationserver.sph.umich.edu/share/results/c300222cf04020a1655e71283f0394b8743974e5b69965701443c11bc27926eb/chr_19.log
wget https://imputationserver.sph.umich.edu/share/results/c0faf5c9522ab7f2349181da6db54b4c7b34e8212318d14ff2e7bf6535c33e50/chr_2.log
wget https://imputationserver.sph.umich.edu/share/results/30e0a842f1d52d711703de92b3b14788f675eeb31aeab2caa347f6edefa10e6c/chr_20.log
wget https://imputationserver.sph.umich.edu/share/results/c9feb51c0af55a1a2e3c6d5d00a55509507288aa5781223c23d3e23cd2c04294/chr_21.log
wget https://imputationserver.sph.umich.edu/share/results/0effbfab0e4b34f48b75602633d84cacd3e70f3d18064283955b2c379f0838df/chr_22.log
wget https://imputationserver.sph.umich.edu/share/results/869da129878c683a8c86b21a98b254d028f9e9746149e55c698fb2a1bee70efd/chr_3.log
wget https://imputationserver.sph.umich.edu/share/results/4b5bf5099f2e797b92bef71aee1e8e0f2aaaf168decae7e48bbb0452e3555cab/chr_4.log
wget https://imputationserver.sph.umich.edu/share/results/ec0c259a0bc5dd9f8718acde8bc1f90a147f76bf9fee70e7e6ccee0442642eb3/chr_5.log
wget https://imputationserver.sph.umich.edu/share/results/b020d39ec691d49c9fd02df60308285201e737e991b06d5ff72c5a090407ed18/chr_6.log
wget https://imputationserver.sph.umich.edu/share/results/8f04f0b6262fcd9dae6a680a024a95a586441c7540cf597e0e11a3b12f878e1f/chr_7.log
wget https://imputationserver.sph.umich.edu/share/results/b7d9532e7046bc3405aed577170546c81a40f51bde103db5e09c7895ad94ced7/chr_8.log
wget https://imputationserver.sph.umich.edu/share/results/62228ebd1edba5c6b60681dad4d841f93f7097a7de5a8d82a16341e9e12d1fb2/chr_9.log


# passcode: Sz7bgjIH0bPEkK

