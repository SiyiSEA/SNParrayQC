#!/bin/sh
#SBATCH --export=ALL # export all environment variables to the batch job.
#SBATCH -p mrcq # submit to the serial queue
#SBATCH --time=24:00:00 # Maximum wall time for the job.
#SBATCH -A Research_Project-MRC190311 # research project to submit under. 
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion 
#SBATCH --output=/lustre/home/sww208/QC/OrganizedSNParray/5_JobReports/05a_downloadMichigan.o
#SBATCH --error=/lustre/home/sww208/QC/OrganizedSNParray/5_JobReports/05a_downloadMichigan.e
#SBATCH --job-name=downloadMichigan

source ./config

cd ${IMPUTEDIR}

# mkdir -p ImputationOutputMichigan

# cd ImputationOutputMichigan || exit

# data for HRC
# wget https://imputationserver.sph.umich.edu/share/results/0638a99edaa602e79d6aaab2d73e4963e2492a55dae76b024d2233295f90bc4a/chr_1.zip
# wget https://imputationserver.sph.umich.edu/share/results/f1af468bb753a08942e74a01dc49ea02e83555a6e684af6bdd165d4179d5460c/chr_10.zip
# wget https://imputationserver.sph.umich.edu/share/results/e19e8259c34de11a49778e3bacc647b178d62e8c725454657d3af084e513f339/chr_11.zip
# wget https://imputationserver.sph.umich.edu/share/results/0bf7beb2b032591bcbf64f189e20ef2b34afc51dc0ff4ce6c82b9128a537e7ce/chr_12.zip
# wget https://imputationserver.sph.umich.edu/share/results/27da8c2b272ddce586cc1d2cc4b5d7d5428e84d3527177315de82d449c838752/chr_13.zip
# wget https://imputationserver.sph.umich.edu/share/results/cf8fb8d441d35e5039ecfd4f98babf7031ef616c1654db576be30e1bfe8b8960/chr_14.zip
# wget https://imputationserver.sph.umich.edu/share/results/438639ffe82cf310e26321bc283e396b363ff0ecd707ab56e135b74a1546a57f/chr_15.zip
# wget https://imputationserver.sph.umich.edu/share/results/f5d141f87e4b6495b1816bec03186234c702717f00b87949b636e35697ccc7b5/chr_16.zip
# wget https://imputationserver.sph.umich.edu/share/results/ec09e4c064a405c3b507bf32673ddb9880caf871fb44ace4c68753d71d426241/chr_17.zip
# wget https://imputationserver.sph.umich.edu/share/results/af4858003a46399f5dc6d3da4d3113b0d6449abd2a11dee9e811f63709d5b2cf/chr_18.zip
# wget https://imputationserver.sph.umich.edu/share/results/41978b7ff5f49ff35f0d53602b1fac8cb2894fb1fa46766555d31791164c7f9d/chr_19.zip
# wget https://imputationserver.sph.umich.edu/share/results/e1d08c440c4c69c81bd8d55cbecd36945bed5a8542876a455ff4b203649aee1f/chr_2.zip
# wget https://imputationserver.sph.umich.edu/share/results/a2148a49ffd92332dbb691d602f9ff7b3676345fdf50c9b4afef4415ee3cf18d/chr_20.zip
# wget https://imputationserver.sph.umich.edu/share/results/e92f9707324e52c2764647cb264f7f200ab22a0cae1f517559bcfe154420b77b/chr_21.zip
# wget https://imputationserver.sph.umich.edu/share/results/055d22d27ec1e571c228e698968d6253ce37b4299761ecc281ae61a8b65bf9e9/chr_22.zip
# wget https://imputationserver.sph.umich.edu/share/results/67323f5ecb7bd72a0b2ce008cefc67d9115a61e942b14991cd1e913326a00d1e/chr_3.zip
# wget https://imputationserver.sph.umich.edu/share/results/4e1ecd8c848f4c49d7a10c6c96d47ae1c4c0ce42949bf23b177eb4adf659e502/chr_4.zip
# wget https://imputationserver.sph.umich.edu/share/results/a87dcea79a44ae0f11b69d224f7a9cc99d81eb41ae0c49fe8f560ef2c150ce68/chr_5.zip
# wget https://imputationserver.sph.umich.edu/share/results/e54245cbf9f8d8a2318ca4826ffadcd076ecc856128c09876427bc25ffc6e612/chr_6.zip
# wget https://imputationserver.sph.umich.edu/share/results/c18b56912e2e654bf931746f350f44a37579d6a8e064b680874e80e7980707ce/chr_7.zip
# wget https://imputationserver.sph.umich.edu/share/results/eab3283eb5c35380f96d795fe7707deb7f45daa61ddf1e11cc7ed54c48883816/chr_8.zip
# wget https://imputationserver.sph.umich.edu/share/results/767287be5e1ebf8daddc8971ff20d47177ea8bcc4c631904cf0b1871c6205485/chr_9.zip
# wget https://imputationserver.sph.umich.edu/share/results/7f320cfc81546f3ef167515cdd489c1aef83fe992079c57eabd310d91452664e/chr_X.zip
# wget https://imputationserver.sph.umich.edu/share/results/01912dbba25adf85b8f30cf1b3a98684094e4f59db9f043a46c60645eca8aadb/results.md5

# wget https://imputationserver.sph.umich.edu/share/results/c531e45dcc19411546d4376995fe6d1f72480dde9d65ca6715d6ced3788aaba1/qcreport.html
# wget https://imputationserver.sph.umich.edu/share/results/f35cde9d64737484cc813f15e3c0824d8b16772e35f7a5fd2c5ab7efb8f596cb/chunks-excluded.txt
# wget https://imputationserver.sph.umich.edu/share/results/61bf07c8720229f845bf20d514d69572c3092e4f58ad507995cd5b3bad530d34/snps-excluded.txt

# passcode: KwrP350XfATza

# 2024 May
mkdir -p ImputationOutputMichigan1000G

cd ImputationOutputMichigan1000G || exit

# wget https://imputationserver.sph.umich.edu/share/results/318b1c08c6457a2ec324902e77cd0ab77332797a574f05f6fe42f8ad11714767/qcreport.html
# wget https://imputationserver.sph.umich.edu/share/results/50b64013ba06fb87c492254d2417eaf9b151476a23b9df5750f9bd23a1e7a28d/chunks-excluded.txt


# wget https://imputationserver.sph.umich.edu/share/results/24823c176ba554a88b0d97dd0d9b0201c0223be9ef25e80e7a905c14c62f8ac6/chr_1.zip
# wget https://imputationserver.sph.umich.edu/share/results/87aab286807aaa25e7c2eb207cead169f9d8853e71fc00ea44f21fbfcb410a42/chr_10.zip
# wget https://imputationserver.sph.umich.edu/share/results/797cc9ae1e925d8e8f9d08cadb28daf1944a3b9212722943a767126052be5093/chr_11.zip
# wget https://imputationserver.sph.umich.edu/share/results/d1cf7368cf119fb6675f1abf7f9d6c67048223e58034cb4b947a6568137136c8/chr_12.zip
# wget https://imputationserver.sph.umich.edu/share/results/32b9ca349b08bbffaf7c1de386aeefbbbdea7764f24c6ffed6985beba7b675b8/chr_13.zip
# wget https://imputationserver.sph.umich.edu/share/results/87b4eece1e179077019c53e93d928d1cbc0bafd331b47d0b8585427a57308ccc/chr_14.zip
# wget https://imputationserver.sph.umich.edu/share/results/974078e8d2d9fbad89a30b46dda1e4843ed9a1adda27ff4075322fa76a3df199/chr_15.zip
# wget https://imputationserver.sph.umich.edu/share/results/fc96300ca7b94fb669004d32efca738b5f89db1f85f14cffc34835bb64ce4e68/chr_16.zip
# wget https://imputationserver.sph.umich.edu/share/results/0ee1ccd26b5ec60c1709b0cf8bea6fcf4bd0da2620d939e9538018b785386d41/chr_17.zip
# wget https://imputationserver.sph.umich.edu/share/results/3ca5d9bcf31bedab74e361b2deb0785011c657f25a4828856b31dc7c39caf6e9/chr_18.zip
# wget https://imputationserver.sph.umich.edu/share/results/d3977ce809ad826f97f359f0a19fbc7d26443232512c69d15e3a6538c7277da5/chr_19.zip
# wget https://imputationserver.sph.umich.edu/share/results/74787751acf203e396f2db083532837359060f926f1c306a701448e1002fd653/chr_2.zip
# wget https://imputationserver.sph.umich.edu/share/results/d4c5d1765e2c3e77d1ba9c0532abab02eb2ff71dddc1f1bca49b498f9b090fc1/chr_20.zip
# wget https://imputationserver.sph.umich.edu/share/results/3859da9d03140bcea9459a5002458a6f42823b6d6e55b07bc73c97ad8a183f51/chr_21.zip
# wget https://imputationserver.sph.umich.edu/share/results/54c852267a7ee3e3efe3ecb2bf540721ef116672939cd6b0602be5339a203775/chr_22.zip
# wget https://imputationserver.sph.umich.edu/share/results/83ef3fc4576548bed60a40843f5879048f7241923b50158ca5428e75c99fc11e/chr_3.zip
# wget https://imputationserver.sph.umich.edu/share/results/b9de74abd832c1d8067d5c5535d336b8734332465d367f8dbd53731e5907d0b6/chr_4.zip
# wget https://imputationserver.sph.umich.edu/share/results/df6664c1582ef38ea47e67236948c459362b01b4d333b9e03278b1ae3ab4fc6d/chr_5.zip
# wget https://imputationserver.sph.umich.edu/share/results/a110d579828534308e389bb8d7eb568815e97d9c2c0bab827335373cbb4cd2cf/chr_6.zip
# wget https://imputationserver.sph.umich.edu/share/results/e389b344e9420f5be2ae99d7addb3b4c71663f94f8a179105683702d66a5413f/chr_7.zip
# wget https://imputationserver.sph.umich.edu/share/results/82ad27526bd95f4738613a5c2a2bf2668c2a8713d2298a548a0e628a9f55e8a2/chr_8.zip
# wget https://imputationserver.sph.umich.edu/share/results/272162651f2b8cfdba858817210bb7e3629ea26efdd7ab4b02fee85d6530caa8/chr_9.zip
# wget https://imputationserver.sph.umich.edu/share/results/afbd73157c8a17b7bef127563ec211667c156e829b354e21dde55bd2fffec329/chr_X.zip
# wget https://imputationserver.sph.umich.edu/share/results/41c9d8437fcae742754484d693e22c3b98c68919f3abbc38e9f4f556a2504fd1/results.md5

# wget https://imputationserver.sph.umich.edu/share/results/0c78758335c9b108a282cf8053a88eeece4139b48722f107b0cfda42dc1304e9/chr_1.log
# wget https://imputationserver.sph.umich.edu/share/results/f6fc3db750bde4eaa43c2b059744fd4117ec39cf21782139c0486976cc7cdca4/chr_10.log
# wget https://imputationserver.sph.umich.edu/share/results/41c0b219484d96f0fcd0b7cb4db2ce3fea890f2eebb48ca10e4d5430f3d73f0b/chr_11.log
# wget https://imputationserver.sph.umich.edu/share/results/44ae6b093663f3e557a6db3a102508b13578e26fa64eb8103668842823fb7038/chr_12.log
# wget https://imputationserver.sph.umich.edu/share/results/341c199864e832950e296b55fbe7df35d56cef12e003b7ec792fbae996cf5a48/chr_13.log
# wget https://imputationserver.sph.umich.edu/share/results/965ecf48e1db98e48b36d64863b1638d0ce362e1b9118ceca36ba6161b716256/chr_14.log
# wget https://imputationserver.sph.umich.edu/share/results/10e7df924528df2cf14705b6081ea4b1e7fd625039a1d08e796211c4af93af98/chr_15.log
# wget https://imputationserver.sph.umich.edu/share/results/15cc784e7d4e5c075a8b2ef89fc3fe9914c7abbdfbde37c40c47c8c1e0796e2b/chr_16.log
# wget https://imputationserver.sph.umich.edu/share/results/aefb2983a1e4a8e18bb6dc9f2e550def7fa092d4b3c3fefe30b5440670fefac3/chr_17.log
# wget https://imputationserver.sph.umich.edu/share/results/b8f56b51b44a51a5aec29386fa4f6b12963d6d0f8d7c1534bea33b4541296fc3/chr_18.log
# wget https://imputationserver.sph.umich.edu/share/results/50dcefe1cd2a246eea79852850262cdf277fd50a555d52cd3ba4cd6956082bc5/chr_19.log
# wget https://imputationserver.sph.umich.edu/share/results/9a51739d693c83134d575b356b0ab9bac87ec4c56e6e46a11d04a3e7c430660b/chr_2.log
# wget https://imputationserver.sph.umich.edu/share/results/245f42f82653c182584fe1baba5803ab247083bbc0a8a47abd7a5d9cc94038aa/chr_20.log
# wget https://imputationserver.sph.umich.edu/share/results/d60ca71b03e82c93e850af99a868de6911bd406d8da5214cdb4e03fd2ac83f2c/chr_21.log
# wget https://imputationserver.sph.umich.edu/share/results/b4fa89d832074e484ff43b4da492f530b52cd073bad7121e500a4f4ee3e60d2c/chr_22.log
# wget https://imputationserver.sph.umich.edu/share/results/39f6c5786aeb7875b0795071f8adba92d501f6ea60a0f16d2f35f85857e6edf4/chr_3.log
# wget https://imputationserver.sph.umich.edu/share/results/9514f4282c03023f779b2730c7f660b0ca2f4651c5e3cbc8b8be9d7c75b49976/chr_4.log
# wget https://imputationserver.sph.umich.edu/share/results/2e0cbe5b02aabc7ae079228463acacee392d2dc58b77c352bf65681ff573d692/chr_5.log
# wget https://imputationserver.sph.umich.edu/share/results/e775b9c36c40d4e54c176814be1385d335e07eb5ceb308ded654808c5044dbe1/chr_6.log
# wget https://imputationserver.sph.umich.edu/share/results/b47da7dcff0727650f1d12be1ae8f3804ef5e2478e71d4d2f037cf7f3f0ea603/chr_7.log
# wget https://imputationserver.sph.umich.edu/share/results/630b658d20e894dacac5f2145b8eea22fdceca680be53e9794568b84d5e0d985/chr_8.log
# wget https://imputationserver.sph.umich.edu/share/results/f61e612ececb59cafd549c0bb8e2f4ddb8319237ad952e03ce10153715df67ad/chr_9.log
# wget https://imputationserver.sph.umich.edu/share/results/2cf1cea30e029b7f21faf38e7ad3a8ced8186dbc912bd29c0686da6c28004f49/chr_X.nonPAR.log
wget https://imputationserver.sph.umich.edu/share/results/74787751acf203e396f2db083532837359060f926f1c306a701448e1002fd653/chr_2.zip
# the password for the imputation results is: 7.epaVdlLA9GuS