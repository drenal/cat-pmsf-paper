# Datasets and scripts for the paper Compositionally constrained sites drive long branch attraction

In [this repository](https://github.com/drenal/cat-pmsf-paper) all the necessary code and datasets have been uploaded to reproduce the results described in the paper Compositionally constrained sites drive long branch attraction and to run the introduced CAT-PMSF pipeline on arbitary datasets.

Structure of the repository:
- [datasets](datasets/) : empirical datasets analysed in the paper
- [datasets/simulation](datasets/simulation/) : simulation dataset
- [scripts](scripts/) : scripts necessary to perform CAT-PMSF pipeline on a dataset
- [step1_iqtree_lg](step1_iqtree_lg/) : results of CAT-PMSF's 1st step applied on the empirical datasets
- [step1_iqtree_lg/simulation](step1_iqtree_lg/simulation/) : correct (good) and incorrect (bad) topology for simulation pipeline
- [step2_pb](step2_pb/) : results of CAT-PMSF's 2nd step applied on the empirical datasets
- [step2_pb/simulation](step2_pb/simulation/) : results of CAT-PMSF's 2nd step applied on the simulated trees
- [step3_iqtree](step3_iqtree/): : results of CAT-PMSF's 3rd step applied on the empirical datasets
- [step3_iqtree/simulation](step3_iqtree/simulation/): : results of CAT-PMSF's 3rd step applied on the simulated trees
