# MV-RTSC

The source code for our paper "Multi-View Robust Tensor-based Subspace Clustering".

Available [here](https://ieeexplore.ieee.org/document/9999185).


# Required toolboxes:
The MATLAB code for MV-RTSC needs some toolboxes:

[1] Tensor decomposition, Xie, Y., Tao, D., Zhang, W., Liu, Y., Zhang, L., & Qu, Y. (2018). 
"On unifying multi-view self-representations for clustering by tensor multi-rank minimization."

[2] Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>

# Data set:
All datasets used in our paper are available at [Google Drive](https://drive.google.com/drive/folders/1T8AITuQZCbCB52PELprIVb75d4-6hi20). 
Each dataset is a mat file containing 2 variables fea or data (i.e., a MATLAB cell of features) and gt (i.e., ground truth label), 
except the file flower17.mat which contains a cell of distance matrices and ground truth since features are unavailable. 

# Citation:

```
@article{9999185,  
author={Al-Sharoa, Esraa M. and Al-Wardat, Mohammad A.},  
journal={IEEE Access},   
title={Multi-View Robust Tensor-Based Subspace Clustering},   
year={2022},  volume={10},  number={},  pages={134292-134306},  
doi={10.1109/ACCESS.2022.3232285}
}
```
