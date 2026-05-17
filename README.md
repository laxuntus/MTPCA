# MTPCA
This repository holds the R-code related to Mixed-type principal component analysis (MTPCA) proposed in the article Model-based sparse mixed-type PCA by Heinonen and Virta

Mixed-type data consists of continuous, binary, integer-valued and positive continuous variables. The data is assumed to come from a probability model, where the parameters of the exponential family distributions come from a mixtures of gaussian latent variables. The proposed method, MTPCA, is based on estimating the covariance matrix of these latent mixtures through the method of moments. A way to sparsify the component loadings is implemented using [SICS](https://github.com/laxuntus/SICS).

## Reference

Lauri Heinonen & Joni Virta (2026). Model-based sparse mixed-type PCA. Preprint.
