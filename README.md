# Cost-sensitive multi-label feature selection

The repository contains the implementation of cost-constrained multi-label feature selection method. The aim of cost-constrained multi-label feature selection is to select a feature subset relevant to multiple labels while satisfying a user-specific maximal admissible budget. This approach allows for building a model with high predictive power, for which the cost of making a prediction for a single instance does not exceed the user-specified budget. The method is based on a novel criterion that combines the relevance and cost of the candidate feature. The relevance measure is derived using the lower bound of the mutual information between the feature subset and label vector. The method includes effective way of determining the cost-factor parameter that controls the trade-off between relevancy and cost.

**How to use?**

Just follow the instructions in the script `experiments.Rmd` which You can find in the root directory of this repository.
