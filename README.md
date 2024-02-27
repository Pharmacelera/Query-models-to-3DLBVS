# Query-models-to-3DLBVS

This folder includes :

-> **An ECPF4 filtered version of the DUD-E+ dataset used in VS benchmarking is available.**

-> **Alignment-independent conformational clustering protocol (clustering_step.py)**

    1.  For each sampled conformer, retrieve all rotatable bonds through a SMARTS pattern.

    2.  Calculate a dihedral angle for each rotatable bond and encode them into this trigonometric tuple: {cos(θi), sin(θi)} where*i* denotes the *i-th* rotatable bond. The encoding step enables to consider the periodic/circular distribution of angles, where an angle of –179 degrees is close to angle of +179 degrees.

    3.  The dihedral vectors of all conformers are projected into their PCA components to reduce their dimensionality, assuring that at least 75% of the variance is retained.

    4.  K-means algorithm runs on the projected dihedral space, spanned by the PCA axes, with a fixed number of clusters. Then the lowest energy conformation of each cluster is selected, this procedure attempts to locate all local minima in the conformational space, although it is not guaranteed. In our case, 5 representative conformations are extracted.
