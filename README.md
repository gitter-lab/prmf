# Pathway-Regularized Matrix Factorization (PRMF)
[![Build Status](https://travis-ci.com/gitter-lab/prmf.svg?branch=master)](https://travis-ci.com/gitter-lab/prmf)

![Graphical abstract of Pathway-Regularized Matrix Factorization](doc/abstract.png)

[Slides from ISMB2018](https://figshare.com/articles/Pathway-Regularized_Matrix_Factorization_Slides/6845648)

## Installation
```
git clone https://github.com/gitter-lab/prmf.git
cd prmf
pip install -r requirements.txt
python setup.py install
```

## Pathway Data Installation
A script to download pathways from KEGG and transform them to the expected graph format is ```script/download_pathways.sh```.
This script runs ```R``` and has separate (installation) dependencies from the core PRMF package.
Its usage is below.
```
# install the installation dependencies
Rscript script/install_pathway_dependencies.R # TODO

# from the prmf root directory
mkdir pathways
script/download_pathways.sh pathways
```

## recount2 Example
Follow these steps if you would like to run PRMF to learn patterns from recount2 and apply those patterns to a rare disease dataset such as Neurofibromatosis.
Skip to the next section if you would like to run PRMF on your own data.

0. Install R dependencies and download a recent harmonized dataset from the NF data portal at https://www.synapse.org/#!Synapse:syn18137070 and put the file 2019-02-01-bioBank_glioma_cNF_pnftidiedData.csv.gz in the PRMF directory. Decompress the file.
```
script/prep_instance.R
# visit the above URL in your browser or use the synapse client below
synapse get syn18137070
gunzip 2019-02-01-bioBank_glioma_cNF_pnftidiedData.csv.gz
```

1. Download the recount2 data archive, unzip the archive, and make the files available in ~/data/recount2_PLIER_data by running prep_recount.sh
```
script/prep_recount.sh
```

2. Transform the data with some light pre-processing for PRMF
```
script/recount_rds_to_tsv.R
```

3. Run PRMF on the transformed data to discover the latent patterns hidden in the data.
```
export OUTDIR=~/data/prmf-results ; script/prmf.py -k 10 --data ~/data/recount2_PLIER_data/recount_rpkm.tsv --manifolds pathways/kegg_ensg/*.graphml --outdir $OUTDIR > $OUTDIR/prmf.out 2> $OUTDIR/prmf.err ; echo "Exit Status: $?" >> $OUTDIR/prmf.out
```

4. PRMF decomposes the sample x gene data into a sample x latent matrix and a gene x latent matrix. These matrices and the decomposition is denoted ```X = U V^T``` and PRMF produces files U.csv and V.csv in the output directory. We next apply the gene x latent matrix in a transfer learning context: we use the latent gene patterns to quantify latent variable activity in NF.
```
transfer_learning.R -y 2019-02-01-bioBank_glioma_cNF_pnftidiedData.csv -z $OUTDIR/V.csv --outdir $OUTDIR
```

5. TODO kolmogorov_smirnov.R

## Running
With the software and necessary data installed, you can now run PRMF.
```
export OUTDIR=prmf_out
prmf_runner.py -k 10 --data ~/data/recount2_PLIER_data/recount_rpkm.tsv --normalize --delimiter '\t' --manifolds pathways/kegg_ensg/*.graphml --outdir $OUTDIR > $OUTDIR/prmf.out 2> $OUTDIR/prmf.err
echo "Exit Status: $?" >> $OUTDIR/prmf.out
```

## Graph Format
The pathways used in PRMF are assumed to be undirected graphs in the .graphml file format where the node attribute.
A sample is given below.
```
<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
         xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
         xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
<!-- Created by igraph -->
  <key id="v_name" for="node" attr.name="name" attr.type="string"/>
  <key id="e_weight" for="edge" attr.name="weight" attr.type="double"/>
  <graph id="G" edgedefault="directed">

    <!-- Node data sample -->
    <node id="n0">
      <data key="v_name">ENSG00000141837</data>
    </node>
    <node id="n1">
      <data key="v_name">ENSG00000148408</data>
    </node>
    ...

    <!-- Edge data sample -->
    <edge source="n2" target="n13">
      <data key="e_weight">1</data>
    </edge>
    <edge source="n2" target="n12">
      <data key="e_weight">1</data>
    </edge>
    ...
```

## Matrix Data Format
PRMF expects CSV-like data which describes gene expression measurements for a set of samples.
A generated dataset which illustrates this is available in the test directory at test/script/PRMF/test_cv/data_header.tsv.
Samples can also be named.
The gene names must correspond to node names as described in the Graph Format section.
```
ENSP0	ENSP1	ENSP2	ENSP3	...
sample1	198.45537672404348	130.95172604578474	94.30322859309669	120.43851578854357	...
sample2	170.9483727406653	103.88265492990088	88.86178431789227	116.51811746441824	...
...
```

## Tests
```
# from the prmf root directory
python test/run_tests.py
```

## Dependencies
- Python (3.7.0)
  - numpy
  - scipy
  - matplotlib
  - sklearn
  - networkx

## Installation Dependencies
- R (3.5.1)
  - KEGGREST
  - KEGGgraph
  - graph
  - igraph
  - argparse
