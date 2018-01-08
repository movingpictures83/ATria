# ATria
# Language: C++
# Input: CSV (network)
# Output: NOA (list of central nodes and centrality values)

PluMA plugin that runs the Ablatio Triadum (ATria) centrality algorithm (Cickovski et al, 2015, 2017).
ATria can run on signed and weighted networks and produces a list of central nodes as both screen
output and as a NOde Attribute (NOA) file for Cytoscape.  The NOA file can subsequently be imported
into Cytoscape resulting in centrality values becoming node attributes and enabling further analysis
and visualization based on these values.

The input network should be specified in CSV format with nodes as rows and columns and entry (i, j)
representing the weight of the edge from node i to node j.

The output is the NOA file, with both centrality value and rank as attributes.
Larger magnitude values indicate higher centrality for both centrality and rank.  This is typically
more convenient for visualization, etc.  Please see attached example.

Note there is also a faster plugin GPUATria implemented in CUDA for NVIDIA graphics cards.

