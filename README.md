# DynamicLouvain

Implementation of a Dynamic Community Detection Algorithm based on the Static Louvain Community Detection Algorithm:

Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008). Fast unfolding of communities in large networks. Journal of Statistical Mechanics: Theory and Experiment, 2008(10). http://doi.org/10.1088/1742-5468/2008/10/P10008

This source code supports the paper published in the following paper:

Cordeiro, M., Sarmento, R. P., & Gama, J. (2016). Dynamic community detection in evolving networks using locality modularity optimization. Social Network Analysis and Mining, 6(1). http://doi.org/10.1007/s13278-016-0325-1

### Build:

```
mkdir build
cd build
cmake ..\src
cmake -G "Visual Studio 11 2012" ..\src
cmake --build . --config Release

```
