# GasSteadyHierarchicalNetworkPartitioningSim.jl 

This package implements an algorithm that first  sets up an hierarchical partition of the underlying gas network into subnetworks and then uses the hierarchical structure to solve the steady state gas flow problem for each pipeline subnetwork to obtain the solution for the full network. The algorithm allows use of both  ideal and non-ideal equations of state.

## Citation
The details of the algorithm are described in the following paper:
([doi link](https://dx.doi.org/10.1109/LCSYS.2025.3533383)): 

```bibtex

@Article{hierarchical-network-partitioning,
  author           = {Shriram Srinivasan and Kaarthik Sundar},
  date             = {2025},
  journal     = {IEEE Control System Letters},
  title            = {Hierarchical Network Partitioning for Solution of Potential-Driven, Steady-State, Nonlinear Network Flow Equations},
  doi              = {10.1109/LCSYS.2025.3533383},
  volume = {8},
  pages = {3368--3373},
}
```


``GasSteadyHierarchicalNetworkPartitioningSim.jl`` is not a registered Julia package. Hence installation of the package should be done as follows:

```julia 
using Pkg
Pkg.add("https://github.com/kaarthiksundar/GasSteadyHierarchicalNetworkPartitioningSim.jl.git")
```

For the API usage, users are referred to the ``test/`` and the ``examples/`` directories.