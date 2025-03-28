# Genetic Algorithm Optimization Framework

A sophisticated multi-level genetic algorithm implementation for optimization problem (here, it is choosing best routes in a connection graph), featuring multiple algorithmic layers and parallel processing capabilities.

## Project Overview

This project implements a hierarchical genetic algorithm framework with multiple levels of optimization strategies:

1. **Base Level**: Standard Genetic Algorithm (`GeneticAlgorithm`)
2. **Island Level**: Parallel Island Model (`PZ_algorithm`)
3. **Meta Level**: Meta-optimization with restart capabilities (`PZ_Meta_Algorithm`)
4. **Parameter-Less Population Pyramid**: implementaiton of sophisticated P3 algorithm (https://www.researchgate.net/publication/266656431_Parameter-less_population_pyramid)
5. **Optimized Problem**: choosing best routes in a connection graph, implemented by Phd. Michał Przewoźniczek (https://www.researchgate.net/profile/Michal-Przewozniczek)

## Key Features

This implementation has several layers:
- generic GA - base layer
- Island model - generic GA exploit some local space, and islands provides more exploration capabilities
- meta algorithm - restarts Island model when stuck

It works more or less like that:
!(images/GA_scheme.png)

The problem works like that
!(images/problem1.png)
We have a graph and few possible routes from one node to the other

!(images/problem2)
There are many connections. Our task is to choose best combination of routes indices:

!(images/problem3.png)
there are e.g. about 2000 nodes that we have to connect, for each one we choose an index of a predefined route. It is then evaluated using some fancy methods based on failure resistance and similar metrics.

The problem itself is was implemented by Phd. Michał Przewoźniczek.

## Project Structure

```
AGChallenge/
├── Core Algorithm Components
│   ├── GeneticAlgorithm.*    # Base genetic algorithm implementation
│   ├── PZ_algorithm.*        # Island model implementation
│   ├── PZ_Meta_Algorithm.*   # Meta-optimization layer
│   └── AbstractEvolutionaryAlgorithm.*  # Base abstract class
├── Supporting Components
│   ├── Individual.*         # Individual representation
│   ├── LinkageTree.*        # Linkage tree implementation
│   ├── ThreadPool.*         # Thread management
│   └── RandomValuesHolder.* # Random number generation
├── Evaluation
│   ├── Evaluator.*          # Fitness evaluation
│   └── NETsimulator.*       # Network simulation
└── Utilities
    ├── Timer.*              # Time tracking
    └── Tools.*              # Helper functions
```

## Key Parameters

Available in PZ_Meta_Algorithm.cpp:
```C++
// Size of the population in each genetic algorithm instance
int param_population_size = 100;

// Probability of performing crossover between two selected parents
float param_cross_prob = 1.0;

// Probability of mutating each gene in an individual
float param_mut_prob = 0.001;

// Probability of choosing the better parent's gene during crossover
// When 0.0, genes are chosen randomly from either parent
// When 0.5, there's a 50% chance to choose the better parent's gene
float param_choose_better_in_cross_prob = 0.0; // previously it was 0.5 both values are ok

// Size of the population segment used for building the linkage tree
// Larger values provide more statistical significance but require more computation
int param_linkage_tree_separation_size = 100;

// Minimum number of genes that can form a cluster in the linkage tree
// Prevents creation of too small, potentially meaningless clusters
int param_linkage_tree_min_cluster = 3;

// Maximum number of genes that can form a cluster in the linkage tree
// Prevents creation of too large clusters that might reduce algorithm flexibility
int param_linkage_tree_max_cluster = 99;

// Whether to use a generic tree structure for linkage learning
// When true, uses a more flexible but potentially slower tree structure
bool param_use_generic_tree = true;

// Whether to output detailed logging information during execution
// Useful for debugging and understanding algorithm behavior
bool param_verbose = false;
```

## Building and Running

The project is built using Visual Studio 2017/2022. To build and run:

1. Open `AGChallenge.sln` in Visual Studio
2. Build the solution
3. Run the executable

## Performance Monitoring

The implementation includes various monitoring features:
- Population entropy tracking
- Fitness quantile analysis
- Best solution tracking across all levels

## Notes

- The Parameter-less Population Pyramid implementation is included but not used in the final version
- The implementation uses a sophisticated linkage tree approach for genetic operations
- The meta-algorithm can restart entire instances when stuck in local optima
- Multi-threading is implemented at the island level for parallel optimization

## Dependencies

- C++17 or later
- Visual Studio 2017/2022
- Standard C++ libraries 