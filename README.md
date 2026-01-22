
Given a population where a small subset of individuals is infected, testing each individual separately is inefficient. Group testing improves efficiency by testing subsets of individuals at once, but classical binary splitting can still incur unnecessary tests.

This project introduces a **subgroup splitting strategy** that dynamically divides the population into uneven or structured subgroups to reduce redundancy and improve efficiency.

The notebook `main.ipynb` evaluates the performance of this approach and compares it directly against **regular binary splitting**.

## Algorithm Description

At a high level, the algorithm works as follows:

1. Begin with a population set and a pooled infection test.
2. If a group tests negative, all individuals in the group are declared uninfected.
3. If a group tests positive, the algorithm splits the group into subgroups using a predefined splitting strategy.
4. The process recurses on positive subgroups until all infected individuals are identified.

Unlike standard binary splitting (which always divides groups evenly), this method uses **subgroup-aware splits** that reduce expected testing cost under sparse infection assumptions.

The code implements four algorithms:

Q1: Quantitative testing with exact infection counts (sequential)

Q2: Quantitative testing with infection range information

Q1-comm: Community-aware version of Q1

Q2-comm: Community-aware version of Q2
