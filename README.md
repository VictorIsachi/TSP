# Algorithms for the TSP developed in the Operations Research 2 course

## Program use and command line arguments
To run the program you need to specify either:
- the file that contains the nodes of the symmetic TSP instance
- some parameters that are used to generate a random dataset.

In the first case you need to specify the command line argument `-file_name <instance.tsp>` where `<instance.tsp>` is the name of a file following the TSPLIB format. `<instannce.tsp>` must have the `TYPE` and `EDGE_WEIGHT_TYPE` fields specified as follows: `TYPE : TSP`, `EDGE_WEIGHT_TYPE : EUC_2D`.

In the second case you need to specify `-x_bound <a> -y_bound <b> -num_nodes <c>`, where `<a>`, `<b>` and `<c>` are **integers > 1**. The TSP instance will be generated as a set of `<c>` nodes in the `[0, <a>]x[0, <b>]` interval. 

You can also specify the following command line arguments:
- `-time_limit` **double**: specifies the maximum number of seconds the program is allowed to run; by default set to DBL_INF.
- `-random_seed` **int**:	specifies the random seed; by default set with time(NULL).
- `-proc_flag` **1 | 2 | 3**: specifies the procedure used to obtain the solution. 1 = sequential, 2 = greedy, 3 = extra-mileage. By default set to 1.	
- `-start_index` **int**: specifies the starting index of the procedure. Not used by the sequential method. Must be **< #nodes** for the greedy method; if set to -1 will pick a starting index at random, if set to -2 will try all the starting indices and keep the best solution; by default set to 0. For the extra-mileage method any integer > -1 will result in the 2 nodes at maximum distance, -1 will result in 2 nodes picked at random, -2 will try all the pairs of nodes and keep the best solution; by default set to 0.  
- `-prob_ign_opt` **double**: specifies the probability of selecting the 2nd best move at a particular iteration (GRASP). Must be between 0.0 and 1.0; by default set to 0.0.
- `-ref_flag` **1 | 2**: specifies the refinement algorithm. 1 = no refinement, 2 = 2-opt move; by default set to 1.
- `-meta_heur_flag` **1 | 2 | 3 | 4 | 5**: specifies the meta-heuristic method used to improve the 2-opt refinement procedure. 1 = no meta-heuristic, 2 = tabu search, 3 = VNS (5-kick), 4 = simulated annealing, 5 = genetic algorithm. By default set to 1.
- `-min_tenure` **int**: specifies the minimum tenure, must be **>= 0**. If not specified will be picked automatically based on the tsp instance size. If specified must also specify `-max_tenure`.
- `-max_tenure` **int**: specifies the maximum tenure, must be **>= min_tenure**. If not specified will be picked automatically based on the tsp instance size. If specified must also specify `-min_tenure`.
- `-min_temp` **int**: specifies the minimum temperature, must be **>= 0**. If not specified will be initialized to 10. If specified must also specify `-max_temperature`.
- `-max_temp` **int**: specifies the maximum temperature, must be **>= min_temperature**. If not specified will be initialized to 100. If specified must also specify `-min_temperature`.
- `-move_weight` **unsigned int**: specifies the weight used in the computation of the probability of the simulated annealing move. By default set to 25.
- `-pop_size` **unsigned int**: specifies the population size used by the genetic algorithm. By default set to 100.

## List of implemented features
```
-command line parser (partial)
-input file TSPLIB format reader (partial)
-points GNUPLOT plotter
-costs (EUC_2D) precomputation (look-up table) function
-sequential tour method function
-tour GNUPLOT plotter
-deterministic greedy heuristic
-deterministic extra-mileage heuristic
-greedy heuristic with GRASP
-extra-mileage heuristic with GRASP
-the 2-opt refinement procedure
-tabu search meta-heuristic
-VNS (5-kick) meta-heuristic
-simulated annealing
-genetic algorithm
```
