#ifndef TSP_H
#define TSP_H

#include <float.h>
#include <string.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <cplex.h>

#define MAX_TIME DBL_MAX
#define DBL_INFY DBL_MAX
#define INT_NINFY ((INT_MIN + 1) / 2)

#define SEQUENTIAL 1
#define GREEDY 2
#define EXTRA_MILEAGE 3

#define NO_REF 1
#define TWO_OPT 2

#define NO_MH 1
#define TABU 2
#define VNS 3
#define SIM_ANNEAL 4
#define GEN 5

#define BENDERS 1
#define CALLBACK 2
#define HARD_FIXING 3
#define LOCAL_BRANCHING 4

#define MIN_RAND_RUNS 1000

#define DIST_INDEX(i, j, n) ((j) - 1 + (i) * ((n) - 1) - ((i) * ((i) + 1)) / 2)

typedef struct point_2d {
	double x_coord;
	double y_coord;
} point_2d_t;

typedef struct tsp_instance {
	//instance input data
	int num_nodes;
	point_2d_t* nodes;

	//instance hyper-parameters
	double time_limit;
	char* input_file_name;
	int x_bound;
	int y_bound;
	unsigned int random_seed;
	unsigned int sol_procedure_flag;
	int starting_index;
	double prob_ign_opt;
	unsigned int refine_flag;
	unsigned int metaheur_flag;
	int min_tenure;
	int max_tenure;
	int min_temperature;
	int max_temperature;
	unsigned int move_weight;
	unsigned int pop_size;
	unsigned int cplex_solver_flag;
	double instance_time_limit;
	double prob_fixing;

	//instance local data
	unsigned int* best_sol;	//indices of the nodes[] in the tour order, it might be better to put this variable into a critical region to guarantee atomicity
	double best_sol_cost;
	double* costs;
	double time_left;
	unsigned int* tabu_list;
	unsigned int num_cycles;
	unsigned int cycle_delimiter;
} tsp_instance_t;

/**
 * @brief read the command line and parse it to determine the tsp instance's parameters (e.g. num_nodes, time_limit, input_file_name, etc.)
 * @param argc number of command line arguments
 * @param argv pointer to the array of command line arguments
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int parse_command_line(const int argc, const char* argv[], tsp_instance_t* instance);

/**
 * @brief read the input file and store the node list if the node file is provided by the user, otherwise generate a random set of points as instance
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int get_data(tsp_instance_t* instance);

/**
 * @brief find the optimal solution to the tsp instance, based on the sol_procedure_flag, the best solution found is stored in the best_sol field of the tsp instance
 * @param instance the instance whose solution must be found
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int tsp_opt(tsp_instance_t* instance);

/**
 * @brief deallocates all the heap memory of the tsp instance
 * @param instance the tsp instance whose memory needs to be cleaned
 */
void free_tsp_instance(tsp_instance_t* instance);

/**
 * @brief plots the points associated with the tsp instance
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int plot_points(tsp_instance_t* instance);

/**
 * @brief computes the Euclidean distance (i.e. cost) between the 2 nodes of index i and j of the tsp instance
 * @param instance the tsp instance
 * @param i index of the first node
 * @param j index of the second node
 * @return Euclidean distance between instance->nodes[i] and instance->nodes[j]
 */
double dist(int i, int j, tsp_instance_t* instance);

/**
 * @brief computes the distance (i.e. cost) between any two nodes of the tsp instance
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int precompute_costs(tsp_instance_t* instance);

/**
 * @brief returns the precomputed distance (i.e. cost) between the nodes at the two specified indices
 * @param i index of the first node
 * @param j index of the second node
 * @param instance the tsp instance
 * @return the distance between the nodes of index i and j of the tsp instance
 */
double lookup_cost(int i, int j, tsp_instance_t* instance);

/**
 * @brief rounds the entries in the cost look-up table to the nearest int
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int int_round_clut(tsp_instance_t* instance);

/**
 * @brief generate a tour by visiting the nodes sequentially (0->1->2->3->...->n-1)
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int tsp_seq_sol(tsp_instance_t* instance);

/**
 * @brief generate a tour by applying the greedy heuristic, the tour is maintained only if it's better than the best known tour
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int tsp_gdy_sol(tsp_instance_t* instance);

/**
 * @brief generate a tour by applying the extra-mileage heuristic, the tour is maintained only if it's better than the best known tour
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int tsp_exm_sol(tsp_instance_t* instance);

/**
 * @brief plots the instance tour
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int plot_tour(tsp_instance_t* instance);

/**
 * @brief given an instance of tsp with an already present solution, refines said solution using the appropriate method defined by the refine_flag
 * @param instance the instance whose solution will be refined
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int ref_sol(tsp_instance_t* instance);

/**
 * @brief refine the tsp solution by using the 2-opt move and the appropriate meta-heuristic specified by metaheur_flag
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int two_opt_ref(tsp_instance_t* instance);

/**
 * @brief produce the tsp solution by using the genetic algorithm meta-heuristic
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int genetic(tsp_instance_t* instance);

/**
 * @brief builds the basinc MILP TSP model (cost, degree constraints and variable type) based on the nodes found in instance->nodes[]
 * @param instance the tsp instance based on whose nodes the model will be created
 * @param env CPLEX environment
 * @param lp CPLEX linear program
 */
void build_model(tsp_instance_t* instance, CPXENVptr env, CPXLPptr lp);

/**
 * @brief plots the cycles found in instance->best_sol[]
 * @param instance the tsp instance whose cycles will be ploted
 * @return 0 if no error was detected, a non-0 value otherwise 
 */
int plot_cycles(tsp_instance_t* instance);

/**
 * @brief iterate over the found solutions until a feasible solution is found (only one connected component)
 * @param instance the tsp instance, must already contain a solution
 * @param env CPLEX environment
 * @param lp CPLEX linear program
 * @param xstar pointer to the solution found using CPLEX
 * @param succ successor list of the solution
 * @param comp component list of the solution
 * @param visited_nodes list of visited nodes of the solution
 * @param num_cycles number of connected components
 */
void benders_method(tsp_instance_t* instance, CPXENVptr env, CPXLPptr lp, double** xstar, int* succ, int* comp, bool* visited_nodes, int* num_cycles);

/**
 * @brief build and optimize the tsp model based on the CPLEX library
 * @param instance the tsp instance whose model will be built and solved
 * @return 0 if no error was detected, a non-0 value otherwise 
 */
int TSPopt(tsp_instance_t* instance);

#endif //TSP_H