#ifndef TSP_H
#define TSP_H

#include <float.h>

#define MAX_TIME DBL_MAX
#define DBL_INFY DBL_MAX

#define SEQUENTIAL 1
#define GREEDY 2
#define EXTRA_MILEAGE 3

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
	unsigned int x_bound;
	unsigned int y_bound;
	unsigned int random_seed;
	unsigned int sol_procedure_flag;
	int starting_index;
	double prob_ign_opt;

	//instance local data
	unsigned int* best_sol;	//indices of the nodes[] in the tour order, it might be better to put this variable into a critical region to guarantee atomicity
	double best_sol_cost;
	double* costs;
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
 * @brief find the optimal solution to the tsp instance, the best solution found is stored in the best_sol field of the tsp instance
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
 * @brief generate a tour by applying the greedy heuristic
 * @param instance the tsp instance
 * @return 0 if no error was detected, a non-0 value otherwise
 */
int tsp_gdy_sol(tsp_instance_t* instance);

/**
 * @brief generate a tour by applying the extra-mileage heuristic
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

#endif //TSP_H