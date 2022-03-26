#include <stdio.h>
#include <stdlib.h>
#include "tsp_utils.h"
#include "tsp.h"

#define TOTAL_TIME_LIMIT 300.0
#define NUM_RERUNS 100
#define NUM_POINTS 1000
#define POINT_BOUND 1000
#define RAND_SEED 1

#define PERF_PROF_FILENAME "perf_prof.csv"

static void set_up_tsp_instance(int argc, char const* argv[], tsp_instance_t* instance) {
	//set-up the tsp instance
	if (parse_command_line(argc, argv, instance)) { fprintf(stderr, "Parsing failed\n"); exit(1); }
	if (get_data(instance)) { free_tsp_instance(instance); fprintf(stderr, "Input file reading failed\n"); exit(1); }
	if (plot_points(instance)) { free_tsp_instance(instance); fprintf(stderr, "Point plotting failed\n"); exit(1); }
	if (precompute_costs(instance)) { free_tsp_instance(instance); fprintf(stderr, "Cost precomputation failed\n"); exit(1); }
	int_round_clut(instance);
}

static void first_try_sol(tsp_instance_t* instance) {
	//first try solution
	if (tsp_opt(instance)) { free_tsp_instance(instance); fprintf(stderr, "Optimization algorithm failed\n"); exit(1); }
	printf("Solution cost: %f\n", instance->best_sol_cost);
	if (plot_tour(instance)) { free_tsp_instance(instance); fprintf(stderr, "Tour plotting failed\n"); exit(1); }
	if (ref_sol(instance)) { free_tsp_instance(instance); fprintf(stderr, "Refinement algorithm failed\n"); exit(1); }
	printf("Solution cost: %f\n", instance->best_sol_cost);
	if (plot_tour(instance)) { free_tsp_instance(instance); fprintf(stderr, "Tour plotting failed\n"); exit(1); }
}

static void cont_search(tsp_instance_t* instance) {
	//continuous search
	while (instance->time_left > 0) {
		double elapsed_timer_start = seconds();
		if (tsp_opt(instance)) { free_tsp_instance(instance); fprintf(stderr, "Optimization algorithm failed\n"); exit(1); }
		if (ref_sol(instance)) { free_tsp_instance(instance); fprintf(stderr, "Refinement algorithm failed\n"); exit(1); }
		double elapsed_timer_stop = seconds();
		instance->time_left -= (elapsed_timer_stop - elapsed_timer_start);
	}
	printf("Time limit reached...\n");
	printf("Solution cost: %f\n", instance->best_sol_cost);
	if (plot_tour(instance)) { free_tsp_instance(instance); fprintf(stderr, "Tour plotting failed\n"); exit(1); }
}

//chap. 2.3 
static void comp_methods_1() {

#if VERBOSE > 1
	{ printf("Comparing greedy_random_start, greedy_random_start + GRASP0.1, greedy_random_start + GRASP0.2\n \
			  extra-mileage_random_start, extra-mileage_random_start + GRASP0.1, extra-mileage_random_start + GRASP0.2...\n"); }
#endif

	tsp_instance_t greedy_inst_1;	//greedy random start
	greedy_inst_1.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	greedy_inst_1.input_file_name = NULL;
	greedy_inst_1.x_bound = -1;
	greedy_inst_1.y_bound = -1;
	greedy_inst_1.num_nodes = -1;
	greedy_inst_1.random_seed = time(NULL);
	greedy_inst_1.sol_procedure_flag = GREEDY;
	greedy_inst_1.starting_index = -1;
	greedy_inst_1.prob_ign_opt = 0.0;
	greedy_inst_1.refine_flag = NO_REF;
	greedy_inst_1.best_sol = NULL;
	greedy_inst_1.best_sol_cost = DBL_INFY;
	greedy_inst_1.costs = NULL;
	greedy_inst_1.time_left = greedy_inst_1.time_limit;

	tsp_instance_t greedy_inst_2;	//greedy random start + GRASP0.1
	greedy_inst_2.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	greedy_inst_2.input_file_name = NULL;
	greedy_inst_2.x_bound = -1;
	greedy_inst_2.y_bound = -1;
	greedy_inst_2.num_nodes = -1;
	greedy_inst_2.random_seed = time(NULL);
	greedy_inst_2.sol_procedure_flag = GREEDY;
	greedy_inst_2.starting_index = -1;
	greedy_inst_2.prob_ign_opt = 0.1;
	greedy_inst_2.refine_flag = NO_REF;
	greedy_inst_2.best_sol = NULL;
	greedy_inst_2.best_sol_cost = DBL_INFY;
	greedy_inst_2.costs = NULL;
	greedy_inst_2.time_left = greedy_inst_2.time_limit;

	tsp_instance_t greedy_inst_3;	//greedy random start + GRASP0.2
	greedy_inst_3.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	greedy_inst_3.input_file_name = NULL;
	greedy_inst_3.x_bound = -1;
	greedy_inst_3.y_bound = -1;
	greedy_inst_3.num_nodes = -1;
	greedy_inst_3.random_seed = time(NULL);
	greedy_inst_3.sol_procedure_flag = GREEDY;
	greedy_inst_3.starting_index = -1;
	greedy_inst_3.prob_ign_opt = 0.2;
	greedy_inst_3.refine_flag = NO_REF;
	greedy_inst_3.best_sol = NULL;
	greedy_inst_3.best_sol_cost = DBL_INFY;
	greedy_inst_3.costs = NULL;
	greedy_inst_3.time_left = greedy_inst_3.time_limit;

	tsp_instance_t extra_mileage_inst_1;	//extra-mileage random start
	extra_mileage_inst_1.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	extra_mileage_inst_1.input_file_name = NULL;
	extra_mileage_inst_1.x_bound = -1;
	extra_mileage_inst_1.y_bound = -1;
	extra_mileage_inst_1.num_nodes = -1;
	extra_mileage_inst_1.random_seed = time(NULL);
	extra_mileage_inst_1.sol_procedure_flag = EXTRA_MILEAGE;
	extra_mileage_inst_1.starting_index = -1;
	extra_mileage_inst_1.prob_ign_opt = 0.0;
	extra_mileage_inst_1.refine_flag = NO_REF;
	extra_mileage_inst_1.best_sol = NULL;
	extra_mileage_inst_1.best_sol_cost = DBL_INFY;
	extra_mileage_inst_1.costs = NULL;
	extra_mileage_inst_1.time_left = extra_mileage_inst_1.time_limit;

	tsp_instance_t extra_mileage_inst_2;	//extra-mileage random start + GRASP0.1
	extra_mileage_inst_2.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	extra_mileage_inst_2.input_file_name = NULL;
	extra_mileage_inst_2.x_bound = -1;
	extra_mileage_inst_2.y_bound = -1;
	extra_mileage_inst_2.num_nodes = -1;
	extra_mileage_inst_2.random_seed = time(NULL);
	extra_mileage_inst_2.sol_procedure_flag = EXTRA_MILEAGE;
	extra_mileage_inst_2.starting_index = -1;
	extra_mileage_inst_2.prob_ign_opt = 0.1;
	extra_mileage_inst_2.refine_flag = NO_REF;
	extra_mileage_inst_2.best_sol = NULL;
	extra_mileage_inst_2.best_sol_cost = DBL_INFY;
	extra_mileage_inst_2.costs = NULL;
	extra_mileage_inst_2.time_left = extra_mileage_inst_2.time_limit;

	tsp_instance_t extra_mileage_inst_3;	//extra-mileage random start + GRASP0.2
	extra_mileage_inst_3.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	extra_mileage_inst_3.input_file_name = NULL;
	extra_mileage_inst_3.x_bound = -1;
	extra_mileage_inst_3.y_bound = -1;
	extra_mileage_inst_3.num_nodes = -1;
	extra_mileage_inst_3.random_seed = time(NULL);
	extra_mileage_inst_3.sol_procedure_flag = EXTRA_MILEAGE;
	extra_mileage_inst_3.starting_index = -1;
	extra_mileage_inst_3.prob_ign_opt = 0.2;
	extra_mileage_inst_3.refine_flag = NO_REF;
	extra_mileage_inst_3.best_sol = NULL;
	extra_mileage_inst_3.best_sol_cost = DBL_INFY;
	extra_mileage_inst_3.costs = NULL;
	extra_mileage_inst_3.time_left = extra_mileage_inst_3.time_limit;

#if VERBOSE > 1
	{ printf("Total time limit: %f...\nTime limit per instance %f...\n", TOTAL_TIME_LIMIT, TOTAL_TIME_LIMIT / NUM_RERUNS); }
#endif

	FILE* perf_prof_file = fopen(PERF_PROF_FILENAME, "w");
	if (perf_prof_file == NULL) { fprintf(stderr, "cannot open the perf_prof_file\n"); exit(1); }

	fprintf(perf_prof_file, "6, greedy_rand_start, greedy_rand_start+GRASP0.1, greedy_rand_start+GRASP0.2, extra-mileage_rand_start, extra-mileage_rand_start+GRASP0.1, extra-mileage_rand_start+GRASP0.2\n");

	srand(RAND_SEED); for (int i = 0; i < MIN_RAND_RUNS + log(1 + RAND_SEED); i++) rand();
	point_2d_t* nodes = (point_2d_t*)malloc(NUM_POINTS * sizeof(point_2d_t));
	int num_costs = NUM_POINTS * (NUM_POINTS - 1) / 2;
	double* costs = (double*)malloc(num_costs * sizeof(double));
	unsigned int random_seed;
	for (int i = 0; i < NUM_RERUNS; i++) {

		//generating the random nodes
		for (int i = 0; i < NUM_POINTS; i++) {
			nodes[i].x_coord = ((double)rand() / RAND_MAX) * POINT_BOUND;
			nodes[i].y_coord = ((double)rand() / RAND_MAX) * POINT_BOUND;
		}
		greedy_inst_1.num_nodes = NUM_POINTS;
		greedy_inst_2.num_nodes = NUM_POINTS;
		greedy_inst_3.num_nodes = NUM_POINTS;
		greedy_inst_1.nodes = nodes;
		greedy_inst_2.nodes = nodes;
		greedy_inst_3.nodes = nodes;
		extra_mileage_inst_1.num_nodes = NUM_POINTS;
		extra_mileage_inst_2.num_nodes = NUM_POINTS;
		extra_mileage_inst_3.num_nodes = NUM_POINTS;
		extra_mileage_inst_1.nodes = nodes;
		extra_mileage_inst_2.nodes = nodes;
		extra_mileage_inst_3.nodes = nodes;

		//precomputing the costs
		for (int i = 0; i < NUM_POINTS - 1; i++) {
			for (int j = i + 1; j < NUM_POINTS; j++) {
				costs[DIST_INDEX(i, j, NUM_POINTS)] = (double)((int)(dist(i, j, &greedy_inst_1) + 0.5));
			}
		}
		greedy_inst_1.costs = costs;
		greedy_inst_2.costs = costs;
		greedy_inst_3.costs = costs;
		extra_mileage_inst_1.costs = costs;
		extra_mileage_inst_2.costs = costs;
		extra_mileage_inst_3.costs = costs;

		//setting-up the random seeds
		random_seed = rand();
		greedy_inst_1.random_seed = random_seed;
		greedy_inst_2.random_seed = random_seed;
		greedy_inst_3.random_seed = random_seed;
		extra_mileage_inst_1.random_seed = random_seed;
		extra_mileage_inst_2.random_seed = random_seed;
		extra_mileage_inst_3.random_seed = random_seed;

		//finding the best solutions
		double greedy_best_sol_cost_1 = DBL_INFY;
		greedy_inst_1.time_left = greedy_inst_1.time_limit;
		while (greedy_inst_1.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&greedy_inst_1);
			ref_sol(&greedy_inst_1);
			if (greedy_inst_1.best_sol_cost < greedy_best_sol_cost_1)
				greedy_best_sol_cost_1 = greedy_inst_1.best_sol_cost;
			free(greedy_inst_1.best_sol);
			greedy_inst_1.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			greedy_inst_1.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double greedy_best_sol_cost_2 = DBL_INFY;
		greedy_inst_2.time_left = greedy_inst_2.time_limit;
		while (greedy_inst_2.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&greedy_inst_2);
			ref_sol(&greedy_inst_2);
			if (greedy_inst_2.best_sol_cost < greedy_best_sol_cost_2)
				greedy_best_sol_cost_2 = greedy_inst_2.best_sol_cost;
			free(greedy_inst_2.best_sol);
			greedy_inst_2.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			greedy_inst_2.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double greedy_best_sol_cost_3 = DBL_INFY;
		greedy_inst_3.time_left = greedy_inst_3.time_limit;
		while (greedy_inst_3.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&greedy_inst_3);
			ref_sol(&greedy_inst_3);
			if (greedy_inst_3.best_sol_cost < greedy_best_sol_cost_3)
				greedy_best_sol_cost_3 = greedy_inst_3.best_sol_cost;
			free(greedy_inst_3.best_sol);
			greedy_inst_3.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			greedy_inst_3.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double extra_mileage_best_sol_cost_1 = DBL_INFY;
		extra_mileage_inst_1.time_left = extra_mileage_inst_1.time_limit;
		while (extra_mileage_inst_1.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&extra_mileage_inst_1);
			ref_sol(&extra_mileage_inst_1);
			if (extra_mileage_inst_1.best_sol_cost < extra_mileage_best_sol_cost_1)
				extra_mileage_best_sol_cost_1 = extra_mileage_inst_1.best_sol_cost;
			free(extra_mileage_inst_1.best_sol);
			extra_mileage_inst_1.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			extra_mileage_inst_1.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double extra_mileage_best_sol_cost_2 = DBL_INFY;
		extra_mileage_inst_2.time_left = extra_mileage_inst_2.time_limit;
		while (extra_mileage_inst_2.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&extra_mileage_inst_2);
			ref_sol(&extra_mileage_inst_2);
			if (extra_mileage_inst_2.best_sol_cost < extra_mileage_best_sol_cost_2)
				extra_mileage_best_sol_cost_2 = extra_mileage_inst_2.best_sol_cost;
			free(extra_mileage_inst_2.best_sol);
			extra_mileage_inst_2.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			extra_mileage_inst_2.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double extra_mileage_best_sol_cost_3 = DBL_INFY;
		extra_mileage_inst_3.time_left = extra_mileage_inst_3.time_limit;
		while (extra_mileage_inst_3.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&extra_mileage_inst_3);
			ref_sol(&extra_mileage_inst_3);
			if (extra_mileage_inst_3.best_sol_cost < extra_mileage_best_sol_cost_3)
				extra_mileage_best_sol_cost_3 = extra_mileage_inst_3.best_sol_cost;
			free(extra_mileage_inst_3.best_sol);
			extra_mileage_inst_3.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			extra_mileage_inst_3.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}

		fprintf(perf_prof_file, "random_seed_%d, %f, %f, %f, %f, %f, %f\n", random_seed, greedy_best_sol_cost_1, greedy_best_sol_cost_2, greedy_best_sol_cost_3, \
			extra_mileage_best_sol_cost_1, extra_mileage_best_sol_cost_2, extra_mileage_best_sol_cost_3);
	}
	free(nodes);
	free(costs);

	fclose(perf_prof_file);
}

//chap. 3.2 fig. 1
static void comp_methods_2() {

#if VERBOSE > 1
	{ printf("Comparing greedy_random_start + 2-opt, greedy_random_start + GRASP0.1 + 2-opt, greedy_random_start + GRASP0.2 + 2-opt\n \
			  extra-mileage_random_start + 2-opt, extra-mileage_random_start + GRASP0.1 + 2-opt, extra-mileage_random_start + GRASP0.2 + 2-opt...\n"); }
#endif

	tsp_instance_t greedy_inst_1;	//greedy random start
	greedy_inst_1.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	greedy_inst_1.input_file_name = NULL;
	greedy_inst_1.x_bound = -1;
	greedy_inst_1.y_bound = -1;
	greedy_inst_1.num_nodes = -1;
	greedy_inst_1.random_seed = time(NULL);
	greedy_inst_1.sol_procedure_flag = GREEDY;
	greedy_inst_1.starting_index = -1;
	greedy_inst_1.prob_ign_opt = 0.0;
	greedy_inst_1.refine_flag = TWO_OPT;
	greedy_inst_1.best_sol = NULL;
	greedy_inst_1.best_sol_cost = DBL_INFY;
	greedy_inst_1.costs = NULL;
	greedy_inst_1.time_left = greedy_inst_1.time_limit;

	tsp_instance_t greedy_inst_2;	//greedy random start + GRASP0.1
	greedy_inst_2.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	greedy_inst_2.input_file_name = NULL;
	greedy_inst_2.x_bound = -1;
	greedy_inst_2.y_bound = -1;
	greedy_inst_2.num_nodes = -1;
	greedy_inst_2.random_seed = time(NULL);
	greedy_inst_2.sol_procedure_flag = GREEDY;
	greedy_inst_2.starting_index = -1;
	greedy_inst_2.prob_ign_opt = 0.1;
	greedy_inst_2.refine_flag = TWO_OPT;
	greedy_inst_2.best_sol = NULL;
	greedy_inst_2.best_sol_cost = DBL_INFY;
	greedy_inst_2.costs = NULL;
	greedy_inst_2.time_left = greedy_inst_2.time_limit;

	tsp_instance_t greedy_inst_3;	//greedy random start + GRASP0.2
	greedy_inst_3.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	greedy_inst_3.input_file_name = NULL;
	greedy_inst_3.x_bound = -1;
	greedy_inst_3.y_bound = -1;
	greedy_inst_3.num_nodes = -1;
	greedy_inst_3.random_seed = time(NULL);
	greedy_inst_3.sol_procedure_flag = GREEDY;
	greedy_inst_3.starting_index = -1;
	greedy_inst_3.prob_ign_opt = 0.2;
	greedy_inst_3.refine_flag = TWO_OPT;
	greedy_inst_3.best_sol = NULL;
	greedy_inst_3.best_sol_cost = DBL_INFY;
	greedy_inst_3.costs = NULL;
	greedy_inst_3.time_left = greedy_inst_3.time_limit;

	tsp_instance_t extra_mileage_inst_1;	//extra-mileage random start
	extra_mileage_inst_1.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	extra_mileage_inst_1.input_file_name = NULL;
	extra_mileage_inst_1.x_bound = -1;
	extra_mileage_inst_1.y_bound = -1;
	extra_mileage_inst_1.num_nodes = -1;
	extra_mileage_inst_1.random_seed = time(NULL);
	extra_mileage_inst_1.sol_procedure_flag = EXTRA_MILEAGE;
	extra_mileage_inst_1.starting_index = -1;
	extra_mileage_inst_1.prob_ign_opt = 0.0;
	extra_mileage_inst_1.refine_flag = TWO_OPT;
	extra_mileage_inst_1.best_sol = NULL;
	extra_mileage_inst_1.best_sol_cost = DBL_INFY;
	extra_mileage_inst_1.costs = NULL;
	extra_mileage_inst_1.time_left = extra_mileage_inst_1.time_limit;

	tsp_instance_t extra_mileage_inst_2;	//extra-mileage random start + GRASP0.1
	extra_mileage_inst_2.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	extra_mileage_inst_2.input_file_name = NULL;
	extra_mileage_inst_2.x_bound = -1;
	extra_mileage_inst_2.y_bound = -1;
	extra_mileage_inst_2.num_nodes = -1;
	extra_mileage_inst_2.random_seed = time(NULL);
	extra_mileage_inst_2.sol_procedure_flag = EXTRA_MILEAGE;
	extra_mileage_inst_2.starting_index = -1;
	extra_mileage_inst_2.prob_ign_opt = 0.1;
	extra_mileage_inst_2.refine_flag = TWO_OPT;
	extra_mileage_inst_2.best_sol = NULL;
	extra_mileage_inst_2.best_sol_cost = DBL_INFY;
	extra_mileage_inst_2.costs = NULL;
	extra_mileage_inst_2.time_left = extra_mileage_inst_2.time_limit;

	tsp_instance_t extra_mileage_inst_3;	//extra-mileage random start + GRASP0.2
	extra_mileage_inst_3.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	extra_mileage_inst_3.input_file_name = NULL;
	extra_mileage_inst_3.x_bound = -1;
	extra_mileage_inst_3.y_bound = -1;
	extra_mileage_inst_3.num_nodes = -1;
	extra_mileage_inst_3.random_seed = time(NULL);
	extra_mileage_inst_3.sol_procedure_flag = EXTRA_MILEAGE;
	extra_mileage_inst_3.starting_index = -1;
	extra_mileage_inst_3.prob_ign_opt = 0.2;
	extra_mileage_inst_3.refine_flag = TWO_OPT;
	extra_mileage_inst_3.best_sol = NULL;
	extra_mileage_inst_3.best_sol_cost = DBL_INFY;
	extra_mileage_inst_3.costs = NULL;
	extra_mileage_inst_3.time_left = extra_mileage_inst_3.time_limit;

#if VERBOSE > 1
	{ printf("Total time limit: %f...\nTime limit per instance %f...\n", TOTAL_TIME_LIMIT, TOTAL_TIME_LIMIT / NUM_RERUNS); }
#endif

	FILE* perf_prof_file = fopen(PERF_PROF_FILENAME, "w");
	if (perf_prof_file == NULL) { fprintf(stderr, "cannot open the perf_prof_file\n"); exit(1); }

	fprintf(perf_prof_file, "6, greedy_rand_start+2-opt, greedy_rand_start+GRASP0.1+2-opt, greedy_rand_start+GRASP0.2+2-opt, \
							extra-mileage_rand_start+2-opt, extra-mileage_rand_start+GRASP0.1+2-opt, extra-mileage_rand_start+GRASP0.2+2-opt\n");

	srand(RAND_SEED); for (int i = 0; i < MIN_RAND_RUNS + log(1 + RAND_SEED); i++) rand();
	point_2d_t* nodes = (point_2d_t*)malloc(NUM_POINTS * sizeof(point_2d_t));
	int num_costs = NUM_POINTS * (NUM_POINTS - 1) / 2;
	double* costs = (double*)malloc(num_costs * sizeof(double));
	unsigned int random_seed;
	for (int i = 0; i < NUM_RERUNS; i++) {

		//generating the random nodes
		for (int i = 0; i < NUM_POINTS; i++) {
			nodes[i].x_coord = ((double)rand() / RAND_MAX) * POINT_BOUND;
			nodes[i].y_coord = ((double)rand() / RAND_MAX) * POINT_BOUND;
		}
		greedy_inst_1.num_nodes = NUM_POINTS;
		greedy_inst_2.num_nodes = NUM_POINTS;
		greedy_inst_3.num_nodes = NUM_POINTS;
		greedy_inst_1.nodes = nodes;
		greedy_inst_2.nodes = nodes;
		greedy_inst_3.nodes = nodes;
		extra_mileage_inst_1.num_nodes = NUM_POINTS;
		extra_mileage_inst_2.num_nodes = NUM_POINTS;
		extra_mileage_inst_3.num_nodes = NUM_POINTS;
		extra_mileage_inst_1.nodes = nodes;
		extra_mileage_inst_2.nodes = nodes;
		extra_mileage_inst_3.nodes = nodes;

		//precomputing the costs
		for (int i = 0; i < NUM_POINTS - 1; i++) {
			for (int j = i + 1; j < NUM_POINTS; j++) {
				costs[DIST_INDEX(i, j, NUM_POINTS)] = (double)((int)(dist(i, j, &greedy_inst_1) + 0.5));
			}
		}
		greedy_inst_1.costs = costs;
		greedy_inst_2.costs = costs;
		greedy_inst_3.costs = costs;
		extra_mileage_inst_1.costs = costs;
		extra_mileage_inst_2.costs = costs;
		extra_mileage_inst_3.costs = costs;

		//setting-up the random seeds
		random_seed = rand();
		greedy_inst_1.random_seed = random_seed;
		greedy_inst_2.random_seed = random_seed;
		greedy_inst_3.random_seed = random_seed;
		extra_mileage_inst_1.random_seed = random_seed;
		extra_mileage_inst_2.random_seed = random_seed;
		extra_mileage_inst_3.random_seed = random_seed;

		//finding the best solutions
		double greedy_best_sol_cost_1 = DBL_INFY;
		greedy_inst_1.time_left = greedy_inst_1.time_limit;
		while (greedy_inst_1.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&greedy_inst_1);
			ref_sol(&greedy_inst_1);
			if (greedy_inst_1.best_sol_cost < greedy_best_sol_cost_1)
				greedy_best_sol_cost_1 = greedy_inst_1.best_sol_cost;
			free(greedy_inst_1.best_sol);
			greedy_inst_1.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			greedy_inst_1.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double greedy_best_sol_cost_2 = DBL_INFY;
		greedy_inst_2.time_left = greedy_inst_2.time_limit;
		while (greedy_inst_2.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&greedy_inst_2);
			ref_sol(&greedy_inst_2);
			if (greedy_inst_2.best_sol_cost < greedy_best_sol_cost_2)
				greedy_best_sol_cost_2 = greedy_inst_2.best_sol_cost;
			free(greedy_inst_2.best_sol);
			greedy_inst_2.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			greedy_inst_2.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double greedy_best_sol_cost_3 = DBL_INFY;
		greedy_inst_3.time_left = greedy_inst_3.time_limit;
		while (greedy_inst_3.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&greedy_inst_3);
			ref_sol(&greedy_inst_3);
			if (greedy_inst_3.best_sol_cost < greedy_best_sol_cost_3)
				greedy_best_sol_cost_3 = greedy_inst_3.best_sol_cost;
			free(greedy_inst_3.best_sol);
			greedy_inst_3.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			greedy_inst_3.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double extra_mileage_best_sol_cost_1 = DBL_INFY;
		extra_mileage_inst_1.time_left = extra_mileage_inst_1.time_limit;
		while (extra_mileage_inst_1.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&extra_mileage_inst_1);
			ref_sol(&extra_mileage_inst_1);
			if (extra_mileage_inst_1.best_sol_cost < extra_mileage_best_sol_cost_1)
				extra_mileage_best_sol_cost_1 = extra_mileage_inst_1.best_sol_cost;
			free(extra_mileage_inst_1.best_sol);
			extra_mileage_inst_1.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			extra_mileage_inst_1.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double extra_mileage_best_sol_cost_2 = DBL_INFY;
		extra_mileage_inst_2.time_left = extra_mileage_inst_2.time_limit;
		while (extra_mileage_inst_2.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&extra_mileage_inst_2);
			ref_sol(&extra_mileage_inst_2);
			if (extra_mileage_inst_2.best_sol_cost < extra_mileage_best_sol_cost_2)
				extra_mileage_best_sol_cost_2 = extra_mileage_inst_2.best_sol_cost;
			free(extra_mileage_inst_2.best_sol);
			extra_mileage_inst_2.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			extra_mileage_inst_2.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double extra_mileage_best_sol_cost_3 = DBL_INFY;
		extra_mileage_inst_3.time_left = extra_mileage_inst_3.time_limit;
		while (extra_mileage_inst_3.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&extra_mileage_inst_3);
			ref_sol(&extra_mileage_inst_3);
			if (extra_mileage_inst_3.best_sol_cost < extra_mileage_best_sol_cost_3)
				extra_mileage_best_sol_cost_3 = extra_mileage_inst_3.best_sol_cost;
			free(extra_mileage_inst_3.best_sol);
			extra_mileage_inst_3.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			extra_mileage_inst_3.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}

		fprintf(perf_prof_file, "random_seed_%d, %f, %f, %f, %f, %f, %f\n", random_seed, greedy_best_sol_cost_1, greedy_best_sol_cost_2, greedy_best_sol_cost_3, \
			extra_mileage_best_sol_cost_1, extra_mileage_best_sol_cost_2, extra_mileage_best_sol_cost_3);
	}
	free(nodes);
	free(costs);

	fclose(perf_prof_file);
}

//chap. 3.2 fig. 2
static void comp_methods_3() {

#if VERBOSE > 1
	{ printf("Comparing greedy_random_start, greedy_random_start + 2-opt\n \
			  extra-mileage_random_start, extra-mileage_random_start + 2-opt...\n");
}
#endif

	tsp_instance_t greedy_inst_1;	//greedy random start
	greedy_inst_1.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	greedy_inst_1.input_file_name = NULL;
	greedy_inst_1.x_bound = -1;
	greedy_inst_1.y_bound = -1;
	greedy_inst_1.num_nodes = -1;
	greedy_inst_1.random_seed = time(NULL);
	greedy_inst_1.sol_procedure_flag = GREEDY;
	greedy_inst_1.starting_index = -1;
	greedy_inst_1.prob_ign_opt = 0.0;
	greedy_inst_1.refine_flag = NO_REF;
	greedy_inst_1.best_sol = NULL;
	greedy_inst_1.best_sol_cost = DBL_INFY;
	greedy_inst_1.costs = NULL;
	greedy_inst_1.time_left = greedy_inst_1.time_limit;

	tsp_instance_t greedy_inst_2;	//greedy random start + GRASP0.1
	greedy_inst_2.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	greedy_inst_2.input_file_name = NULL;
	greedy_inst_2.x_bound = -1;
	greedy_inst_2.y_bound = -1;
	greedy_inst_2.num_nodes = -1;
	greedy_inst_2.random_seed = time(NULL);
	greedy_inst_2.sol_procedure_flag = GREEDY;
	greedy_inst_2.starting_index = -1;
	greedy_inst_2.prob_ign_opt = 0.0;
	greedy_inst_2.refine_flag = TWO_OPT;
	greedy_inst_2.best_sol = NULL;
	greedy_inst_2.best_sol_cost = DBL_INFY;
	greedy_inst_2.costs = NULL;
	greedy_inst_2.time_left = greedy_inst_2.time_limit;

	tsp_instance_t extra_mileage_inst_1;	//extra-mileage random start
	extra_mileage_inst_1.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	extra_mileage_inst_1.input_file_name = NULL;
	extra_mileage_inst_1.x_bound = -1;
	extra_mileage_inst_1.y_bound = -1;
	extra_mileage_inst_1.num_nodes = -1;
	extra_mileage_inst_1.random_seed = time(NULL);
	extra_mileage_inst_1.sol_procedure_flag = EXTRA_MILEAGE;
	extra_mileage_inst_1.starting_index = -1;
	extra_mileage_inst_1.prob_ign_opt = 0.0;
	extra_mileage_inst_1.refine_flag = NO_REF;
	extra_mileage_inst_1.best_sol = NULL;
	extra_mileage_inst_1.best_sol_cost = DBL_INFY;
	extra_mileage_inst_1.costs = NULL;
	extra_mileage_inst_1.time_left = extra_mileage_inst_1.time_limit;

	tsp_instance_t extra_mileage_inst_2;	//extra-mileage random start + GRASP0.1
	extra_mileage_inst_2.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	extra_mileage_inst_2.input_file_name = NULL;
	extra_mileage_inst_2.x_bound = -1;
	extra_mileage_inst_2.y_bound = -1;
	extra_mileage_inst_2.num_nodes = -1;
	extra_mileage_inst_2.random_seed = time(NULL);
	extra_mileage_inst_2.sol_procedure_flag = EXTRA_MILEAGE;
	extra_mileage_inst_2.starting_index = -1;
	extra_mileage_inst_2.prob_ign_opt = 0.0;
	extra_mileage_inst_2.refine_flag = TWO_OPT;
	extra_mileage_inst_2.best_sol = NULL;
	extra_mileage_inst_2.best_sol_cost = DBL_INFY;
	extra_mileage_inst_2.costs = NULL;
	extra_mileage_inst_2.time_left = extra_mileage_inst_2.time_limit;

#if VERBOSE > 1
	{ printf("Total time limit: %f...\nTime limit per instance %f...\n", TOTAL_TIME_LIMIT, TOTAL_TIME_LIMIT / NUM_RERUNS); }
#endif

	FILE* perf_prof_file = fopen(PERF_PROF_FILENAME, "w");
	if (perf_prof_file == NULL) { fprintf(stderr, "cannot open the perf_prof_file\n"); exit(1); }

	fprintf(perf_prof_file, "4, greedy_rand_start, greedy_rand_start+2-opt, extra-mileage_rand_start, extra-mileage_rand_start+2-opt\n");

	srand(RAND_SEED); for (int i = 0; i < MIN_RAND_RUNS + log(1 + RAND_SEED); i++) rand();
	point_2d_t* nodes = (point_2d_t*)malloc(NUM_POINTS * sizeof(point_2d_t));
	int num_costs = NUM_POINTS * (NUM_POINTS - 1) / 2;
	double* costs = (double*)malloc(num_costs * sizeof(double));
	unsigned int random_seed;
	for (int i = 0; i < NUM_RERUNS; i++) {

		//generating the random nodes
		for (int i = 0; i < NUM_POINTS; i++) {
			nodes[i].x_coord = ((double)rand() / RAND_MAX) * POINT_BOUND;
			nodes[i].y_coord = ((double)rand() / RAND_MAX) * POINT_BOUND;
		}
		greedy_inst_1.num_nodes = NUM_POINTS;
		greedy_inst_2.num_nodes = NUM_POINTS;
		greedy_inst_1.nodes = nodes;
		greedy_inst_2.nodes = nodes;
		extra_mileage_inst_1.num_nodes = NUM_POINTS;
		extra_mileage_inst_2.num_nodes = NUM_POINTS;
		extra_mileage_inst_1.nodes = nodes;
		extra_mileage_inst_2.nodes = nodes;

		//precomputing the costs
		for (int i = 0; i < NUM_POINTS - 1; i++) {
			for (int j = i + 1; j < NUM_POINTS; j++) {
				costs[DIST_INDEX(i, j, NUM_POINTS)] = (double)((int)(dist(i, j, &greedy_inst_1) + 0.5));
			}
		}
		greedy_inst_1.costs = costs;
		greedy_inst_2.costs = costs;
		extra_mileage_inst_1.costs = costs;
		extra_mileage_inst_2.costs = costs;

		//setting-up the random seeds
		random_seed = rand();
		greedy_inst_1.random_seed = random_seed;
		greedy_inst_2.random_seed = random_seed;
		extra_mileage_inst_1.random_seed = random_seed;
		extra_mileage_inst_2.random_seed = random_seed;

		//finding the best solutions
		double greedy_best_sol_cost_1 = DBL_INFY;
		greedy_inst_1.time_left = greedy_inst_1.time_limit;
		while (greedy_inst_1.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&greedy_inst_1);
			ref_sol(&greedy_inst_1);
			if (greedy_inst_1.best_sol_cost < greedy_best_sol_cost_1)
				greedy_best_sol_cost_1 = greedy_inst_1.best_sol_cost;
			free(greedy_inst_1.best_sol);
			greedy_inst_1.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			greedy_inst_1.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double greedy_best_sol_cost_2 = DBL_INFY;
		greedy_inst_2.time_left = greedy_inst_2.time_limit;
		while (greedy_inst_2.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&greedy_inst_2);
			ref_sol(&greedy_inst_2);
			if (greedy_inst_2.best_sol_cost < greedy_best_sol_cost_2)
				greedy_best_sol_cost_2 = greedy_inst_2.best_sol_cost;
			free(greedy_inst_2.best_sol);
			greedy_inst_2.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			greedy_inst_2.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double extra_mileage_best_sol_cost_1 = DBL_INFY;
		extra_mileage_inst_1.time_left = extra_mileage_inst_1.time_limit;
		while (extra_mileage_inst_1.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&extra_mileage_inst_1);
			ref_sol(&extra_mileage_inst_1);
			if (extra_mileage_inst_1.best_sol_cost < extra_mileage_best_sol_cost_1)
				extra_mileage_best_sol_cost_1 = extra_mileage_inst_1.best_sol_cost;
			free(extra_mileage_inst_1.best_sol);
			extra_mileage_inst_1.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			extra_mileage_inst_1.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double extra_mileage_best_sol_cost_2 = DBL_INFY;
		extra_mileage_inst_2.time_left = extra_mileage_inst_2.time_limit;
		while (extra_mileage_inst_2.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&extra_mileage_inst_2);
			ref_sol(&extra_mileage_inst_2);
			if (extra_mileage_inst_2.best_sol_cost < extra_mileage_best_sol_cost_2)
				extra_mileage_best_sol_cost_2 = extra_mileage_inst_2.best_sol_cost;
			free(extra_mileage_inst_2.best_sol);
			extra_mileage_inst_2.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			extra_mileage_inst_2.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}

		fprintf(perf_prof_file, "random_seed_%d, %f, %f, %f, %f\n", random_seed, greedy_best_sol_cost_1, greedy_best_sol_cost_2, \
			extra_mileage_best_sol_cost_1, extra_mileage_best_sol_cost_2);
	}
	free(nodes);
	free(costs);

	fclose(perf_prof_file);
}

static void comp_methods_0() {
	//compare the solution methods

#if VERBOSE > 1
	{ printf("Comparing greedy_random_start + GRASP0.1 + 2-opt and extra-mileage_random_start + GRASP0.1 + 2-opt...\n"); }
#endif

	tsp_instance_t greedy_inst;	//greedy random start GRASP0.1 and 2-opt
	greedy_inst.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	greedy_inst.input_file_name = NULL;
	greedy_inst.x_bound = -1;
	greedy_inst.y_bound = -1;
	greedy_inst.num_nodes = -1;
	greedy_inst.random_seed = time(NULL);
	greedy_inst.sol_procedure_flag = GREEDY;
	greedy_inst.starting_index = -1;
	greedy_inst.prob_ign_opt = 0.1;
	greedy_inst.refine_flag = TWO_OPT;
	greedy_inst.best_sol = NULL;
	greedy_inst.best_sol_cost = DBL_INFY;
	greedy_inst.costs = NULL;
	greedy_inst.time_left = greedy_inst.time_limit;

	tsp_instance_t extra_mileage_inst;	//extra-mileage random start GRASP0.1 and 2-opt
	extra_mileage_inst.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	extra_mileage_inst.input_file_name = NULL;
	extra_mileage_inst.x_bound = -1;
	extra_mileage_inst.y_bound = -1;
	extra_mileage_inst.num_nodes = -1;
	extra_mileage_inst.random_seed = time(NULL);
	extra_mileage_inst.sol_procedure_flag = EXTRA_MILEAGE;
	extra_mileage_inst.starting_index = -1;
	extra_mileage_inst.prob_ign_opt = 0.1;
	extra_mileage_inst.refine_flag = TWO_OPT;
	extra_mileage_inst.best_sol = NULL;
	extra_mileage_inst.best_sol_cost = DBL_INFY;
	extra_mileage_inst.costs = NULL;
	extra_mileage_inst.time_left = extra_mileage_inst.time_limit;

#if VERBOSE > 1
	{ printf("Total time limit: %f...\nTime limit per instance %f...\n", TOTAL_TIME_LIMIT, TOTAL_TIME_LIMIT / NUM_RERUNS); }
#endif

	FILE* perf_prof_file = fopen(PERF_PROF_FILENAME, "w");
	if (perf_prof_file == NULL) { fprintf(stderr, "cannot open the perf_prof_file\n"); exit(1); }

	fprintf(perf_prof_file, "2, greedy_rand_start+GRASP0.1+2-opt, extra-mileage_rand_start+GRASP0.1+2-opt\n");

	srand(RAND_SEED); for (int i = 0; i < MIN_RAND_RUNS + log(1 + RAND_SEED); i++) rand();
	point_2d_t* nodes = (point_2d_t*)malloc(NUM_POINTS * sizeof(point_2d_t));
	int num_costs = NUM_POINTS * (NUM_POINTS - 1) / 2;
	double* costs = (double*)malloc(num_costs * sizeof(double));
	unsigned int random_seed;
	for (int i = 0; i < NUM_RERUNS; i++) {

		//generating the random nodes
		for (int i = 0; i < NUM_POINTS; i++) {
			nodes[i].x_coord = ((double)rand() / RAND_MAX) * POINT_BOUND;
			nodes[i].y_coord = ((double)rand() / RAND_MAX) * POINT_BOUND;
		}
		greedy_inst.num_nodes = NUM_POINTS;
		extra_mileage_inst.num_nodes = NUM_POINTS;
		greedy_inst.nodes = nodes;
		extra_mileage_inst.nodes = nodes;

		//precomputing the costs
		for (int i = 0; i < NUM_POINTS - 1; i++) {
			for (int j = i + 1; j < NUM_POINTS; j++) {
				costs[DIST_INDEX(i, j, NUM_POINTS)] = (double)((int)(dist(i, j, &greedy_inst) + 0.5));
			}
		}
		greedy_inst.costs = costs;
		extra_mileage_inst.costs = costs;

		//setting-up the random seeds
		random_seed = rand();
		greedy_inst.random_seed = random_seed;
		extra_mileage_inst.random_seed = random_seed;

		//finding the best solutions
		double greedy_best_sol_cost = DBL_INFY;
		greedy_inst.time_left = greedy_inst.time_limit;
		while (greedy_inst.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&greedy_inst);
			ref_sol(&greedy_inst);
			if (greedy_inst.best_sol_cost < greedy_best_sol_cost)
				greedy_best_sol_cost = greedy_inst.best_sol_cost;
			free(greedy_inst.best_sol);
			greedy_inst.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			greedy_inst.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}

		double extra_mileage_best_sol_cost = DBL_INFY;
		extra_mileage_inst.time_left = extra_mileage_inst.time_limit;
		while (extra_mileage_inst.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&extra_mileage_inst);
			ref_sol(&extra_mileage_inst);
			if (extra_mileage_inst.best_sol_cost < extra_mileage_best_sol_cost)
				extra_mileage_best_sol_cost = extra_mileage_inst.best_sol_cost;
			free(extra_mileage_inst.best_sol);
			extra_mileage_inst.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			extra_mileage_inst.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		
		fprintf(perf_prof_file, "random_seed_%d, %f, %f\n", random_seed, greedy_best_sol_cost, extra_mileage_best_sol_cost);
	}
	free(nodes);
	free(costs);
	
	fclose(perf_prof_file);
}

int main(int argc, char const* argv[]) {

#if VERBOSE > 0 
	{
		printf("Input command:\n");
		for (int i = 0; i < argc - 1; i++)
			printf("%s ", argv[i]);
		printf("%s\n", argv[argc - 1]);
	}
#endif

	double beg_time = seconds();

	////uncomment for best solution at first try////
	tsp_instance_t instance;
	set_up_tsp_instance(argc, argv, &instance);
	first_try_sol(&instance);
	free_tsp_instance(&instance);

	////uncomment for searching until time runs out////
	/*tsp_instance_t instance;
	set_up_tsp_instance(argc, argv, &instance);
	cont_search(&instance);
	free_tsp_instance(&instance);*/

	////uncomment for perf prof of chap 2.3////
	//comp_methods_1();

	////uncomment for perf prof of chap 3.2 fig 1////
	//comp_methods_2();

	////uncomment for perf prof of chap 3.2 fig 2////
	//comp_methods_3();

	//performance comparisons among solution strategies greedy_random_start+GRASP0.1+2opt and extra-mileage_random_start+GRASP0.1+2opt
	//comp_methods_0();

	double end_time = seconds();

#if VERBOSE > 0
	{
		printf("Execution time: %fs\n", end_time - beg_time);
	}
#endif

	return 0;
}