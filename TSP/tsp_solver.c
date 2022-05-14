#include <stdio.h>
#include <stdlib.h>
#include "tsp_utils.h"
#include "tsp.h"

#define TOTAL_TIME_LIMIT 12000.0
#define NUM_RERUNS 100
#define NUM_POINTS 300
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

static void genetic_procedure(tsp_instance_t* instance) {
	//genetic procedure
	if (tsp_opt(instance)) { free_tsp_instance(instance); fprintf(stderr, "Optimization algorithm failed\n"); exit(1); }
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

static void cplex_sol(tsp_instance_t* instance) {
	//cplex solution
	if (TSPopt(instance)) { free_tsp_instance(instance); fprintf(stderr, "Optimization algorithm failed\n"); exit(1); }
	//printf("Solution cost: %f\n", instance->best_sol_cost);
	if (plot_cycles(instance)) { free_tsp_instance(instance); fprintf(stderr, "Cycles plotting failed\n"); exit(1); }
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
	greedy_inst_1.metaheur_flag = NO_MH;
	greedy_inst_1.best_sol = NULL;
	greedy_inst_1.best_sol_cost = DBL_INFY;
	greedy_inst_1.costs = NULL;
	greedy_inst_1.time_left = greedy_inst_1.time_limit;
	greedy_inst_1.min_tenure = -1;
	greedy_inst_1.max_tenure = -1;

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
	greedy_inst_2.metaheur_flag = NO_MH;
	greedy_inst_2.best_sol = NULL;
	greedy_inst_2.best_sol_cost = DBL_INFY;
	greedy_inst_2.costs = NULL;
	greedy_inst_2.time_left = greedy_inst_2.time_limit;
	greedy_inst_2.min_tenure = -1;
	greedy_inst_2.max_tenure = -1;

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
	greedy_inst_3.metaheur_flag = NO_MH;
	greedy_inst_3.best_sol = NULL;
	greedy_inst_3.best_sol_cost = DBL_INFY;
	greedy_inst_3.costs = NULL;
	greedy_inst_3.time_left = greedy_inst_3.time_limit;
	greedy_inst_3.min_tenure = -1;
	greedy_inst_3.max_tenure = -1;

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
	extra_mileage_inst_1.metaheur_flag = NO_MH;
	extra_mileage_inst_1.best_sol = NULL;
	extra_mileage_inst_1.best_sol_cost = DBL_INFY;
	extra_mileage_inst_1.costs = NULL;
	extra_mileage_inst_1.time_left = extra_mileage_inst_1.time_limit;
	extra_mileage_inst_1.min_tenure = -1;
	extra_mileage_inst_1.max_tenure = -1;

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
	extra_mileage_inst_2.metaheur_flag = NO_MH;
	extra_mileage_inst_2.best_sol = NULL;
	extra_mileage_inst_2.best_sol_cost = DBL_INFY;
	extra_mileage_inst_2.costs = NULL;
	extra_mileage_inst_2.time_left = extra_mileage_inst_2.time_limit;
	extra_mileage_inst_2.min_tenure = -1;
	extra_mileage_inst_2.max_tenure = -1;

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
	extra_mileage_inst_3.metaheur_flag = NO_MH;
	extra_mileage_inst_3.best_sol = NULL;
	extra_mileage_inst_3.best_sol_cost = DBL_INFY;
	extra_mileage_inst_3.costs = NULL;
	extra_mileage_inst_3.time_left = extra_mileage_inst_3.time_limit;
	extra_mileage_inst_3.min_tenure = -1;
	extra_mileage_inst_3.max_tenure = -1;

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
	greedy_inst_1.metaheur_flag = NO_MH;
	greedy_inst_1.best_sol = NULL;
	greedy_inst_1.best_sol_cost = DBL_INFY;
	greedy_inst_1.costs = NULL;
	greedy_inst_1.time_left = greedy_inst_1.time_limit;
	greedy_inst_1.min_tenure = -1;
	greedy_inst_1.max_tenure = -1;

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
	greedy_inst_2.metaheur_flag = NO_MH;
	greedy_inst_2.best_sol = NULL;
	greedy_inst_2.best_sol_cost = DBL_INFY;
	greedy_inst_2.costs = NULL;
	greedy_inst_2.time_left = greedy_inst_2.time_limit;
	greedy_inst_2.min_tenure = -1;
	greedy_inst_2.max_tenure = -1;

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
	greedy_inst_3.metaheur_flag = NO_MH;
	greedy_inst_3.best_sol = NULL;
	greedy_inst_3.best_sol_cost = DBL_INFY;
	greedy_inst_3.costs = NULL;
	greedy_inst_3.time_left = greedy_inst_3.time_limit;
	greedy_inst_3.min_tenure = -1;
	greedy_inst_3.max_tenure = -1;

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
	extra_mileage_inst_1.metaheur_flag = NO_MH;
	extra_mileage_inst_1.best_sol = NULL;
	extra_mileage_inst_1.best_sol_cost = DBL_INFY;
	extra_mileage_inst_1.costs = NULL;
	extra_mileage_inst_1.time_left = extra_mileage_inst_1.time_limit;
	extra_mileage_inst_1.min_tenure = -1;
	extra_mileage_inst_1.max_tenure = -1;

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
	extra_mileage_inst_2.metaheur_flag = NO_MH;
	extra_mileage_inst_2.best_sol = NULL;
	extra_mileage_inst_2.best_sol_cost = DBL_INFY;
	extra_mileage_inst_2.costs = NULL;
	extra_mileage_inst_2.time_left = extra_mileage_inst_2.time_limit;
	extra_mileage_inst_2.min_tenure = -1;
	extra_mileage_inst_2.max_tenure = -1;

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
	extra_mileage_inst_3.metaheur_flag = NO_MH;
	extra_mileage_inst_3.best_sol = NULL;
	extra_mileage_inst_3.best_sol_cost = DBL_INFY;
	extra_mileage_inst_3.costs = NULL;
	extra_mileage_inst_3.time_left = extra_mileage_inst_3.time_limit;
	extra_mileage_inst_3.min_tenure = -1;
	extra_mileage_inst_3.max_tenure = -1;

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
	greedy_inst_1.metaheur_flag = NO_MH;
	greedy_inst_1.best_sol = NULL;
	greedy_inst_1.best_sol_cost = DBL_INFY;
	greedy_inst_1.costs = NULL;
	greedy_inst_1.time_left = greedy_inst_1.time_limit;
	greedy_inst_1.min_tenure = -1;
	greedy_inst_1.max_tenure = -1;

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
	greedy_inst_2.metaheur_flag = NO_MH;
	greedy_inst_2.best_sol = NULL;
	greedy_inst_2.best_sol_cost = DBL_INFY;
	greedy_inst_2.costs = NULL;
	greedy_inst_2.time_left = greedy_inst_2.time_limit;
	greedy_inst_2.min_tenure = -1;
	greedy_inst_2.max_tenure = -1;

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
	extra_mileage_inst_1.metaheur_flag = NO_MH;
	extra_mileage_inst_1.best_sol = NULL;
	extra_mileage_inst_1.best_sol_cost = DBL_INFY;
	extra_mileage_inst_1.costs = NULL;
	extra_mileage_inst_1.time_left = extra_mileage_inst_1.time_limit;
	extra_mileage_inst_1.min_tenure = -1;
	extra_mileage_inst_1.max_tenure = -1;

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
	extra_mileage_inst_2.metaheur_flag = NO_MH;
	extra_mileage_inst_2.best_sol = NULL;
	extra_mileage_inst_2.best_sol_cost = DBL_INFY;
	extra_mileage_inst_2.costs = NULL;
	extra_mileage_inst_2.time_left = extra_mileage_inst_2.time_limit;
	extra_mileage_inst_2.min_tenure = -1;
	extra_mileage_inst_2.max_tenure = -1;

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

//chap. 4.5 fig. 1
static void comp_methods_4() {

#if VERBOSE > 1
	{ printf("Comparing greedy_random_start + 2-opt + tabu search, greedy_random_start + 2-opt + VNS\n \
			  greedy_random_start + 2-opt + simulated annealing, genetic algorithm (greedy_random_start + 2-opt)...\n"); }
#endif

	tsp_instance_t tabu_search;	//tabu search
	tabu_search.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	tabu_search.input_file_name = NULL;
	tabu_search.x_bound = -1;
	tabu_search.y_bound = -1;
	tabu_search.num_nodes = -1;
	tabu_search.nodes = NULL;
	tabu_search.random_seed = time(NULL);
	tabu_search.sol_procedure_flag = GREEDY;
	tabu_search.starting_index = -1;
	tabu_search.prob_ign_opt = 0.0;
	tabu_search.refine_flag = TWO_OPT;
	tabu_search.metaheur_flag = TABU;
	tabu_search.min_tenure = -1;
	tabu_search.max_tenure = -1;
	tabu_search.min_temperature = -1;
	tabu_search.max_temperature = -1;
	tabu_search.move_weight = 25;
	tabu_search.pop_size = 100;
	tabu_search.best_sol = NULL;
	tabu_search.best_sol_cost = DBL_INFY;
	tabu_search.costs = NULL;
	tabu_search.time_left = tabu_search.time_limit;
	tabu_search.tabu_list = NULL;

	tsp_instance_t vns;	//VNS
	vns.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	vns.input_file_name = NULL;
	vns.x_bound = -1;
	vns.y_bound = -1;
	vns.num_nodes = -1;
	vns.nodes = NULL;
	vns.random_seed = time(NULL);
	vns.sol_procedure_flag = GREEDY;
	vns.starting_index = -1;
	vns.prob_ign_opt = 0.0;
	vns.refine_flag = TWO_OPT;
	vns.metaheur_flag = VNS;
	vns.min_tenure = -1;
	vns.max_tenure = -1;
	vns.min_temperature = -1;
	vns.max_temperature = -1;
	vns.move_weight = 25;
	vns.pop_size = 100;
	vns.best_sol = NULL;
	vns.best_sol_cost = DBL_INFY;
	vns.costs = NULL;
	vns.time_left = vns.time_limit;
	vns.tabu_list = NULL;

	tsp_instance_t sim_ann;	//simulated annealing
	sim_ann.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	sim_ann.input_file_name = NULL;
	sim_ann.x_bound = -1;
	sim_ann.y_bound = -1;
	sim_ann.num_nodes = -1;
	sim_ann.nodes = NULL;
	sim_ann.random_seed = time(NULL);
	sim_ann.sol_procedure_flag = GREEDY;
	sim_ann.starting_index = -1;
	sim_ann.prob_ign_opt = 0.0;
	sim_ann.refine_flag = TWO_OPT;
	sim_ann.metaheur_flag = SIM_ANNEAL;
	sim_ann.min_tenure = -1;
	sim_ann.max_tenure = -1;
	sim_ann.min_temperature = -1;
	sim_ann.max_temperature = -1;
	sim_ann.move_weight = 25;
	sim_ann.pop_size = 100;
	sim_ann.best_sol = NULL;
	sim_ann.best_sol_cost = DBL_INFY;
	sim_ann.costs = NULL;
	sim_ann.time_left = sim_ann.time_limit;
	sim_ann.tabu_list = NULL;

	tsp_instance_t gen_alg;	//genetic algorithm
	gen_alg.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	gen_alg.input_file_name = NULL;
	gen_alg.x_bound = -1;
	gen_alg.y_bound = -1;
	gen_alg.num_nodes = -1;
	gen_alg.nodes = NULL;
	gen_alg.random_seed = time(NULL);
	gen_alg.sol_procedure_flag = GREEDY;
	gen_alg.starting_index = -1;
	gen_alg.prob_ign_opt = 0.0;
	gen_alg.refine_flag = TWO_OPT;
	gen_alg.metaheur_flag = GEN;
	gen_alg.min_tenure = -1;
	gen_alg.max_tenure = -1;
	gen_alg.min_temperature = -1;
	gen_alg.max_temperature = -1;
	gen_alg.move_weight = 25;
	gen_alg.pop_size = 100;
	gen_alg.best_sol = NULL;
	gen_alg.best_sol_cost = DBL_INFY;
	gen_alg.costs = NULL;
	gen_alg.time_left = gen_alg.time_limit;
	gen_alg.tabu_list = NULL;

#if VERBOSE > 1
	{ printf("Total time limit: %f...\nTime limit per instance %f...\n", TOTAL_TIME_LIMIT, TOTAL_TIME_LIMIT / NUM_RERUNS); }
#endif

	FILE* perf_prof_file = fopen(PERF_PROF_FILENAME, "w");
	if (perf_prof_file == NULL) { fprintf(stderr, "cannot open the perf_prof_file\n"); exit(1); }

	fprintf(perf_prof_file, "4, greedy_rand_start+2-opt+tabu, greedy_rand_start+2-opt+vns, greedy_rand_start+2-opt+sim_anneal, greedy_rand_start+2-opt+genetic\n");

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
		tabu_search.num_nodes = NUM_POINTS;
		tabu_search.nodes = nodes;
		vns.num_nodes = NUM_POINTS;
		vns.nodes = nodes;
		sim_ann.num_nodes = NUM_POINTS;
		sim_ann.nodes = nodes;
		gen_alg.num_nodes = NUM_POINTS;
		gen_alg.nodes = nodes;

		//precomputing the costs
		for (int i = 0; i < NUM_POINTS - 1; i++) {
			for (int j = i + 1; j < NUM_POINTS; j++) {
				costs[DIST_INDEX(i, j, NUM_POINTS)] = (double)((int)(dist(i, j, &tabu_search) + 0.5));
			}
		}
		tabu_search.costs = costs;
		vns.costs = costs;
		sim_ann.costs = costs;
		gen_alg.costs = costs;

		//setting-up the random seeds
		random_seed = rand();
		tabu_search.random_seed = random_seed;
		vns.random_seed = random_seed;
		sim_ann.random_seed = random_seed;
		gen_alg.random_seed = random_seed;

		//finding the best solutions
		double tabu_search_best_sol_cost = DBL_INFY;
		tabu_search.time_left = tabu_search.time_limit;
		while (tabu_search.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&tabu_search);
			ref_sol(&tabu_search);
			if (tabu_search.best_sol_cost < tabu_search_best_sol_cost)
				tabu_search_best_sol_cost = tabu_search.best_sol_cost;
			free(tabu_search.best_sol);
			tabu_search.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			tabu_search.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double vns_best_sol_cost = DBL_INFY;
		vns.time_left = vns.time_limit;
		while (vns.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&vns);
			ref_sol(&vns);
			if (vns.best_sol_cost < vns_best_sol_cost)
				vns_best_sol_cost = vns.best_sol_cost;
			free(vns.best_sol);
			vns.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			vns.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double sim_ann_best_sol_cost = DBL_INFY;
		sim_ann.time_left = sim_ann.time_limit;
		while (sim_ann.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&sim_ann);
			ref_sol(&sim_ann);
			if (sim_ann.best_sol_cost < sim_ann_best_sol_cost)
				sim_ann_best_sol_cost = sim_ann.best_sol_cost;
			free(sim_ann.best_sol);
			sim_ann.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			sim_ann.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double gen_alg_best_sol_cost = DBL_INFY;
		gen_alg.time_left = gen_alg.time_limit;
		while (gen_alg.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&gen_alg);
			if (gen_alg.best_sol_cost < gen_alg_best_sol_cost)
				gen_alg_best_sol_cost = gen_alg.best_sol_cost;
			free(gen_alg.best_sol);
			gen_alg.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			gen_alg.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}

		fprintf(perf_prof_file, "random_seed_%d, %f, %f, %f, %f\n", random_seed, tabu_search_best_sol_cost, vns_best_sol_cost, sim_ann_best_sol_cost, gen_alg_best_sol_cost);
	}
	free(nodes);
	free(costs);

	fclose(perf_prof_file);
}

//chap. 4.5 fig. 2
static void comp_methods_5() {

#if VERBOSE > 1
	{ printf("Comparing greedy_random_start + 2-opt + tabu search, greedy_random_start + 2-opt + VNS, greedy_random_start + 2-opt + simulated annealing...\n"); }
#endif

	tsp_instance_t tabu_search;	//tabu search
	tabu_search.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	tabu_search.input_file_name = NULL;
	tabu_search.x_bound = -1;
	tabu_search.y_bound = -1;
	tabu_search.num_nodes = -1;
	tabu_search.nodes = NULL;
	tabu_search.random_seed = time(NULL);
	tabu_search.sol_procedure_flag = GREEDY;
	tabu_search.starting_index = -1;
	tabu_search.prob_ign_opt = 0.0;
	tabu_search.refine_flag = TWO_OPT;
	tabu_search.metaheur_flag = TABU;
	tabu_search.min_tenure = -1;
	tabu_search.max_tenure = -1;
	tabu_search.min_temperature = -1;
	tabu_search.max_temperature = -1;
	tabu_search.move_weight = 25;
	tabu_search.pop_size = 100;
	tabu_search.best_sol = NULL;
	tabu_search.best_sol_cost = DBL_INFY;
	tabu_search.costs = NULL;
	tabu_search.time_left = tabu_search.time_limit;
	tabu_search.tabu_list = NULL;

	tsp_instance_t vns;	//VNS
	vns.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	vns.input_file_name = NULL;
	vns.x_bound = -1;
	vns.y_bound = -1;
	vns.num_nodes = -1;
	vns.nodes = NULL;
	vns.random_seed = time(NULL);
	vns.sol_procedure_flag = GREEDY;
	vns.starting_index = -1;
	vns.prob_ign_opt = 0.0;
	vns.refine_flag = TWO_OPT;
	vns.metaheur_flag = VNS;
	vns.min_tenure = -1;
	vns.max_tenure = -1;
	vns.min_temperature = -1;
	vns.max_temperature = -1;
	vns.move_weight = 25;
	vns.pop_size = 100;
	vns.best_sol = NULL;
	vns.best_sol_cost = DBL_INFY;
	vns.costs = NULL;
	vns.time_left = vns.time_limit;
	vns.tabu_list = NULL;

	tsp_instance_t sim_ann;	//simulated annealing
	sim_ann.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
	sim_ann.input_file_name = NULL;
	sim_ann.x_bound = -1;
	sim_ann.y_bound = -1;
	sim_ann.num_nodes = -1;
	sim_ann.nodes = NULL;
	sim_ann.random_seed = time(NULL);
	sim_ann.sol_procedure_flag = GREEDY;
	sim_ann.starting_index = -1;
	sim_ann.prob_ign_opt = 0.0;
	sim_ann.refine_flag = TWO_OPT;
	sim_ann.metaheur_flag = SIM_ANNEAL;
	sim_ann.min_tenure = -1;
	sim_ann.max_tenure = -1;
	sim_ann.min_temperature = -1;
	sim_ann.max_temperature = -1;
	sim_ann.move_weight = 25;
	sim_ann.pop_size = 100;
	sim_ann.best_sol = NULL;
	sim_ann.best_sol_cost = DBL_INFY;
	sim_ann.costs = NULL;
	sim_ann.time_left = sim_ann.time_limit;
	sim_ann.tabu_list = NULL;

#if VERBOSE > 1
	{ printf("Total time limit: %f...\nTime limit per instance %f...\n", TOTAL_TIME_LIMIT, TOTAL_TIME_LIMIT / NUM_RERUNS); }
#endif

	FILE* perf_prof_file = fopen(PERF_PROF_FILENAME, "w");
	if (perf_prof_file == NULL) { fprintf(stderr, "cannot open the perf_prof_file\n"); exit(1); }

	fprintf(perf_prof_file, "3, greedy_rand_start+2-opt+tabu, greedy_rand_start+2-opt+vns, greedy_rand_start+2-opt+sim_anneal\n");

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
		tabu_search.num_nodes = NUM_POINTS;
		tabu_search.nodes = nodes;
		vns.num_nodes = NUM_POINTS;
		vns.nodes = nodes;
		sim_ann.num_nodes = NUM_POINTS;
		sim_ann.nodes = nodes;

		//precomputing the costs
		for (int i = 0; i < NUM_POINTS - 1; i++) {
			for (int j = i + 1; j < NUM_POINTS; j++) {
				costs[DIST_INDEX(i, j, NUM_POINTS)] = (double)((int)(dist(i, j, &tabu_search) + 0.5));
			}
		}
		tabu_search.costs = costs;
		vns.costs = costs;
		sim_ann.costs = costs;

		//setting-up the random seeds
		random_seed = rand();
		tabu_search.random_seed = random_seed;
		vns.random_seed = random_seed;
		sim_ann.random_seed = random_seed;

		//finding the best solutions
		double tabu_search_best_sol_cost = DBL_INFY;
		tabu_search.time_left = tabu_search.time_limit;
		while (tabu_search.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&tabu_search);
			ref_sol(&tabu_search);
			if (tabu_search.best_sol_cost < tabu_search_best_sol_cost)
				tabu_search_best_sol_cost = tabu_search.best_sol_cost;
			free(tabu_search.best_sol);
			tabu_search.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			tabu_search.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double vns_best_sol_cost = DBL_INFY;
		vns.time_left = vns.time_limit;
		while (vns.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&vns);
			ref_sol(&vns);
			if (vns.best_sol_cost < vns_best_sol_cost)
				vns_best_sol_cost = vns.best_sol_cost;
			free(vns.best_sol);
			vns.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			vns.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}
		double sim_ann_best_sol_cost = DBL_INFY;
		sim_ann.time_left = sim_ann.time_limit;
		while (sim_ann.time_left > 0) {
			double elapsed_timer_start = seconds();
			tsp_opt(&sim_ann);
			ref_sol(&sim_ann);
			if (sim_ann.best_sol_cost < sim_ann_best_sol_cost)
				sim_ann_best_sol_cost = sim_ann.best_sol_cost;
			free(sim_ann.best_sol);
			sim_ann.best_sol = NULL;
			double elapsed_timer_stop = seconds();
			sim_ann.time_left -= (elapsed_timer_stop - elapsed_timer_start);
		}

		fprintf(perf_prof_file, "random_seed_%d, %f, %f, %f\n", random_seed, tabu_search_best_sol_cost, vns_best_sol_cost, sim_ann_best_sol_cost);
	}
	free(nodes);
	free(costs);

	fclose(perf_prof_file);
}

//chap. 5.4
static void comp_methods_6() {

#if VERBOSE > 1
	{ printf("Comparing Benders' method, Branch and cut (callbacks)...\n"); 
	  printf("Number of nodes: %d...\nNumber of reruns %d...\n", NUM_POINTS, NUM_RERUNS);
	  printf("Total time limit: %f...\nTime limit per instance %f...\n", TOTAL_TIME_LIMIT, TOTAL_TIME_LIMIT / NUM_RERUNS); }
#endif

	FILE* perf_prof_file = fopen(PERF_PROF_FILENAME, "w");
	if (perf_prof_file == NULL) { fprintf(stderr, "cannot open the perf_prof_file\n"); exit(1); }

	fprintf(perf_prof_file, "2, benders, callback\n");

	srand(RAND_SEED); for (int i = 0; i < MIN_RAND_RUNS + log(1 + RAND_SEED); i++) rand();
	point_2d_t* nodes = (point_2d_t*)malloc(NUM_POINTS * sizeof(point_2d_t));
	int num_costs = NUM_POINTS * (NUM_POINTS - 1) / 2;
	double* costs = (double*)malloc(num_costs * sizeof(double));
	unsigned int random_seed;
	for (int i = 0; i < NUM_RERUNS; i++) {

		tsp_instance_t benders;	//benders' method
		benders.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
		benders.input_file_name = NULL;
		benders.x_bound = -1;
		benders.y_bound = -1;
		benders.sol_procedure_flag = GREEDY;
		benders.starting_index = -1;
		benders.prob_ign_opt = 0.0;
		benders.refine_flag = TWO_OPT;
		benders.metaheur_flag = NO_MH;
		benders.min_tenure = -1;
		benders.max_tenure = -1;
		benders.min_temperature = -1;
		benders.max_temperature = -1;
		benders.move_weight = 25;
		benders.pop_size = 100;
		benders.cplex_solver_flag = BENDERS;
		benders.best_sol = NULL;
		benders.best_sol_cost = DBL_INFY;
		benders.time_left = benders.time_limit;
		benders.tabu_list = NULL;

		tsp_instance_t bnc;	//Branch and cut
		bnc.time_limit = TOTAL_TIME_LIMIT / NUM_RERUNS;
		bnc.input_file_name = NULL;
		bnc.x_bound = -1;
		bnc.y_bound = -1;
		bnc.sol_procedure_flag = GREEDY;
		bnc.starting_index = -1;
		bnc.prob_ign_opt = 0.0;
		bnc.refine_flag = TWO_OPT;
		bnc.metaheur_flag = NO_MH;
		bnc.min_tenure = -1;
		bnc.max_tenure = -1;
		bnc.min_temperature = -1;
		bnc.max_temperature = -1;
		bnc.move_weight = 25;
		bnc.pop_size = 100;
		bnc.cplex_solver_flag = CALLBACK;
		bnc.best_sol = NULL;
		bnc.best_sol_cost = DBL_INFY;
		bnc.time_left = bnc.time_limit;
		bnc.tabu_list = NULL;

		//generating the random nodes
		for (int i = 0; i < NUM_POINTS; i++) {
			nodes[i].x_coord = ((double)rand() / RAND_MAX) * POINT_BOUND;
			nodes[i].y_coord = ((double)rand() / RAND_MAX) * POINT_BOUND;
		}
		benders.num_nodes = NUM_POINTS;
		benders.nodes = nodes;
		benders.num_cycles = benders.num_nodes;
		benders.cycle_delimiter = benders.num_nodes;
		bnc.num_nodes = NUM_POINTS;
		bnc.nodes = nodes;
		bnc.num_cycles = bnc.num_nodes;
		bnc.cycle_delimiter = bnc.num_nodes;

		//precomputing the costs
		for (int i = 0; i < NUM_POINTS - 1; i++) {
			for (int j = i + 1; j < NUM_POINTS; j++) {
				costs[DIST_INDEX(i, j, NUM_POINTS)] = (double)((int)(dist(i, j, &benders) + 0.5));
			}
		}
		benders.costs = costs;
		bnc.costs = costs;

		//setting-up the random seeds
		random_seed = rand();
		benders.random_seed = random_seed;
		bnc.random_seed = random_seed;

		//finding the best solutions
		double elapsed_timer_start, elapsed_timer_stop;

		double benders_sol_cost = DBL_INFY;
		double benders_computation_time = 0;
		elapsed_timer_start = seconds();
		if (TSPopt(&benders)) { free_tsp_instance(&benders); fprintf(stderr, "Optimization algorithm failed\n"); exit(1); }
		elapsed_timer_stop = seconds();
		benders_sol_cost = benders.best_sol_cost;
		free(benders.best_sol);
		benders.best_sol = NULL;
		benders.best_sol_cost = DBL_INFY;
		benders_computation_time = (elapsed_timer_stop - elapsed_timer_start);

		double bnc_sol_cost = DBL_INFY;
		double bnc_computation_time = 0;
		elapsed_timer_start = seconds();
		if (TSPopt(&bnc)) { free_tsp_instance(&bnc); fprintf(stderr, "Optimization algorithm failed\n"); exit(1); }
		elapsed_timer_stop = seconds();
		bnc_sol_cost = bnc.best_sol_cost;
		free(bnc.best_sol);
		bnc.best_sol = NULL;
		bnc.best_sol_cost = DBL_INFY;
		bnc_computation_time = (elapsed_timer_stop - elapsed_timer_start);

#if VERBOSE > 1
		{ printf("Benders' method: cost = %f, time = %fs; Branch and cut: cost = %f, time = %fs\n", benders_sol_cost, benders_computation_time, bnc_sol_cost, bnc_computation_time); }
#endif

		fprintf(perf_prof_file, "random_seed_%d, %f, %f\n", random_seed, benders_computation_time, bnc_computation_time);
	}
	free(nodes);
	free(costs);

	fclose(perf_prof_file);
}

int main(int argc, char const* argv[]) {

#if VERBOSE > 0 
	{ printf("Input command:\n"); for (int i = 0; i < argc - 1; i++) printf("%s ", argv[i]); printf("%s\n", argv[argc - 1]); }
#endif

	double beg_time = seconds();

	//////uncomment for best solution at first try////
	//tsp_instance_t instance;
	//set_up_tsp_instance(argc, argv, &instance);
	//first_try_sol(&instance);
	//free_tsp_instance(&instance);

	//////uncomment for searching until time runs out////
	//tsp_instance_t instance;
	//set_up_tsp_instance(argc, argv, &instance);
	//cont_search(&instance);
	//free_tsp_instance(&instance);

	//////uncomment for the genetic algorithm instance (does not refine)////
	//tsp_instance_t instance;
	//set_up_tsp_instance(argc, argv, &instance);
	//genetic_procedure(&instance);
	//free_tsp_instance(&instance);

	////uncomment for method using cplex////
	tsp_instance_t instance;
	set_up_tsp_instance(argc, argv, &instance);
	cplex_sol(&instance);
	free_tsp_instance(&instance);

	////uncomment for perf prof of chap 2.3////
	//comp_methods_1();

	////uncomment for perf prof of chap 3.2 fig 1////
	//comp_methods_2();

	////uncomment for perf prof of chap 3.2 fig 2////
	//comp_methods_3();

	////uncomment for perf prof of chap 4.5 fig 1////
	//comp_methods_4();

	////uncomment for perf prof of chap 4.5 fig 2////
	//comp_methods_5();

	////uncomment for perf prof of chap 4.5 fig 2////
	//comp_methods_6();

	double end_time = seconds();

#if VERBOSE > 0
	{ printf("Execution time: %fs\n", end_time - beg_time); }
#endif

	return 0;
}