#include "tsp_utils.h"
#include "tsp.h"
#include <string.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#define INSTANCE_POINTS_PLOT "points.dat"
#define INSTANCE_TOUR_PLOT "tour.dat"
#define GNUPLOT_COMMAND_POINTS "gnuplot commands_points.txt"
#define GNUPLOT_COMMAND_TOUR "gnuplot commands_tour.txt"

#define TABU_LIST_MAX 100
#define TABU_LIST_MAX_DIV 10
#define TABU_LIST_MIN 75
#define TABU_LIST_MIN_DIV 15

int parse_command_line(const int argc, const char* argv[], tsp_instance_t* instance) {

#if VERBOSE > 2
	{ printf("Parsing command line args...\n"); }
#endif

	assert(instance != NULL);

	//default hyperparameter values
	instance->time_limit = MAX_TIME;
	instance->input_file_name = NULL;
	instance->x_bound = -1;
	instance->y_bound = -1;
	instance->num_nodes = -1;
	instance->random_seed = time(NULL);
	instance->sol_procedure_flag = SEQUENTIAL;
	instance->starting_index = 0;
	instance->prob_ign_opt = 0.0;
	instance->refine_flag = NO_REF;
	instance->metaheur_flag = NO_MH;
	instance->min_tenure = -1;
	instance->max_tenure = -1;

	int help = 0;
	if (argc < 1) help = 1;
	// #pragma omp parallel for num_threads(OMP_NUM_THREADS)
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-file_name") == 0) { instance->input_file_name = argv[++i]; continue; } 					// input file
		if (strcmp(argv[i], "-time_limit") == 0) { instance->time_limit = atof(argv[++i]); continue; }					// total time limit
		if (strcmp(argv[i], "-x_bound") == 0) { instance->x_bound = atoi(argv[++i]); continue; }						// random nodes hor. limit
		if (strcmp(argv[i], "-y_bound") == 0) { instance->y_bound = atoi(argv[++i]); continue; }						// random nodes ver. limit
		if (strcmp(argv[i], "-num_nodes") == 0) { instance->num_nodes = atoi(argv[++i]); continue; }					// number of random nodes
		if (strcmp(argv[i], "-random_seed") == 0) { instance->random_seed = atoi(argv[++i]); continue; }				// random seed
		if (strcmp(argv[i], "-proc_flag") == 0) { instance->sol_procedure_flag = atoi(argv[++i]); continue; }			// flag indicating the type of sol. procedure
		if (strcmp(argv[i], "-start_index") == 0) { instance->starting_index = atoi(argv[++i]); continue; }				// index of the initial node
		if (strcmp(argv[i], "-prob_ign_opt") == 0) { instance->prob_ign_opt = atof(argv[++i]); continue; }				// prob. to go to the 2nd best node
		if (strcmp(argv[i], "-ref_flag") == 0) { instance->refine_flag = atoi(argv[++i]); continue; }					// flag indicating the refinement procedure
		if (strcmp(argv[i], "-metaheur_flag") == 0) { instance->metaheur_flag = atoi(argv[++i]); continue; }			// flag indicating the metaheuristic method
		if (strcmp(argv[i], "-min_tenure") == 0) { instance->min_tenure = atoi(argv[++i]); continue; }					// the minimum size of the tenure (tabu search)
		if (strcmp(argv[i], "-max_tenure") == 0) { instance->max_tenure = atoi(argv[++i]); continue; }					// the maximum size of the tenure (tabu search)
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; } 													    // help
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; } 													// help
		help = 1;
	}

	//finalize parameter initialization
	instance->best_sol = NULL;
	instance->best_sol_cost = DBL_INFY;
	instance->costs = NULL;
	instance->time_left = instance->time_limit;
	instance->tabu_list = NULL;

	if (help || (VERBOSE > 1)) {
		printf("\n\navailable parameters -------------------------------------------------------------------------\n");
		printf("-file_name %s (.tsp file containing the input data.\n if no file_name is provided a random instance" \
			" will be generated)\n", instance->input_file_name);
		printf("-time_limit %E (maximum number of seconds the program is allowed to run)\n", instance->time_limit);
		printf("-x_bound %d (used when generating a random instance)\n", instance->x_bound);
		printf("-y_bound %d (used when generating a random instance)\n", instance->y_bound);
		printf("-num_nodes %d (used when generating a random instance)\n", instance->num_nodes);
		printf("-random_seed %d\n", instance->random_seed);
		printf("-proc_flag %d (1: SEQUENTIAL; 2: GREEDY; 3: EXTRA_MILEAGE)\n", instance->sol_procedure_flag);
		printf("-start_index %d (0: default; -1: random; -2: exhaustive)\n", instance->starting_index);
		printf("-prob_ign_opt %f (probability of selecting the second best node at each iteration)\n", instance->prob_ign_opt);
		printf("-ref_flag %d (1: NO_REF; 2: TWO_OPT)\n", instance->refine_flag);
		printf("-metaheur_flag %d (1: NO_MH; 2: TABU; 3: VNS)\n", instance->metaheur_flag);
		printf("-min_tenure %d (the minumum size of the tenure (tabu search))\n", instance->min_tenure);
		printf("-max_tenure %d (the maximum size of the tenure (tabu search))\n", instance->max_tenure);
		printf("\nenter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}

	return help;
}

static int read_input_file(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Reading input file...\n"); }
#endif

	assert(instance != NULL);

	//input file containing the tsp data
	FILE* input_file = fopen(instance->input_file_name, "r");
	if (input_file == NULL) { fprintf(stderr, "could not open input file\n"); return 1; }

	//default values for the data
	// instance->num_nodes = -1;	//ALREADY DONE DURING COMMAND LINE PARSING

	//used for parsing each line
	char line[180];
	char* par_name;
	char* token1;
	char* token2;

	int active_section = 0; // =1 NODE_COORD_SECTION

	while (fgets(line, sizeof(line), input_file) != NULL) {

#if VERBOSE > 4 
		{printf("%s", line); fflush(NULL); }
#endif

		if (strlen(line) <= 1) continue; // skip empty lines

		par_name = strtok(line, " :");
#if VERBOSE > 7 
		{printf("parameter \"%s\" \n", par_name); fflush(NULL); }
#endif

		if (strncmp(par_name, "NAME", 4) == 0) {
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "COMMENT", 7) == 0) {
			active_section = 0;
#if VERBOSE > 7 
			{
				token1 = strtok(NULL, "");
				printf("solving instance of %s\n\n", token1);
			}
#endif
			continue;
		}

		if (strncmp(par_name, "TYPE", 4) == 0) {
			token1 = strtok(NULL, " :");
			if (strncmp(token1, "TSP", 3) != 0) { fprintf(stderr, "only TYPE : TSP is accepted\n"); return 1; }
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "DIMENSION", 9) == 0) {
			if (instance->num_nodes >= 0) { fprintf(stderr, "redefining the number of nodes in the DIMENSION section\n"); return 1; }
			token1 = strtok(NULL, " :");
			instance->num_nodes = atoi(token1);
#if VERBOSE > 7 
			{printf("num_nodes %d\n", instance->num_nodes); }
#endif
			instance->nodes = (point_2d_t*)malloc(instance->num_nodes * sizeof(point_2d_t));
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0) {
			token1 = strtok(NULL, " :");
			if (strncmp(token1, "EUC_2D", 6) != 0) { fprintf(stderr, "only EDGE_WEIGHT_TYPE : EUC_2D is accepted\n"); return 1; }
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "NODE_COORD_SECTION", 18) == 0) {
			if (instance->num_nodes <= 0) { fprintf(stderr, "DIMENSION should appear before the NODE_COORD_SECTION\n"); return 1; }
			active_section = 1;
			continue;
		}

		if (strncmp(par_name, "EOF", 3) == 0) {
			active_section = 0;
			break;
		}

		if (active_section == 1) { // within NODE_COORD_SECTION
			int i = atoi(par_name) - 1;
			if (i < 0 || i >= instance->num_nodes) { fprintf(stderr, "node index out of bounds\n"); return 1; }
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			instance->nodes[i].x_coord = atof(token1);
			instance->nodes[i].y_coord = atof(token2);
#if VERBOSE > 7
			{printf("node %4d at coordinates (%15.7lf, %15.7lf)\n", i + 1, instance->nodes[i].x_coord, instance->nodes[i].y_coord); }
#endif
			continue;
		}

		fprintf(stderr, "failed to parse line, last acrive_section: %d\n", active_section);
		return 1;
	}

	fclose(input_file);

#if VERBOSE > 3
	{
		printf("read %d nodes:\n", instance->num_nodes);
		for (int i = 0; i < instance->num_nodes; i++) {
			printf("%d: (%f, %f)\n", i + 1, instance->nodes[i].x_coord, instance->nodes[i].y_coord);
		}
	}
#endif

	return 0;
}

static int generate_random_nodes(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Generating random instance...\n"); }
#endif

	assert(instance != NULL);

#if VERBOSE > 3 
	{printf("num_nodes %d\n", instance->num_nodes); }
#endif
	if (instance->num_nodes < 1) { fprintf(stderr, "random instance with less than one node\n"); return 1; }

	instance->nodes = (point_2d_t*)malloc(instance->num_nodes * sizeof(point_2d_t));

#if VERBOSE > 3
	{printf("random bounds [(0, %d), (0, %d)]\n", instance->x_bound, instance->y_bound); }
#endif
	if (instance->x_bound < 1 || instance->y_bound < 1) { fprintf(stderr, "random bounds less than one\n"); return 1; }
	srand(instance->random_seed); for (int i = 0; i < MIN_RAND_RUNS + log(1 + instance->random_seed); i++) rand();
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
	for (int i = 0; i < instance->num_nodes; i++) {
		instance->nodes[i].x_coord = ((double)rand() / RAND_MAX) * instance->x_bound;
		instance->nodes[i].y_coord = ((double)rand() / RAND_MAX) * instance->y_bound;
#if VERBOSE > 7
		{printf("generated random node (%f, %f)\n", instance->nodes[i].x_coord, instance->nodes[i].x_coord); }
#endif
	}

#if VERBOSE > 3
	{
		printf("generated %d random nodes:\n", instance->num_nodes);
		for (int i = 0; i < instance->num_nodes; i++) {
			printf("%d: (%f, %f)\n", i + 1, instance->nodes[i].x_coord, instance->nodes[i].y_coord);
		}
	}
#endif

	return 0;
}

int get_data(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Preparing the instance nodes...\n"); }
#endif

	assert(instance != NULL);
	if (instance->input_file_name != NULL)
		return read_input_file(instance);
	else
		return generate_random_nodes(instance);
}

int plot_points(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Plotting points...\n"); }
#endif

	assert(instance != NULL);
	FILE* points_plot_file = fopen(INSTANCE_POINTS_PLOT, "w");
	if (points_plot_file == NULL) { fprintf(stderr, "cannot open the points_plot_file\n"); return 1; }

	//save the points in the file using the gnuplot format
	for (int i = 0; i < instance->num_nodes; i++) {
		fprintf(points_plot_file, "%f %f\n\n", instance->nodes[i].x_coord, instance->nodes[i].y_coord);
	}

	fclose(points_plot_file);

	//execute the commands found in the commands.txt file. these commands will plot the points using gnuplot
	system(GNUPLOT_COMMAND_POINTS);

	return 0;
}

double dist(int i, int j, tsp_instance_t* instance) {

	assert(instance != NULL);
	assert(0 <= i && i <= instance->num_nodes);
	assert(0 <= i && i <= instance->num_nodes);

	double dx = fabs(instance->nodes[i].x_coord - instance->nodes[j].x_coord);
	double dy = fabs(instance->nodes[i].y_coord - instance->nodes[j].y_coord);

	double dist = sqrt(dx * dx + dy * dy);

#if VERBOSE > 8
	{
		printf("dx = %f, dy = %f, dist = %f\n", dx, dy, dist);
	}
#endif

	return dist;
}

int precompute_costs(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Filling in the cost look-up table...\n"); }
#endif

	assert(instance != NULL);

	//creating the array that contains the costs
	int num_costs = (instance->num_nodes * (instance->num_nodes - 1)) / 2;
#if VERBOSE > 2
	{
		printf("Allocateing %d-length array (%dbytes) for the cost look-up table\n", num_costs, num_costs * sizeof(double));
	}
#endif

	if (instance->costs == NULL) {
		instance->costs = (double*)malloc(num_costs * sizeof(double));
		if (instance->costs == NULL) { fprintf(stderr, "could not allocate memory for the cost look-up table\n"); return 1; }
	}
	else { fprintf(stderr, "already precomputed the cost look-up table\n"); return 1; }

	for (int i = 0; i < instance->num_nodes - 1; i++) {
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
		for (int j = i + 1; j < instance->num_nodes; j++) {
			instance->costs[DIST_INDEX(i, j, instance->num_nodes)] = dist(i, j, instance);
#if VERBOSE > 3
			{
				printf("The distance between %d and %d is %f, allocated at index %d\n", i, j, instance->costs[DIST_INDEX(i, j, instance->num_nodes)], \
					DIST_INDEX(i, j, instance->num_nodes));
			}
#endif
		}
	}

	return 0;
}

double lookup_cost(int i, int j, tsp_instance_t* instance) {

	assert(instance != NULL);
	assert(0 <= i && i <= instance->num_nodes);
	assert(0 <= i && i <= instance->num_nodes);

	if (i == j) return 0.0;
	return i < j ? instance->costs[DIST_INDEX(i, j, instance->num_nodes)] : instance->costs[DIST_INDEX(j, i, instance->num_nodes)];
}

int int_round_clut(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Rounding to the nearest int the cost look-up table...\n"); }
#endif

	assert(instance != NULL);

	for (int i = 0; i < instance->num_nodes - 1; i++) {
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
		for (int j = i + 1; j < instance->num_nodes; j++) {
			instance->costs[DIST_INDEX(i, j, instance->num_nodes)] = (double)((int)(instance->costs[DIST_INDEX(i, j, instance->num_nodes)] + 0.5));
#if VERBOSE > 3
			{
				printf("The distance between %d and %d is %f, allocated at index %d\n", i, j, instance->costs[DIST_INDEX(i, j, instance->num_nodes)], \
					DIST_INDEX(i, j, instance->num_nodes));
			}
#endif
		}
	}
}

int tsp_opt(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Solving the tsp instance...\n"); }
#endif

	assert(instance != NULL);

	if (instance->num_nodes < 2) {
		printf("No tour possible with %d nodes\n", instance->num_nodes);
	}
	else if (instance->num_nodes < 4) {
		instance->sol_procedure_flag = SEQUENTIAL;
		instance->refine_flag = NO_REF;
		instance->metaheur_flag = NO_MH;
		return tsp_seq_sol(instance);
	}
	else {
		switch (instance->sol_procedure_flag) {
		case SEQUENTIAL:
			return tsp_seq_sol(instance);
		case GREEDY:
			return tsp_gdy_sol(instance);
		case EXTRA_MILEAGE:
			return tsp_exm_sol(instance);
		default:
			return 1;
		}
	}
}

int tsp_seq_sol(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Applying the sequential method...\n"); }
#endif

	if (instance->best_sol == NULL) {
		instance->best_sol = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));
		if (instance->best_sol == NULL) { fprintf(stderr, "could not allocate memory for the best solution array\n"); return 1; }
	}
	instance->best_sol_cost = 0.0;

	for (int i = 0; i < instance->num_nodes; i++) {
		instance->best_sol[i] = (unsigned int)i;
		instance->best_sol_cost += lookup_cost(i, (i + 1) % instance->num_nodes, instance);
	}

	return 0;
}

static void tsp_gdy_sol_si(unsigned int starting_index, tsp_instance_t* instance) {

	srand(instance->random_seed); for (int i = 0; i < MIN_RAND_RUNS + log(1 + instance->random_seed); i++) rand();

#if VERBOSE > 2
	{ printf("Applying the greedy method with initial idex %d...\n", starting_index); }
#endif

	unsigned int* nodes2visit = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));
	for (unsigned int i = 0; i < instance->num_nodes; i++) {
		nodes2visit[i] = i;
	}

#if VERBOSE > 8
	{ printf("nodes2visit: [ "); for (int i = 0; i < instance->num_nodes; i++) printf("(%d, %d) ", i, nodes2visit[i]); printf("]\n"); }
#endif

	unsigned int current_node_index = starting_index;	//node under examination
	unsigned int current_tour_node_index = 0;	//node position in the tour
	swap_uint_array(current_node_index, current_tour_node_index, nodes2visit);	//put the starting node index at the beginning of the tour

#if VERBOSE > 9
	{ printf("Swaped %d with %d. Now node2visit[%d] = %d and node2visit[%d] = %d\n", current_node_index, current_tour_node_index,
		current_node_index, nodes2visit[current_node_index], current_tour_node_index, nodes2visit[current_tour_node_index]); }
#endif
#if VERBOSE > 8
	{ printf("nodes2visit: [ "); for (int i = 0; i < instance->num_nodes; i++) printf("(%d, %d) ", i, nodes2visit[i]); printf("]\n"); }
#endif

	current_node_index = current_tour_node_index;	//update the current node index
	unsigned int next_node_index_1 = ++current_tour_node_index;	//indicate that we are looking for the next tour position
	unsigned int next_node_index_2 = next_node_index_1;	//index of the second best choice

	for (int i = current_tour_node_index; i < instance->num_nodes; i++) {	//iterate over the tour positions

#if VERBOSE > 8
		{ printf("Looking for element %d of the tour\n", current_tour_node_index); }
#endif

		double nearest_neighbor_dist_1 = lookup_cost(nodes2visit[current_node_index], nodes2visit[next_node_index_1], instance);
		double nearest_neighbor_dist_2 = nearest_neighbor_dist_1;
		for (int j = current_tour_node_index + 1; j < instance->num_nodes; j++) {	//look for the nodes closest to the current node under exam.
			if (lookup_cost(nodes2visit[current_node_index], nodes2visit[j], instance) < nearest_neighbor_dist_2) {
				if (lookup_cost(nodes2visit[current_node_index], nodes2visit[j], instance) < nearest_neighbor_dist_1) {	//best node
					next_node_index_2 = next_node_index_1;
					nearest_neighbor_dist_2 = nearest_neighbor_dist_1;
					next_node_index_1 = j;
					nearest_neighbor_dist_1 = lookup_cost(nodes2visit[current_node_index], nodes2visit[next_node_index_1], instance);
				}
				else {	//second best node
					next_node_index_2 = j;
					nearest_neighbor_dist_2 = lookup_cost(nodes2visit[current_node_index], nodes2visit[next_node_index_2], instance);
				}
			}
		}

#if VERBOSE > 8
		{ printf("Nearest neighbors %d, %d of %d at distance %f, %f\n", nodes2visit[next_node_index_1], nodes2visit[next_node_index_2], nodes2visit[current_node_index],
			nearest_neighbor_dist_1, nearest_neighbor_dist_2); }
#endif

		current_node_index = ((double)rand() / (RAND_MAX + 1)) < instance->prob_ign_opt ? next_node_index_2 : next_node_index_1;	//choose best or second best
		swap_uint_array(current_node_index, current_tour_node_index, nodes2visit);

#if VERBOSE > 9
		{ printf("Swaped %d with %d. Now node2visit[%d] = %d and node2visit[%d] = %d\n", current_node_index, current_tour_node_index,
			current_node_index, nodes2visit[current_node_index], current_tour_node_index, nodes2visit[current_tour_node_index]); }
#endif
#if VERBOSE > 8
		{ printf("nodes2visit: [ "); for (int i = 0; i < instance->num_nodes; i++) printf("(%d, %d) ", i, nodes2visit[i]); printf("]\n"); }
#endif
		current_node_index = current_tour_node_index;
		next_node_index_1 = ++current_tour_node_index;
		next_node_index_2 = next_node_index_1;

	}

	double current_cost = 0.0;
	for (int i = 0; i < instance->num_nodes; i++) {
		current_cost += lookup_cost(nodes2visit[i], nodes2visit[(i + 1) % instance->num_nodes], instance);

#if VERBOSE > 8
		{ printf("Cost of edge [(%f, %f), (%f, %f)]: %f; total path cost: %f\n", instance->nodes[nodes2visit[i]].x_coord, instance->nodes[nodes2visit[i]].y_coord,
			instance->nodes[nodes2visit[(i + 1) % instance->num_nodes]].x_coord, instance->nodes[nodes2visit[(i + 1) % instance->num_nodes]].y_coord,
			lookup_cost(nodes2visit[i], nodes2visit[(i + 1) % instance->num_nodes], instance), current_cost); }
#endif

	}
#if VERBOSE > 4
	{ printf("Cost: %f\n", current_cost); }
#endif

	if (current_cost < instance->best_sol_cost) {
		instance->best_sol_cost = current_cost;
		for (int i = 0; i < instance->num_nodes; i++) {
			instance->best_sol[i] = nodes2visit[i];

#if VERBOSE > 2
			{ printf("Element of index %d of the tour is the node of index %d [%f, %f]\n", i, instance->best_sol[i], instance->nodes[instance->best_sol[i]].x_coord,
				instance->nodes[instance->best_sol[i]].y_coord); }
#endif
		}
	}

	free(nodes2visit);
}

int tsp_gdy_sol(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Applying the greedy heuristic...\n"); }
#endif

	if (instance->best_sol == NULL) {
		instance->best_sol = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));
		if (instance->best_sol == NULL) { fprintf(stderr, "could not allocate memory for the best solution array\n"); return 1; }
	}
	instance->best_sol_cost = DBL_INFY;

	int starting_index;
	if (instance->starting_index >= 0)
		starting_index = instance->starting_index;
	else if (instance->starting_index == -1) {
		srand(instance->random_seed); for (int i = 0; i < MIN_RAND_RUNS + log(1 + instance->random_seed); i++) rand();
		starting_index = (int)(((double)rand() / RAND_MAX) * (instance->num_nodes - 1) + 0.5);
	}
	else
		starting_index = -2;

#if VERBOSE > 2
	{ printf("Initial node index: %d\n", starting_index); }
#endif

	if (starting_index >= 0) {
#if VERBOSE > 3
		{ printf("Starting point %d (%f, %f)\n", starting_index, instance->nodes[starting_index].x_coord, instance->nodes[starting_index].y_coord); }
#endif
		tsp_gdy_sol_si(starting_index, instance);
	}
	else {
		for (int i = 0; i < instance->num_nodes; i++) {
#if VERBOSE > 3
			{ printf("Starting point %d (%f, %f)\n", i, instance->nodes[i].x_coord, instance->nodes[i].y_coord); }
#endif
			tsp_gdy_sol_si(i, instance);
		}
	}

	return 0;
}

static void tsp_exm_sol_se(unsigned int start_index, unsigned int end_index, tsp_instance_t* instance) {

	srand(instance->random_seed); for (int i = 0; i < MIN_RAND_RUNS + log(1 + instance->random_seed); i++) rand();

#if VERBOSE > 2
	{ printf("Applying the extra-mileage method with starting pair (%d, %d)...\n", start_index, end_index); }
#endif

	unsigned int num_nodes2visit = instance->num_nodes;
	unsigned int* nodes2visit = (unsigned int*)malloc(num_nodes2visit * sizeof(unsigned int));
	for (unsigned int i = 0; i < num_nodes2visit; i++) {
		nodes2visit[i] = i;
	}

#if VERBOSE > 8
	{ printf("nodes2visit: [ "); for (int i = 0; i < num_nodes2visit; i++) printf("(%d, %d) ", i, nodes2visit[i]); printf("]\n"); }
#endif

	int* successor = (int*)malloc(instance->num_nodes * sizeof(int));
	for (int i = 0; i < instance->num_nodes; i++) {
		successor[i] = -1;	//-1 indicates no successor
	}

#if VERBOSE > 8
	{ printf("successor: [ "); for (int i = 0; i < instance->num_nodes; i++) printf("(%d, %d) ", i, successor[i]); printf("]\n"); }
#endif

	successor[start_index] = end_index;
	successor[end_index] = start_index;
	nodes2visit[end_index] = nodes2visit[--num_nodes2visit];
	nodes2visit[start_index] = nodes2visit[--num_nodes2visit];

#if VERBOSE > 8
	{ printf("nodes2visit after initialization: [ "); for (int i = 0; i < num_nodes2visit; i++) printf("(%d, %d) ", i, nodes2visit[i]); printf("]\n"); }
#endif
#if VERBOSE > 8
	{ printf("successor after initialization: [ "); for (int i = 0; i < instance->num_nodes; i++) printf("(%d, %d) ", i, successor[i]); printf("]\n"); }
#endif

	while (num_nodes2visit > 0) {
		unsigned int node2add_index_1 = 0;	//index of the nodes2visit (not the node index)
		unsigned int node2add_index_2 = node2add_index_1;	//index of the second best node found
		unsigned int starting_node_1 = start_index;
		unsigned int starting_node_2 = starting_node_1;
		double min_cost_1 = lookup_cost(starting_node_1, nodes2visit[node2add_index_1], instance) + lookup_cost(nodes2visit[node2add_index_1], successor[starting_node_1], instance) - \
			lookup_cost(starting_node_1, successor[starting_node_1], instance);
		double min_cost_2 = min_cost_1;

#if VERBOSE > 8
		{ printf("Number of nodes not in the tour: %d\n", num_nodes2visit); }
#endif

		for (int i = 0; i < instance->num_nodes; i++) if (successor[i] != -1) {	//iterate over the already visited nodes
			for (int j = 1; j < num_nodes2visit; j++) {	//iterate over the unseen nodes
				if (lookup_cost(i, nodes2visit[j], instance) + lookup_cost(nodes2visit[j], successor[i], instance) - \
					lookup_cost(i, successor[i], instance) < min_cost_2) {
					if (lookup_cost(i, nodes2visit[j], instance) + lookup_cost(nodes2visit[j], successor[i], instance) - \
						lookup_cost(i, successor[i], instance) < min_cost_1) {	//best node
						node2add_index_2 = node2add_index_1;
						starting_node_2 = starting_node_1;
						min_cost_2 = min_cost_1;
						node2add_index_1 = j;
						starting_node_1 = i;
						min_cost_1 = lookup_cost(starting_node_1, nodes2visit[node2add_index_1], instance) + lookup_cost(nodes2visit[node2add_index_1], successor[starting_node_1], instance) - \
							lookup_cost(starting_node_1, successor[starting_node_1], instance);
					}
					else {	//second best node
						node2add_index_2 = j;
						starting_node_2 = i;
						min_cost_2 = lookup_cost(starting_node_2, nodes2visit[node2add_index_2], instance) + lookup_cost(nodes2visit[node2add_index_2], successor[starting_node_2], instance) - \
							lookup_cost(starting_node_2, successor[starting_node_2], instance);
					}
				}
			}
		}
		if (((double)rand() / (RAND_MAX + 1)) < instance->prob_ign_opt) {	//second best
			successor[nodes2visit[node2add_index_2]] = successor[starting_node_2];
			successor[starting_node_2] = nodes2visit[node2add_index_2];
			nodes2visit[node2add_index_2] = nodes2visit[--num_nodes2visit];

#if VERBOSE > 8
			{ printf("Substituted (%f, %f) -> (%f, %f) with (%f, %f) -> (%f, %f) -> (%f, %f)\n", instance->nodes[starting_node_2].x_coord, instance->nodes[starting_node_2].y_coord,
				instance->nodes[successor[successor[starting_node_2]]].x_coord, instance->nodes[successor[successor[starting_node_2]]].y_coord, instance->nodes[starting_node_2].x_coord,
				instance->nodes[starting_node_2].y_coord, instance->nodes[successor[starting_node_2]].x_coord, instance->nodes[successor[starting_node_2]].y_coord,
				instance->nodes[successor[successor[starting_node_2]]].x_coord, instance->nodes[successor[successor[starting_node_2]]].y_coord); }
#endif
		}
		else {	//best
			successor[nodes2visit[node2add_index_1]] = successor[starting_node_1];
			successor[starting_node_1] = nodes2visit[node2add_index_1];
			nodes2visit[node2add_index_1] = nodes2visit[--num_nodes2visit];

#if VERBOSE > 8
			{ printf("Substituted (%f, %f) -> (%f, %f) with (%f, %f) -> (%f, %f) -> (%f, %f)\n", instance->nodes[starting_node_1].x_coord, instance->nodes[starting_node_1].y_coord,
				instance->nodes[successor[successor[starting_node_1]]].x_coord, instance->nodes[successor[successor[starting_node_1]]].y_coord, instance->nodes[starting_node_1].x_coord,
				instance->nodes[starting_node_1].y_coord, instance->nodes[successor[starting_node_1]].x_coord, instance->nodes[successor[starting_node_1]].y_coord,
				instance->nodes[successor[successor[starting_node_1]]].x_coord, instance->nodes[successor[successor[starting_node_1]]].y_coord); }
#endif
		}
	}

	double current_cost = 0.0;
	unsigned int current_node = start_index;
	for (int i = 0; i < instance->num_nodes; i++) {
		current_cost += lookup_cost(current_node, successor[current_node], instance);

#if VERBOSE > 8
		{ printf("Cost of edge [(%f, %f), (%f, %f)]: %f; total path cost: %f\n", instance->nodes[current_node].x_coord, instance->nodes[current_node].y_coord,
			instance->nodes[successor[current_node]].x_coord, instance->nodes[successor[current_node]].y_coord, lookup_cost(current_node, successor[current_node], instance),
			current_cost); }
#endif
		current_node = successor[current_node];
	}

#if VERBOSE > 4
	{ printf("Cost: %f\n", current_cost); }
#endif

	if (current_cost < instance->best_sol_cost) {
		instance->best_sol_cost = current_cost;
		for (int i = 0; i < instance->num_nodes; i++) {
			instance->best_sol[i] = current_node;
			current_node = successor[current_node];

#if VERBOSE > 2
			{ printf("Element of index %d of the tour is the node of index %d [%f, %f]\n", i, instance->best_sol[i], instance->nodes[instance->best_sol[i]].x_coord,
				instance->nodes[instance->best_sol[i]].y_coord); }
#endif
		}
	}

	free(nodes2visit);
	free(successor);
}

int tsp_exm_sol(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Applying the extra-mileage heuristic...\n"); }
#endif

	if (instance->best_sol == NULL) {
		instance->best_sol = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));
		if (instance->best_sol == NULL) { fprintf(stderr, "could not allocate memory for the best solution array\n"); return 1; }
	}
	instance->best_sol_cost = DBL_INFY;

	int starting_index;
	if (instance->starting_index >= 0)
		starting_index = 0;
	else if (instance->starting_index == -1) {
		srand(instance->random_seed); for (int i = 0; i < MIN_RAND_RUNS + log(1 + instance->random_seed); i++) rand();
		int i = (int)(((double)rand() / RAND_MAX) * (instance->num_nodes - 1) + 0.5);
		int j = (int)(((double)rand() / RAND_MAX) * (instance->num_nodes - 1) + 0.5);
		starting_index = i * instance->num_nodes + j + 1; //the +1 is needed to make starting index go from 1 to n^2 (instead of 0 to n^2 - 1), to avoid collision with the default
	}
	else
		starting_index = -2;

#if VERBOSE > 2
	{ printf("Initial edge index: %d\n", starting_index); }
#endif

	if (starting_index == -2) {
		for (int i = 0; i < instance->num_nodes - 1; i++) {
			for (int j = i + 1; j < instance->num_nodes; j++) {
#if VERBOSE > 3
				{ printf("Starting pair (%d, %d) [(%f, %f), (%f, %f)]\n", i, j, instance->nodes[i].x_coord, instance->nodes[i].y_coord,
					instance->nodes[j].x_coord, instance->nodes[j].y_coord); }
#endif
				tsp_exm_sol_se(i, j, instance);
			}
		}
	}
	else if (starting_index == 0) {
		int start_max = 0, end_max = 0;
		for (int i = 0; i < instance->num_nodes - 1; i++) {
			for (int j = i + 1; j < instance->num_nodes; j++) {
				if (lookup_cost(i, j, instance) > lookup_cost(start_max, end_max, instance)) {
					start_max = i;
					end_max = j;
				}
			}
		}
#if VERBOSE > 3
		{ printf("Starting pair (%d, %d) [(%f, %f), (%f, %f)]\n", start_max, end_max, instance->nodes[start_max].x_coord, instance->nodes[start_max].y_coord,
			instance->nodes[end_max].x_coord, instance->nodes[end_max].y_coord); }
#endif
		tsp_exm_sol_se(start_max, end_max, instance);
	}
	else {
		int j = (starting_index - 1) % instance->num_nodes;
		int i = (starting_index - 1 - j) / instance->num_nodes;
#if VERBOSE > 3
		{ printf("Starting pair (%d, %d) [(%f, %f), (%f, %f)]\n", i, j, instance->nodes[i].x_coord, instance->nodes[i].y_coord,
			instance->nodes[j].x_coord, instance->nodes[j].y_coord); }
#endif
		tsp_exm_sol_se(i, j, instance);
	}

	return 0;
}

int ref_sol(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Refining the tsp solution...\n"); }
#endif

	assert(instance != NULL);
	assert(instance->best_sol != NULL);

	switch (instance->refine_flag) {
	case NO_REF:
		return 0;
	case TWO_OPT:
		return two_opt_ref(instance);
	default:
		return 1;
	}
}

int two_opt_ref(tsp_instance_t* instance) {

	srand(instance->random_seed); for (int i = 0; i < MIN_RAND_RUNS + log(1 + instance->random_seed); i++) rand();

#if VERBOSE > 1
	{ printf("Applying the 2-opt move refinement procedure...\n"); }
#endif

	if (instance->metaheur_flag == TABU) {

#if VERBOSE > 1
		{ printf("Applying the tabu search meta-heuristic...\n"); }
#endif

		instance->tabu_list = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));
		if (instance->tabu_list == NULL) { fprintf(stderr, "could not allocate memory for the tabu list\n"); return 1; }
		for (int i = 0; i < instance->num_nodes; i++) {
			instance->tabu_list[i] = INT_NINFY;
		}

		if (instance->min_tenure < 0) { instance->min_tenure = min(TABU_LIST_MIN, instance->num_nodes / TABU_LIST_MIN_DIV); }
		if (instance->max_tenure < 0) { instance->max_tenure = min(TABU_LIST_MAX, instance->num_nodes / TABU_LIST_MAX_DIV); }

#if VERBOSE > 2
		{ printf("Tenure range [%d, %d]\n", instance->min_tenure, instance->max_tenure); }
#endif
	}

	unsigned int current_iter = 0;
	unsigned int* current_sol = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));
	memcpy(current_sol, instance->best_sol, instance->num_nodes * sizeof(unsigned int));
	double current_sol_cost = instance->best_sol_cost;
	bool done = false;
	while (!done && instance->time_left > 0) {

		double start_time = seconds();

		unsigned int cross_node_a_index = -1;
		unsigned int cross_node_b_index = -1;
		double move_cost = DBL_INFY;

		for (int i = 0; i < instance->num_nodes - 2; i++) {
			for (int j = i + 2; j < instance->num_nodes; j++) {

				if ((lookup_cost(current_sol[i], current_sol[j], instance) + \
					lookup_cost(current_sol[(i + 1) % instance->num_nodes], current_sol[(j + 1) % instance->num_nodes], instance)) - \
					(lookup_cost(current_sol[i], current_sol[(i + 1) % instance->num_nodes], instance) + \
						lookup_cost(current_sol[j], current_sol[(j + 1) % instance->num_nodes], instance)) < move_cost) {	//best 2-opt move so far
					if (instance->metaheur_flag == TABU) {	//tabu search
						if (current_iter - instance->tabu_list[current_sol[i]] <= instance->min_tenure + (current_iter % (instance->max_tenure + 1 - instance->min_tenure))) {	//in tabu list
#if VERBOSE > 5
							{ printf("Node %d detected to be in the tabu list\n", current_sol[i]); }
#endif	
							if (current_sol_cost + move_cost < instance->best_sol_cost) {	//aspiration criterion (improving best sol)
								cross_node_a_index = i;
								cross_node_b_index = j;
								move_cost = (lookup_cost(current_sol[cross_node_a_index], current_sol[cross_node_b_index], instance) + \
									lookup_cost(current_sol[(cross_node_a_index + 1) % instance->num_nodes], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance)) - \
									(lookup_cost(current_sol[cross_node_a_index], current_sol[(cross_node_a_index + 1) % instance->num_nodes], instance) + \
										lookup_cost(current_sol[cross_node_b_index], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance));
#if VERBOSE > 5
								{ printf("Aspiration criterion! Best cost: %f; new cost: %f\n", instance->best_sol_cost, current_sol_cost + move_cost); }
#endif
							}
						}
						else {	//not in tabu list
#if VERBOSE > 5
							{ printf("Node %d not detected to be in the tabu list\n", current_sol[i]); }
#endif
							cross_node_a_index = i;
							cross_node_b_index = j;
							move_cost = (lookup_cost(current_sol[cross_node_a_index], current_sol[cross_node_b_index], instance) + \
								lookup_cost(current_sol[(cross_node_a_index + 1) % instance->num_nodes], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance)) - \
								(lookup_cost(current_sol[cross_node_a_index], current_sol[(cross_node_a_index + 1) % instance->num_nodes], instance) + \
									lookup_cost(current_sol[cross_node_b_index], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance));
						}
					}
					else {	//no tabu search
						cross_node_a_index = i;
						cross_node_b_index = j;
						move_cost = (lookup_cost(current_sol[cross_node_a_index], current_sol[cross_node_b_index], instance) + \
							lookup_cost(current_sol[(cross_node_a_index + 1) % instance->num_nodes], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance)) - \
							(lookup_cost(current_sol[cross_node_a_index], current_sol[(cross_node_a_index + 1) % instance->num_nodes], instance) + \
								lookup_cost(current_sol[cross_node_b_index], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance));
					}

#if VERBOSE > 8
					{ printf("The cost of going  %d -> %d is %f; cost of going %d -> %d is %f\n", current_sol[cross_node_a_index],
						current_sol[(cross_node_a_index + 1) % instance->num_nodes],
						lookup_cost(current_sol[cross_node_a_index], current_sol[(cross_node_a_index + 1) % instance->num_nodes], instance), current_sol[cross_node_b_index],
						current_sol[(cross_node_b_index + 1) % instance->num_nodes],
						lookup_cost(current_sol[cross_node_b_index], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance));
					printf("Substituted these edges with %d -> %d, cost %f and %d -> %d, cost %f\n", current_sol[cross_node_a_index], current_sol[cross_node_b_index],
						lookup_cost(current_sol[cross_node_a_index], current_sol[cross_node_b_index], instance), current_sol[(cross_node_a_index + 1) % instance->num_nodes],
						current_sol[(cross_node_b_index + 1) % instance->num_nodes],
						lookup_cost(current_sol[(cross_node_a_index + 1) % instance->num_nodes], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance)); }
#endif
				}
			}
		}

		if (cross_node_a_index < 0 || cross_node_b_index < 0) { fprintf(stderr, "2-opt found no feasible move: aborting\n"); return 1; }

#if VERBOSE > 2
		{ printf("Best move found: node %d, node %d, cost %f\n", current_sol[cross_node_a_index], current_sol[cross_node_b_index], move_cost); }
#endif

		if (move_cost >= 0 && instance->metaheur_flag == NO_MH) done = true;
		else if (move_cost >= 0 && instance->metaheur_flag == VNS) {

#if VERBOSE > 1
			{ printf("Applying the VNS meta-heuristic kick...\n"); }
#endif
			
			if (instance->num_nodes < 5) { fprintf(stderr, "Cannot apply VNS (5-kick) with less than 5 nodes\n"); return 1; }

			//indices of the 5-kick nodes
			int node_indices[] = {-1, -1, -1, -1, -1};
			int unitialized_indices = 5;

			while (unitialized_indices > 0) {
				int temp = (int)floor(((double)rand() / (RAND_MAX + 1)) * instance->num_nodes);	//[0, instance->num_nodes - 1]
				bool duplicate = false;
				for (int i = 0; i < 5 - unitialized_indices; i++)
					if (node_indices[i] == temp)
						duplicate = true;
				if (!duplicate)
					node_indices[5 - unitialized_indices--] = temp;
			}
			qsort(node_indices, 5, sizeof(int), cmp_int);

#if VERBOSE > 2
			{ printf("5-kick nodes(indices): %d(%d), %d(%d), %d(%d), %d(%d), %d(%d)\n", current_sol[node_indices[0]], node_indices[0], current_sol[node_indices[1]], node_indices[1],
				current_sol[node_indices[2]], node_indices[2], current_sol[node_indices[3]], node_indices[3], current_sol[node_indices[4]], node_indices[4]); }
#endif
			
			unsigned int* temp = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));
			memcpy(temp, current_sol, instance->num_nodes * sizeof(unsigned int));
			int tour_index = 0;
			for (int i = 0; i <= node_indices[0]; i++)	//0->a
				current_sol[tour_index++] = temp[i];
			for (int i = node_indices[1] + 1; i <= node_indices[2]; i++)	//b'->c
				current_sol[tour_index++] = temp[i];
			for (int i = node_indices[3] + 1; i <= node_indices[4]; i++)	//d'->e
				current_sol[tour_index++] = temp[i];
			for (int i = node_indices[0] + 1; i <= node_indices[1]; i++)	//a'->b
				current_sol[tour_index++] = temp[i];
			for (int i = node_indices[2] + 1; i <= node_indices[3]; i++)	//c'->d
				current_sol[tour_index++] = temp[i];
			for (int i = node_indices[4] + 1; i < instance->num_nodes; i++)	//e'->0
				current_sol[tour_index++] = temp[i];
			current_sol_cost -= lookup_cost(temp[node_indices[0]], temp[(node_indices[0] + 1) % instance->num_nodes], instance);	//removing a->a'
			current_sol_cost -= lookup_cost(temp[node_indices[1]], temp[(node_indices[1] + 1) % instance->num_nodes], instance);	//removing b->b'
			current_sol_cost -= lookup_cost(temp[node_indices[2]], temp[(node_indices[2] + 1) % instance->num_nodes], instance);	//removing c->c'
			current_sol_cost -= lookup_cost(temp[node_indices[3]], temp[(node_indices[3] + 1) % instance->num_nodes], instance);	//removing d->d'
			current_sol_cost -= lookup_cost(temp[node_indices[4]], temp[(node_indices[4] + 1) % instance->num_nodes], instance);	//removing e->e'
			current_sol_cost += lookup_cost(temp[node_indices[0]], temp[(node_indices[1] + 1) % instance->num_nodes], instance);	//adding   a->b'
			current_sol_cost += lookup_cost(temp[node_indices[2]], temp[(node_indices[3] + 1) % instance->num_nodes], instance);	//adding   c->d'
			current_sol_cost += lookup_cost(temp[node_indices[4]], temp[(node_indices[0] + 1) % instance->num_nodes], instance);	//adding   e->a'
			current_sol_cost += lookup_cost(temp[node_indices[1]], temp[(node_indices[2] + 1) % instance->num_nodes], instance);	//adding   b->c'
			current_sol_cost += lookup_cost(temp[node_indices[3]], temp[(node_indices[4] + 1) % instance->num_nodes], instance);	//adding   d->e'
			free(temp);

#if VERBOSE > 2
			{ printf("New tour has %d nodes and cost %f\n", tour_index, current_sol_cost); }
#endif

#if VERBOSE > 6
			{ for (int i = 0; i < instance->num_nodes; i++) printf("Element of index %d of the tour is the node of index %d [%f, %f]\n", i, current_sol[i],
					instance->nodes[current_sol[i]].x_coord, instance->nodes[current_sol[i]].y_coord); }
#endif

			if (current_sol_cost < instance->best_sol_cost) {
				memcpy(instance->best_sol, current_sol, instance->num_nodes * sizeof(unsigned int));
				instance->best_sol_cost = current_sol_cost;

#if VERBOSE > 2
				{ printf("Lucky kick (random kick improved best known sol)!\nSolution cost %f\n", instance->best_sol_cost); for (int i = 0; i < instance->num_nodes; i++)
					printf("Element of index %d of the tour is the node of index %d [%f, %f]\n", i, instance->best_sol[i],
						instance->nodes[instance->best_sol[i]].x_coord, instance->nodes[instance->best_sol[i]].y_coord); }
#endif
			}
		}
		else {
			//perform the 2-opt move (reverse current_sol from element of index cross_node_a_index + 1 to element of index cross_node_b_index)
			for (int i = 0; i < floor((cross_node_b_index - cross_node_a_index) / 2); i++) {	
				unsigned int temp = current_sol[cross_node_a_index + 1 + i];					
				current_sol[cross_node_a_index + 1 + i] = current_sol[cross_node_b_index - i];
				current_sol[cross_node_b_index - i] = temp;
			}
			current_sol_cost += move_cost;

#if VERBOSE > 4
			{ printf("New solution cost: %f\n", current_sol_cost); }
#endif

			if (instance->metaheur_flag == TABU) { 
				instance->tabu_list[current_sol[cross_node_a_index]] = current_iter; 

#if VERBOSE > 5
				{ printf("Added node %d to the tabu list at iteration %d\nTenure: %d\nTabu list:\n", current_sol[cross_node_a_index], current_iter,
					instance->min_tenure + (current_iter % (instance->max_tenure + 1 - instance->min_tenure)));
				for (int i = 0; i < instance->num_nodes; i++) 
					if (current_iter - instance->tabu_list[current_sol[i]] <= instance->min_tenure + (current_iter % (instance->max_tenure + 1 - instance->min_tenure))) 
						printf("Node %d, added at iteration %d\n", current_sol[i], instance->tabu_list[current_sol[i]]); printf("\n"); }
#endif
			}

			if (current_sol_cost < instance->best_sol_cost) {
				memcpy(instance->best_sol, current_sol, instance->num_nodes * sizeof(unsigned int));
				instance->best_sol_cost = current_sol_cost;

#if VERBOSE > 2
				{ printf("Solution cost %f\n", instance->best_sol_cost); for (int i = 0; i < instance->num_nodes; i++) 
					printf("Element of index %d of the tour is the node of index %d [%f, %f]\n", i, instance->best_sol[i],
						instance->nodes[instance->best_sol[i]].x_coord, instance->nodes[instance->best_sol[i]].y_coord); }
#endif
			}
		}
		double end_time = seconds();
		instance->time_left -= (end_time - start_time);

#if VERBOSE > 2
		{ if (instance->time_left <= 0) printf("Time limit reached during 2-opt refinement...\n"); }
#endif

		current_iter++;
	}

	if (instance->metaheur_flag == TABU) { free(instance->tabu_list); instance->tabu_list = NULL; }
	free(current_sol);

	return 0;
}

int plot_tour(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Plotting tour...\n"); }
#endif

	assert(instance != NULL);
	FILE* tour_plot_file = fopen(INSTANCE_TOUR_PLOT, "w");
	if (tour_plot_file == NULL) { fprintf(stderr, "cannot open the tour_plot_file\n"); return 1; }

	//save the tour in the file using the gnuplot format
	for (int i = 0; i < instance->num_nodes; i++) {
		fprintf(tour_plot_file, "%f %f\n", instance->nodes[instance->best_sol[i]].x_coord, instance->nodes[instance->best_sol[i]].y_coord);
	}
	fprintf(tour_plot_file, "%f %f\n", instance->nodes[instance->best_sol[0]].x_coord, instance->nodes[instance->best_sol[0]].y_coord);

	fclose(tour_plot_file);

	//execute the commands found in the commands.txt file. these commands will plot the points using gnuplot
	system(GNUPLOT_COMMAND_TOUR);

	return 0;
}

void free_tsp_instance(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Freeing memory...\n"); }
#endif

	assert(instance != NULL);

	assert(instance->nodes != NULL);
	free(instance->nodes);
	assert(instance->costs != NULL);
	free(instance->costs);
	assert(instance->best_sol != NULL);
	free(instance->best_sol);
}