#include "tsp_utils.h"
#include "tsp.h"
#include <string.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define MIN_RAND_RUNS 1000
#define INSTANCE_POINTS_PLOT "points.dat"
#define INSTANCE_TOUR_PLOT "tour.dat"
#define GNUPLOT_COMMAND_POINTS "gnuplot commands_points.txt"
#define GNUPLOT_COMMAND_TOUR "gnuplot commands_tour.txt"

#define DIST_INDEX(i, j, n) ((j) - 1 + (i) * ((n) - 1) - ((i) * ((i) + 1)) / 2)

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
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; } 													    // help
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; } 													// help
		help = 1;
	}

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
		printf("\nenter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}

	return help;
}

static int read_input_file(tsp_instance_t* instance) {

#if VERBOSE > 2
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

#if VERBOSE > 2
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
	srand(instance->random_seed); for (int i = 0; i < MIN_RAND_RUNS + log(1+instance->random_seed) ; i++) rand();
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

#if VERBOSE > 2
	{ printf("Preparing the instance nodes...\n"); }
#endif

	assert(instance != NULL);
	if (instance->input_file_name != NULL)
		return read_input_file(instance);
	else
		return generate_random_nodes(instance);
}

int plot_points(tsp_instance_t* instance) {

#if VERBOSE > 2
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

static double dist(int i, int j, tsp_instance_t* instance) {

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

#if VERBOSE > 2
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
	instance->costs = (double*)malloc(num_costs * sizeof(double));
	if (instance->costs == NULL) { fprintf(stderr, "could not allocate memory for the cost look-up table\n"); return 1; }

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

#if VERBOSE > 2
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

#if VERBOSE > 2
	{ printf("Solving the tsp instance...\n"); }
#endif

	assert(instance != NULL);

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

int tsp_seq_sol(tsp_instance_t* instance) {

#if VERBOSE > 2
	{ printf("Applying the sequential method...\n"); }
#endif

	assert(instance != NULL);

	instance->best_sol = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));
	if (instance->best_sol == NULL) { fprintf(stderr, "could not allocate memory for the best solution array\n"); return 1; }
	instance->best_sol_cost = 0.0;

	for (int i = 0; i < instance->num_nodes; i++) {
		instance->best_sol[i] = (unsigned int)i;
		instance->best_sol_cost += lookup_cost(i, (i + 1) % instance->num_nodes, instance);
	}

	return 0;
}

static void tsp_gdy_sol_si(unsigned int starting_index, tsp_instance_t* instance) {

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
	unsigned int next_node_index = ++current_tour_node_index;	//indicate that we are looking for the next tour position

	for (int i = current_tour_node_index; i < instance->num_nodes; i++) {	//iterate over the tour positions
#if VERBOSE > 8
		{ printf("Looking for element %d of the tour\n", current_tour_node_index); }
#endif
		double nearest_neighbor_dist = lookup_cost(nodes2visit[current_node_index], nodes2visit[next_node_index], instance);
		for (int j = current_tour_node_index + 1; j < instance->num_nodes; j++) {	//look for the node closest to the current node under exam.
			if (lookup_cost(nodes2visit[current_node_index], nodes2visit[j], instance) < nearest_neighbor_dist) {
				next_node_index = j;
				nearest_neighbor_dist = lookup_cost(nodes2visit[current_node_index], nodes2visit[next_node_index], instance);
			}
		}
#if VERBOSE > 8
		{ printf("Nearest neighbor %d of %d at distance %f\n", nodes2visit[next_node_index], nodes2visit[current_node_index], nearest_neighbor_dist); }
#endif
		current_node_index = next_node_index;
		swap_uint_array(current_node_index, current_tour_node_index, nodes2visit);
#if VERBOSE > 9
		{ printf("Swaped %d with %d. Now node2visit[%d] = %d and node2visit[%d] = %d\n", current_node_index, current_tour_node_index,
			current_node_index, nodes2visit[current_node_index], current_tour_node_index, nodes2visit[current_tour_node_index]); }
#endif
#if VERBOSE > 8
		{ printf("nodes2visit: [ "); for (int i = 0; i < instance->num_nodes; i++) printf("(%d, %d) ", i, nodes2visit[i]); printf("]\n"); }
#endif
		current_node_index = current_tour_node_index;
		next_node_index = ++current_tour_node_index;

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

#if VERBOSE > 2
	{ printf("Applying the greedy heuristic...\n"); }
#endif

	assert(instance != NULL);

	instance->best_sol = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));
	if (instance->best_sol == NULL) { fprintf(stderr, "could not allocate memory for the best solution array\n"); return 1; }
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

}

int tsp_exm_sol(tsp_instance_t* instance) {

#if VERBOSE > 2
	{ printf("Applying the extra-mileage heuristic...\n"); }
#endif

	assert(instance != NULL);

	instance->best_sol = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));
	if (instance->best_sol == NULL) { fprintf(stderr, "could not allocate memory for the best solution array\n"); return 1; }
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
			for (int j = i + 1; j < instance->nodes; j++) {
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

int plot_tour(tsp_instance_t* instance) {

#if VERBOSE > 2
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

#if VERBOSE > 2
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