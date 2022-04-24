#include "tsp_utils.h"
#include "tsp.h"

#define INSTANCE_POINTS_PLOT "points.dat"
#define INSTANCE_TOUR_PLOT "tour.dat"
#define GNUPLOT_COMMAND_POINTS "gnuplot commands_points.txt"
#define GNUPLOT_COMMAND_TOUR "gnuplot commands_tour.txt"

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
	instance->min_temperature = -1;
	instance->max_temperature = -1;
	instance->move_weight = 25;
	instance->pop_size = 100;

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
		if (strcmp(argv[i], "-min_temp") == 0) { instance->min_temperature = atoi(argv[++i]); continue; }				// the minimum temperature (simulated annealing)
		if (strcmp(argv[i], "-max_temp") == 0) { instance->max_temperature = atoi(argv[++i]); continue; }				// the maximum temperature (simulated annealing)
		if (strcmp(argv[i], "-move_weight") == 0) { instance->move_weight = atoi(argv[++i]); continue; }				// used to determine p. acceptance (simulated annealing)
		if (strcmp(argv[i], "-pop_size") == 0) { instance->pop_size = atoi(argv[++i]); continue; }						// population size (genetic algorithm)
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; } 													    // help
		if (strcmp(argv[i], "--help") == 0) { help = 1; continue; } 													// help
		help = 1;
	}

	//finalize parameter initialization
	instance->nodes = NULL;
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
		printf("-metaheur_flag %d (1: NO_MH; 2: TABU; 3: VNS; 4: SIM_ANNEAL; 5: GEN)\n", instance->metaheur_flag);
		printf("-min_tenure %d (the minumum size of the tenure (tabu search))\n", instance->min_tenure);
		printf("-max_tenure %d (the maximum size of the tenure (tabu search))\n", instance->max_tenure);
		printf("-min_temp %d (the minumum temperature (simulated annealing))\n", instance->min_temperature);
		printf("-max_temp %d (the maximum temperature (simulated annealing))\n", instance->max_temperature);
		printf("-move_weight %d (scaling factor used to determine the probability of acceptance (simulated annealing))\n", instance->move_weight);
		printf("-pop_size %d (population size (genetic algorithm))\n", instance->pop_size);
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
	printf("read %d nodes:\n", instance->num_nodes);
	for (int i = 0; i < instance->num_nodes; i++)
		printf("%d: (%f, %f)\n", i + 1, instance->nodes[i].x_coord, instance->nodes[i].y_coord);
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
	assert(0 <= j && j <= instance->num_nodes);

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

	if (instance->nodes != NULL) free(instance->nodes);
	if (instance->costs != NULL) free(instance->costs);
	if (instance->best_sol != NULL) free(instance->best_sol);
}