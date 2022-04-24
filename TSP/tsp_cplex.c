#include "tsp_utils.h"
#include "tsp.h"

#define LP_MODEL_FILENAME "model.lp"

#define INSTANCE_TOUR_PLOT "tour.dat"
#define GNUPLOT_COMMAND_TOUR "gnuplot commands_tour.txt"

static int xpos(int i, int j, tsp_instance_t* instance) {
	if (i == j) print_error("i == j in xpos()", 0);
	if (i > j) return xpos(j, i, instance);
	int pos = i * instance->num_nodes + j - ((i + 1) * (i + 2)) / 2;
	return pos;
}

void build_model(tsp_instance_t* instance, CPXENVptr env, CPXLPptr lp) {
	int error;
	//double zero = 0.0;//NOT USED, WHAT IS THIS FOR?
	char binary = 'B';

	char** cname = (char**)calloc(1, sizeof(char*));	// (char **) required by cplex...
	cname[0] = (char*)calloc(100, sizeof(char));

	//add binary var.s x(i,j) for i < j  
	for (int i = 0; i < instance->num_nodes; i++){
		for (int j = i + 1; j < instance->num_nodes; j++){
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);	// ... x(1,2), x(1,3) ....
			double obj = lookup_cost(i, j, instance);	// cost == distance   
			double lb = 0.0;
			double ub = 1.0;
			error = CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname);
			if (error) print_error("Wrong CPXnewcols() on x var.s", error);
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, instance)) print_error("Wrong position for x var.s (xpos error)", 0);
		}
	}

	// add the degree constraints 
	int* index = (int*)calloc(instance->num_nodes, sizeof(int));
	double* value = (double*)calloc(instance->num_nodes, sizeof(double));
	for (int h = 0; h < instance->num_nodes; h++) {	// add the degree constraint on node h
		double rhs = 2.0;
		char sense = 'E';	// 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h + 1);
		int num_non_zeros = 0;
		for (int i = 0; i < instance->num_nodes; i++) {
			if (i == h) continue;
			index[num_non_zeros] = xpos(i, h, instance);
			value[num_non_zeros] = 1.0;
			num_non_zeros++;
		}
		int izero = 0;
		error = CPXaddrows(env, lp, 0, 1, num_non_zeros, &rhs, &sense, &izero, index, value, NULL, &cname[0]);
		if (error) print_error("CPXaddrows() error", error);
	}
	free(value);
	free(index);

	free(cname[0]);
	free(cname);

#if VERBOSE > 0
	{ CPXwriteprob(env, lp, LP_MODEL_FILENAME, NULL); }
#endif
}

int plot_cycles(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Plotting cycles...\n"); }
#endif

	assert(instance != NULL);
	FILE* tour_plot_file = fopen(INSTANCE_TOUR_PLOT, "w");
	if (tour_plot_file == NULL) { fprintf(stderr, "cannot open the tour_plot_file\n"); return 1; }

	//save the cycles in the file using the gnuplot format
	int cycle_head = -1;
	for (int i = 0; i < instance->num_nodes + instance->num_cycles; i++) {
		if (instance->best_sol[i] != instance->cycle_delimiter) {
			if (cycle_head == -1) cycle_head = instance->best_sol[i];
			fprintf(tour_plot_file, "%f %f\n", instance->nodes[instance->best_sol[i]].x_coord, instance->nodes[instance->best_sol[i]].y_coord);
		}
		else {
			if (cycle_head != -1) fprintf(tour_plot_file, "%f %f\n\n", instance->nodes[cycle_head].x_coord, instance->nodes[cycle_head].y_coord);
			cycle_head = -1;
		}	
	}

	fclose(tour_plot_file);

	//execute the commands found in the commands.txt file. these commands will plot the points using gnuplot
	system(GNUPLOT_COMMAND_TOUR);

	return 0;
}

static void build_model_ds(const tsp_instance_t* instance, const double* xstar, int* succ, int* comp, bool* visited_nodes, int* num_cycles) {
	for (int i = 0; i < instance->num_nodes; i++)
		succ[i] = -1;
	for (int i = 0; i < instance->num_nodes; i++)
		comp[i] = -1;
	for (int i = 0; i < instance->num_nodes; i++)
		visited_nodes[i] = false;
	*num_cycles = 0;
	for (int i = 0; i < instance->num_nodes; i++) {
		int current_node = i;
		while (!visited_nodes[current_node]) {
			visited_nodes[current_node] = true;
			for (int j = 0; j < instance->num_nodes; j++) {
				if ((j != current_node) && (xstar[xpos(current_node, j, instance)] > 0.5) && (!visited_nodes[j])) {
					succ[current_node] = j;
					comp[current_node] = *num_cycles + 1;
					current_node = j;
					break;
				}
			}
		}
		if (succ[current_node] == -1) {
			succ[current_node] = i;
			comp[current_node] = ++(*num_cycles);
		}
	}
}

static void print_cycles(const tsp_instance_t* instance, const int* succ, const int* comp, const int num_cycles) {
	int next_node2visit;
	bool* visited_nodes = (bool*)malloc(instance->num_nodes * sizeof(bool));
	if (visited_nodes == NULL) print_error("print_cycles() could not allocate memory for visited_nodes[]", 0);
	for (int i = 0; i < instance->num_nodes; i++)
		visited_nodes[i] = false;

	printf("Solution cycles:\n");
	for (int i = 0; i < instance->num_nodes; i++) {
		next_node2visit = succ[i];
		if (!visited_nodes[i]) {
			visited_nodes[i] = true;
			printf("(%d)%d", comp[i], i + 1);
			while (!visited_nodes[next_node2visit]) {
				visited_nodes[next_node2visit] = true;
				printf("->(%d)%d", comp[next_node2visit], next_node2visit + 1);
				next_node2visit = succ[next_node2visit];
			}
		}
		if (next_node2visit == i) {
			visited_nodes[next_node2visit] = true;
			printf("->(%d)%d\n", comp[next_node2visit], next_node2visit + 1);
		}
	}
	printf("Number of cycles: %d\n", num_cycles);

	free(visited_nodes);
}

static void update_tsp_best_sol(tsp_instance_t* instance, const int* succ, const unsigned int num_cycles, const unsigned int cycle_delimiter) {
	unsigned int* current_sol = (unsigned int*)malloc((instance->num_nodes + num_cycles) * sizeof(unsigned int));
	if (current_sol == NULL) print_error("update_tsp_best_sol() could not allocate memory for current_sol[]", 0);
	bool* visited_nodes = (bool*)malloc(instance->num_nodes * sizeof(bool));
	if (visited_nodes == NULL) print_error("update_tsp_best_sol() could not allocate memory for visited_nodes[]", 0);
	for (int i = 0; i < instance->num_nodes; i++)
		visited_nodes[i] = false;
	double sol_cost = 0.0;
	int best_sol_node_index = 0;
	int next_node2visit = 0;
	for (int i = 0; i < instance->num_nodes; i++) {
		next_node2visit = succ[i];
		if (!visited_nodes[i]) {
			visited_nodes[i] = true;
			current_sol[best_sol_node_index++] = i;
			//printf("%d starter of a new cycle\n", i + 1);
			while (!visited_nodes[next_node2visit]) {
				visited_nodes[next_node2visit] = true;
				current_sol[best_sol_node_index++] = next_node2visit;
				//printf("%d added to the current cycle\n", next_node2visit + 1);
				sol_cost += lookup_cost(current_sol[best_sol_node_index - 2], current_sol[best_sol_node_index - 1], instance);
#if VERBOSE > 9
				{ printf("Cost of going from node %d to node %d is %f\n", current_sol[best_sol_node_index - 2] + 1, current_sol[best_sol_node_index - 1] + 1,
					lookup_cost(current_sol[best_sol_node_index - 2], current_sol[best_sol_node_index - 1], instance)); }
#endif
				next_node2visit = succ[next_node2visit];
			}
		}
		if (next_node2visit == i) {
			visited_nodes[next_node2visit] = true;
			current_sol[best_sol_node_index++] = cycle_delimiter;
			//printf("cycle node sequence finished\n");
			sol_cost += lookup_cost(current_sol[best_sol_node_index - 2], i, instance);
#if VERBOSE > 9
			{ printf("Cost of going from node %d to node %d is %f\n", current_sol[best_sol_node_index - 2] + 1, i + 1,
				lookup_cost(current_sol[best_sol_node_index - 2], i, instance)); }
#endif
		}
	}
	if (best_sol_node_index != (instance->num_nodes + num_cycles)) print_error("update_tsp_best_sol()'s best_sol[] size and nodes do not match", 0);
	if (true) {
		instance->best_sol_cost = sol_cost;
		if (instance->best_sol != NULL) {
			free(instance->best_sol);
			instance->best_sol = NULL;
		}
		instance->best_sol = (unsigned int*)malloc((instance->num_nodes + num_cycles) * sizeof(unsigned int));
		if (instance->best_sol == NULL) print_error("update_tsp_best_sol() could not allocate memory for instance->best_sol[]", 0);
		memcpy(instance->best_sol, current_sol, (instance->num_nodes + num_cycles) * sizeof(unsigned int));
		instance->num_cycles = num_cycles;
		instance->cycle_delimiter = cycle_delimiter;
	}
	free(current_sol);
	free(visited_nodes);

#if VERBOSE > 1
	printf("best_sol[] after update_tsp_best_sol():\n");
	for (int i = 0; i < instance->num_nodes + num_cycles; i++) {
		if (instance->best_sol[i] != cycle_delimiter)
			printf(" %d ", instance->best_sol[i] + 1);
		else
			printf(" - ");
	}
	printf("\nSolution cost: %f\n", instance->best_sol_cost);
#endif
}

void benders_method(tsp_instance_t* instance, CPXENVptr env, CPXLPptr lp, double** xstar, int* succ, int* comp, bool* visited_nodes, int* num_cycles) {
	int error;
	int iteration = 0;
	while (*num_cycles > 1) {
		iteration++;
		int num_constraints;
		if (*num_cycles == 2)
			num_constraints = 1;
		else
			num_constraints = *num_cycles;
		int ncols = CPXgetnumcols(env, lp);
		int* index = (int*)calloc(ncols, sizeof(int));
		double* coeff = (double*)calloc(ncols, sizeof(double));
		char** cname = (char**)calloc(1, sizeof(char*));
		cname[0] = (char*)calloc(100, sizeof(char));
		for (int k = 1; k <= num_constraints; k++) {
			int num_non_zeros = 0;
			double rhs = -1.0;
			char sense = 'L';
			int izero = 0;
			sprintf(cname[0], "iteration(%d).component(%d)", iteration, k);
			for (int i = 0; i < instance->num_nodes; i++) {
				if (comp[i] == k) {
					rhs += 1;
					for (int j = i + 1; j < instance->num_nodes; j++) {
						if (comp[j] == k) {
							index[num_non_zeros] = xpos(i, j, instance);
							coeff[num_non_zeros] = 1.0;
							num_non_zeros++;
						}
					}
				}
			}
			error = CPXaddrows(env, lp, 0, 1, num_non_zeros, &rhs, &sense, &izero, index, coeff, NULL, &cname[0]);
			if (error) print_error("CPXaddrows() error", error);
		}
#if VERBOSE > 0
		{ CPXwriteprob(env, lp, LP_MODEL_FILENAME, NULL); }
#endif
		free(index);
		free(coeff);
		free(cname[0]);
		free(cname);
		error = CPXmipopt(env, lp);
		if (error) print_error("CPXmipopt() error", error);
		free(*xstar);
		*xstar = (double*)calloc(ncols, sizeof(double));
		error = CPXgetx(env, lp, *xstar, 0, ncols - 1);
		if (error) print_error("CPXgetx() error", error);
		build_model_ds(instance, *xstar, succ, comp, visited_nodes, num_cycles);
		printf("---New solution found has %d cycles---\n", *num_cycles);
		update_tsp_best_sol(instance, succ, *num_cycles, instance->num_nodes);
	}
}

int TSPopt(tsp_instance_t* instance) {
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error", error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model");
	if (error) print_error("CPXcreateprob() error", error);

	build_model(instance, env, lp);

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
#if VERBOSE > 1
	{ CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); }	// Cplex output on screen
#endif
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, instance->random_seed);
	CPXsetdblparam(env, CPX_PARAM_TILIM, instance->time_limit);

	error = CPXmipopt(env, lp);
	if (error) print_error("CPXmipopt() error", error);

	//optimal solution
	int ncols = CPXgetnumcols(env, lp);
	double* xstar = (double*)calloc(ncols, sizeof(double));
	error = CPXgetx(env, lp, xstar, 0, ncols - 1);
	if (error) print_error("CPXgetx() error", error);

#if VERBOSE > 4	//print the solution edges
	for (int i = 0; i < instance->num_nodes; i++)
		for (int j = i + 1; j < instance->num_nodes; j++)
			if (xstar[xpos(i, j, instance)] > 0.5)
				printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1);
#endif
	
	//compute the successor list, component list, visited nodes list and number of cycles
	int* succ = (int*)malloc(instance->num_nodes * sizeof(int));
	if (succ == NULL) print_error("Could not allocate succ[]", 0);
	int* comp = (int*)malloc(instance->num_nodes * sizeof(int));
	if (comp == NULL) print_error("Could not allocate comp[]", 0);
	bool* visited_nodes = (bool*)malloc(instance->num_nodes * sizeof(bool));
	if (visited_nodes == NULL) print_error("Could not allocate visited_nodes[]", 0);
	int num_cycles;
	build_model_ds(instance, xstar, succ, comp, visited_nodes, &num_cycles);
	
	//print the found cycles
	print_cycles(instance, succ, comp, num_cycles);
	
	//save the cycles in the instance best sol 
	update_tsp_best_sol(instance, succ, num_cycles, instance->num_nodes);
	
	//plot the best sol cycles
	error = plot_cycles(instance);
	if (error) print_error("plot_cycles() error", 0);

	//find the optimal feasible solution using benders' method
	benders_method(instance, env, lp, &xstar, succ, comp, visited_nodes, &num_cycles);
	print_cycles(instance, succ, comp, num_cycles);
	error = plot_cycles(instance);
	if (error) print_error("plot_cycles() error", 0);
	
	free(xstar);
	free(succ);
	free(comp);
	free(visited_nodes);

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return 0;
}