#include "tsp_utils.h"
#include "tsp.h"

#define LP_MODEL_FILENAME "model.lp"

#define INSTANCE_TOUR_PLOT "tour.dat"
#define GNUPLOT_COMMAND_TOUR "gnuplot commands_tour.txt"

#define EPSILON_PERCENT (0.0001)	/*0.01%*/

#define NODE_LIM (100)

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

static void get_local_sol(const tsp_instance_t* instance, const int* succ, const unsigned int num_cycles, const unsigned int cycle_delimiter, 
	unsigned int* current_sol, double* sol_cost) {
	bool* visited_nodes = (bool*)malloc(instance->num_nodes * sizeof(bool)); if (visited_nodes == NULL) print_error("refine_and_update_tsp_best_sol() could not allocate memory for visited_nodes[]", 0);
	for (int i = 0; i < instance->num_nodes; i++)
		visited_nodes[i] = false;

	//makes current_sol[] contain the list of nodes of the tsp solution (sequence of cycles)
	//makes sol_cost contain the cost of the current_sol[] (sum of the cycles' costs)
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
				*sol_cost += lookup_cost(current_sol[best_sol_node_index - 2], current_sol[best_sol_node_index - 1], instance);
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
			*sol_cost += lookup_cost(current_sol[best_sol_node_index - 2], i, instance);
#if VERBOSE > 9
			{ printf("Cost of going from node %d to node %d is %f\n", current_sol[best_sol_node_index - 2] + 1, i + 1,
				lookup_cost(current_sol[best_sol_node_index - 2], i, instance)); }
#endif
		}
	}
	if (best_sol_node_index != (instance->num_nodes + num_cycles)) print_error("get_local_sol()'s best_sol[] size and nodes do not match", 0);

	free(visited_nodes);
}

static void get_new_opt(tsp_instance_t* instance, const double new_sol_cost, unsigned int* new_sol, const unsigned int num_cycles, const unsigned int cycle_delimiter) {
	instance->best_sol_cost = new_sol_cost;
	if (instance->best_sol != NULL) {
		free(instance->best_sol);
		instance->best_sol = NULL;
	}
	instance->best_sol = (unsigned int*)malloc((instance->num_nodes + num_cycles) * sizeof(unsigned int));
	if (instance->best_sol == NULL) print_error("get_new_opt() could not allocate memory for instance->best_sol[]", 0);
	memcpy(instance->best_sol, new_sol, (instance->num_nodes + num_cycles) * sizeof(unsigned int));
	instance->num_cycles = num_cycles;
	instance->cycle_delimiter = cycle_delimiter;
}

static void print_best_sol(const tsp_instance_t* instance) {
	printf("best_sol[]:\n");
	for (int i = 0; i < instance->num_nodes + instance->num_cycles; i++) {
		if (instance->best_sol[i] != instance->cycle_delimiter)
			printf(" %d ", instance->best_sol[i] + 1);
		else
			printf(" - ");
	}
	printf("\nSolution cost: %f\n", instance->best_sol_cost);
}

static void update_tsp_best_sol(tsp_instance_t* instance, const int* succ, const unsigned int num_cycles, const unsigned int cycle_delimiter) {
	unsigned int* current_sol = (unsigned int*)malloc((instance->num_nodes + num_cycles) * sizeof(unsigned int));
	if (current_sol == NULL) print_error("update_tsp_best_sol() could not allocate memory for current_sol[]", 0);
	double sol_cost = 0.0;
	get_local_sol(instance, succ, num_cycles, cycle_delimiter, current_sol, &sol_cost);

	//TO BE DELETED
	/*
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
	*/

	if (sol_cost < instance->best_sol_cost) {
		get_new_opt(instance, sol_cost, current_sol, num_cycles, cycle_delimiter);

		//TO BE DELETED
		/*
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
		*/
	}
	free(current_sol);

	//TO BE DELTED
	/*
	free(visited_nodes);
	*/

#if VERBOSE > 1
	print_best_sol(instance);

	//TO BE DELETED
	/*
	printf("best_sol[] after update_tsp_best_sol():\n");
	for (int i = 0; i < instance->num_nodes + num_cycles; i++) {
		if (instance->best_sol[i] != cycle_delimiter)
			printf(" %d ", instance->best_sol[i] + 1);
		else
			printf(" - ");
	}
	printf("\nSolution cost: %f\n", instance->best_sol_cost);
	*/
#endif
}

static void refine_and_update_tsp_best_sol(tsp_instance_t* instance, const int* succ, const unsigned int num_cycles, const unsigned int cycle_delimiter) {
	if (num_cycles != 1) print_error("refine_and_update_tsp_best_sol() num_cycles != 1", num_cycles);

	unsigned int* current_sol = (unsigned int*)malloc((instance->num_nodes + num_cycles) * sizeof(unsigned int));
	if (current_sol == NULL) print_error("refine_and_update_tsp_best_sol() could not allocate memory for current_sol[]", 0);
	double sol_cost = 0.0;
	get_local_sol(instance, succ, num_cycles, cycle_delimiter, current_sol, &sol_cost);

	//TO BE DELETED
	/*
	bool* visited_nodes = (bool*)malloc(instance->num_nodes * sizeof(bool)); if (visited_nodes == NULL) print_error("refine_and_update_tsp_best_sol() could not allocate memory for visited_nodes[]", 0);
	for (int i = 0; i < instance->num_nodes; i++)
		visited_nodes[i] = false;

	//makes current_sol[] contain the list of nodes of the tsp solution (sequence of cycles)
	//makes sol_cost contain the cost of the current_sol[] (sum of the cycles' costs)
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
	if (best_sol_node_index != (instance->num_nodes + num_cycles)) print_error("refine_and_update_tsp_best_sol()'s best_sol[] size and nodes do not match", 0);
	*/

	//dummy instance to be refined
	tsp_instance_t curr_instance;
	memcpy(&curr_instance, instance, sizeof(tsp_instance_t));
	curr_instance.best_sol = (unsigned int*)malloc((instance->num_nodes + num_cycles) * sizeof(unsigned int)); if (curr_instance.best_sol == NULL) print_error("refine_and_update_tsp_best_sol() could not allocate memory for curr_instance.best_sol[]", 0);
	memcpy(curr_instance.best_sol, current_sol, (instance->num_nodes + num_cycles) * sizeof(unsigned int));
	curr_instance.best_sol_cost = sol_cost;
	curr_instance.num_cycles = num_cycles;
	curr_instance.cycle_delimiter = cycle_delimiter;
	curr_instance.metaheur_flag = NO_MH;
	two_opt_ref(&curr_instance);

	//update the best solution
	if (curr_instance.best_sol_cost < instance->best_sol_cost) {
		get_new_opt(instance, curr_instance.best_sol_cost, curr_instance.best_sol, curr_instance.num_cycles, curr_instance.cycle_delimiter);

		//TO BE DELETED
		/*
		instance->best_sol_cost = curr_instance.best_sol_cost;
		if (instance->best_sol != NULL) {
			free(instance->best_sol);
			instance->best_sol = NULL;
		}
		instance->best_sol = (unsigned int*)malloc((instance->num_nodes + num_cycles) * sizeof(unsigned int));
		if (instance->best_sol == NULL) print_error("refine_and_update_tsp_best_sol() could not allocate memory for instance->best_sol[]", 0);
		memcpy(instance->best_sol, curr_instance.best_sol, (instance->num_nodes + num_cycles) * sizeof(unsigned int));
		instance->num_cycles = num_cycles;
		instance->cycle_delimiter = cycle_delimiter;
		*/
	}
	free(current_sol);

	//TO BE DELETED
	/*
	free(visited_nodes);
	*/

	free(curr_instance.best_sol);

#if VERBOSE > 1
	print_best_sol(instance);

	//TO BE DELETED
	/*
	printf("best_sol[] after refine_and_update_tsp_best_sol():\n");
	for (int i = 0; i < instance->num_nodes + num_cycles; i++) {
		if (instance->best_sol[i] != cycle_delimiter)
			printf(" %d ", instance->best_sol[i] + 1);
		else
			printf(" - ");
	}
	printf("\nSolution cost: %f\n", instance->best_sol_cost);
	*/
#endif
}

//TO BE DELETED
/*
static void patch_and_update_tsp_best_sol(tsp_instance_t* instance, const int* succ, const int* comp, int num_cycles, const unsigned int cycle_delimiter, const bool add_SECs) {
	printf("Patching the cycles...\n");
	int* local_succ = (int*)malloc(instance->num_nodes * sizeof(int)); if (local_succ == NULL) print_error("patch_and_update_tsp_best_sol() could not allocate memory for local_succ[]", 0);
	memcpy(local_succ, succ, instance->num_nodes * sizeof(int));
	int* local_comp = (int*)malloc(instance->num_nodes * sizeof(int)); if (local_comp == NULL) print_error("patch_and_update_tsp_best_sol() could not allocate memory for local_comp[]", 0);
	memcpy(local_comp, comp, instance->num_nodes * sizeof(int));

	while (num_cycles > 1) {
		int node_a = -1;
		int node_b = -1;
		double move_cost = DBL_INFY;
		for (int i = 0; i < instance->num_nodes; i++) {
			for (int j = 0; j < instance->num_nodes; j++) {
				if (local_comp[i] < local_comp[j]) {
					if ((lookup_cost(i, local_succ[j], instance) + lookup_cost(j, local_succ[i], instance)) -
						(lookup_cost(i, local_succ[i], instance) + lookup_cost(j, local_succ[j], instance)) < move_cost) {
						node_a = i;
						node_b = j;
						move_cost = (lookup_cost(i, local_succ[j], instance) + lookup_cost(j, local_succ[i], instance)) -
							(lookup_cost(i, local_succ[i], instance) + lookup_cost(j, local_succ[j], instance));
					}
				}
			}
		}
		if (node_a < 0 || node_a >(instance->num_nodes)) printf("\n\nnode_a %d is out of bound\n\n\n", node_a);
		if (node_b < 0 || node_b >(instance->num_nodes)) printf("\n\nnode_b %d is out of bound\n\n\n", node_b);
#if VERBOSE > 4
		{ printf("Merging components: node_a %d (comp %d), node_b %d (comp %d), cost %f\n", node_a + 1, local_comp[node_a], node_b + 1, local_comp[node_b], move_cost); }
#endif
		int node_a_succ = local_succ[node_a];
		int node_b_succ = local_succ[node_b];
#if VERBOSE > 7
		{ printf("Before update: local_succ[%d] = %d, local_succ[%d] = %d\n", node_a + 1, local_succ[node_a] + 1, node_b + 1, local_succ[node_b] + 1); }
#endif
		local_succ[node_a] = node_b_succ;
		local_succ[node_b] = node_a_succ;
#if VERBOSE > 7
		{ printf("After update: local_succ[%d] = %d, local_succ[%d] = %d\n", node_a + 1, local_succ[node_a] + 1, node_b + 1, local_succ[node_b] + 1); }
#endif
		for (int i = 0; i < instance->num_nodes; i++)
			if ((i != node_b) && (local_comp[i] == local_comp[node_b])) {
#if VERBOSE > 7
				{ printf("local_comp[%d] went from %d to ", i + 1, local_comp[i]); }
#endif
				local_comp[i] = local_comp[node_a];
#if VERBOSE > 7
				{ printf("%d\n", local_comp[i]); }
#endif
			}
#if VERBOSE > 7
		{ printf("local_comp[%d] went from %d to ", node_b + 1, local_comp[node_b]); }
#endif
		local_comp[node_b] = local_comp[node_a];
#if VERBOSE > 7
		{ printf("%d\n", local_comp[node_b]); }
#endif
		num_cycles--;
#if VERBOSE > 4
		print_cycles(instance, local_succ, local_comp, num_cycles);
#endif
	}

	refine_and_update_tsp_best_sol(instance, local_succ, num_cycles, cycle_delimiter);

	free(local_succ);
	free(local_comp);
}
*/

static void add_sec(const tsp_instance_t* instance, CPXENVptr env, CPXLPptr lp, CPXCALLBACKCONTEXTptr context, int* index, double* coeff, const int* comp, char* cname,
	int iteration, int component) {
	int error;
	int num_non_zeros = 0;
	double rhs = -1.0;
	char sense = 'L';
	int izero = 0;
	sprintf(cname, "iteration(%d).component(%d)", iteration, component);
	for (int i = 0; i < instance->num_nodes; i++) {
		if (comp[i] == component) {
			rhs += 1;
			for (int j = i + 1; j < instance->num_nodes; j++) {
				if (comp[j] == component) {
					index[num_non_zeros] = xpos(i, j, instance);
					coeff[num_non_zeros] = 1.0;
					num_non_zeros++;
				}
			}
		}
	}
	if (env != NULL && lp != NULL && context == NULL) {
		error = CPXaddrows(env, lp, 0, 1, num_non_zeros, &rhs, &sense, &izero, index, coeff, NULL, &cname);
		if (error) print_error("CPXaddrows() error", error);
	}
	else if (env == NULL && lp == NULL && context != NULL) {
		error = CPXcallbackrejectcandidate(context, 1, num_non_zeros, &rhs, &sense, &izero, index, coeff);	//reject the solution and adds one cut 
		if (error) print_error("CPXcallbackrejectcandidate() error", 0);
	}
	else
		print_error("Could not resolve how to add SEC", 0);
}

static void add_single_cycle_secs(const tsp_instance_t* instance, CPXENVptr env, CPXLPptr lp, CPXCALLBACKCONTEXTptr context, const int* comp,
	int iteration, int num_cycles, int sol_cardinality) {
	int num_constraints;
	if (num_cycles < 2)
		num_constraints = 0;
	else if (num_cycles == 2)
		num_constraints = 1;
	else
		num_constraints = num_cycles;
	int* index = (int*)calloc(sol_cardinality, sizeof(int));
	double* coeff = (double*)calloc(sol_cardinality, sizeof(double));
	char* cname = (char*)calloc(100, sizeof(char));
	for (int k = 1; k <= num_constraints; k++)
		add_sec(instance, env, lp, context, index, coeff, comp, cname, iteration, k);
	free(index);
	free(coeff);
	free(cname);
}

static void patch_and_update_tsp_best_sol(tsp_instance_t* instance, CPXENVptr env, CPXLPptr lp, int ncols, const int* succ,
	const int* comp, int num_cycles, const unsigned int cycle_delimiter, int iteration, const bool add_SECs) {
	//used to maintain the local (i.e. pached) solution
	int* local_succ = (int*)malloc(instance->num_nodes * sizeof(int)); if (local_succ == NULL) print_error("add_sec_and_patch_update_tsp_best_sol() could not allocate memory for local_succ[]", 0);
	memcpy(local_succ, succ, instance->num_nodes * sizeof(int));
	int* local_comp = (int*)malloc(instance->num_nodes * sizeof(int)); if (local_comp == NULL) print_error("add_sec_and_patch_update_tsp_best_sol() could not allocate memory for local_comp[]", 0);
	memcpy(local_comp, comp, instance->num_nodes * sizeof(int));

	//patch until only one cycle left, if add_SECs each patch adds a new SEC
	int* index = (int*)calloc(ncols, sizeof(int));
	double* coeff = (double*)calloc(ncols, sizeof(double));
	char* cname = (char*)calloc(100, sizeof(char));
	while (num_cycles > 1) {
		int node_a = -1;
		int node_b = -1;
		double move_cost = DBL_INFY;
		for (int i = 0; i < instance->num_nodes; i++) {
			for (int j = 0; j < instance->num_nodes; j++) {
				if (local_comp[i] < local_comp[j]) {
					if ((lookup_cost(i, local_succ[j], instance) + lookup_cost(j, local_succ[i], instance)) -
						(lookup_cost(i, local_succ[i], instance) + lookup_cost(j, local_succ[j], instance)) < move_cost) {
						node_a = i;
						node_b = j;
						move_cost = (lookup_cost(i, local_succ[j], instance) + lookup_cost(j, local_succ[i], instance)) -
							(lookup_cost(i, local_succ[i], instance) + lookup_cost(j, local_succ[j], instance));
					}
				}
			}
		}
		if (node_a < 0 || node_a >(instance->num_nodes - 1)) printf("\n\nnode_a %d is out of bound\n\n\n", node_a);
		if (node_b < 0 || node_b >(instance->num_nodes - 1)) printf("\n\nnode_b %d is out of bound\n\n\n", node_b);
#if VERBOSE > 4
		{ printf("Merging components: node_a %d (comp %d), node_b %d (comp %d), cost %f\n", node_a + 1, local_comp[node_a], node_b + 1, local_comp[node_b], move_cost); }
#endif
		int node_a_succ = local_succ[node_a];
		int node_b_succ = local_succ[node_b];
#if VERBOSE > 7
		{ printf("Before update: local_succ[%d] = %d, local_succ[%d] = %d\n", node_a + 1, local_succ[node_a] + 1, node_b + 1, local_succ[node_b] + 1); }
#endif
		local_succ[node_a] = node_b_succ;
		local_succ[node_b] = node_a_succ;
#if VERBOSE > 7
		{ printf("After update: local_succ[%d] = %d, local_succ[%d] = %d\n", node_a + 1, local_succ[node_a] + 1, node_b + 1, local_succ[node_b] + 1); }
#endif
		for (int i = 0; i < instance->num_nodes; i++)
			if ((i != node_b) && (local_comp[i] == local_comp[node_b])) {
#if VERBOSE > 7
				{ printf("local_comp[%d] went from %d to ", i + 1, local_comp[i]); }
#endif
				local_comp[i] = local_comp[node_a];
#if VERBOSE > 7
				{ printf("%d\n", local_comp[i]); }
#endif
			}
#if VERBOSE > 7
		{ printf("local_comp[%d] went from %d to ", node_b + 1, local_comp[node_b]); }
#endif
		local_comp[node_b] = local_comp[node_a];
#if VERBOSE > 7
		{ printf("%d\n", local_comp[node_b]); }
#endif
		num_cycles--;
#if VERBOSE > 4
		print_cycles(instance, local_succ, local_comp, num_cycles);
#endif
		if (add_SECs) {
			add_sec(instance, env, lp, NULL, index, coeff, comp, cname, iteration, local_comp[node_a]);
#if VERBOSE > 0
			{ CPXwriteprob(env, lp, LP_MODEL_FILENAME, NULL); }
#endif
		}
	}
	free(index);
	free(coeff);
	free(cname);

	refine_and_update_tsp_best_sol(instance, local_succ, num_cycles, cycle_delimiter);

	free(local_succ);
	free(local_comp);
}

static void add_sec_and_patch_update_tsp_best_sol(tsp_instance_t* instance, CPXENVptr env, CPXLPptr lp, int ncols, const int* succ, 
	const int* comp, int num_cycles, const unsigned int cycle_delimiter, int iteration) {
	printf("Patching the cycles, adding SEC constraints and updating best solution...\n");
	
	//adding the single cycle SECs
	if (num_cycles > 1)
		add_single_cycle_secs(instance, env, lp, NULL, comp, iteration, num_cycles, ncols);
	
	//TO BE DELETED
	/*
	int num_constraints;
	if (num_cycles < 2)
		num_constraints = 0;
	else if (num_cycles == 2)
		num_constraints = 1;
	else
		num_constraints = num_cycles;
	int* index = (int*)calloc(ncols, sizeof(int));
	double* coeff = (double*)calloc(ncols, sizeof(double));
	char* cname = (char*)calloc(100, sizeof(char));
	for (int k = 1; k <= num_constraints; k++)
		add_sec(instance, env, lp, NULL, index, coeff, comp, cname, iteration, k);
	*/

	patch_and_update_tsp_best_sol(instance, env, lp, ncols,succ, comp, num_cycles, cycle_delimiter, iteration, true);

	//TO BE DELETED
	/*
	//used to maintain the local (i.e. pached) solution
	int* local_succ = (int*)malloc(instance->num_nodes * sizeof(int)); if (local_succ == NULL) print_error("add_sec_and_patch_update_tsp_best_sol() could not allocate memory for local_succ[]", 0);
	memcpy(local_succ, succ, instance->num_nodes * sizeof(int));
	int* local_comp = (int*)malloc(instance->num_nodes * sizeof(int)); if (local_comp == NULL) print_error("add_sec_and_patch_update_tsp_best_sol() could not allocate memory for local_comp[]", 0);
	memcpy(local_comp, comp, instance->num_nodes * sizeof(int));

	//patch until only one cycle left, each patch adds a new SEC
	int* index = (int*)calloc(ncols, sizeof(int));
	double* coeff = (double*)calloc(ncols, sizeof(double));
	char* cname = (char*)calloc(100, sizeof(char));
	while (num_cycles > 1) {
		int node_a = -1;
		int node_b = -1;
		double move_cost = DBL_INFY;
		for (int i = 0; i < instance->num_nodes; i++) {
			for (int j = 0; j < instance->num_nodes; j++) {
				if (local_comp[i] < local_comp[j]) {
					if ((lookup_cost(i, local_succ[j], instance) + lookup_cost(j, local_succ[i], instance)) -
						(lookup_cost(i, local_succ[i], instance) + lookup_cost(j, local_succ[j], instance)) < move_cost) {
						node_a = i;
						node_b = j;
						move_cost = (lookup_cost(i, local_succ[j], instance) + lookup_cost(j, local_succ[i], instance)) -
							(lookup_cost(i, local_succ[i], instance) + lookup_cost(j, local_succ[j], instance));
					}
				}
			}
		}
		if (node_a < 0 || node_a > (instance->num_nodes - 1)) printf("\n\nnode_a %d is out of bound\n\n\n", node_a);
		if (node_b < 0 || node_b > (instance->num_nodes - 1)) printf("\n\nnode_b %d is out of bound\n\n\n", node_b);
#if VERBOSE > 4
		{ printf("Merging components: node_a %d (comp %d), node_b %d (comp %d), cost %f\n", node_a + 1, local_comp[node_a], node_b + 1, local_comp[node_b], move_cost); }
#endif
		int node_a_succ = local_succ[node_a];
		int node_b_succ = local_succ[node_b];
#if VERBOSE > 7
		{ printf("Before update: local_succ[%d] = %d, local_succ[%d] = %d\n", node_a + 1, local_succ[node_a] + 1, node_b + 1, local_succ[node_b] + 1); }
#endif
		local_succ[node_a] = node_b_succ;
		local_succ[node_b] = node_a_succ;
#if VERBOSE > 7
		{ printf("After update: local_succ[%d] = %d, local_succ[%d] = %d\n", node_a + 1, local_succ[node_a] + 1, node_b + 1, local_succ[node_b] + 1); }
#endif
		for (int i = 0; i < instance->num_nodes; i++)
			if ((i != node_b) && (local_comp[i] == local_comp[node_b])) {
#if VERBOSE > 7
				{ printf("local_comp[%d] went from %d to ", i + 1, local_comp[i]); }
#endif
				local_comp[i] = local_comp[node_a];
#if VERBOSE > 7
				{ printf("%d\n", local_comp[i]); }
#endif
			}
#if VERBOSE > 7
		{ printf("local_comp[%d] went from %d to ", node_b + 1, local_comp[node_b]); }
#endif
		local_comp[node_b] = local_comp[node_a];
#if VERBOSE > 7
		{ printf("%d\n", local_comp[node_b]); }
#endif
		num_cycles--;
#if VERBOSE > 4
		print_cycles(instance, local_succ, local_comp, num_cycles);
#endif
		add_sec(instance, env, lp, NULL, index, coeff, comp, cname, iteration, local_comp[node_a]);
#if VERBOSE > 0
		{ CPXwriteprob(env, lp, LP_MODEL_FILENAME, NULL); }
#endif
	}
	free(index);
	free(coeff);
	free(cname);

	refine_and_update_tsp_best_sol(instance, local_succ, num_cycles, cycle_delimiter);

	free(local_succ);
	free(local_comp);
	*/
}

static int CPXPUBLIC cplex_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle) {
	tsp_instance_t* instance = (tsp_instance_t*)userhandle;
	int sol_cardinality = ((instance->num_nodes) * (instance->num_nodes - 1)) / 2;
	double* xstar = (double*)malloc(sol_cardinality * sizeof(double));
	double objval = CPX_INFBOUND;
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, sol_cardinality - 1, &objval)) print_error("CPXcallbackgetcandidatepoint() error", 0);

	// get some random information at the node
	int mythread = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &mythread);
	int mynode = -1; CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &mynode);
	double incumbent = CPX_INFBOUND; CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
#if VERBOSE > 4
	{ printf("callback at node %d, thread %d, incumbent %.3lf\n", mynode, mythread, incumbent); }
#endif

	int* succ = (int*)malloc(instance->num_nodes * sizeof(int)); if (succ == NULL) print_error("Could not allocate succ[]", 0);
	int* comp = (int*)malloc(instance->num_nodes * sizeof(int)); if (comp == NULL) print_error("Could not allocate comp[]", 0);
	bool* visited_nodes = (bool*)malloc(instance->num_nodes * sizeof(bool)); if (visited_nodes == NULL) print_error("Could not allocate visited_nodes[]", 0);
	int num_cycles;
	build_model_ds(instance, xstar, succ, comp, visited_nodes, &num_cycles);
	//adding the single cycle SECs
	if (num_cycles > 1)
		add_single_cycle_secs(instance, NULL, NULL, context, comp, -1, num_cycles, sol_cardinality);

	//TO BE DELETED
	/*
	int num_constraints;
	if (num_cycles < 2)
		num_constraints = 0;
	else if (num_cycles == 2)
		num_constraints = 1;
	else
		num_constraints = num_cycles;
	int* index = (int*)calloc(sol_cardinality, sizeof(int));
	double* coeff = (double*)calloc(sol_cardinality, sizeof(double));
	char* cname = (char*)calloc(100, sizeof(char));
	for (int k = 1; k <= num_constraints; k++)
		add_sec(instance, NULL, NULL, context, index, coeff, comp, cname, sol_cardinality, k);
	free(index);
	free(coeff);
	free(cname);
	*/

	free(xstar);
	free(succ);
	free(comp);
	free(visited_nodes);

	return 0;

/*
	...

		int nnz = 0;
	... if xstart is infeasible, find a violated cutand store it in the usual Cplex's data structute (rhs, sense, nnz, index and value)

		if (nnz > 0) // means that the solution is infeasible and a violated cut has been found
		{
			int izero = 0;
			if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value)) print_error("CPXcallbackrejectcandidate() error"); // reject the solution and adds one cut 

			//if ( CPXcallbackrejectcandidate(context, 0, NULL, NULL, NULL, NULL, NULL, NULL) ) print_error("CPXcallbackrejectcandidate() error"); // just reject the solution without adding cuts (less effective)
		}

	free(xstar);
	return 0;
*/
}

void benders_method(tsp_instance_t* instance, CPXENVptr env, CPXLPptr lp, double** xstar, int* succ, int* comp, bool* visited_nodes, int* num_cycles) {
	//initialize the incumbent
	instance->sol_procedure_flag = GREEDY;
	instance->starting_index = -1;
	instance->prob_ign_opt = 0.0;
	instance->refine_flag = TWO_OPT;
	instance->metaheur_flag = NO_MH;
	if (tsp_opt(instance)) { free_tsp_instance(instance); fprintf(stderr, "Optimization algorithm failed\n"); exit(1); }
	if (ref_sol(instance)) { free_tsp_instance(instance); fprintf(stderr, "Refinement algorithm failed\n"); exit(1); }

	int error;
	int iteration = 0;
	double lower_bound = -DBL_INFY;
	double upper_bound = instance->best_sol_cost;
	while (((upper_bound - lower_bound) / upper_bound > EPSILON_PERCENT) && (instance->time_left > 0)) {
		double start_time = seconds();
		iteration++;

		//limit the number of nodes and add the upper bound to cplex (commenting these out might improve performance)
		CPXsetdblparam(env, CPX_PARAM_CUTUP, upper_bound);
		CPXsetintparam(env, CPX_PARAM_NODELIM, NODE_LIM);

		//solve the model
		error = CPXmipopt(env, lp);
		if (error) print_error("CPXmipopt() error", error);

		//get the solution
		int ncols = CPXgetnumcols(env, lp);
		//error = CPXgetx(env, lp, *xstar, 0, ncols - 1);
		//if (error) print_error("CPXgetx() error", error);
		if (CPXgetx(env, lp, *xstar, 0, ncols - 1)) continue;	//in case no solution was fount (cost higher than upper bound) retry

		//get the model data structures
		build_model_ds(instance, *xstar, succ, comp, visited_nodes, num_cycles);
		printf("---New solution found has %d cycles---\n", *num_cycles);

		//update the lower bound
		double temp;
		error = CPXgetbestobjval(env, lp, &temp);
		if (error) print_error("CPXgetbestobjval() error", error);
		lower_bound = max(lower_bound, temp);

		//patch the solution and update the best_sol
		add_sec_and_patch_update_tsp_best_sol(instance, env, lp, ncols, succ, comp, *num_cycles, instance->num_nodes, iteration);

		//update the upper bound
		upper_bound = instance->best_sol_cost;
		printf("Iteration %d: -Lower bound: %f; -Upper bound: %f; -Percent difference: %f%%\n\n", iteration, lower_bound, upper_bound, ((upper_bound - lower_bound) / upper_bound)*100);
		
		double end_time = seconds();
		instance->time_left -= (end_time - start_time);
		if (instance->time_left <= 0)
			printf("Benders' method has exhausted the time...\n");
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

	//OBSOLETE
	/*
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
	patch_and_update_tsp_best_sol(instance, succ, comp, num_cycles, instance->num_nodes);
	
	//plot the best sol cycles
	error = plot_cycles(instance);
	if (error) print_error("plot_cycles() error", 0);
	*/

	//necessary data structures
	int ncols = CPXgetnumcols(env, lp);
	double* xstar = (double*)calloc(ncols, sizeof(double)); if (xstar == NULL) print_error("Could not allocate xstar[]", 0);
	int* succ = (int*)malloc(instance->num_nodes * sizeof(int)); if (succ == NULL) print_error("Could not allocate succ[]", 0);
	int* comp = (int*)malloc(instance->num_nodes * sizeof(int)); if (comp == NULL) print_error("Could not allocate comp[]", 0);
	bool* visited_nodes = (bool*)malloc(instance->num_nodes * sizeof(bool)); if (visited_nodes == NULL) print_error("Could not allocate visited_nodes[]", 0);
	int num_cycles;

	if (instance->cplex_solver_flag == BENDERS) {	//find the optimal feasible solution using benders' method
		benders_method(instance, env, lp, &xstar, succ, comp, visited_nodes, &num_cycles);
		error = plot_cycles(instance); if (error) print_error("plot_cycles() error", 0);
	}
	else if (instance->cplex_solver_flag == CALLBACK) {	//find the optimal feasible solution using the callback function
		CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE;	//lazyconstraints
		error = CPXcallbacksetfunc(env, lp, contextid, cplex_callback, instance); if (error) print_error("CPXcallbacksetfunc() error", error);
		error = CPXmipopt(env, lp); if (error) print_error("CPXmipopt() error", error);	// with the callback installed
		error = CPXgetx(env, lp, xstar, 0, ncols - 1); if (error) print_error("CPXgetx() error", error);	//optimal solution if time limit not reached, a feasible one otherwise
		build_model_ds(instance, xstar, succ, comp, visited_nodes, &num_cycles);	//compute the successor list, component list, visited nodes list and number of cycles
		//print_cycles(instance, succ, comp, num_cycles);	//print the found cycles
		update_tsp_best_sol(instance, succ, num_cycles, instance->num_nodes);	//save the cycles in the instance best sol 
		error = plot_cycles(instance); if (error) print_error("plot_cycles() error", 0);	//plot the best sol cycles
	}
	else
		print_error("No cplex_solver defined", 0);

	free(xstar);
	free(succ);
	free(comp);
	free(visited_nodes);

	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	return 0;
}