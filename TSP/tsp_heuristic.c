#include "tsp_utils.h"
#include "tsp.h"

#define INSTANCE_POINTS_PLOT "points.dat"
#define INSTANCE_TOUR_PLOT "tour.dat"
#define GNUPLOT_COMMAND_POINTS "gnuplot commands_points.txt"
#define GNUPLOT_COMMAND_TOUR "gnuplot commands_tour.txt"

#define TABU_LIST_MAX 100
#define TABU_LIST_MAX_DIV 10
#define TABU_LIST_MIN 75
#define TABU_LIST_MIN_DIV 15

#define DEFAULT_MAX_TEMP 100
#define DEFAULT_MIN_TEMP 10

#define NUM_CHILDREN ((instance->pop_size) / 2)

int tsp_opt(tsp_instance_t* instance) {

#if VERBOSE > 1
	{ printf("Solving the tsp instance...\n"); }
#endif

	assert(instance != NULL);
	assert(instance->sol_procedure_flag != NULL);
	assert(instance->metaheur_flag != NULL);

	if (instance->num_nodes < 2) {
		printf("No tour possible with %d nodes\n", instance->num_nodes);
	}
	else if (instance->num_nodes < 4) {
		instance->sol_procedure_flag = SEQUENTIAL;
		instance->refine_flag = NO_REF;
		instance->metaheur_flag = NO_MH;
		return tsp_seq_sol(instance);
	}
	else if (instance->num_nodes == 4) {
		instance->sol_procedure_flag = GREEDY;
		instance->refine_flag = TWO_OPT;
		instance->metaheur_flag = NO_MH;
		return tsp_gdy_sol(instance);
	}
	else if (instance->metaheur_flag == GEN) {
		instance->sol_procedure_flag = SEQUENTIAL;
		instance->refine_flag = NO_REF;
		tsp_seq_sol(instance);
		return genetic(instance);
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

	//initializing the tabu list data
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
	bool computed_avg_cost = false;
	double avg_move_cost = 0.0;
	while (!done && instance->time_left > 0) {

		double start_time = seconds();

		int cross_node_a_index = -1;
		int cross_node_b_index = -1;
		double move_cost = DBL_INFY;

		//random 2-opt move
		if (instance->metaheur_flag == SIM_ANNEAL) {
			cross_node_a_index = (unsigned int)floor(((double)rand() / (RAND_MAX + 1)) * instance->num_nodes);	//[0, instance->num_nodes - 1]
			do {
				cross_node_b_index = (unsigned int)floor(((double)rand() / (RAND_MAX + 1)) * instance->num_nodes);	//[0, instance->num_nodes - 1]
			} while (cross_node_a_index == cross_node_b_index - 1 || cross_node_a_index == cross_node_b_index || cross_node_a_index == cross_node_b_index + 1);	//while the 2 nodes are cosecutive
			if (cross_node_a_index > cross_node_b_index) {	//enforce cross_node_a_index < cross_node_b_index
				unsigned int temp = cross_node_a_index;
				cross_node_a_index = cross_node_b_index;
				cross_node_b_index = temp;
			}
			move_cost = (lookup_cost(current_sol[cross_node_a_index], current_sol[cross_node_b_index], instance) + \
				lookup_cost(current_sol[(cross_node_a_index + 1) % instance->num_nodes], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance)) - \
				(lookup_cost(current_sol[cross_node_a_index], current_sol[(cross_node_a_index + 1) % instance->num_nodes], instance) + \
					lookup_cost(current_sol[cross_node_b_index], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance));

#if VERBOSE > 2
			{ printf("Random move: node %d, node %d, cost %f\n", current_sol[cross_node_a_index], current_sol[cross_node_b_index], move_cost); }
#endif
		}
		//looking for the min cost 2-opt move
		else {
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
						{ if (cross_node_a_index < 0 || cross_node_b_index < 0) { printf("Infeasible 2-opt move\n"); }
						else { printf("The cost of going  %d -> %d is %f; cost of going %d -> %d is %f\n", current_sol[cross_node_a_index],
							current_sol[(cross_node_a_index + 1) % instance->num_nodes],
							lookup_cost(current_sol[cross_node_a_index], current_sol[(cross_node_a_index + 1) % instance->num_nodes], instance), current_sol[cross_node_b_index],
							current_sol[(cross_node_b_index + 1) % instance->num_nodes],
							lookup_cost(current_sol[cross_node_b_index], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance));
						printf("Substituted these edges with %d -> %d, cost %f and %d -> %d, cost %f\n", current_sol[cross_node_a_index], current_sol[cross_node_b_index],
							lookup_cost(current_sol[cross_node_a_index], current_sol[cross_node_b_index], instance), current_sol[(cross_node_a_index + 1) % instance->num_nodes],
							current_sol[(cross_node_b_index + 1) % instance->num_nodes],
							lookup_cost(current_sol[(cross_node_a_index + 1) % instance->num_nodes], current_sol[(cross_node_b_index + 1) % instance->num_nodes], instance));} 
						}
#endif
					}
				}
			}

			if (cross_node_a_index < 0 || cross_node_b_index < 0) { fprintf(stderr, "2-opt found no feasible move: aborting\n"); return 1; }

#if VERBOSE > 2
			{ printf("Best move found: node %d, node %d, cost %f\n", current_sol[cross_node_a_index], current_sol[cross_node_b_index], move_cost); }
#endif
		}

		//no improveing 2-opt move and no meta-heuristic: finish refinement
		if (move_cost >= 0 && instance->metaheur_flag == NO_MH) done = true;
		//5-kick
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
		//simulated annealing stochastic move
		else if (move_cost >= 0 && instance->metaheur_flag == SIM_ANNEAL) {
#if VERBOSE > 3
			{ printf("Deciding whether to accept the move (simulated annealing)...\n"); }
#endif
			//compute the average cost of a refined solution
			if (!computed_avg_cost) {
				for (int i = 0; i < instance->num_nodes - 1; i++)
					for (int j = i + 1; j < instance->num_nodes; j++)
						avg_move_cost += lookup_cost(i, j, instance);
				avg_move_cost /= ((instance->num_nodes * (instance->num_nodes - 1)) / 2);
				computed_avg_cost = true;
			}

			if (instance->min_temperature < 0) { instance->min_temperature = DEFAULT_MIN_TEMP; }
			if (instance->max_temperature < 0) { instance->max_temperature = DEFAULT_MAX_TEMP; }

			int temperature = instance->max_temperature - (current_iter % (instance->max_temperature + 1 - instance->min_temperature));
			double prob_acceptance = exp(-(instance->move_weight * move_cost * instance->max_temperature) / (temperature * avg_move_cost));

#if VERBOSE > 2
			{ printf("Temperature range [%d, %d], average move cost %f, temperature %d, prob_acceptance %f\n", 
				instance->min_temperature, instance->max_temperature, avg_move_cost, temperature, prob_acceptance); }
#endif

			//accept move
			if ((double)rand() / RAND_MAX < prob_acceptance) {
#if VERBOSE > 3
				{ printf("Move accepted...\n"); }
#endif
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

				if (current_sol_cost < instance->best_sol_cost) {
					memcpy(instance->best_sol, current_sol, instance->num_nodes * sizeof(unsigned int));
					instance->best_sol_cost = current_sol_cost;

#if VERBOSE > 2
					{ printf("Lucky swap (random 2-opt move improved best known sol)!\nSolution cost %f\n", instance->best_sol_cost); for (int i = 0; i < instance->num_nodes; i++)
						printf("Element of index %d of the tour is the node of index %d [%f, %f]\n", i, instance->best_sol[i],
							instance->nodes[instance->best_sol[i]].x_coord, instance->nodes[instance->best_sol[i]].y_coord); }
#endif
				}
			}
		}
		//applying the 2-opt move
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

#if VERBOSE > 3
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

static tsp_instance_t* mate(tsp_instance_t* population, tsp_instance_t* instance, double pop_cprob) {

#if VERBOSE > 3
	{ printf("Mateing...\n"); }
#endif

	//generate random number that will determine the indices od the 2 parents
	double parent_1_pos = ((double)rand() / (RAND_MAX + 1)) * pop_cprob;
	double parent_2_pos = ((double)rand() / (RAND_MAX + 1)) * pop_cprob;

	//population indices of the 2 parents
	int parent_1_index = -1;
	int parent_2_index = -1;

	//accumulates the probability distribution
	double accumulator = 0;

	//find the indices of the parents
	for (int i = 0; i < instance->pop_size; i++) {
		accumulator += (1 / population[i].best_sol_cost);
		if ((accumulator > parent_1_pos) && (parent_1_index == -1))
			parent_1_index = i;
		if ((accumulator > parent_2_pos) && (parent_2_index == -1))
			parent_2_index = i;
	}
	
#if VERBOSE > 2
	if ((parent_1_index < 0) || (parent_1_index > instance->pop_size - 1)) { fprintf(stderr, "ERROR: Parent 1 for mating not found\n"); return NULL; }
	if ((parent_2_index < 0) || (parent_2_index > instance->pop_size - 1)) { fprintf(stderr, "ERROR: Parent 2 for mating not found\n"); return NULL; }
	printf("Mates: %d and %d\n", parent_1_index, parent_2_index);
#endif

	//child
	tsp_instance_t* child = (tsp_instance_t*)malloc(sizeof(tsp_instance_t));
	if (child == NULL) { fprintf(stderr, "ERROR: Could not allocate child\n"); return NULL; }

	//make child identical to parent 1 (except for the optimal solution)
	memcpy(child, &population[parent_1_index], sizeof(tsp_instance_t));
	child->best_sol = (unsigned int*)malloc(instance->num_nodes * sizeof(unsigned int));

	//updated half of the child's best sol with the nodes of parent 1
	for (int i = 0; i < (instance->num_nodes / 2); i++)
		child->best_sol[i] = population[parent_1_index].best_sol[i];

	//updated half of the child's best sol with the nodes of parent 2
	for (int i = (instance->num_nodes / 2); i < instance->num_nodes; i++)
		child->best_sol[i] = population[parent_2_index].best_sol[i];

	//short-cut child solution (eliminate loops) 
	bool* visited_nodes = (bool*)malloc(instance->num_nodes * sizeof(bool));
	for (int i = 0; i < instance->num_nodes; i++)
		visited_nodes[i] = false;
	int num_tour_nodes = 0;
	for (int i = 0; i < instance->num_nodes; i++) {
		if (!visited_nodes[child->best_sol[i]]) {
			child->best_sol[num_tour_nodes++] = child->best_sol[i];
			visited_nodes[child->best_sol[i]] = true;
		}
	}

	//fill up the rest of the solution array with the unvisited nodes (used for the in-place greedy alg.)
	int num_unvisited_nodes = 0;
	for (int i = 0; i < instance->num_nodes; i++)
		if (!visited_nodes[i])
			child->best_sol[num_tour_nodes + num_unvisited_nodes++] = i;
	if (num_tour_nodes + num_unvisited_nodes != instance->num_nodes) { fprintf(stderr, "ERROR: Total number of nodes not matching\n"); return NULL; }

	//visit the unvisited nodes using the greedy in-place method
	while (num_unvisited_nodes > 0) {
		int next_tour_elem_ind = -1;
		double next_move_cost = DBL_INFY;
		for (int i = num_tour_nodes; i < instance->num_nodes; i++) {
			if (lookup_cost(child->best_sol[(num_tour_nodes - 1) % instance->num_nodes], child->best_sol[i], instance) < next_move_cost) {
				next_tour_elem_ind = i;
				next_move_cost = lookup_cost(child->best_sol[num_tour_nodes - 1], child->best_sol[next_tour_elem_ind], instance);
			}
		}
		if (next_tour_elem_ind < num_tour_nodes || next_tour_elem_ind >= instance->num_nodes) { fprintf(stderr, "ERROR: Invalid next element index\n"); return NULL; }

#if VERBOSE > 4
		{ printf("Swapping child->best_sol[%d] = %d with child->best_sol[%d] = %d\n", num_tour_nodes, child->best_sol[num_tour_nodes],
			next_tour_elem_ind, child->best_sol[next_tour_elem_ind]); }
#endif
		swap_uint_array(num_tour_nodes, next_tour_elem_ind, child->best_sol);
#if VERBOSE > 4
		{ printf("After swapping child->best_sol[%d] = %d with child->best_sol[%d] = %d\n", num_tour_nodes, child->best_sol[num_tour_nodes],
			next_tour_elem_ind, child->best_sol[next_tour_elem_ind]); }
#endif
		num_tour_nodes++;
		num_unvisited_nodes--;
	}

	//update the solution cost
	double sol_cost = 0.0;
	for (int i = 0; i < instance->num_nodes; i++) {
		sol_cost += lookup_cost(child->best_sol[i], child->best_sol[(i + 1) % instance->num_nodes], instance);
	}
	child->best_sol_cost = sol_cost;

#if VERBOSE > 8	//check for solution integrity
	printf("Child cost %f\n", child->best_sol_cost);
	int* node_deg = (int*)malloc(instance->num_nodes * sizeof(int));
	for (int i = 0; i < instance->num_nodes; i++) { node_deg[i] = 0; }
	for (int i = 0; i < instance->num_nodes; i++) { node_deg[child->best_sol[i]]++; node_deg[child->best_sol[(i + 1) % instance->num_nodes]]++; }
	bool degree_violation = false;
	for (int i = 0; i < instance->num_nodes; i++) { if (node_deg[i] != 2) degree_violation = true; }
	if (degree_violation) printf("Degree property violated: MALFORMED CHILD\n"); else printf("Degree property asserted: FEASIBLE CHILD\n");
	free(node_deg);
#endif

	free(visited_nodes);

	return child;
}

static int procreate(tsp_instance_t* population, tsp_instance_t* instance) {

#if VERBOSE > 2
	{ printf("Generating new population...\n"); }
#endif

	tsp_instance_t* new_population = (tsp_instance_t*)malloc((instance->pop_size + NUM_CHILDREN) * sizeof(tsp_instance_t));
	if (new_population == NULL) { printf("Could not allocate memory for the new generation\n"); return 1; }
	memcpy(new_population, population, instance->pop_size * sizeof(tsp_instance_t));

#if VERBOSE > 8
	printf("Beginning of procreation\n");
	printf("population[i]:\n");
	for (int i = 0; i < instance->pop_size; i++) {
		printf("0x%p i = %d [ ", &population[i], i);
		for (int j = 0; j < instance->num_nodes; j++)
			printf("%d ", population[i].best_sol[j]);
		printf("]\n");
	}
	printf("new_population[i]:\n");
	for (int i = 0; i < instance->pop_size; i++) {
		printf("0x%p i = %d [ ", &new_population[i], i);
		for (int j = 0; j < instance->num_nodes; j++)
			printf("%d ", new_population[i].best_sol[j]);
		printf("]\n");
	}
#endif

	//cumulative probability of the population (before mating)
	double pop_cprob = 0;
	for (int i = 0; i < instance->pop_size; i++)
		pop_cprob += (1 / new_population[i].best_sol_cost);

	//generated the new population
	for (int i = 0; i < NUM_CHILDREN; i++) {
		tsp_instance_t* child = mate(population, instance, pop_cprob);
		if (child == NULL) { printf("No child returned\n"); return 1; }
		new_population[instance->pop_size + i] = *(child);
		free(child);

		if (ref_sol(&new_population[instance->pop_size + i])) { free_tsp_instance(&new_population[instance->pop_size + i]); fprintf(stderr, "Refinement algorithm failed\n"); exit(1); }
#if VERBOSE > 8
		{ printf("Solution cost after refinement: %f\n", new_population[instance->pop_size + i].best_sol_cost); }
#endif
		//if (plot_tour(&new_population[instance->pop_size + i])) { free_tsp_instance(&new_population[instance->pop_size + i]); fprintf(stderr, "Tour plotting failed\n"); exit(1); }

#if VERBOSE > 6
		printf("Newly added child, at position %d, has cost %f and is [ ", instance->pop_size + i, new_population[instance->pop_size + i].best_sol_cost);
		for (int j = 0; j < instance->num_nodes; j++)
			printf("%d ", new_population[instance->pop_size + i].best_sol[j]);
		printf("]\n");
		//plot_tour(&new_population[instance->pop_size + i]);
#endif
	}

#if VERBOSE > 8
	printf("End of procreation, before pruning\n");
	printf("population[i]:\n");
	for (int i = 0; i < instance->pop_size; i++) {
		printf("0x%p i = %d [ ", &population[i], i);
		for (int j = 0; j < instance->num_nodes; j++)
			printf("%d ", population[i].best_sol[j]);
		printf("]\n");
	}
	printf("new_population[i]:\n");
	for (int i = 0; i < instance->pop_size + NUM_CHILDREN; i++) {
		printf("0x%p i = %d [ ", &new_population[i], i);
		for (int j = 0; j < instance->num_nodes; j++)
			printf("%d ", new_population[i].best_sol[j]);
		printf("]\n");
	}
#endif

	//prune the population
	for (int i = 0; i < instance->pop_size; i++) {

		//cumulative probability of the population (after mating, during pruning)
		double pop_cprob = 0;
		for (int j = 0; j < instance->pop_size + NUM_CHILDREN - i; j++)
			pop_cprob += new_population[j].best_sol_cost;

		//determine the index of the next element to be saved for the next generation
		double curr_tsp_inst_pos = ((double)rand() / (RAND_MAX + 1)) * pop_cprob;
		int curr_tsp_inst_ind = -1;
		double accumulator = 0;
		for (int j = 0; j < instance->pop_size + NUM_CHILDREN - i; j++) {
			accumulator += new_population[j].best_sol_cost;
			if (accumulator > curr_tsp_inst_pos) {
				curr_tsp_inst_ind = j;
				break;
			}
		}
		if ((curr_tsp_inst_ind < 0) || (curr_tsp_inst_ind > instance->pop_size + NUM_CHILDREN - 1 - i)) { printf("Pruning population out of bound\n"); return 1; }

		//add the instance to be saved to the population and update the new population array
		population[i] = new_population[curr_tsp_inst_ind];
		new_population[curr_tsp_inst_ind] = new_population[instance->pop_size + NUM_CHILDREN - 1 - i];

#if VERBOSE > 8
		{ printf("Pruning step: population[%d] = new_population[%d], new_population[%d] = new_population[%d]\n", i, curr_tsp_inst_ind, curr_tsp_inst_ind, 
			instance->pop_size + NUM_CHILDREN - 1 - i); }
#endif
	}

#if VERBOSE > 8
	printf("End of procreation, after pruning\n");
	printf("population[i]:\n");
	for (int i = 0; i < instance->pop_size; i++) {
		printf("0x%p i = %d [ ", &population[i], i);
		for (int j = 0; j < instance->num_nodes; j++)
			printf("%d ", population[i].best_sol[j]);
		printf("]\n");
	}
	printf("new_population[i]:\n");
	for (int i = 0; i < instance->pop_size + NUM_CHILDREN; i++) {
		printf("0x%p i = %d [ ", &new_population[i], i);
		for (int j = 0; j < instance->num_nodes; j++)
			printf("%d ", new_population[i].best_sol[j]);
		printf("]\n");
	}
#endif

	//eliminate the instances that were not selected
	for (int i = 0; i < NUM_CHILDREN; i++)
		free(new_population[i].best_sol);
	free(new_population);

	return 0;
}

static void update_champion(tsp_instance_t* population, tsp_instance_t* instance) {

#if VERBOSE > 0
	{ printf("Current champion cost: %f", instance->best_sol_cost); }
#endif

	for (int i = 0; i < instance->pop_size; i++) {
		if (population[i].best_sol_cost < instance->best_sol_cost) {
			memcpy(instance->best_sol, population[i].best_sol, instance->num_nodes * sizeof(unsigned int));
			instance->best_sol_cost = population[i].best_sol_cost;
		}
	}
	//makes sure that the champion from the previous iteration is present in the current iteration
	memcpy(population[0].best_sol, instance->best_sol, instance->num_nodes * sizeof(unsigned int));
	population[0].best_sol_cost = instance->best_sol_cost;

#if VERBOSE > 0
	{ printf("; new champion cost: %f\n", instance->best_sol_cost); }
#endif
}

int genetic(tsp_instance_t* instance) {

	srand(instance->random_seed); for (int i = 0; i < MIN_RAND_RUNS + log(1 + instance->random_seed); i++) rand();

#if VERBOSE > 1
	{ printf("Applying the genetic algorithm metaheuristic...\n"); }
#endif

	tsp_instance_t* population = (tsp_instance_t*)malloc(instance->pop_size * sizeof(tsp_instance_t));
	if (population == NULL) { fprintf(stderr, "Could not allocate memory for the population\n"); return 1; }

	//initialize the population
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
	for (int i = 0; i < instance->pop_size; i++) {
		population[i].time_limit = MAX_TIME;
		population[i].input_file_name = NULL;
		population[i].x_bound = -1;
		population[i].y_bound = -1;
		population[i].num_nodes = instance->num_nodes;
		population[i].nodes = instance->nodes;
		population[i].random_seed = ((double)rand() / RAND_MAX) * UINT_MAX;	//the random seed is by def. an unsigned int
		population[i].sol_procedure_flag = GREEDY;
		population[i].starting_index = -1;
		population[i].prob_ign_opt = 0.5;
		population[i].refine_flag = TWO_OPT;
		population[i].metaheur_flag = NO_MH;
		population[i].min_tenure = -1;
		population[i].max_tenure = -1;
		population[i].min_temperature = -1;
		population[i].max_temperature = -1;
		population[i].move_weight = 0;
		population[i].pop_size = 0;
		population[i].best_sol = NULL;
		population[i].best_sol_cost = DBL_INFY;
		population[i].costs = instance->costs;
		population[i].time_left = population[i].time_limit;
		population[i].tabu_list = NULL;

		if (tsp_opt(&population[i])) { free_tsp_instance(&population[i]); fprintf(stderr, "Optimization algorithm failed\n"); exit(1); }
#if VERBOSE > 8
		{ printf("Solution cost greedy+GRASP0.5: %f\n", population[i].best_sol_cost); }
#endif
		//if (plot_tour(&population[i])) { free_tsp_instance(&population[i]); fprintf(stderr, "Tour plotting failed\n"); exit(1); }
		if (ref_sol(&population[i])) { free_tsp_instance(&population[i]); fprintf(stderr, "Refinement algorithm failed\n"); exit(1); }
#if VERBOSE > 8
		{ printf("Solution cost after refinement: %f\n", population[i].best_sol_cost); }
#endif
		//if (plot_tour(instance)) { free_tsp_instance(instance); fprintf(stderr, "Tour plotting failed\n"); exit(1); }
	}

#if VERBOSE > 8
	printf("Population initialized\n");
	printf("population[i]:\n");
	for (int i = 0; i < instance->pop_size; i++) {
		printf("0x%p i = %d [ ", &population[i], i);
		for (int j = 0; j < instance->num_nodes; j++)
			printf("%d ", population[i].best_sol[j]);
		printf("]\n");
	}
#endif

	update_champion(population, instance);

#if VERBOSE > 8
	printf("Population after champion update\n");
	printf("population[i]:\n");
	for (int i = 0; i < instance->pop_size; i++) {
		printf("0x%p i = %d [ ", &population[i], i);
		for (int j = 0; j < instance->num_nodes; j++)
			printf("%d ", population[i].best_sol[j]);
		printf("]\n");
	}
#endif

	int generations = 0;
	while (instance->time_left > 0) {
		
		double start_time = seconds();

		if (procreate(population, instance)) { printf("Procreation failed, aborting...\n"); return 1; }
		update_champion(population, instance);
		generations++;

#if VERBOSE > 8
		printf("Population after champion update\n");
		printf("population[i]:\n");
		for (int i = 0; i < instance->pop_size; i++) {
			printf("0x%p i = %d [ ", &population[i], i);
			for (int j = 0; j < instance->num_nodes; j++)
				printf("%d ", population[i].best_sol[j]);
			printf("]\n");
		}
#endif

		double end_time = seconds();
		instance->time_left -= (end_time - start_time);

#if VERBOSE > 0
		{ if (instance->time_left <= 0) printf("Time limit reached in the generic algorithm, %d generations executed...\n", generations); }
#endif
	}

	for (int i = 0; i < instance->pop_size; i++)
		free(population[i].best_sol);
	free(population);

	return 0;
}