#include <stdio.h>
#include <stdlib.h>
#include "tsp_utils.h"
#include "tsp.h"

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

	tsp_instance_t instance;
	if (parse_command_line(argc, argv, &instance)) { fprintf(stderr, "Parsing failed\n"); exit(1); }
	if (get_data(&instance)) { free_tsp_instance(&instance); fprintf(stderr, "Input file reading failed\n"); exit(1); }
	if (plot_points(&instance)) { free_tsp_instance(&instance); fprintf(stderr, "Point plotting failed\n"); exit(1); }
	if (precompute_costs(&instance)) { free_tsp_instance(&instance); fprintf(stderr, "Cost precomputation failed\n"); exit(1); }
	int_round_clut(&instance);
	if (tsp_opt(&instance)) { free_tsp_instance(&instance); fprintf(stderr, "Optimization algorithm failed\n"); exit(1); }
	printf("Solution cost: %f\n", instance.best_sol_cost);
	if (plot_tour(&instance)) { free_tsp_instance(&instance); fprintf(stderr, "Tour plotting failed\n"); exit(1); }
	if (ref_sol(&instance)) { free_tsp_instance(&instance); fprintf(stderr, "Refinement algorithm failed\n"); exit(1); }
	printf("Solution cost: %f\n", instance.best_sol_cost);
	if (plot_tour(&instance)) { free_tsp_instance(&instance); fprintf(stderr, "Tour plotting failed\n"); exit(1); }

	//TESTS
//	printf("Distance between node 99 and 98: %f\n", lookup_cost(99, 98, &instance));

	double end_time = seconds();

#if VERBOSE > 0
	{
		printf("Execution time: %fs\n", end_time - beg_time);
	}
#endif

	free_tsp_instance(&instance);
	return 0;
}