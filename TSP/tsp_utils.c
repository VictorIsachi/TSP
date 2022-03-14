#include <time.h>
#include "tsp_utils.h"

double seconds() {
	return (double)clock() / CLOCKS_PER_SEC;
}

void swap_uint_array(unsigned int a, unsigned int b, unsigned int array[]) {
	unsigned int temp = array[a];
	array[a] = array[b];
	array[b] = temp;
}