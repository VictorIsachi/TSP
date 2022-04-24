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

int cmp_int(const void* a, const void* b) {
	return (*(int*)a - *(int*)b);
}

void print_error(const char* error_message, const int error_code) {
	if (error_code)
		printf("\n\n ERROR(code %d): %s \n\n", error_code, error_message);
	else
		printf("\n\n ERROR: %s \n\n", error_message);
	fflush(NULL);
	exit(1);
}