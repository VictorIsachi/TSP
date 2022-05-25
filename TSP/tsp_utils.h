#ifndef TSP_UTILS_H
#define TSP_UTILS_H

//uncomment to disable assertions
// #define NDEBUG

#include <assert.h>

//need to figure out why -UVERBOSE does not undefine the macro
#define VERBOSE 2

#define OMP_NUM_THREADS 4

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

/**
 * @brief count the number of seconds elapsed from the beginning of the program
 * @return the number of seconds elapesed
 */
double seconds();

/**
 * @brief swaps the a-th and the b-th elements of the array
 * @param a index of the first element
 * @param b index of the second element
 * @param array of unsigned int
 */
void swap_uint_array(unsigned int a, unsigned int b, unsigned int array[]);

/**
 * @brief function that can be used in qsort() to compare integer values
 * @return a - b
 */
int cmp_int(const void* a, const void* b);

/**
 * @brief prints an error message and the error code (if any), then terminates execution
 * @param error_message the message to be displayed
 * @param error_code the code error, if 0 then no error code will be displayed
 */
void print_error(const char* error_message, const int error_code);

#endif //TSP_UTILS_H