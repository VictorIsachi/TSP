#ifndef TSP_UTILS_H
#define TSP_UTILS_H

//uncomment to disable assertions
// #define NDEBUG

//might be usefull in the future
// #define LOOP for(;;)
// #define GENERIC_FUNC(out_type, in_type, func_name)       \
// out_type in_type##_##func_name(in_type x, in_type y){    \
//  /*func body*/                                       \
// }                                                        \
// 
// inline int max_int(int a, int b) { return a < b ? b : a; }

#include <assert.h>

//need to figure out why -UVERBOSE does not undefine the macro
#define VERBOSE 3

#define OMP_NUM_THREADS 4

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

#endif //TSP_UTILS_H