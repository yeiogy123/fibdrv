#ifndef _BIG_NUM_H_
#define _BIG_NUM_H_
#include<stddef.h>

/*
 * big_num data structure
 * num[0] stores for lsb
 * num[size-1] stores for msb
 * sign represent for positive or negative
 * */
typedef struct _big_num{
	unsigned int *num;
	unsigned int size;
	int sign;
}big_num;
/*
 * to transform input to string
 * */
char *big_num_to_string(const big_num *src);
/*
 * allocate a big_num data structure for a given size
 * */
big_num *big_num_allocate(size_t size);
/*
 * free the big_num data structure
 * */
int big_num_free(big_num *src);
/*
 * copy the src to dest
 * */
int big_num_resize(big_num *src, size_t size);
int big_num_cpy(big_num *dest, big_num *src);
/*
 * swap the given to big_num
 * */
void big_num_swap(big_num *a, big_num *b);
/*
 * compare two given big_num
 * */
int big_num_cmp(const big_num *src1, const big_num *src2);
/* 
 * left shift on big_num
 *  */
static int big_num_clz(const big_num *src);
void big_num_left_shift(big_num *src, size_t shift);
/*
 * target = src1 + src2
 * */
void big_num_add(const big_num *src1, const big_num *src2, big_num *target);
/* 
 * target = src1 - src2
 * */
void big_num_sub(const big_num *src1, const big_num *src2, big_num *target);
/* 
 * target = src1 * src2
 * */
void big_num_mul(const big_num *src1, const big_num *src2, big_num *target);
/* 
 * calculate nth Fibonacci number
 * */
void big_num_fib_fast_doubling(big_num *dest, unsigned int n);
void big_num_fib(big_num *dest, unsigned int n);
#endif
