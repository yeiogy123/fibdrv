#include"big_num.h"
#include<linux/slab.h>
#define DIV_ROUNDUP(x, len) (((x) + (len) -1) / (len))
#define MAX(a, b) ((a)>(b)? (a):(b))
#ifndef SWAP
#define SWAP(a, b)			\
	do{				\
		typeof(a)_tmp = a;	\
		a = b;			\
		b = _tmp;		\
	}while(0)
#endif
big_num *big_num_allocate(size_t size){
	big_num *new = (big_num*)kmalloc(sizeof(big_num), GFP_KERNEL);
	new->num = (int*)kmalloc(sizeof(int) * size, GFP_KERNEL);
	memset(new->num, 0, sizeof(int) * size);
	new->size = size;
	new->sign = 0;
	return new;
}

char *big_num_to_string(const big_num *src){
    // log10(x) = log2(x) / log2(10) ~= log2(x) / 3.322
    size_t len = (8 * sizeof(int) * src->size) / 3 + 2 + src->sign;
    char *s = kmalloc(len, GFP_KERNEL);
    char *p = s;

    memset(s, '0', len - 1);
    s[len - 1] = '\0';

    /* src.number[0] contains least significant bits */
    for (int i = src->size - 1; i >= 0; i--) {
        /* walk through every bit of bn */
        for (unsigned int d = 1U << 31; d; d >>= 1) {
            /* binary -> decimal string */
            int carry = !!(d & src->num[i]);
            for (int j = len - 2; j >= 0; j--) {
                s[j] += s[j] - '0' + carry;
                carry = (s[j] > '9');
                if (carry)
                    s[j] -= 10;
            }
        }
    }
    // skip leading zero
    while (p[0] == '0' && p[1] != '\0') {
        p++;
    }
    if (src->sign)
        *(--p) = '-';
    memmove(s, p, strlen(p) + 1);
    return s;
}
static int big_num_msb(const big_num *src){
    return src->size * 32 - big_num_clz(src);
}

/*
 * return -1 if it doesn't free successfully, else return 0
 * */
int big_num_free(big_num *src){
	if(!src) return -1;
	kfree(src->num);
	kfree(src);
	return 0;
}
int big_num_resize(big_num *src, size_t size){
	if(!src) return -1; //no src
	if(size == src->size) return 0; // same size 
	src->num = (int*)krealloc(src->num, sizeof(int) * size, GFP_KERNEL);
	if(!src->num) return-1; 
	return 0;
}
int big_num_cpy(big_num *dest, big_num *src){
	if(big_num_resize(dest, src->size) < 0) return -1;
	dest->sign = src->sign;
	memcpy(dest->num, src->num, src->size * sizeof(int));
	return 0;
}
void big_num_swap(big_num *a, big_num *b){
	big_num *tmp = a;
	*a = *b;
	*b = *tmp;
}
int big_num_cmp(const big_num *src1, const big_num *src2){
	if(src1->size > src2->size) return 1;
	else if(src1->size < src2->size) return -1;
	else{
		int position = src1->size - 1;
		while(position >= 0){
			if(src1->num > src2->num) return 1;
			else if (src1->num < src2->num) return -1;
			position--;
		}
		return 0;
	}
}
static int big_num_clz(const big_num *src){
    int cnt = 0;
    for (int i = src->size - 1; i >= 0; i--) {
        if (src->num[i]) {
            cnt += __builtin_clz(src->num[i]);
            return cnt;
        } else {
            cnt += 32;
        }
    }
    return cnt;
}
void big_num_left_shift(big_num *src, size_t shift){
	size_t z = big_num_clz(src);
	shift %= 32;
	if(!shift) return;
	if(shift > z) big_num_resize(src, src->size + 1);
	else big_num_resize(src, src->size);
	for(int i = src->size - 1; i > 0; i--)
		src->num[i] = src->num[i] << shift | src->num[i - 1] >> (32 - shift);
	src->num[0] <<= shift;
}
static void big_num_mult_add(big_num *c, int offset, unsigned long long int x){
    unsigned long long int carry = 0;
    for (int i = offset; i < c->size; i++) {
        carry += c->num[i] + (x & 0xFFFFFFFF);
        c->num[i] = carry;
        carry >>= 32;
        x >>= 32;
        if (!x && !carry)  // done
            return;
    }
}
static void big_num_do_add(const big_num *a, const big_num *b, big_num *c){
    // max digits = max(sizeof(a) + sizeof(b)) + 1
    int d = MAX(big_num_msb(a), big_num_msb(b)) + 1;
    d = DIV_ROUNDUP(d, 32) + !d;
    big_num_resize(c, d);  // round up, min size = 1

    unsigned long long int carry = 0;
    for (int i = 0; i < c->size; i++) {
        unsigned int tmp1 = (i < a->size) ? a->num[i] : 0;
        unsigned int tmp2 = (i < b->size) ? b->num[i] : 0;
        carry += (unsigned long long int) tmp1 + tmp2;
        c->num[i] = carry;
        carry >>= 32;
    }

    if (!c->num[c->size - 1] && c->size > 1)
        big_num_resize(c, c->size - 1);
}

static void big_num_do_sub(const big_num *a, const big_num *b, big_num *c){
    // max digits = max(sizeof(a) + sizeof(b))
    int d = MAX(a->size, b->size);
    big_num_resize(c, d);
    long long int carry = 0;
    for (int i = 0; i < c->size; i++) {
        unsigned int tmp1 = (i < a->size) ? a->num[i] : 0;
        unsigned int tmp2 = (i < b->size) ? b->num[i] : 0;
        carry = (long long int) tmp1 - tmp2 - carry;
        if (carry < 0) {
            c->num[i] = carry + (1LL << 32);
            carry = 1;
        } else {
            c->num[i] = carry;
            carry = 0;
        }
    }

    d = big_num_clz(c) / 32;
    if (d == c->size)
        --d;
    big_num_resize(c, c->size - d);
}
void big_num_add(const big_num *src1, const big_num *src2, big_num *target){
    if(src1->sign == src2->sign) {  // both positive or negative
        big_num_do_add(src1, src2, target);
        target->sign = src1->sign;
    } else {          // different sign
        if (src1->sign)  // let a > 0, b < 0
            SWAP(src1, src2);
        int cmp = big_num_cmp(src1, src2);
        if (cmp > 0) {
            /* |a| > |b| and b < 0, hence c = a - |b| */
            big_num_do_sub(src1, src2, target);
            target->sign = 0;
        } else if (cmp < 0) {
            /* |a| < |b| and b < 0, hence c = -(|b| - |a|) */
            big_num_do_sub(src2, src1, target);
            target->sign = 1;
        } else {
            /* |a| == |b| */
            big_num_resize(target, 1);
            target->num[0] = 0;
            target->sign = 0;
        }
    }
}
void big_num_sub(const big_num *src1, const big_num *src2, big_num *target){
	big_num tmp = *src1;
	tmp.sign ^= 1;
	big_num_add(target, &tmp, &*src2);
}
void big_num_mul(const big_num *src1, const big_num *src2, big_num *target){
    int d = big_num_msb(src1) + big_num_msb(src2);
    d = DIV_ROUNDUP(d, 32) + !d;  // round up, min size = 1
    big_num *tmp;
    if (target == src1 || target == src2) {
        tmp = target;  // save c
        target = big_num_allocate(d);
    } else {
        tmp = NULL;
        for (int i = 0; i < target->size; i++)
            target->num[i] = 0;
        big_num_resize(target, d);
    }

    for (int i = 0; i < src1->size; i++) {
        for (int j = 0; j < src2->size; j++) {
            unsigned long long int carry = 0;
            carry = (unsigned long long int) src1->num[i] * src2->num[j];
            big_num_mult_add(target, i + j, carry);
        }
    }
    target->sign = src1->sign ^ src2->sign;

    if (tmp) {
        big_num_swap(tmp, target);  // restore c
        big_num_free(target);
    }
}
void big_num_fib_fast_doubling(big_num *dest, unsigned int n){
    big_num_resize(dest, 1);
    if (n < 2) {
        dest->num[0] = n;
        return;
    }
    big_num *f1 = dest, *f2 = big_num_allocate(1);
    f1->num[0] = 0;
    f2->num[0] = 1;
    big_num *t1 = big_num_allocate(1), *t2 = big_num_allocate(1);
    for (unsigned int i = 1U << 31; i; i >>= 1) {
        /* F(2t) = F(t) * [ 2 * F(t+1) – F(t) ] */
        big_num_cpy(t1, f2);
        big_num_left_shift(t1, 1);
        big_num_sub(t1, f1, t1);
        big_num_mul(t1, f1, t1);
        /* F(2k+1) = F(k)^2 + F(k+1)^2 */
        big_num_mul(f1, f1, f1);
        big_num_mul(f2, f2, f2);
        big_num_add(t2, f1, f2);
        if (n & i) {
            big_num_cpy(f1, t2);
            big_num_cpy(f2, t1);
            big_num_add(f2, t2, f2);
        } else {
            big_num_cpy(f1, t1);
            big_num_cpy(f2, t2);
        }
    }
    big_num_free(f2);
    big_num_free(t1);
    big_num_free(t2);
}
void big_num_fib(big_num *dest, unsigned int n){
    big_num_resize(dest, 1);
    if (n < 2) {  // Fib(0) = 0, Fib(1) = 1
        dest->num[0] = n;
        return;
    }
    big_num *a = big_num_allocate(1);
    big_num *b = big_num_allocate(1);
    dest->num[0] = 1;
    for (unsigned int i = 1; i < n; i++) {
        big_num_cpy(b, dest);
        big_num_add(dest, a, dest);
        big_num_swap(a, b);
    }
    big_num_free(a);
    big_num_free(b);
}
