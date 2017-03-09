# ifndef LFM_H
# define LFM_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <stdint.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include "murmurhash.h"



using namespace std;

class LFM{
public:
	LFM();
	void lfm_cnt_offer(const void *buf, uint32_t len);
	int64_t generate_permutation(int n);
	uint32_t lfm_cnt_card();
	void lfm_cnt_fini();
	vector< vector<int> > permutation;
private:
	uint8_t num_of_substr;
	uint32_t num_of_permutations;
	uint8_t hash_len;
	uint8_t substr_len;
	uint32_t *bitmap;

	

	void initialize_permutation(int num_of_substr);
	void generate_permutation(vector<int> prefix, vector<int> remain);

	void initialization();
	
	

	int16_t num_of_trail_zeros(uint32_t i, int16_t substr_len );
	//uint64_t murmurhash(void *buf, uint32_t len, uint32_t seed);
	uint32_t harmonic_mean(vector<int> nums);
	int factorial(int n);
};


# endif