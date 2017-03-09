# include "LFM.h"

LFM::LFM(){
	initialization();
}

void LFM::initialization(){
	num_of_substr = 8;
	hash_len = 32;

	substr_len = hash_len / num_of_substr;
	/*num_of_permutations = 1;
	for (int i = 1; i <= num_of_substr; i++){
		num_of_permutations *= i;
	}*/
	num_of_permutations = 258;
	initialize_permutation(num_of_substr);
	vector< vector<int> > tmp;
	//for (auto iter=permutation.begin(); iter!=permutation.end(); iter++)
	int fac = factorial(num_of_substr);
	int maxSize = permutation.size();

	/*for (auto iter1=permutation.cbegin(); iter1!=permutation.cend(); iter1++){
		for (auto iter2=(*iter1).cbegin(); iter2!=(*iter1).cend(); iter2++){
			cout << *iter2 << " ";
		}
		cout << endl;
	}*/

	for (uint32_t i=0; i<num_of_permutations; i++){//get 1024 permutations randomly
		int randN = rand()%maxSize;
		tmp.push_back(permutation[randN]);
	}
	permutation = tmp;
	num_of_permutations = permutation.size();

	bitmap = new uint32_t[num_of_permutations];
	memset(bitmap, 0, num_of_permutations*4);
}

void LFM::lfm_cnt_fini(){
	delete [] bitmap;
}

//generate permulations based on the number of substring.
//e.g. if num_of_substr = 3, there will be 6 permutations:
// 012 021 102 120 210 201
void LFM::initialize_permutation(int num_of_substr){
	vector<int> prefix, remain;
	for (int i=0; i<num_of_substr; i++){
		remain.push_back(i);
	}
	generate_permutation(prefix, remain);
}

//generate all permutations
void LFM::generate_permutation(vector<int> prefix, vector<int> remain){
	
	for (vector<int>::iterator it=remain.begin(); it!=remain.end(); it++){
		vector<int> next_prefix(prefix), next_remain(remain);
		next_prefix.push_back(*it);
		next_remain.erase(remove(next_remain.begin(),next_remain.end(),*it),next_remain.end()); 

		if (next_remain.size() == 0){
			permutation.push_back(next_prefix);
		}else{
			generate_permutation(next_prefix, next_remain);
		}
	}
}

// search the index of the leftmost 1 bit in each substring
int16_t LFM::num_of_trail_zeros(uint32_t x, int16_t substr_len){
	//if(!leftmostOne)
	//	i = ~i;// search the index of the leftmost 0

	if (x == 0)
		return -1;

    /*int16_t offset = 1;
	while ((i = i>>1) != 0){
		offset ++;
	}
	return substr_len - offset;*/

	x |= (x >> 1);
	x |= (x >> 2);
	x |= (x >> 4);
	x |= (x >> 8);
	x |= (x >> 16);

	x -= (x >> 1) & 0x55555555;
	x = ((x >> 2) & 0x33333333) + (x & 0x33333333);
	x = ((x >> 4) + x) & 0x0F0F0F0F;
	x += (x >> 8);
	x += (x >> 16);

	return (substr_len - x & 0x0000003F);

}

void LFM::lfm_cnt_offer(const void *buf, uint32_t len){
	//hash value 
    uint64_t x;
	uint32_t seed = -1;
    x = murmurhash3((const char*)buf, (uint32_t) strlen((const char *)buf), seed);
	// if substr_len is 8, then the mask will be ...000 1111 1111
	uint64_t mask = (1 << substr_len) - 1;

	// get the offset of each substr
	vector<int16_t> substr_offsets;
	uint64_t substr;
	int16_t tmp;
	for (uint8_t i = 0; i < num_of_substr; i ++){
		substr = x & mask;
		tmp = num_of_trail_zeros(substr, substr_len);// 64 -> 32 bit
		substr_offsets.push_back(tmp);
		x = x >> substr_len;
	}
	for (uint32_t i = 0; i < num_of_permutations; i ++){
		int offset = 0;
		for (int j = 0; j < num_of_substr; j ++){
			int id = permutation[i][j];
			if (substr_offsets[id] == -1)
				offset += substr_len;
			else {
				offset += substr_offsets[id];
				break;
			}
		}
		bitmap[i] = bitmap[i] | (1 << (31 - offset));
	}
}

uint32_t LFM::lfm_cnt_card(){
	vector<int> result;
	uint32_t anti;
	for(uint32_t i=0; i<num_of_permutations; i++){
		anti = ~bitmap[i];
		result.push_back(num_of_trail_zeros(anti, 32));
	}
	/*cout << "Results:" << endl;
	for(unsigned int i = 0; i<result.size(); i++)
		cout << i << ":"<<result[i] << "\t";*/
	uint32_t tmp = 0;
	/*tmp = 1 << harmonic_mean(result);
	cout << "(Alternative) Power of harmonic:" << tmp << endl;*/

	for(unsigned int i = 0; i<result.size(); i++){
		result[i] = 1 << result[i];
		//cout << "Result" << i << ":"<<result[i] << "\t";
	}
	tmp = harmonic_mean(result)*1.5;
	//cout << "Harmonic of power:" << tmp << endl;
	return tmp;//1.5 is a deviation
}


uint32_t LFM::harmonic_mean(vector<int> nums){
	int len = nums.size();
	vector<int>::iterator iter;
	uint32_t result = 0;
	double tmp = 0;
	for (iter=nums.begin();iter!=nums.end();iter++){
		tmp += 1 / (double)(*iter);
	}
	result = (uint32_t)( len / tmp);
	return result;
}

int LFM::factorial(int n){
	int m = 1;
	for(int i=1; i<=n; i++){
		m*=i;
	}
	return m;
}