#include "hyperloglog.h"
#include "LFM.h"
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>

using namespace std;

#define BUFFERSIZE 256



int main(int argc, char *argv[]){
	int methodType = -1;
	if( argc==2 && argv[1][0]=='0'){
		methodType = 0;
	}
	
	if( argc==2 && argv[1][0]=='1'){
		methodType = 1;
	}
	
	string filename = "zipf.txt";
	char buffer[256] = {};

	while(methodType != 0 && methodType != 1){
		cout << "Choose a Method:" << endl;
		cout << "LFM:0" << endl;
		cout << "hyperloglog:1" << endl;
		cin >> methodType;
		if( methodType == 1 || methodType == 0)
			break;
	}
	clock_t start, end;
	start = clock();
	ifstream in(filename);
	if(!in.is_open()){
		cout << "Open file failed." << endl;
		return -1;
	}
	
	if(methodType){
		Hyperloglog hll;
		//cout << "Begin processing date" << endl;
		int count = 0, n = 0;
		while (!in.eof() )  {
			count++;
			if (count == 1000){
				count = 0;
				n++;
				//cout << n << "000 date processed." << endl;
			}
			in.getline (buffer,BUFFERSIZE);
			hll.hll_cnt_offer(buffer, strlen(buffer));
			memset(buffer, 0, BUFFERSIZE);
		}
		cout << endl;
		int64_t result = hll.hll_cnt_card();
		cout << "HLL cardinality: " << result << endl;
		hll.hll_cnt_fini();
	}else{
		LFM lfm;
		/*for (auto iter1=lfm.permutation.cbegin(); iter1!=lfm.permutation.cend(); iter1++){
			for (auto iter2=(*iter1).cbegin(); iter2!=(*iter1).cend(); iter2++){
				cout << *iter2 << " ";
			}
			cout << endl;
		}
		cout << "Total number of permutations:" << lfm.permutation.size() << endl;*/
		//middle = clock();
		//cout << "Begin processing date" << endl;
		int count = 0, n = 0;
		while (!in.eof() )  {
			count++;
			//cout << count;
			if (count == 1000){
				count = 0;
				n++;
				//cout << n << "000 date processed." << endl;
			}
			in.getline (buffer,BUFFERSIZE);
			lfm.lfm_cnt_offer(buffer, strlen(buffer));
			memset(buffer, 0, BUFFERSIZE);
		}
		cout << endl;
		//cout << "Process over" << endl;
		uint32_t result = lfm.lfm_cnt_card();
		cout << "LFM carinality: " << result << endl;
		lfm.lfm_cnt_fini();
		/*double time=(middle-start)/(double)CLOCKS_PER_SEC;
		cout << time << "seconds used for generating permutations." << endl;*/
	}
	end = clock();
	double time=(end-start)/(double)CLOCKS_PER_SEC;
	cout << time << "seconds used totally." << endl;


	//system("pause");
	return 0;
}
