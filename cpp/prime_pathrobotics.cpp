/*
   Author: V Kotteda 
   Date:   Apr 2022
   c++   -pthread -g   prime_path.cpp  -lgtest -lgtest_main -lgmock -lgmock_main
*/

#include <iostream>
#include<vector>
using namespace std;
#include "gtest/gtest.h"

void checkprime(int n, int m, vector<int>& vec){

	if(n==1 && m==1) {
		vec.push_back(1);
		return;
	}

	for (int j=n; j<=m; j++) {
		bool prime=true;
		if(j<2) {
			cout << "value" << j << endl;
		}
		else if(j >1) {
			for(int i=2; i<j-1; i++) {
				if(j%i==0)  {
					if(m == n) vec.push_back(0);
					prime= false;
					break;
				}
			}
		if(prime==true) vec.push_back(j);
		}
	
  	 }
   return;
}

TEST(test1, prime1)
{
        vector<int> vec;
        checkprime(1,1,vec);
        EXPECT_EQ(1, vec[0]);
}

TEST(test2, prime4)
{
        vector<int> vec;
        checkprime(4,4,vec);
        EXPECT_EQ(0, vec[0]);
}

int main(int argc, char **argv) {
	int n=0, m= 10;
	vector<int> vec1;
	checkprime(n,m,vec1);
	for(auto x:vec1) cout << x << endl; 

	testing::InitGoogleTest(&argc, argv);
	  return RUN_ALL_TESTS();

	return 0;
}
