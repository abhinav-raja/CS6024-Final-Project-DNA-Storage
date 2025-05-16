#include <bits/stdc++.h>
#include "Clone.cpp"
#include "LCS2.cpp"
using namespace std;

//length of strands and number of pairs to create
const int len = 3000;
const int numPairs = 50;

int main(){

    //initialise the rng
    unsigned sd = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 rng(sd);

    // //create random strands for half the dataset
    vector<pair<string, string>> testData(numPairs);
    // for(int i = 0; i < numPairs; i++){
    //     testData[i].first = MakeStrand(len, rng);
    //     testData[i].second = MakeStrand(len, rng);
    // }

    //create copied strands for the other half
    for(int i = 0; i < numPairs; i++){
        testData[i].first = MakeStrand(len, rng);
        testData[i].second = CopyStrand(testData[i].first, rng, 0.1, 0.1, 0.1);
    }

    vector<int> resultsOrig(numPairs), resultsModified(numPairs);
    //run the original LCS algorithm
    unsigned startTimeOrig = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    for(int i = 0; i < numPairs; i++){
        auto x = ComputeLCS2(testData[i].first, testData[i].second);
        resultsOrig[i] = x.Len();
    }
    unsigned endTimeOrig = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    //run the modified algorithm
    unsigned startTimeModified = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();
    for(int i = 0; i < numPairs; i++){
        auto x = FindLCSFast(testData[i].first, testData[i].second);
        resultsModified[i] = x.Len();
    }
    unsigned endTimeModified = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now().time_since_epoch()).count();

    //print timings
    cout << "Original: " << (endTimeOrig - startTimeOrig)/1000.0 << "s\n";
    cout << "Modified: " << (endTimeModified - startTimeModified)/1000.0 << "s\n";


    //verify results:
    for(int i = 0; i < numPairs; i++){
        if(resultsOrig[i] != resultsModified[i]){
            cout << "Mismatch!!\n";
        }
    }
    return 0;
}