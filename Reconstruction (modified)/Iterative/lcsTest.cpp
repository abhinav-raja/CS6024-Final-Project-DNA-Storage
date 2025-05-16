#include "LCS2.cpp"
#include <bits/stdc++.h>
using namespace std;

string y = "AGCATGGTTA";
string x = "ATGGTCCATG";


int main(){
    LCS2 l1 = ComputeLCS2(x, y);
    cout << l1 << '\n';
    LCS2 l2 = FindLCSFast(x, y);
    cout << l2 << '\n';

    // for(const int el : l1.String1LCSIndexes()){
    //     cout << el << " ";
    // }
    // cout << '\n';

    cout << l1.Len() << '\n';
    cout << l2.Len() << '\n';
    return 0;
}