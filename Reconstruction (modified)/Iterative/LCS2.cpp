#include "LCS2.hpp"
#include "CommonSubstring2.hpp"
#include "CommonSubstring2.cpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <unordered_map>
using namespace std;

const int NA = -1;
const int hashMod = 1e9 + 7;
const int hashPrime = 31;

LCS2::LCS2(const string& string1, const string& string2) :
		string1(string1), string2(string2), commonSubstrings(), commonSubstringsIndex(-1), commonActive(false) {

}

void LCS2::AddLetter(const char& letter, const int firstIndex, const int secondIndex) {
	if (letter == '_') {
		commonActive = false;
	}
	else { // letter is from GTAC
		if (commonActive) { // There is an active common substring. just add to it.
			commonSubstrings[commonSubstringsIndex].AddLetter(letter);
		}
		else { // no active common substrings. Open a new one.
			commonSubstrings.push_back(CommonSubstring2(letter, firstIndex, secondIndex));
			commonSubstringsIndex++;
			commonActive = true;
		}
	}
}

int LCS2::Len() const {
	int totalLen = 0;
	for (vector<CommonSubstring2>::const_iterator it = commonSubstrings.begin(); it != commonSubstrings.end(); it++) {
		totalLen += it->Len();
	}
	return totalLen;
}

vector<vector<int> > LCS2::Edges() const {
	vector<vector<int> > edges;
	for (vector<CommonSubstring2>::const_iterator it = commonSubstrings.begin(); it != commonSubstrings.end(); it++) {
		vector<int> current(3);
		int start1, end1, start2;
		it->StartEnd(&start1, &end1, &start2);
		current[0] = start1;
		current[1] = end1;
		current[2] = start2;
		edges.push_back(current);
	}
	return edges;
}

vector<int> LCS2::String1LCSIndexes() const {
	vector<int> string1Indexes(string1.size());
	for (vector<CommonSubstring2>::const_iterator it = commonSubstrings.begin(); it != commonSubstrings.end(); it++) {
		pair<int, int> range = it->Range1();
		for (int index = range.first; index < range.second; index++) {
			string1Indexes[index] = 1;
		}
	}
	return string1Indexes;
}

vector<int> LCS2::String1Mirrors() const {
	int start1, end1, start2;
	int vectorSize = string1.size();
	vector<vector<int> > edges = Edges();
	vector<int> mirrors(vectorSize, NA);
	for (vector<vector<int> >::iterator it = edges.begin(); it != edges.end(); it++) {
		start1 = (*it)[0];
		end1 = (*it)[1];
		start2 = (*it)[2];
		for (int index = start1; index < end1; index++) {
			mirrors[index] = start2 + index - start1;
		}
	}
	return mirrors;
}

vector<string> LCS2::DelStringGaps(const vector<int>& string1Mirrors) const {
	int vectorSize = string1.size() + 1;
	vector<string> string1MirrorsGaps(vectorSize);
	int lastIndex = -1;
	for (unsigned int index = 0; index < string1Mirrors.size(); index++) {
		if (string1Mirrors[index] != NA) {
			int startPos = lastIndex + 1;
			int substrLen = string1Mirrors[index] - startPos;
			string1MirrorsGaps[index] = string2.substr(startPos, substrLen);
			lastIndex = string1Mirrors[index];
		}
	}

	// End Gap
	int startPos = lastIndex + 1;
	string1MirrorsGaps[vectorSize - 1] = string2.substr(startPos);

	return string1MirrorsGaps;
}

string LCS2::RepStringWithLowercaseGap(const int patternLen, const int index, const vector<int>& string1Mirrors,
		const vector<string>& string1MirrorsGaps) const {
	assert(patternLen > 1);
	string posStr;
	// handle first letter in string1
	if (index == 0) {
		posStr += 'S';
		for (int j = 0; j < patternLen - 1; j++) {
			if (string1Mirrors[index + j] != NA) {
				string gap2 = string1MirrorsGaps[index + j];
				transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
				posStr += gap2;
				posStr += string1[index + j];
			}
			else {
				posStr += NOT_IN_LCS;
			}
		}
		return posStr;
	}

	int lettersToEnd = string1.size() - index;
	int maxPatternLen = ((patternLen - 1) < lettersToEnd) ? patternLen - 1 : lettersToEnd;
	for (int j = -1; j < maxPatternLen; j++) {
		if (string1Mirrors[index + j] != NA) {
			if (j > -1) {
				string gap2 = string1MirrorsGaps[index + j];
				transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
				posStr += gap2;
			}
			posStr += string1[index + j];
		}
		else {
			posStr += NOT_IN_LCS;
		}
	}

	if (patternLen - 1 > lettersToEnd) {
		string gap2 = string1MirrorsGaps[string1.size()];
		transform(gap2.begin(), gap2.end(), gap2.begin(), ::tolower);
		posStr += gap2;
		posStr += 'S';
	}
	return posStr;
}

string LCS2::RepStringWithNumberGap(const int patternLen, const int index, const vector<int>& string1Mirrors,
		const vector<string>& string1MirrorsGaps) const {
	assert(patternLen > 1);
	string posStr;
	// handle first letter in string1
	if (index == 0) {
		posStr += 'S';
		for (int j = 0; j < patternLen - 1; j++) {
			if (string1Mirrors[index + j] != NA) {
				string gap2str = string1MirrorsGaps[index + j];
				string gap2 = to_string(gap2str.size());
				posStr += gap2;
				posStr += string1[index + j];
			}
			else {
				posStr += NOT_IN_LCS;
			}
		}
		return posStr;
	}

	int lettersToEnd = string1.size() - index;
	int maxPatternLen = ((patternLen - 1) < lettersToEnd) ? patternLen - 1 : lettersToEnd;
	for (int j = -1; j < maxPatternLen; j++) {
		if (string1Mirrors[index + j] != NA) {
			if (j > -1) {
				string gap2str = string1MirrorsGaps[index + j];
				string gap2 = to_string(gap2str.size());
				posStr += gap2;
			}
			posStr += string1[index + j];
		}
		else {
			posStr += NOT_IN_LCS;
		}
	}

	if (patternLen - 1 > lettersToEnd) {
		string gap2str = string1MirrorsGaps[string1.size()];
		string gap2 = to_string(gap2str.size());
		posStr += gap2;
		posStr += 'S';
	}
	return posStr;
}

string LCS2::RepStringNoGap(const int patternLen, const int index, const vector<int>& string1Mirrors) const {
	assert(patternLen > 1);
	string posStr;
	// handle first letter in string1
	if (index == 0) {
		posStr += 'S';
		for (int j = 0; j < patternLen - 1; j++) {
			if (string1Mirrors[index + j] != NA) {
				posStr += string1[index + j];
			}
			else {
				posStr += NOT_IN_LCS;
			}
		}
		return posStr;
	}

	int lettersToEnd = string1.size() - index;
	int maxPatternLen = ((patternLen - 1) < lettersToEnd) ? patternLen - 1 : lettersToEnd;
	for (int j = -1; j < maxPatternLen; j++) {
		if (string1Mirrors[index + j] != NA) {
			posStr += string1[index + j];
		}
		else {
			posStr += NOT_IN_LCS;
		}
	}

	if (patternLen - 1 > lettersToEnd) {
		posStr += 'S';
	}
	return posStr;
}

vector<string> LCS2::RepStringWithLowercaseGapArray(const int patternLen) const {
	int vectorLen = string1.size() + 1;
	vector<int> string1Mirrors = String1Mirrors();
	vector<string> string1MirrorsGaps = DelStringGaps(string1Mirrors);
	vector<string> repArray;
	for (int index = 0; index < vectorLen; index++) {
		string posStr = RepStringWithLowercaseGap(patternLen, index, string1Mirrors, string1MirrorsGaps);
		repArray.push_back(posStr);
	}
	return repArray;
}

vector<string> LCS2::RepStringWithNumberGapArray(const int patternLen) const {
	int vectorLen = string1.size() + 1;
	vector<int> string1Mirrors = String1Mirrors();
	vector<string> string1MirrorsGaps = DelStringGaps(string1Mirrors);
	vector<string> repArray;
	for (int index = 0; index < vectorLen; index++) {
		string posStr = RepStringWithNumberGap(patternLen, index, string1Mirrors, string1MirrorsGaps);
		repArray.push_back(posStr);
	}
	return repArray;
}

vector<string> LCS2::RepStringNoGapArray(const int patternLen) const {
	int vectorLen = string1.size() + 1;
	vector<int> string1Mirrors = String1Mirrors();
	vector<string> repArray;
	for (int index = 0; index < vectorLen; index++) {
		string posStr = RepStringNoGap(patternLen, index, string1Mirrors);
		repArray.push_back(posStr);
	}
	return repArray;
}

std::ostream& operator<<(std::ostream& os, const LCS2& a) {
	if (a.commonSubstrings.empty()) {
		return os;
	}
	for (unsigned commonIndex = 0; commonIndex < a.commonSubstrings.size() - 1; commonIndex++) {
		os << a.commonSubstrings[commonIndex];
		os << "_";
	}
	os << a.commonSubstrings[a.commonSubstrings.size() - 1];
	return os;
}


// //fast exponentiation under a mod using binary exponentiation
// int fastModExp(int x, int y, int mod){
// 	int result = 1;
// 	int currPower = x;
// 	for(int i = 0; i < 32; i++){
// 		if(y & 1) result = result * currPower % mod;
// 		currPower = currPower * currPower % mod;
// 		y >>= 1;
// 	}
// 	return result;
// }

// //find hashes of all length 'len' substrings of string X, and return a vector of pairs (hash, substring starting index)
// vector<pair<int, int>> FindAllSubstringHashes(const string &X, int len) {
// 	//initialise necessary constant
// 	const int primePowLen = fastModExp(hashPrime, len, hashMod);
	
// 	int currHash = 0;
// 	vector<pair<int, int>> result;
// 	//calculate and store the first hash
// 	for(int i = 0; i < len; i++){
// 		currHash = currHash * hashPrime % hashMod;
// 		currHash = currHash + (X[i] - 'A' + 1) % hashMod;
// 	}
// 	result.push_back({currHash, 0});
// 	//calculate for other substrings
// 	for(int i = len; i < (int)X.size(); i++){
// 		currHash = currHash * hashPrime % hashMod;
// 		currHash = currHash + (X[i] - 'A' + 1) % hashMod;
// 		currHash = currHash + hashMod - (X[i - len] * primePowLen % hashMod) % hashMod;
// 		result.push_back({currHash, i});
// 	}

// 	return result;
// }

// //compares substring hash lists of two different strings and check if they have a common substring hash
// pair<int, int> CheckForCommonHash(vector<pair<int, int>> hashList1, vector<pair<int, int>> hashList2){
// 	//concatenate the lists, mark the pairs from the second list, and sort the combined list
// 	int len1 = hashList1.size();
// 	for(pair<int, int> &el : hashList2){
// 		el.first += len1;
// 		hashList1.push_back(el);
// 	}
// 	sort(hashList1.begin(), hashList1.end());

// 	//iterate over consecutive pairs and check if there are elements with the same hash and from different lists
// 	for(int i = 1; i < (int)hashList1.size(); i++){
// 		if(hashList1[i-1].first == hashList1[i].first){
// 			if(hashList1[i-1].second < len1 && hashList1[i].second >= len1){
// 				return {hashList1[i-1].second, hashList1[i].second - len1};
// 			}
// 		}
// 	}
// 	//return (-1,-1) if no matching pair found
// 	return pair<int, int>(-1, -1);
// }

// //find LCS of two strings using binary search
// LCS2 FindLCSFast(const string &X, const string &Y){
// 	int n = X.size();
// 	int m = Y.size();
// 	int N = min(n, m);

// 	//variables to store result:
// 	pair<int, int> matchPosns = {-1, -1};
// 	int matchLength = 0;
// 	//initialise binary search pointers
// 	int l = -1, r = N + 1;
// 	//binary search
// 	while(r - l > 1){
// 		int mid = (l + r) / 2;
// 		auto xHashes = FindAllSubstringHashes(X, mid);
// 		auto yHashes = FindAllSubstringHashes(Y, mid);
// 		auto commonPosns = CheckForCommonHash(xHashes, yHashes);
// 		//if there's no common position, check for smaller lengths:
// 		if(commonPosns.first == -1 && commonPosns.second == -1){
// 			r = mid;
// 		} else {
// 			l = mid;
// 		}
// 		//save the result
// 		matchPosns = commonPosns;
// 		matchLength = mid;
// 	}

// 	//format the output
// 	LCS2 result = LCS2(X, Y);
// 	for(int i = 0; i < matchLength; i++){
// 		result.AddLetter(X[matchPosns.first + i], matchPosns.first + i, matchPosns.second + i);
// 	}
// 	return result;
// }

//match node required for Kuo-Cross
struct matchNode {
	int i; 
	int j;
	matchNode* ptr;

	matchNode(int i_, int j_, matchNode* ptr_){
		i = i_;
		j = j_;
		ptr = ptr_;
	}
};

//find LCS via Kuo-Cross algorithm:
LCS2 FindLCSFast(const string& X, const string& Y){
	int n = X.size();
	int m = Y.size();

	auto basePtr = new matchNode(-1, -1, nullptr);
	vector<matchNode*> links(n+1, basePtr);
	vector<int> T(m+1, m+1);

	unordered_map<char, vector<int>> XCharPosns, YCharPosns;
	for(int i = 0; i < n; i++) XCharPosns[X[i]].push_back(i);
	for(int i = 0; i < m; i++) YCharPosns[Y[i]].push_back(i);

	for(int i = 0; i < n; i++){
		int temp = -1;
		int k = 0;
		int r = 0;
		auto c = links[0];
		for(int posn : YCharPosns[X[i]]){
			if(posn > temp){
				k++;
				while (posn > T[k]) k++;
				temp = T[k];
				T[k] = posn;
				auto prev = k > 0 ? links[k-1] : basePtr;
				links[r] = c;
				r = k;
				c = new matchNode(i, posn, prev);
			}
		}
		links[r] = c;
	}

	//find LCS length:
	int lcsLength = 0;
	for(int i = 0; i < m + 1; i++){
		if(T[i] != m + 1) lcsLength = i;
	}
	
	//retrieve the LCS:
	auto currNode = links[lcsLength];
	vector<int> lcsPosnsX, lcsPosnsY;
	while(currNode->ptr){
		lcsPosnsX.push_back(currNode->i);
		lcsPosnsY.push_back(currNode->j);
		currNode = currNode->ptr;
	}
	reverse(lcsPosnsX.begin(), lcsPosnsX.end());
	reverse(lcsPosnsY.begin(), lcsPosnsY.end());

	LCS2 res(X, Y);
	for(int i = 0; i < lcsLength; i++){
		res.AddLetter(X[lcsPosnsX[i]], lcsPosnsX[i], lcsPosnsY[i]);
		//res.AddLetter('A', i, i);
	}
	return res;
}


// find one LCS
LCS2 FindLCS(const string& X, const string& Y, int m, int n, const vector<vector<int> >& L) {

// If we reach end of either string, return a vector with empty LCS2
	if (m == 0 || n == 0) {
		return LCS2(X, Y);
	}

// If the last characters of X and Y are the same
	if (X[m - 1] == Y[n - 1]) {
		// recurse for X[0..m-2] and Y[0..n-2] in the matrix
		LCS2 lcs2 = FindLCS(X, Y, m - 1, n - 1, L);

		// append current character to LCS of substring X[0..m-2] and Y[0..n-2].
		lcs2.AddLetter(X[m - 1], m - 1, n - 1);
		return lcs2;
	}

// If the last characters of X and Y are not same
	else {

		// If LCS can be constructed from top side of the matrix, recurse for X[0..m-2] and Y[0..n-1]
		if (L[m - 1][n] >= L[m][n - 1]) {
			LCS2 lcs2 = FindLCS(X, Y, m - 1, n, L);
			lcs2.AddLetter('_', m - 1, n);
			return lcs2;
		}
		// If LCS can be constructed from left side of the matrix, recurse for X[0..m-1] and Y[0..n-2]
		else {
			LCS2 lcs2 = FindLCS(X, Y, m, n - 1, L);
			lcs2.AddLetter('_', m, n - 1);
			return lcs2;
		}
	}
}

/* Returns length of LCS for X[0..m-1], Y[0..n-1] */
int LCSLength(const string& X, const string& Y, int m, int n, vector<vector<int> >& L) {
// Build L[m+1][n+1] in bottom up fashion
	for (int i = 0; i <= m; i++) {
		for (int j = 0; j <= n; j++) {
			if (i == 0 || j == 0)
				L[i][j] = 0;
			else if (X[i - 1] == Y[j - 1])
				L[i][j] = L[i - 1][j - 1] + 1;
			else
				L[i][j] = max(L[i - 1][j], L[i][j - 1]);
		}
	}
	return L[m][n];
}

LCS2 ComputeLCS2(const string& X, const string& Y) {

	int m = X.length(), n = Y.length();
	vector<vector<int> > L(m + 1, vector<int>(n + 1));

	LCSLength(X, Y, m, n, L);
	return FindLCS(X, Y, m, n, L);
}

//modified version:
// LCS2 ComputeLCS2(const string& X, const string& Y) {
// 	return FindLCSFast(X, Y);
// }

int ComputeEditDistWithoutReplace(const string& X, const string& Y) {

	int m = X.length(), n = Y.length();
	vector<vector<int> > L(m + 1, vector<int>(n + 1));

	int lcsLen = LCSLength(X, Y, m, n, L);
	return X.size() - lcsLen + Y.size() - lcsLen;
}

