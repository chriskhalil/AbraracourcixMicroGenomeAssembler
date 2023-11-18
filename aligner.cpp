#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <set>
#include <chrono>
#include <fstream>
#include <sstream>
#include <utility>
#include <omp.h>
#include <regex>
#include <string_view>
#include "poopy-aligner.cpp"

using std::string;
using std::vector;
using std::string_view;
using std::cout;
using std::unordered_map;
using std::endl;


class SuffixArray{
    private:
    public:
    vector<long> sortCharacters(string_view S){
        size_t l = S.size();
        vector<long> order(l,0); //init vector of size l to 0
        unordered_map<char, long> count;

        for(size_t i{0};i<l;++i){
            count[S[i]]= count.find(S[i]) != count.end() ? count[S[i]] + 1 : 1;
        }
        std::vector<char> charList;
        for (const auto& kv : count) {
            charList.push_back(kv.first);
        }
        std::sort(charList.begin(), charList.end());

        char prevChar = charList[0];
        for (size_t i = 1; i < charList.size(); ++i) {
            count[charList[i]] += count[prevChar];
            prevChar = charList[i];
        }

        for (long i = l - 1; i >= 0; --i) {
            char c = S[i];
            count[c] = count[c] - 1;
            order[count[c]] = i;
        }

        return order;

    }
    std::vector<long> computeCharClasses(string_view S, const std::vector<long>& order) {
    long l = S.size();
    std::vector<long> charClass(l, 0);
    charClass[order[0]] = 0;

    for (long i = 1; i < l; ++i) {
        if (S[order[i]] != S[order[i - 1]]) {
            charClass[order[i]] = charClass[order[i - 1]] + 1;
        } else {
            charClass[order[i]] = charClass[order[i - 1]];
            }
        }
        return charClass;
    }

    std::vector<long> sortDoubled(string_view S, long L, const std::vector<long>& order, const std::vector<long>& _class) {
    long sLen = S.length();
    std::vector<long> count(sLen, 0);
    std::vector<long> newOrder(sLen, 0);

    for (long i = 0; i < sLen; ++i) {
        count[_class[i]] += 1;
    }

    for (long j = 1; j < sLen; ++j) {
        count[j] += count[j - 1];
    }

    for (long i = sLen - 1; i >= 0; --i) {
        long start = (order[i] - L + sLen) % sLen;
        long cl = _class[start];
        count[cl] -= 1;
        newOrder[count[cl]] = start;
    }

    return newOrder;
}

std::vector<long> updateClasses(const std::vector<long>& newOrder, const std::vector<long>& _class, long L) {
    long n = newOrder.size();
    std::vector<long> newClass(n, 0);
    newClass[newOrder[0]] = 0;

    for (long i = 1; i < n; ++i) {
        long curr = newOrder[i];
        long prev = newOrder[i - 1];
        long mid = curr + L;
        long midPrev = (prev + L) % n;

        if (_class[curr] != _class[prev] || _class[mid] != _class[midPrev]) {
            newClass[curr] = newClass[prev] + 1;
        } else {
            newClass[curr] = newClass[prev];
        }
    }

    return newClass;
   }
    vector<long> buildSuffixArray(string_view S){
       size_t sLen= S.size();
       auto order= sortCharacters(S);
       auto _class=computeCharClasses(S,order);
       size_t L=1;
       while(L < sLen){
        order =sortDoubled(S,L,order,_class);
        _class=updateClasses(order,_class,L);
        L=2*L;
       }
       return order;

    }
};
struct Occurrence {
    long position;
    int score;
    string pattern;

    Occurrence(long pos, int scr, const string& pat) : position(pos), score(scr), pattern(pat) {}
};
struct PatPos{
    long position;
    string pattern;

    PatPos(long pos, const string& pat) : position(pos), pattern(pat) {}
};
struct Scoring{
    unsigned int match;
    unsigned int mismatch; 
    unsigned int gap;
    unsigned int extend;
    Scoring(const unsigned int mat, const unsigned int mis, const unsigned int gay, const unsigned int ext) : match(mat), mismatch(mis), gap(gay), extend(ext) {}
};
/// @brief ////////////////////////////////////
class ApproxPatternMatching {
public:
    void run(string_view text,vector<string> patterns,int d, const unsigned int threads = 1, const unsigned int match = 1, const unsigned int mismatch = 4, const unsigned int gap = 6, const unsigned int extend = 1) {
        const Scoring score = Scoring(match, mismatch, gap, extend);

        std::cout << "+-----------------------+" << endl;
        std::cout << "| Reference Genome size |" << text.size() << endl;
        std::cout << "| Patterns count        |" << patterns.size() << endl;
        std::cout << "| Allowed mismatches    |" << d << endl;
        std::cout << "| Assigned threads      |" << threads << endl;
        std::cout << "+-----------------------+" << endl;

        auto startT = std::chrono::high_resolution_clock::now();

        // auto comparator = [](const std::vector<long>& a, const std::vector<long>& b) {return a.size() < b.size();};
        unordered_map<string, vector<long>> perfectScores; // dictionary to store the indexes of patterns that match perfectly

        // get the index and pattern occurences (score computed inside)
        vector<Occurrence> occs = findOccurrences(text, patterns, d, perfectScores, threads, score);
        cout<<"poopy occur"<<endl;
        auto endT = std::chrono::high_resolution_clock::now();

        cout<<"Num of non-perfect occurences:"<<occs.size()<<endl;
        cout<<endl;
        cout << "Elapsed Time to Find All Occurences: " << std::chrono::duration_cast<std::chrono::milliseconds>(endT - startT).count() << " ms" << endl;

        // build dictionary to keep track of how many copies of a read exist
        unordered_map<string, unsigned int> patternDict = getReadDict(patterns);
        cout<<"\nNumber of Unique Patterns: "<<patternDict.size()<<endl;
        cout<<"poopy dict"<<endl;
        // Aligning reads to reference
        vector<char> finalStringVector = alignPerfect(perfectScores, patternDict, text.length(), threads);
        cout<<"poopy perfect"<<endl;
        alignString(finalStringVector, occs, patternDict, threads);
        cout<<"poopy align"<<endl;
        string finalString(finalStringVector.begin(),finalStringVector.end());
        // alignString(text, occs, patternDict, threads);
        endT = std::chrono::high_resolution_clock::now();
        std::cout << "Total Elapsed Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(endT - startT).count() << " ms" << endl;
        cout<<endl;
        std::ofstream finalFile("output.txt");
        //finalFile<<text<<endl;
        finalFile<<finalString;
        finalFile.close();
        cout<<"Saved string to output.txt"<<endl;
    }

private:
    vector<long> bwtFromSuffixArray(string_view text, const std::vector<long>& order, const std::vector<char>& alphabet,std::unordered_map<char, std::vector<long>>& counts,std::unordered_map<char, long>& starts);
    bool approxMatching(string_view p1, string_view p2, int d);
    vector<Occurrence> findOccurrences(string_view text, const vector<string>& patterns, int d, unordered_map<string, vector<long>>& perfectScores, const unsigned int threads, const Scoring& scorer);
    void alignString(vector<char>& finalStringVector, vector<Occurrence>& orderedOcc, unordered_map<string, unsigned int>& patternDict, const unsigned int threads);
    unordered_map<string, unsigned int> getReadDict(vector<string> patterns);
    vector<char> alignPerfect(unordered_map<string, vector<long>>& perfectDict, unordered_map<string, unsigned int>& patternDict, const unsigned long textSize, const unsigned int threads);
    int computeScore(string_view text, const string &pattern, int d, const Scoring& score);
};

unordered_map<string, unsigned int> ApproxPatternMatching::getReadDict(vector<string> patterns){
    unordered_map<string, unsigned int> readDict;
    for(string str : patterns){
        if (readDict.find(str) == readDict.end()) {
            readDict[str] = 1;
        } else {
            readDict[str]++;
        }
    }
    return readDict;
}
vector<long> ApproxPatternMatching::bwtFromSuffixArray(string_view text, const vector<long>& order, const vector<char>& alphabet,std::unordered_map<char, vector<long>>& counts,std::unordered_map<char, long>& starts) {
    long l = text.size();
    vector<long> bwt(l);
    for (int i = 0; i < l; ++i) {
        bwt[i] = text[(order[i] + l - 1) % l];
    }

    for (char ch : alphabet) {
        counts[ch] = vector<long>(l + 1, 0);
    }

    for (long i = 0; i < l; ++i) {
        char currChar = bwt[i];
        for (const auto& entry : counts) {
            counts[entry.first][i + 1] = counts[entry.first][i];
        }
        counts[currChar][i + 1] += 1;
    }

    long currIndex = 0;
    for (char ch : alphabet) {
        starts[ch] = currIndex;
        currIndex += counts[ch][l];
    }
    return bwt;
}
bool ApproxPatternMatching::approxMatching(string_view p1, string_view p2, const int d) {
    int error = 0;
    for (size_t i = 0; i < p1.size(); ++i) {
        if (p1[i] != p2[i]) {
            error += 1;
            if (error > d) {
                return false;
            }
        }
    }
    return true;
}
int ApproxPatternMatching::computeScore(string_view text, const string &pattern, int d, const Scoring& scorer) {
    int score = 0;

    //+1 for match, -4 for mismatch, -6 for gap open, -1 for gap extension
    for (size_t i = 0; i < pattern.size(); ++i) {
        if (text[i] == pattern[i]) {
            score += scorer.match; // Match
        } else {
            score -= scorer.mismatch; // Mismatch
        }
    }

    int gapCount = 0;
    for (char ch : pattern) {
        if (ch == ' ') {
            gapCount += 1;
        }
    }

    score -= scorer.gap * std::min(gapCount, d); // Gap open
    score -= scorer.extend * std::max(0, gapCount - d); // Gap extension

    return score;
}
string negativeIndexSubstring(const std::string& str, int start, int off) {
    // Calculate positive indices from negative indices
    int len = str.length();
    int startPos = (start >= 0) ? start : len + start; // if negative go backwards
    int endPos = (startPos+off)%len;
    if (startPos + off >= len){ // need to wrap around
        std::string temp = str.substr(startPos,len-1);
        temp += str.substr(0,endPos);
        return temp;
    }
    else{
        return str.substr(startPos, off);
    }
}

size_t NegIndexToPos(const size_t str_size, long long index) {
    // Handle negative indexing
    if (index < 0) {
        index = str_size + index;
    }
    // Check if the index is within bounds
    if (index >= 0 && index < static_cast<long long>(str_size)) {
        return index;
    } else {
        // Handle index out of bounds (you may want to throw an exception or handle it differently)
        std::cerr << "Index out of bounds\n";
        return '\0';  // Return a default value, you can modify this accordingly
    }
}
string specialSubstr(const string& str, long start, long end){
    long len = str.length();
    long startPos = (start >= 0) ? start : len + start;
    long endPos = (end >= 0) ? end : len + end;

    // Perform substring operation
    if (startPos < 0) startPos = 0;
    if (endPos > len) endPos = len;
    if (startPos > endPos) startPos = endPos;

    return str.substr(startPos, endPos - startPos);
}
vector<Occurrence> ApproxPatternMatching::findOccurrences(string_view text, const std::vector<std::string>& patterns, int d, unordered_map<string, vector<long>>& perfectScores, const unsigned int threads, const Scoring& scorer) {
    SuffixArray suffixArray;
    vector<long> order = suffixArray.buildSuffixArray(text);

    vector<char> alphabet = {'$', 'a', 'c', 'g', 't'}; //  check if the alphabet should be capetilizaed or not
    bool keep = false;
    for (unsigned int i=0; i<alphabet.size();i++){
        if(text[0] == alphabet[i]){
            keep = true;
            break;
        }
    }
    if (!keep){
        alphabet = {'$', 'A', 'C', 'G', 'T'};
    }

    std::unordered_map<char, vector<long>> counts;
    std::unordered_map<char, long> starts;

    std::vector<long> bwt = bwtFromSuffixArray(text, order, alphabet,counts,starts);
    std::vector<Occurrence> occs;
    
    #pragma omp parallel for num_threads(threads) // Parallelize the loop using OpenMP
    for (unsigned int p = 0; p < patterns.size(); ++p) {
        const std::string& pattern = patterns[p];
        std::set<long> currOccs;
        long n = pattern.size();
        long k = n / (d + 1);
        vector<std::pair<string, int>> seeds;
        for (long i = 0; i < d; ++i) {
            seeds.emplace_back(pattern.substr(i * k, k), i * k);
        }
        seeds.emplace_back(pattern.substr(d * k, n - d * k), d * k);
        for (const auto& seed : seeds) {
            long top = 0;
            long bottom = bwt.size() - 1;
            long currIndex = seed.first.size() - 1;
            while (top <= bottom) {
                if (currIndex >= 0) {
                    char symbol = seed.first[currIndex];
                    currIndex -= 1;
                    // cout<< symbol<<endl; /////////////////////////////////////////////////////////// for testing ////////////////
                    if (counts[symbol][bottom + 1] - counts[symbol][top] > 0) {
                        top = starts[symbol] + counts[symbol][top];
                        bottom = starts[symbol] + counts[symbol][bottom + 1] - 1;
                    } else {
                        break;
                    }
                } else {
                    for (long i = top; i <= bottom; ++i) {
                        currOccs.insert(order[i] - seed.second);
                    }
                    break;
                }
            }
        }


        for (long occ : currOccs) {
            string_view sub = text.substr(NegIndexToPos(text.size(),occ), n);
            if (approxMatching(sub, pattern, d)) {
                int score = computeScore(sub, pattern, d, scorer); // could be merged with approxMatching for improvement (maybe)
                if (score == pattern.length()){ // perfect match
                    #pragma omp critical(perfectScore) // creates a crititcal data block to ensure no race conditions occur
                    if (perfectScores.count(pattern)){ // check if vector with key already exists
                        perfectScores[pattern].push_back(occ);
                    }
                    else {
                        perfectScores[pattern] = {occ}; // create new vector
                    }
                }else{
                    #pragma omp critical(occsCritical) // creates a crititcal data block to ensure no race conditions occur
                    occs.emplace_back(Occurrence(occ, score, pattern));
                }
                
            }
        }
    }
    return occs;
}
vector<char> ApproxPatternMatching::alignPerfect(unordered_map<string, vector<long>>& perfectDict, unordered_map<string, unsigned int>& patternDict, const unsigned long textSize, const unsigned int threads){
    
    vector<char> conString(textSize-1, '-'); //  the vector that keeps track of the final string (-1 for the $)
    vector<string> sortedKeys; //  vector to store the patterns in descending order of how many copies exist
    for (const auto& entry : perfectDict) {
        sortedKeys.push_back(entry.first); // store the keys
    }
    std::sort(sortedKeys.begin(),sortedKeys.end(),[&patternDict](const std::string& a, const std::string& b){ // sort in descending order of how many copies exist
        return patternDict[a] > patternDict[b];
    });

    vector<PatPos> toWrite;
    // Assign indices to keys
    for (const std::string& key : sortedKeys) {
        vector<long>& indices = perfectDict[key];
        unsigned int count = patternDict[key];
        for (long int& index : indices) {
            if (count > 0 && index < textSize) {
                // Bind the index to the key
                toWrite.emplace_back(PatPos(index,key));
                count--;
            }
        }
        // Update the remaining usage count for the key
        patternDict[key] = count;
    }
    
    // write in parallel
    #pragma omp parallel for num_threads(threads) // Parallelize the loop using OpenMP
    for (unsigned int p = 0; p < toWrite.size(); ++p) {
        for (unsigned int j = 0; j < toWrite[p].pattern.length(); j++){ // go through the characters of the string
            conString[toWrite[p].position+j] = toWrite[p].pattern[j];
        }
    }

    return conString;
}
void ApproxPatternMatching::alignString(vector<char>& finalStringVector, vector<Occurrence>& orderedOcc, unordered_map<string, unsigned int>& patternDict, const unsigned int threads){
    unsigned long dashCount = 0;
    #pragma omp parallel for reduction(+:dashCount) num_threads(threads) // reduction(+:dashCount) clause ensures that each thread has its own private copy of dashCount and the final result is combined after all threads finish execution
    for (unsigned long i = 0; i < finalStringVector.size(); i++) {
        if (finalStringVector[i] == '-') {
            dashCount++;
        }
    }

    // sort the remaining patterns based on index then score
    std::sort(orderedOcc.begin(),orderedOcc.end(),[](const Occurrence& a,const Occurrence& b){ // sort based on index then by score
        if (a.position != b.position){
            return a.position < b.position;
        }else{
            return a.score < b.score;
        }
    });

    string pattern;
    //unsigned long lastDash = 0, dashIndex, dashPatIndex;
    for (unsigned long i = 0; i < orderedOcc.size() && dashCount > 0; i++){ // keep adding till out of empty gaps
        // get soonest dash
        //auto it = std::find(finalStringVector.begin() + lastDash, finalStringVector.end(), '-'); // get first occurrence of '-' that we already know doesn't exist in the first lastDash indices
        //lastDash = dashIndex;
        //dashIndex = std::distance(finalStringVector.begin(), it);

        pattern = orderedOcc[i].pattern;
        //  && orderedOcc[i].position <= dashIndex && (orderedOcc[i].position + pattern.length()-1) >= dashIndex
        if(patternDict[pattern] > 0){ // we are still allowed to use the pattern and this pattern can replace a -
            --patternDict[pattern]; // reduce how many times we can still use this pattern
            //dashPatIndex = dashIndex - orderedOcc[i].position;// get index to reduce redundant checking at the start
            for (unsigned long j = 0; j < pattern.length(); j++){ // try to copy 
                if(finalStringVector[orderedOcc[i].position + j] == '-'){
                    finalStringVector[orderedOcc[i].position + j] = pattern[j];
                    --dashCount;
                }
            }
        }
    }
}

vector<string> loadStrings(const string& filename, const bool& check = false){ // returns the reads then the reference as the last element to be popped out
    vector<string> reads;
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "Error opening file: " << filename<< std::endl;
        return reads;
    }
    std::string line;
    while (std::getline(f,line)){
        if (line.length() > 2){
            if (check){
                reads.push_back(std::regex_replace(line, std::regex("[\r\n\t]+$"), ""));  // check  for \r and remove it then add to vector
            }
            else{
                reads.push_back(line);
            }
        }
    }
    f.close();
    return reads;
}

string buildReference(const vector<string> refVector){
    string reference;
    //cout<<"building reference..."<<std::endl;
    for (unsigned int i=0; i<refVector.size();i++){
        reference += refVector[i];
    }
    reference += "$";
    //cout<<"done building reference"<<std::endl;
    return reference;
}

int main(){
    ApproxPatternMatching p;
    vector<string> reads;
    string reference;
    reads = loadStrings("tests/coronavirus_reads.txt");
    reference = buildReference(loadStrings("tests/coronavirus_genome.txt"));
    p.run(reference,reads,2,16);

    return 0;
};