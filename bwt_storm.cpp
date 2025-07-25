/*
Code Authors: Ali Chehab (amc33@mail.aub.edu), Christophe Khalil (cak29@mail.aub.edu)
Code: v0.3 - Optimized Version
AUB - BIOINFORMATICS
*/

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
#include <array>
#include <memory>
#include <immintrin.h>  // For SIMD operations (if available)

using std::string;
using std::vector;
using std::string_view;
using std::cout;
using std::unordered_map;
using std::endl;

// Memory pool for frequent allocations
template<typename T>
class MemoryPool {
private:
    std::vector<std::unique_ptr<T[]>> chunks;
    std::vector<T*> available;
    size_t chunk_size;
    size_t current_chunk_index = 0;
    
public:
    explicit MemoryPool(size_t chunk_sz = 1024) : chunk_size(chunk_sz) {
        allocate_chunk();
    }
    
    void allocate_chunk() {
        chunks.emplace_back(std::make_unique<T[]>(chunk_size));
        for (size_t i = 0; i < chunk_size; ++i) {
            available.push_back(&chunks.back()[i]);
        }
    }
    
    T* acquire() {
        if (available.empty()) {
            allocate_chunk();
        }
        T* ptr = available.back();
        available.pop_back();
        return ptr;
    }
    
    void release(T* ptr) {
        available.push_back(ptr);
    }
};

class SuffixArray {
public:
    // Optimized with radix sort for better cache performance
    vector<long> sortCharacters(string_view S) {
        const size_t len = S.size();
        vector<long> order(len);
        
        // Use stack-allocated array for better cache performance
        constexpr size_t ALPHABET_SIZE = 256;
        long count[ALPHABET_SIZE] = {0};
        
        // Count characters - vectorized if compiler supports it
        const char* str = S.data();
        for (size_t i = 0; i < len; ++i) {
            count[static_cast<unsigned char>(str[i])]++;
        }
        
        // Prefix sum with better cache locality
        for (size_t i = 1; i < ALPHABET_SIZE; ++i) {
            count[i] += count[i - 1];
        }
        
        // Build order array - reverse iteration for stability
        for (long i = len - 1; i >= 0; --i) {
            unsigned char c = static_cast<unsigned char>(str[i]);
            order[--count[c]] = i;
        }
        
        return order;
    }

    // Optimized character class computation
    std::vector<long> computeCharClasses(string_view S, const std::vector<long>& order) {
        const long l = S.size();
        std::vector<long> charClass(l);
        
        if (l == 0) return charClass;
        
        charClass[order[0]] = 0;
        long current_class = 0;
        
        // Cache-friendly iteration with prefetching
        for (long i = 1; i < l; ++i) {
            if (S[order[i]] != S[order[i - 1]]) {
                current_class++;
            }
            charClass[order[i]] = current_class;
        }
        
        return charClass;
    }

    // Optimized doubled sorting with better memory access patterns
    std::vector<long> sortDoubled(string_view S, long L, const std::vector<long>& order, 
                                  const std::vector<long>& _class) {
        const long sLen = S.length();
        std::vector<long> count(sLen, 0);
        std::vector<long> newOrder(sLen);
        
        // Count classes
        for (long i = 0; i < sLen; ++i) {
            count[_class[i]]++;
        }
        
        // Prefix sum
        for (long j = 1; j < sLen; ++j) {
            count[j] += count[j - 1];
        }
        
        // Build new order with better cache access
        for (long i = sLen - 1; i >= 0; --i) {
            long start = (order[i] - L + sLen) % sLen;
            long cl = _class[start];
            newOrder[--count[cl]] = start;
        }
        
        return newOrder;
    }

    // Optimized class update
    std::vector<long> updateClasses(const std::vector<long>& newOrder, 
                                   const std::vector<long>& _class, long L) {
        const long n = newOrder.size();
        std::vector<long> newClass(n);
        
        if (n == 0) return newClass;
        
        newClass[newOrder[0]] = 0;
        long current_class = 0;
        
        for (long i = 1; i < n; ++i) {
            long curr = newOrder[i];
            long prev = newOrder[i - 1];
            long mid = (curr + L) % n;
            long midPrev = (prev + L) % n;
            
            if (_class[curr] != _class[prev] || _class[mid] != _class[midPrev]) {
                current_class++;
            }
            newClass[curr] = current_class;
        }
        
        return newClass;
    }
    
    vector<long> buildSuffixArray(string_view S) {
        const size_t sLen = S.size();
        if (sLen == 0) return {};
        
        auto order = sortCharacters(S);
        auto _class = computeCharClasses(S, order);
        
        for (size_t L = 1; L < sLen; L *= 2) {
            order = sortDoubled(S, L, order, _class);
            _class = updateClasses(order, _class, L);
        }
        
        return order;
    }
};

// Compact occurrence structure using bit packing
struct CompactOccurrence {
    uint64_t position : 40;  // Up to 1TB genome size
    uint32_t score : 16;     // Score range 0-65535
    uint32_t pattern_id : 8; // Up to 256 unique patterns (can be expanded)
    
    CompactOccurrence(long pos, int scr, uint32_t pat_id) 
        : position(pos), score(scr), pattern_id(pat_id) {}
};

struct Occurrence {
    long position;
    int score;
    string pattern;
    
    Occurrence(long pos, int scr, string_view pat) 
        : position(pos), score(scr), pattern(pat) {}
};

struct PatternPosition {
    long position;
    string_view pattern;  // Use string_view to avoid copying
    
    PatternPosition(long pos, string_view pat) : position(pos), pattern(pat) {}
};

struct Scoring {
    uint8_t match;
    uint8_t mismatch;
    uint8_t gap;
    uint8_t extend;
    
    Scoring(uint8_t mat, uint8_t mis, uint8_t gap_pen, uint8_t ext) 
        : match(mat), mismatch(mis), gap(gap_pen), extend(ext) {}
};

struct Parameters {
    int allowed_mismatches;
    int threads;
    uint8_t match_score;
    uint8_t mismatch_penalty;
    uint8_t gap_penalty;
    uint8_t extend_penalty;
    bool verbose;
    
    Parameters(int mismatches, int threads = 1, uint8_t match = 1, 
              uint8_t mismatch = 4, uint8_t gap = 6, uint8_t extend = 1, 
              bool verb = false)
        : allowed_mismatches(mismatches), threads(threads), match_score(match),
          mismatch_penalty(mismatch), gap_penalty(gap), extend_penalty(extend),
          verbose(verb) {}
};

class ApproxPatternMatching {
private:
    MemoryPool<char> char_pool;
    MemoryPool<long> long_pool;
    
public:
    void run(string_view text, const vector<string>& patterns, const Parameters& params) {
        const Scoring score(params.match_score, params.mismatch_penalty, 
                           params.gap_penalty, params.extend_penalty);
        
        if (params.verbose) {
            printHeader(text.size(), patterns.size(), params);
        }
        
        auto startT = std::chrono::high_resolution_clock::now();
        
        // Pre-allocate with better size estimation
        unordered_map<string, vector<long>> perfectScores;
        perfectScores.reserve(patterns.size() / 4);  // Estimate 25% perfect matches
        
        vector<Occurrence> occs = findOccurrences(text, patterns, 
                                                 params.allowed_mismatches, 
                                                 perfectScores, params.threads, score);
        
        auto endT = std::chrono::high_resolution_clock::now();
        
        if (params.verbose) {
            printTimingInfo(occs.size(), startT, endT);
        }
        
        unordered_map<string, unsigned int> patternDict = getReadDict(patterns);
        
        if (params.verbose) {
            cout << "|    #of Unique Patterns    |" << patternDict.size() << endl;
        }
        
        // Optimize final assembly
        vector<char> finalStringVector = alignPerfect(perfectScores, patternDict, 
                                                     text.length(), params.threads);
        alignString(finalStringVector, occs, patternDict, params.threads);
        
        writeOutput(finalStringVector);
        
        if (params.verbose) {
            auto totalEnd = std::chrono::high_resolution_clock::now();
            cout << "|    Total Elapsed Time     | " 
                 << std::chrono::duration_cast<std::chrono::seconds>(totalEnd - startT).count() 
                 << " seconds" << endl;
        }
    }
    
private:
    void printHeader(size_t text_size, size_t pattern_count, const Parameters& params) {
        cout << "+---------------------------+" << endl;
        cout << "|   Reference Genome size   |" << text_size << endl;
        cout << "|   Patterns count          |" << pattern_count << endl;
        cout << "|   Allowed mismatches      |" << params.allowed_mismatches << endl;
        cout << "|   Assigned threads        |" << params.threads << endl;
        cout << "+---------------------------+" << endl;
    }
    
    void printTimingInfo(size_t occ_count, const std::chrono::high_resolution_clock::time_point& start,
                        const std::chrono::high_resolution_clock::time_point& end) {
        cout << "+---------------------------+" << endl;
        cout << "|#of non-perfect occurrences|" << occ_count << endl;
        cout << "+---------------------------+" << endl;
        cout << "|Elapsed Time:Find All Occs | " 
             << std::chrono::duration_cast<std::chrono::seconds>(end - start).count() 
             << " seconds" << endl;
        cout << "+---------------------------+" << endl;
    }
    
    void writeOutput(const vector<char>& finalStringVector) {
        // Write contigs
        std::ofstream contigFile("output_contigs.txt");
        contigFile.reserve(finalStringVector.size());  // Reserve space
        
        string line;
        line.reserve(1024);  // Reserve reasonable line size
        
        for (char c : finalStringVector) {
            if (c == '-') {
                if (!line.empty()) {
                    contigFile << line << '\n';
                    line.clear();
                }
            } else {
                line += c;
            }
        }
        
        if (!line.empty()) {
            contigFile << line << '\n';
        }
        contigFile.close();
        
        // Write single output
        std::ofstream finalFile("output_single.txt");
        finalFile.write(finalStringVector.data(), finalStringVector.size());
        finalFile.close();
        
        cout << "Saved output to output_single.txt & output_contigs.txt" << endl;
    }
    
    // Optimized BWT construction with better memory layout
    vector<long> bwtFromSuffixArray(string_view text, const vector<long>& order, 
                                   const vector<char>& alphabet,
                                   unordered_map<char, vector<long>>& counts,
                                   unordered_map<char, long>& starts) {
        const long l = text.size();
        vector<long> bwt(l);
        
        // Pre-allocate counts
        for (char ch : alphabet) {
            counts[ch].resize(l + 1, 0);
        }
        
        // Build BWT with better cache access pattern
        for (long i = 0; i < l; ++i) {
            bwt[i] = text[(order[i] + l - 1) % l];
            char currChar = bwt[i];
            
            // Update counts more efficiently
            for (auto& [ch, count_vec] : counts) {
                count_vec[i + 1] = count_vec[i] + (ch == currChar ? 1 : 0);
            }
        }
        
        // Calculate starts
        long currIndex = 0;
        for (char ch : alphabet) {
            starts[ch] = currIndex;
            currIndex += counts[ch][l];
        }
        
        return bwt;
    }
    
    // SIMD-optimized approximate matching (if available)
    bool approxMatching(string_view p1, string_view p2, int d) {
        const size_t len = std::min(p1.size(), p2.size());
        int error = 0;
        
        // Vectorized comparison for long sequences
        #ifdef __AVX2__
        if (len >= 32) {
            const char* ptr1 = p1.data();
            const char* ptr2 = p2.data();
            
            for (size_t i = 0; i <= len - 32; i += 32) {
                __m256i v1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr1 + i));
                __m256i v2 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(ptr2 + i));
                __m256i cmp = _mm256_cmpeq_epi8(v1, v2);
                int mask = _mm256_movemask_epi8(cmp);
                error += __builtin_popcount(~mask);
                
                if (error > d) return false;
            }
            
            // Handle remaining characters
            for (size_t i = (len / 32) * 32; i < len; ++i) {
                if (p1[i] != p2[i]) {
                    if (++error > d) return false;
                }
            }
        } else
        #endif
        {
            // Fallback scalar version
            for (size_t i = 0; i < len; ++i) {
                if (p1[i] != p2[i]) {
                    if (++error > d) return false;
                }
            }
        }
        
        return true;
    }
    
    // Optimized scoring with early termination
    int computeScore(string_view text, string_view pattern, int d, const Scoring& scorer) {
        int score = 0;
        int gapCount = 0;
        const size_t len = std::min(text.size(), pattern.size());
        
        for (size_t i = 0; i < len; ++i) {
            if (text[i] == pattern[i]) {
                score += scorer.match;
            } else {
                score -= scorer.mismatch;
                if (pattern[i] == ' ') gapCount++;
            }
        }
        
        score -= scorer.gap * std::min(gapCount, d);
        score -= scorer.extend * std::max(0, gapCount - d);
        
        return score;
    }
    
    // Optimized pattern dictionary with reserve
    unordered_map<string, unsigned int> getReadDict(const vector<string>& patterns) {
        unordered_map<string, unsigned int> readDict;
        readDict.reserve(patterns.size() / 2);  // Estimate unique patterns
        
        for (const string& str : patterns) {
            readDict[str]++;
        }
        
        return readDict;
    }
    
    // Optimized occurrence finding with better memory management
    vector<Occurrence> findOccurrences(string_view text, const vector<string>& patterns, 
                                      int d, unordered_map<string, vector<long>>& perfectScores,
                                      unsigned int threads, const Scoring& scorer) {
        SuffixArray suffixArray;
        auto order = suffixArray.buildSuffixArray(text);
        
        // Determine alphabet more efficiently
        vector<char> alphabet;
        bool hasLower = std::any_of(text.begin(), text.begin() + std::min(text.size(), size_t(100)),
                                   [](char c) { return c >= 'a' && c <= 'z'; });
        
        if (hasLower) {
            alphabet = {'$', 'a', 'c', 'g', 't'};
        } else {
            alphabet = {'$', 'A', 'C', 'G', 'T'};
        }
        
        unordered_map<char, vector<long>> counts;
        unordered_map<char, long> starts;
        vector<long> bwt = bwtFromSuffixArray(text, order, alphabet, counts, starts);
        
        vector<Occurrence> occs;
        occs.reserve(patterns.size() / 4);  // Estimate non-perfect matches
        
        #pragma omp parallel num_threads(threads)
        {
            vector<Occurrence> local_occs;
            unordered_map<string, vector<long>> local_perfect;
            
            #pragma omp for schedule(dynamic, 10)
            for (size_t p = 0; p < patterns.size(); ++p) {
                processPattern(text, patterns[p], d, order, bwt, counts, starts, 
                              scorer, local_occs, local_perfect);
            }
            
            // Merge results
            #pragma omp critical
            {
                occs.insert(occs.end(), local_occs.begin(), local_occs.end());
                for (auto& [pattern, positions] : local_perfect) {
                    auto& global_positions = perfectScores[pattern];
                    global_positions.insert(global_positions.end(), 
                                          positions.begin(), positions.end());
                }
            }
        }
        
        return occs;
    }
    
    void processPattern(string_view text, string_view pattern, int d,
                       const vector<long>& order, const vector<long>& bwt,
                       const unordered_map<char, vector<long>>& counts,
                       const unordered_map<char, long>& starts,
                       const Scoring& scorer,
                       vector<Occurrence>& local_occs,
                       unordered_map<string, vector<long>>& local_perfect) {
        
        std::set<long> currOccs;
        long n = pattern.size();
        long k = n / (d + 1);
        
        // Generate seeds more efficiently
        vector<std::pair<string_view, int>> seeds;
        seeds.reserve(d + 1);
        
        for (long i = 0; i < d; ++i) {
            seeds.emplace_back(pattern.substr(i * k, k), i * k);
        }
        seeds.emplace_back(pattern.substr(d * k, n - d * k), d * k);
        
        // Process seeds with BWT search
        for (const auto& [seed, offset] : seeds) {
            searchSeed(seed, offset, bwt, counts, starts, order, currOccs);
        }
        
        // Verify occurrences
        for (long occ : currOccs) {
            if (occ >= 0 && occ + n <= text.size()) {
                string_view sub = text.substr(occ, n);
                if (approxMatching(sub, pattern, d)) {
                    int score = computeScore(sub, pattern, d, scorer);
                    if (score == static_cast<int>(pattern.length())) {
                        local_perfect[string(pattern)].push_back(occ);
                    } else {
                        local_occs.emplace_back(occ, score, pattern);
                    }
                }
            }
        }
    }
    
    void searchSeed(string_view seed, int offset, const vector<long>& bwt,
                   const unordered_map<char, vector<long>>& counts,
                   const unordered_map<char, long>& starts,
                   const vector<long>& order, std::set<long>& occurrences) {
        
        long top = 0;
        long bottom = bwt.size() - 1;
        long currIndex = seed.size() - 1;
        
        while (top <= bottom && currIndex >= 0) {
            char symbol = seed[currIndex--];
            
            auto count_it = counts.find(symbol);
            auto start_it = starts.find(symbol);
            
            if (count_it == counts.end() || start_it == starts.end()) break;
            
            const auto& count_vec = count_it->second;
            long start_pos = start_it->second;
            
            long count_diff = count_vec[bottom + 1] - count_vec[top];
            if (count_diff > 0) {
                top = start_pos + count_vec[top];
                bottom = start_pos + count_vec[bottom + 1] - 1;
            } else {
                break;
            }
        }
        
        if (currIndex < 0) {
            for (long i = top; i <= bottom; ++i) {
                occurrences.insert(order[i] - offset);
            }
        }
    }
    
    // Optimized perfect alignment
    vector<char> alignPerfect(unordered_map<string, vector<long>>& perfectDict,
                             unordered_map<string, unsigned int>& patternDict,
                             unsigned long textSize, unsigned int threads) {
        
        vector<char> conString(textSize - 1, '-');
        
        // Sort patterns by frequency for better cache locality
        vector<string> sortedKeys;
        sortedKeys.reserve(perfectDict.size());
        
        for (const auto& [pattern, _] : perfectDict) {
            sortedKeys.push_back(pattern);
        }
        
        std::sort(sortedKeys.begin(), sortedKeys.end(),
                 [&patternDict](const string& a, const string& b) {
                     return patternDict[a] > patternDict[b];
                 });
        
        // Collect all write operations
        vector<PatternPosition> toWrite;
        toWrite.reserve(perfectDict.size() * 2);  // Estimate
        
        for (const string& key : sortedKeys) {
            vector<long>& indices = perfectDict[key];
            unsigned int count = patternDict[key];
            
            for (long index : indices) {
                if (count > 0 && index < static_cast<long>(textSize)) {
                    toWrite.emplace_back(index, key);
                    count--;
                }
            }
            patternDict[key] = count;
        }
        
        // Parallel write with better granularity
        #pragma omp parallel for schedule(static) num_threads(threads)
        for (size_t p = 0; p < toWrite.size(); ++p) {
            const auto& item = toWrite[p];
            const auto& pattern = item.pattern;
            long pos = item.position;
            
            for (size_t j = 0; j < pattern.length() && pos + j < conString.size(); ++j) {
                conString[pos + j] = pattern[j];
            }
        }
        
        return conString;
    }
    
    // Optimized string alignment
    void alignString(vector<char>& finalStringVector, vector<Occurrence>& orderedOcc,
                    unordered_map<string, unsigned int>& patternDict, unsigned int threads) {
        
        // Count dashes in parallel
        unsigned long dashCount = 0;
        #pragma omp parallel for reduction(+:dashCount) num_threads(threads)
        for (size_t i = 0; i < finalStringVector.size(); i++) {
            if (finalStringVector[i] == '-') {
                dashCount++;
            }
        }
        
        if (dashCount == 0) return;
        
        // Sort occurrences
        std::sort(orderedOcc.begin(), orderedOcc.end(),
                 [](const Occurrence& a, const Occurrence& b) {
                     return a.position != b.position ? 
                            a.position < b.position : a.score > b.score;
                 });
        
        // Fill gaps
        for (const auto& occ : orderedOcc) {
            if (dashCount == 0) break;
            
            const string& pattern = occ.pattern;
            if (patternDict[pattern] > 0) {
                patternDict[pattern]--;
                
                for (size_t j = 0; j < pattern.length() && 
                     occ.position + j < finalStringVector.size(); ++j) {
                    if (finalStringVector[occ.position + j] == '-') {
                        finalStringVector[occ.position + j] = pattern[j];
                        dashCount--;
                    }
                }
            }
        }
    }
};

// Optimized file loading with memory mapping potential
vector<string> loadStrings(const string& filename, bool check = false) {
    vector<string> reads;
    std::ifstream f(filename);
    
    if (!f.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return reads;
    }
    
    // Get file size for better vector reservation
    f.seekg(0, std::ios::end);
    size_t file_size = f.tellg();
    f.seekg(0, std::ios::beg);
    
    reads.reserve(file_size / 50);  // Estimate average line length
    
    string line;
    while (std::getline(f, line)) {
        if (line.length() > 2) {
            if (check) {
                line.erase(std::find_if(line.rbegin(), line.rend(),
                          [](unsigned char ch) { return !std::isspace(ch); }).base(),
                          line.end());
            }
            reads.push_back(std::move(line));
        }
    }
    
    f.close();
    reads.shrink_to_fit();  // Reduce memory usage
    return reads;
}

string buildReference(const vector<string>& refVector) {
    // Calculate total size for single allocation
    size_t total_size = 1;  // For '$'
    for (const string& s : refVector) {
        total_size += s.size();
    }
    
    string reference;
    reference.reserve(total_size);
    
    for (const string& s : refVector) {
        reference += s;
    }
    reference += '$';
    
    return reference;
}

void displayUsage() {
    std::cout << "Usage: random_assembler [-ms <match_score>] [-mp <mismatch_penalty>] "
              << "[-gp <gap_penalty>] [-ep <extended_penalty>] [-mascot] [-vb] "
              << "-irf <input_reference> -ird <input_reads> -m <allowed_mismatches> "
              << "-t <num_threads>\n";
}

void PrintMascot(){
       cout<<"\e[48;2;255;231;148m                                                            \e[m\n";
   cout<<"\e[48;2;255;231;148m                                                            \e[m\n";
   cout<<"\e[48;2;255;231;148m           \e[38;2;92;84;54;48;2;151;127;84m▄\e[38;2;1;1;1;48;2;39;36;21m▄\e[38;2;0;0;0;48;2;67;58;36m▄\e[38;2;20;16;10;48;2;174;158;100m▄\e[38;2;120;107;69;48;2;245;221;142m▄\e[38;2;231;209;134;48;2;255;231;148m▄\e[48;2;255;231;148m   \e[38;2;254;229;147;48;2;255;231;148m▄\e[38;2;55;51;38;48;2;229;207;133m▄\e[38;2;166;157;117;48;2;248;224;143m▄\e[38;2;177;163;110;48;2;255;231;148m▄\e[38;2;173;159;106;48;2;255;231;148m▄\e[38;2;182;168;112;48;2;255;231;148m▄\e[38;2;152;139;96;48;2;251;227;146m▄\e[38;2;162;145;98;48;2;253;228;147m▄\e[38;2;179;161;103;48;2;255;231;148m▄\e[38;2;211;183;118;48;2;255;231;148m▄\e[38;2;237;215;137;48;2;255;231;148m▄\e[48;2;255;231;148m                             \e[m\n";
   cout<<"\e[48;2;255;231;148m      \e[38;2;248;224;144;48;2;255;231;148m▄\e[38;2;111;101;65;48;2;251;227;146m▄\e[38;2;19;17;11;48;2;216;194;125m▄\e[38;2;17;15;10;48;2;204;186;119m▄\e[38;2;17;15;10;48;2;206;187;120m▄\e[38;2;16;13;9;48;2;147;129;87m▄\e[38;2;1;1;1;48;2;10;8;6m▄\e[48;2;0;0;0m \e[38;2;69;70;72;48;2;0;0;0m▄\e[38;2;138;140;141;48;2;8;7;5m▄\e[38;2;23;23;22;48;2;84;74;49m▄\e[38;2;69;63;36;48;2;246;223;143m▄\e[38;2;237;215;137;48;2;255;231;148m▄\e[48;2;255;231;148m \e[38;2;222;201;129;48;2;245;221;142m▄\e[38;2;130;134;138;48;2;54;56;57m▄\e[38;2;139;147;154;48;2;151;159;166m▄\e[38;2;46;49;51;48;2;127;133;136m▄\e[38;2;59;63;66;48;2;112;118;120m▄\e[38;2;157;165;173;48;2;102;107;105m▄\e[38;2;167;176;185;48;2;88;93;98m▄\e[38;2;150;158;166;48;2;68;70;70m▄\e[38;2;90;93;95;48;2;65;65;60m▄\e[38;2;13;13;13;48;2;88;79;67m▄\e[38;2;78;79;80;48;2;89;82;57m▄\e[38;2;103;107;107;48;2;176;160;102m▄\e[38;2;120;111;72;48;2;251;227;146m▄\e[38;2;231;210;134;48;2;255;231;148m▄\e[48;2;255;231;148m                          \e[m\n";
   cout<<"\e[48;2;255;231;148m      \e[38;2;235;213;136;48;2;194;176;113m▄\e[38;2;97;88;56;48;2;34;31;20m▄\e[38;2;2;2;2;48;2;0;0;0m▄\e[48;2;0;0;0m  \e[38;2;86;87;88;48;2;12;12;12m▄\e[38;2;244;245;245;48;2;54;54;54m▄\e[38;2;246;246;246;48;2;83;83;83m▄\e[38;2;189;189;189;48;2;100;101;102m▄\e[38;2;226;226;226;48;2;238;239;240m▄\e[38;2;255;255;255;48;2;237;237;237m▄\e[38;2;249;249;249;48;2;165;168;171m▄\e[38;2;171;170;170;48;2;121;115;92m▄\e[38;2;141;134;110;48;2;228;207;134m▄\e[38;2;96;97;91;48;2;110;100;72m▄\e[38;2;79;81;82;48;2;145;152;160m▄\e[38;2;145;153;159;48;2;50;53;55m▄\e[38;2;170;180;189;48;2;125;132;139m▄\e[38;2;146;149;152;48;2;172;182;191m▄\e[38;2;112;107;101;48;2;164;173;181m▄\e[38;2;160;125;89;48;2;151;154;157m▄\e[38;2;8;7;5;48;2;38;40;42m▄\e[38;2;77;60;41;48;2;1;1;1m▄\e[38;2;152;114;76;48;2;13;9;6m▄\e[38;2;185;143;95;48;2;38;35;30m▄\e[38;2;196;149;100;48;2;121;120;118m▄\e[38;2;57;46;31;48;2;2;2;2m▄\e[38;2;6;5;3;48;2;59;53;34m▄\e[38;2;126;102;68;48;2;246;220;141m▄\e[38;2;188;180;156;48;2;237;215;138m▄\e[38;2;158;155;149;48;2;213;194;127m▄\e[38;2;79;80;81;48;2;169;152;100m▄\e[38;2;1;1;1;48;2;153;139;89m▄\e[38;2;61;55;35;48;2;206;187;119m▄\e[38;2;229;207;133;48;2;249;226;144m▄\e[48;2;255;231;148m                   \e[m\n";
   cout<<"\e[48;2;255;231;148m       \e[38;2;255;231;148;48;2;226;205;131m▄\e[38;2;253;228;146;48;2;152;126;83m▄\e[38;2;238;215;138;48;2;41;38;26m▄\e[38;2;193;174;112;48;2;4;4;3m▄\e[38;2;64;56;37;48;2;47;49;50m▄\e[38;2;17;16;16;48;2;210;210;211m▄\e[38;2;46;46;46;48;2;251;251;251m▄\e[38;2;161;161;161;48;2;254;254;254m▄\e[38;2;228;228;228;48;2;254;254;254m▄\e[38;2;253;253;253;48;2;255;255;255m▄\e[38;2;218;220;221;48;2;238;238;238m▄\e[38;2;126;126;126;48;2;124;125;126m▄\e[38;2;176;179;181;48;2;178;178;179m▄\e[38;2;136;144;150;48;2;72;74;76m▄\e[38;2;172;182;191;48;2;155;163;170m▄\e[38;2;137;139;138;48;2;169;179;187m▄\e[38;2;154;85;12;48;2;137;131;124m▄\e[38;2;242;122;0;48;2;171;104;35m▄\e[38;2;169;131;92;48;2;115;89;60m▄\e[38;2;222;171;114;48;2;205;156;103m▄\e[38;2;185;147;99;48;2;95;74;47m▄\e[38;2;255;192;128;48;2;231;174;116m▄\e[48;2;255;192;128m \e[38;2;249;188;127;48;2;255;192;128m▄\e[38;2;218;170;121;48;2;248;187;126m▄\e[38;2;253;191;127;48;2;213;163;109m▄\e[38;2;244;184;122;48;2;77;62;40m▄\e[38;2;52;42;26;48;2;11;10;10m▄\e[38;2;184;134;78;48;2;226;220;213m▄\e[38;2;254;251;247;48;2;248;248;248m▄\e[38;2;247;247;247;48;2;155;156;156m▄\e[38;2;245;245;245;48;2;84;84;84m▄\e[38;2;211;211;211;48;2;118;112;90m▄\e[38;2;50;50;49;48;2;140;126;81m▄\e[38;2;7;6;4;48;2;143;129;83m▄\e[38;2;10;9;6;48;2;169;153;99m▄\e[38;2;31;27;17;48;2;215;195;124m▄\e[38;2;173;155;98;48;2;253;229;146m▄\e[48;2;255;231;148m               \e[m\n";
   cout<<"\e[48;2;255;231;148m          \e[38;2;190;171;111;48;2;187;168;108m▄\e[38;2;40;36;23;48;2;33;29;18m▄\e[38;2;1;1;0;48;2;0;0;0m▄\e[38;2;87;88;90;48;2;116;117;117m▄\e[38;2;162;162;162;48;2;255;255;255m▄\e[38;2;139;138;135;48;2;249;249;249m▄\e[38;2;141;135;116;48;2;235;236;236m▄\e[38;2;127;110;78;48;2;172;172;171m▄\e[38;2;97;89;61;48;2;93;96;97m▄\e[38;2;105;111;116;48;2;64;66;69m▄\e[38;2;165;174;182;48;2;148;156;163m▄\e[38;2;63;46;28;48;2;145;151;158m▄\e[38;2;247;125;0;48;2;125;67;9m▄\e[48;2;255;128;0m \e[38;2;233;117;5;48;2;244;123;1m▄\e[38;2;246;183;120;48;2;236;181;124m▄\e[38;2;236;183;128;48;2;196;157;115m▄\e[38;2;214;164;112;48;2;229;178;125m▄\e[38;2;144;115;79;48;2;224;171;116m▄\e[38;2;231;176;118;48;2;251;190;126m▄\e[38;2;186;144;98;48;2;242;184;124m▄\e[38;2;204;159;110;48;2;207;163;116m▄\e[38;2;202;157;106;48;2;249;188;126m▄\e[38;2;230;177;123;48;2;254;191;128m▄\e[38;2;248;185;122;48;2;238;185;129m▄\e[38;2;85;48;8;48;2;90;55;9m▄\e[38;2;113;105;87;48;2;151;146;141m▄\e[38;2;176;173;163;48;2;253;253;253m▄\e[38;2;205;205;205;48;2;179;182;186m▄\e[38;2;191;191;192;48;2;47;47;47m▄\e[38;2;22;22;18;48;2;13;11;7m▄\e[38;2;88;80;51;48;2;54;48;31m▄\e[38;2;211;183;118;48;2;89;80;52m▄\e[38;2;255;231;148;48;2;120;108;69m▄\e[38;2;255;231;148;48;2;252;228;146m▄\e[48;2;255;231;148m               \e[m\n";
   cout<<"\e[48;2;255;231;148m          \e[38;2;255;231;148;48;2;243;220;141m▄\e[38;2;255;231;148;48;2;196;178;114m▄\e[38;2;255;231;148;48;2;176;159;102m▄\e[38;2;251;228;146;48;2;180;163;107m▄\e[38;2;223;202;129;48;2;172;157;106m▄\e[38;2;212;193;122;48;2;178;161;104m▄\e[38;2;221;200;127;48;2;218;196;125m▄\e[38;2;247;223;144;48;2;250;222;146m▄\e[38;2;129;112;87;48;2;154;142;100m▄\e[38;2;113;119;126;48;2;120;127;134m▄\e[38;2;118;84;40;48;2;116;117;118m▄\e[38;2;255;128;0;48;2;198;99;0m▄\e[48;2;255;128;0m  \e[38;2;135;76;13;48;2;182;97;10m▄\e[38;2;179;140;100;48;2;180;142;99m▄\e[38;2;252;190;127;48;2;228;175;120m▄\e[38;2;248;187;125;48;2;190;147;103m▄\e[38;2;249;188;125;48;2;171;136;90m▄\e[38;2;109;87;63;48;2;103;82;58m▄\e[38;2;255;192;128;48;2;216;165;109m▄\e[38;2;252;190;126;48;2;192;149;100m▄\e[38;2;187;143;99;48;2;116;92;62m▄\e[38;2;106;87;60;48;2;111;89;64m▄\e[38;2;181;146;108;48;2;250;189;125m▄\e[38;2;225;156;84;48;2;161;114;67m▄\e[38;2;141;85;24;48;2;127;105;59m▄\e[38;2;207;186;120;48;2;195;176;115m▄\e[38;2;253;229;147;48;2;198;183;133m▄\e[38;2;245;222;142;48;2;116;111;90m▄\e[38;2;243;220;141;48;2;40;35;26m▄\e[38;2;243;220;141;48;2;34;30;21m▄\e[38;2;251;227;145;48;2;175;153;99m▄\e[48;2;255;231;148m                 \e[m\n";
   cout<<"\e[48;2;255;231;148m            \e[38;2;221;200;127;48;2;255;231;148m▄\e[38;2;137;123;76;48;2;233;210;135m▄\e[38;2;136;111;60;48;2;95;80;53m▄\e[38;2;149;117;60;48;2;47;43;28m▄\e[38;2;104;75;32;48;2;40;30;12m▄\e[38;2;103;54;2;48;2;129;102;52m▄\e[38;2;155;80;0;48;2;129;107;63m▄\e[38;2;70;36;0;48;2;55;51;46m▄\e[38;2;208;105;0;48;2;217;116;14m▄\e[48;2;255;128;0m  \e[38;2;250;126;0;48;2;254;128;0m▄\e[38;2;228;158;87;48;2;169;114;57m▄\e[38;2;253;191;128;48;2;249;188;126m▄\e[48;2;255;192;128m \e[38;2;254;191;128;48;2;233;178;122m▄\e[38;2;240;184;127;48;2;184;148;101m▄\e[38;2;252;191;128;48;2;245;186;124m▄\e[38;2;255;192;128;48;2;240;188;134m▄\e[38;2;219;167;111;48;2;216;166;111m▄\e[38;2;110;89;60;48;2;214;164;111m▄\e[38;2;150;117;81;48;2;239;181;120m▄\e[38;2;172;132;91;48;2;253;191;127m▄\e[38;2;244;185;124;48;2;188;146;103m▄\e[38;2;255;192;128;48;2;95;75;50m▄\e[38;2;227;173;118;48;2;170;148;96m▄\e[38;2;143;109;75;48;2;204;184;118m▄\e[38;2;179;149;97;48;2;236;214;137m▄\e[38;2;179;162;104;48;2;253;229;147m▄\e[38;2;246;223;143;48;2;255;231;148m▄\e[48;2;255;231;148m                  \e[m\n";
   cout<<"\e[48;2;255;231;148m          \e[38;2;254;230;148;48;2;255;231;148m▄\e[38;2;104;95;61;48;2;212;182;117m▄\e[38;2;42;36;24;48;2;77;64;38m▄\e[38;2;96;78;50;48;2;149;77;0m▄\e[38;2;93;50;0;48;2;192;97;0m▄\e[38;2;219;112;0;48;2;209;105;0m▄\e[38;2;220;110;0;48;2;209;105;0m▄\e[38;2;216;108;0;48;2;194;97;0m▄\e[38;2;230;122;6;48;2;212;110;2m▄\e[38;2;151;84;8;48;2;124;64;2m▄\e[38;2;220;111;0;48;2;209;105;0m▄\e[48;2;255;128;0m  \e[38;2;243;125;4;48;2;249;126;2m▄\e[38;2;175;123;68;48;2;251;172;93m▄\e[38;2;167;134;99;48;2;255;192;128m▄\e[48;2;255;192;128m \e[38;2;252;190;126;48;2;252;190;127m▄\e[38;2;244;184;124;48;2;245;186;126m▄\e[38;2;253;191;128;48;2;254;191;128m▄\e[48;2;255;192;128m  \e[38;2;255;192;128;48;2;231;175;117m▄\e[38;2;253;191;128;48;2;186;146;105m▄\e[38;2;229;181;131;48;2;181;140;97m▄\e[38;2;251;189;126;48;2;184;142;98m▄\e[38;2;255;192;128;48;2;250;189;127m▄\e[48;2;255;192;128m  \e[38;2;255;192;128;48;2;241;181;121m▄\e[38;2;251;189;126;48;2;140;110;73m▄\e[38;2;86;66;43;48;2;168;150;95m▄\e[38;2;206;187;120;48;2;251;227;145m▄\e[48;2;255;231;148m                 \e[m\n";
   cout<<"\e[48;2;255;231;148m          \e[38;2;224;203;130;48;2;240;216;139m▄\e[38;2;94;79;45;48;2;157;133;88m▄\e[38;2;143;73;1;48;2;119;101;62m▄\e[38;2;240;121;1;48;2;107;54;0m▄\e[38;2;163;102;37;48;2;153;82;5m▄\e[38;2;219;167;111;48;2;69;50;30m▄\e[38;2;247;186;124;48;2;91;70;42m▄\e[38;2;252;190;128;48;2;77;59;40m▄\e[38;2;193;146;98;48;2;96;71;38m▄\e[38;2;102;62;15;48;2;126;80;22m▄\e[38;2;230;117;0;48;2;223;112;0m▄\e[48;2;255;128;0m \e[38;2;251;126;0;48;2;255;128;0m▄\e[38;2;154;85;17;48;2;171;90;4m▄\e[38;2;98;80;56;48;2;33;23;7m▄\e[38;2;165;129;91;48;2;138;110;82m▄\e[38;2;145;73;0;48;2;196;145;93m▄\e[38;2;225;114;0;48;2;209;148;85m▄\e[38;2;130;68;1;48;2;134;98;58m▄\e[38;2;205;157;107;48;2;237;179;119m▄\e[48;2;255;192;128m           \e[38;2;186;146;100;48;2;146;113;79m▄\e[38;2;22;19;15;48;2;96;86;56m▄\e[48;2;255;231;148m                 \e[m\n";
   cout<<"\e[48;2;255;231;148m         \e[38;2;232;210;135;48;2;249;225;144m▄\e[38;2;110;88;47;48;2;142;128;80m▄\e[38;2;51;25;0;48;2;146;79;6m▄\e[38;2;10;5;0;48;2;243;122;0m▄\e[38;2;121;68;4;48;2;226;115;2m▄\e[38;2;175;123;64;48;2;111;73;29m▄\e[38;2;249;187;125;48;2;89;70;47m▄\e[38;2;221;169;112;48;2;85;65;43m▄\e[38;2;2;2;1;48;2;84;65;44m▄\e[38;2;22;19;13;48;2;140;107;72m▄\e[38;2;105;59;8;48;2;129;77;19m▄\e[38;2;237;120;0;48;2;238;120;0m▄\e[48;2;255;128;0m \e[38;2;197;100;3;48;2;232;117;0m▄\e[38;2;123;93;63;48;2;184;119;52m▄\e[38;2;128;81;34;48;2;127;99;70m▄\e[38;2;188;98;0;48;2;175;98;17m▄\e[38;2;127;64;0;48;2;255;128;0m▄\e[48;2;255;128;0m \e[38;2;255;128;0;48;2;231;118;0m▄\e[38;2;248;125;0;48;2;83;57;28m▄\e[38;2;42;26;10;48;2;239;181;120m▄\e[38;2;163;123;82;48;2;255;192;128m▄\e[38;2;253;191;128;48;2;255;192;128m▄\e[48;2;255;192;128m       \e[38;2;251;189;127;48;2;255;192;128m▄\e[38;2;82;62;42;48;2;125;97;66m▄\e[38;2;63;39;1;48;2;29;28;18m▄\e[38;2;138;125;80;48;2;254;230;147m▄\e[48;2;255;231;148m                \e[m\n";
   cout<<"\e[48;2;255;231;148m         \e[38;2;243;220;141;48;2;223;202;130m▄\e[38;2;213;193;124;48;2;87;78;51m▄\e[38;2;121;98;57;48;2;95;85;54m▄\e[38;2;240;121;2;48;2;101;55;6m▄\e[38;2;196;100;0;48;2;247;124;0m▄\e[38;2;183;95;4;48;2;194;118;39m▄\e[38;2;137;93;47;48;2;225;172;117m▄\e[38;2;208;161;112;48;2;222;170;114m▄\e[38;2;132;102;66;48;2;22;17;10m▄\e[38;2;137;108;73;48;2;192;147;100m▄\e[38;2;131;90;45;48;2;146;99;43m▄\e[38;2;189;97;0;48;2;229;116;0m▄\e[38;2;236;119;0;48;2;254;127;0m▄\e[38;2;75;52;29;48;2;170;95;20m▄\e[38;2;188;132;71;48;2;109;88;63m▄\e[38;2;118;65;6;48;2;69;43;13m▄\e[38;2;166;84;0;48;2;147;76;0m▄\e[48;2;255;128;0m    \e[38;2;255;128;0;48;2;243;123;0m▄\e[38;2;254;128;0;48;2;124;66;0m▄\e[38;2;154;91;27;48;2;73;53;26m▄\e[38;2;120;91;61;48;2;164;126;84m▄\e[38;2;123;93;63;48;2;196;148;99m▄\e[38;2;89;75;56;48;2;202;153;103m▄\e[38;2;85;62;26;48;2;207;158;106m▄\e[38;2;118;73;16;48;2;207;158;106m▄\e[38;2;123;63;2;48;2;193;146;98m▄\e[38;2;132;66;0;48;2;192;144;95m▄\e[38;2;193;97;0;48;2;154;108;61m▄\e[38;2;113;57;0;48;2;107;60;8m▄\e[38;2;168;86;0;48;2;238;120;0m▄\e[38;2;201;107;7;48;2;95;61;19m▄\e[38;2;155;137;85;48;2;231;210;134m▄\e[38;2;253;229;147;48;2;255;231;148m▄\e[48;2;255;231;148m              \e[m\n";
   cout<<"\e[48;2;255;231;148m          \e[38;2;229;207;134;48;2;234;211;136m▄\e[38;2;169;113;37;48;2;151;102;38m▄\e[38;2;116;78;25;48;2;210;106;0m▄\e[38;2;65;34;0;48;2;70;36;0m▄\e[38;2;254;128;0;48;2;242;122;0m▄\e[38;2;251;126;0;48;2;200;102;2m▄\e[38;2;232;118;0;48;2;73;51;27m▄\e[38;2;33;28;11;48;2;210;165;117m▄\e[38;2;130;97;60;48;2;173;134;93m▄\e[38;2;60;31;0;48;2;45;31;16m▄\e[38;2;18;10;0;48;2;60;31;0m▄\e[38;2;96;48;0;48;2;97;49;0m▄\e[38;2;121;84;47;48;2;118;81;34m▄\e[38;2;136;101;59;48;2;154;110;63m▄\e[38;2;82;48;10;48;2;157;79;0m▄\e[38;2;247;124;1;48;2;252;126;0m▄\e[48;2;255;128;0m     \e[38;2;176;94;6;48;2;198;100;0m▄\e[38;2;163;126;82;48;2;170;119;68m▄\e[38;2;255;192;128;48;2;249;188;126m▄\e[38;2;255;192;128;48;2;254;192;128m▄\e[38;2;189;142;95;48;2;182;139;95m▄\e[38;2;142;75;6;48;2;137;75;10m▄\e[38;2;255;128;0;48;2;237;120;0m▄\e[38;2;255;128;0;48;2;253;127;0m▄\e[48;2;255;128;0m \e[38;2;237;119;0;48;2;252;127;0m▄\e[38;2;129;68;0;48;2;221;112;0m▄\e[38;2;40;21;0;48;2;47;24;0m▄\e[38;2;199;102;0;48;2;219;111;0m▄\e[38;2;163;115;42;48;2;111;90;48m▄\e[38;2;230;208;133;48;2;239;217;139m▄\e[48;2;255;231;148m              \e[m\n";
   cout<<"\e[48;2;255;231;148m          \e[38;2;252;228;146;48;2;236;213;136m▄\e[38;2;182;160;101;48;2;155;109;41m▄\e[38;2;158;142;91;48;2;145;126;78m▄\e[38;2;176;90;0;48;2;111;56;0m▄\e[38;2;174;111;34;48;2;168;89;6m▄\e[38;2;170;111;34;48;2;113;64;6m▄\e[38;2;243;123;0;48;2;240;121;0m▄\e[48;2;251;126;0m \e[38;2;182;92;0;48;2;180;92;3m▄\e[38;2;192;97;0;48;2;175;90;4m▄\e[38;2;252;127;0;48;2;116;58;0m▄\e[38;2;119;66;3;48;2;225;114;1m▄\e[38;2;207;160;112;48;2;95;75;53m▄\e[38;2;150;109;60;48;2;113;86;55m▄\e[38;2;206;105;0;48;2;192;95;0m▄\e[48;2;255;128;0m     \e[38;2;190;99;4;48;2;228;115;0m▄\e[38;2;130;98;65;48;2;96;71;45m▄\e[38;2;254;192;128;48;2;241;185;128m▄\e[48;2;255;192;128m \e[38;2;254;191;128;48;2;255;192;128m▄\e[38;2;171;132;91;48;2;189;142;94m▄\e[38;2;135;68;1;48;2;143;74;5m▄\e[38;2;252;127;0;48;2;255;128;0m▄\e[48;2;255;128;0m   \e[38;2;255;128;0;48;2;216;110;0m▄\e[38;2;101;54;0;48;2;24;12;0m▄\e[38;2;127;68;0;48;2;170;86;0m▄\e[38;2;132;69;0;48;2;194;112;22m▄\e[38;2;152;130;76;48;2;192;174;111m▄\e[38;2;249;225;144;48;2;254;230;148m▄\e[48;2;255;231;148m             \e[m\n";
   cout<<"\e[48;2;255;231;148m            \e[38;2;215;193;124;48;2;146;132;85m▄\e[38;2;117;68;5;48;2;191;100;2m▄\e[38;2;125;114;72;48;2;130;105;53m▄\e[38;2;118;97;53;48;2;154;109;45m▄\e[38;2;196;100;0;48;2;236;118;0m▄\e[38;2;26;13;0;48;2;140;71;0m▄\e[38;2;116;61;0;48;2;162;83;0m▄\e[38;2;92;50;0;48;2;172;89;0m▄\e[38;2;166;84;0;48;2;242;122;0m▄\e[38;2;99;79;59;48;2;36;27;15m▄\e[38;2;252;190;127;48;2;229;174;119m▄\e[38;2;141;104;65;48;2;124;89;52m▄\e[38;2;189;95;0;48;2;210;108;3m▄\e[48;2;255;128;0m    \e[38;2;88;44;0;48;2;192;99;0m▄\e[38;2;45;37;29;48;2;62;38;13m▄\e[38;2;255;192;128;48;2;232;175;116m▄\e[48;2;255;192;128m  \e[38;2;241;181;121;48;2;247;186;124m▄\e[38;2;73;56;37;48;2;107;89;62m▄\e[38;2;138;101;62;48;2;161;104;43m▄\e[38;2;200;101;2;48;2;235;119;0m▄\e[48;2;255;128;0m   \e[38;2;254;127;0;48;2;216;110;0m▄\e[38;2;27;14;0;48;2;154;81;0m▄\e[38;2;92;48;0;48;2;43;23;0m▄\e[38;2;201;104;0;48;2;104;54;0m▄\e[38;2;142;96;34;48;2;153;116;55m▄\e[38;2;231;209;134;48;2;240;217;139m▄\e[48;2;255;231;148m             \e[m\n";
   cout<<"\e[48;2;255;231;148m            \e[38;2;255;231;148;48;2;254;230;148m▄\e[38;2;255;231;148;48;2;185;167;106m▄\e[38;2;255;231;148;48;2;220;200;128m▄\e[38;2;243;220;141;48;2;200;181;116m▄\e[38;2;38;35;22;48;2;97;60;14m▄\e[38;2;7;1;0;48;2;1;1;0m▄\e[38;2;83;43;0;48;2;32;18;0m▄\e[38;2;151;78;0;48;2;58;32;0m▄\e[38;2;180;92;0;48;2;68;35;0m▄\e[38;2;80;42;0;48;2;37;24;2m▄\e[38;2;33;28;16;48;2;105;85;61m▄\e[38;2;125;66;1;48;2;168;108;42m▄\e[38;2;228;116;0;48;2;144;77;7m▄\e[38;2;220;112;0;48;2;248;125;1m▄\e[48;2;255;128;0m \e[38;2;254;128;0;48;2;255;128;0m▄\e[38;2;121;61;0;48;2;226;114;0m▄\e[38;2;55;38;14;48;2;101;53;0m▄\e[38;2;227;177;126;48;2;175;140;104m▄\e[38;2;154;124;92;48;2;232;175;118m▄\e[38;2;153;121;85;48;2;255;192;128m▄\e[38;2;127;101;68;48;2;231;174;116m▄\e[38;2;133;101;68;48;2;44;39;23m▄\e[38;2;250;189;126;48;2;154;118;81m▄\e[38;2;171;133;86;48;2;137;104;69m▄\e[38;2;62;40;12;48;2;70;38;3m▄\e[38;2;244;123;0;48;2;251;126;0m▄\e[48;2;255;128;0m  \e[38;2;187;94;0;48;2;255;128;0m▄\e[38;2;86;44;0;48;2;62;31;0m▄\e[38;2;68;38;5;48;2;123;63;0m▄\e[38;2;205;105;0;48;2;255;128;0m▄\e[38;2;98;50;0;48;2;115;66;12m▄\e[38;2;183;164;104;48;2;217;197;126m▄\e[48;2;255;231;148m             \e[m\n";
   cout<<"\e[48;2;255;231;148m           \e[38;2;241;219;140;48;2;255;231;148m▄\e[38;2;98;84;51;48;2;255;231;148m▄\e[38;2;16;10;6;48;2;255;231;148m▄\e[38;2;31;6;3;48;2;253;229;147m▄\e[38;2;166;2;2;48;2;144;121;78m▄\e[38;2;249;0;0;48;2;170;10;6m▄\e[38;2;248;0;0;48;2;186;0;0m▄\e[38;2;175;20;1;48;2;121;36;3m▄\e[38;2;196;99;0;48;2;222;113;0m▄\e[38;2;255;128;0;48;2;249;125;0m▄\e[38;2;149;77;0;48;2;154;81;0m▄\e[38;2;119;0;5;48;2;25;13;0m▄\e[38;2;60;6;0;48;2;120;64;0m▄\e[38;2;81;43;0;48;2;213;109;0m▄\e[38;2;245;123;0;48;2;233;118;0m▄\e[38;2;228;115;0;48;2;243;122;0m▄\e[38;2;99;54;6;48;2;215;108;0m▄\e[38;2;192;140;88;48;2;139;79;17m▄\e[38;2;254;192;128;48;2;189;145;96m▄\e[38;2;255;192;128;48;2;253;191;128m▄\e[38;2;255;192;128;48;2;229;175;120m▄\e[38;2;255;192;128;48;2;209;160;109m▄\e[38;2;255;192;128;48;2;217;166;111m▄\e[38;2;255;192;128;48;2;253;191;128m▄\e[48;2;255;192;128m \e[38;2;255;192;128;48;2;238;181;121m▄\e[38;2;171;129;84;48;2;99;71;32m▄\e[38;2;133;108;62;48;2;150;81;0m▄\e[38;2;182;94;0;48;2;237;119;0m▄\e[38;2;254;128;0;48;2;255;128;0m▄\e[38;2;104;53;0;48;2;110;56;0m▄\e[38;2;98;55;0;48;2;93;51;0m▄\e[38;2;45;16;0;48;2;78;40;2m▄\e[38;2;179;93;4;48;2;172;88;3m▄\e[38;2;87;47;4;48;2;93;52;3m▄\e[38;2;116;96;59;48;2;171;150;91m▄\e[38;2;198;178;115;48;2;255;231;148m▄\e[38;2;249;225;144;48;2;255;231;148m▄\e[38;2;252;229;147;48;2;255;231;148m▄\e[48;2;255;231;148m          \e[m\n";
   cout<<"\e[48;2;255;231;148m         \e[38;2;232;210;135;48;2;254;230;148m▄\e[38;2;112;76;49;48;2;209;189;121m▄\e[38;2;171;1;2;48;2;89;49;32m▄\e[38;2;252;0;0;48;2;121;0;2m▄\e[38;2;255;0;0;48;2;194;0;2m▄\e[38;2;255;0;0;48;2;189;0;4m▄\e[38;2;255;0;0;48;2;187;0;4m▄\e[38;2;255;0;0;48;2;166;0;5m▄\e[38;2;244;0;0;48;2;175;0;4m▄\e[38;2;107;13;5;48;2;145;6;1m▄\e[38;2;99;52;1;48;2;100;53;1m▄\e[38;2;188;96;0;48;2;232;117;0m▄\e[38;2;86;45;0;48;2;56;27;0m▄\e[38;2;124;0;4;48;2;199;0;6m▄\e[38;2;225;1;1;48;2;242;18;0m▄\e[38;2;131;8;0;48;2;78;18;0m▄\e[38;2;39;18;0;48;2;122;64;0m▄\e[38;2;200;93;0;48;2;226;115;0m▄\e[38;2;164;36;4;48;2;80;46;4m▄\e[38;2;119;67;46;48;2;170;128;86m▄\e[38;2;143;110;75;48;2;255;192;128m▄\e[38;2;221;171;119;48;2;255;192;128m▄\e[48;2;255;192;128m    \e[38;2;247;187;126;48;2;255;192;128m▄\e[38;2;188;145;100;48;2;255;192;128m▄\e[38;2;132;100;68;48;2;231;179;126m▄\e[38;2;55;43;32;48;2;178;142;105m▄\e[38;2;102;51;6;48;2;107;59;3m▄\e[38;2;132;66;1;48;2;199;101;0m▄\e[38;2;110;42;0;48;2;90;46;0m▄\e[38;2;129;6;0;48;2;120;64;0m▄\e[38;2;170;2;3;48;2;161;5;8m▄\e[38;2;76;34;1;48;2;115;57;0m▄\e[38;2;107;55;0;48;2;95;49;0m▄\e[38;2;159;11;0;48;2;172;3;1m▄\e[38;2;254;0;0;48;2;122;0;1m▄\e[38;2;220;0;0;48;2;65;14;10m▄\e[38;2;145;2;2;48;2;125;106;69m▄\e[38;2;102;49;33;48;2;214;194;125m▄\e[38;2;168;150;97;48;2;255;231;148m▄\e[38;2;198;179;115;48;2;255;231;148m▄\e[38;2;216;191;123;48;2;255;231;148m▄\e[38;2;250;224;143;48;2;255;231;148m▄\e[48;2;255;231;148m     \e[m\n";
   cout<<"\e[48;2;255;231;148m        \e[38;2;227;205;132;48;2;255;231;148m▄\e[38;2;58;4;3;48;2;150;129;83m▄\e[38;2;246;0;1;48;2;126;7;5m▄\e[38;2;253;0;0;48;2;245;0;0m▄\e[38;2;249;0;0;48;2;255;0;0m▄\e[38;2;246;0;0;48;2;255;0;0m▄\e[38;2;245;0;0;48;2;255;0;0m▄▄▄\e[38;2;243;0;0;48;2;244;0;1m▄\e[38;2;161;5;0;48;2;112;13;5m▄\e[38;2;153;77;0;48;2;169;85;0m▄\e[38;2;132;67;0;48;2;227;115;0m▄\e[38;2;39;18;0;48;2;140;72;0m▄\e[38;2;226;0;7;48;2;157;0;5m▄\e[38;2;252;0;0;48;2;167;0;2m▄\e[38;2;236;0;1;48;2;145;2;4m▄\e[38;2;192;3;9;48;2;123;16;4m▄\e[38;2;50;0;2;48;2;85;19;0m▄\e[38;2;235;0;6;48;2;246;3;0m▄\e[38;2;247;0;0;48;2;210;2;2m▄\e[38;2;252;0;0;48;2;182;8;6m▄\e[38;2;250;0;0;48;2;115;32;27m▄\e[38;2;164;26;0;48;2;56;42;28m▄\e[38;2;159;142;0;48;2;65;50;33m▄\e[38;2;189;169;0;48;2;69;52;35m▄\e[38;2;15;14;0;48;2;64;48;33m▄\e[38;2;63;56;0;48;2;52;40;28m▄\e[38;2;201;180;0;48;2;32;27;16m▄\e[38;2;230;206;0;48;2;35;7;3m▄\e[38;2;119;11;0;48;2;192;0;2m▄\e[38;2;176;12;1;48;2;164;70;0m▄\e[38;2;57;28;0;48;2;85;46;0m▄\e[38;2;187;9;1;48;2;168;25;0m▄\e[38;2;252;0;0;48;2;244;0;0m▄\e[38;2;218;0;0;48;2;220;1;0m▄\e[38;2;103;39;1;48;2;100;47;0m▄\e[38;2;172;89;0;48;2;206;106;0m▄\e[38;2;35;12;0;48;2;120;16;3m▄\e[38;2;94;0;0;48;2;250;0;0m▄\e[38;2;219;0;2;48;2;162;0;2m▄\e[38;2;245;0;1;48;2;87;0;1m▄\e[38;2;255;0;0;48;2;103;0;4m▄\e[38;2;255;0;0;48;2;142;17;15m▄\e[38;2;255;0;0;48;2;170;28;22m▄\e[38;2;254;0;0;48;2;79;33;23m▄\e[38;2;175;0;0;48;2;67;58;37m▄\e[38;2;21;5;3;48;2;207;186;120m▄\e[38;2;188;161;112;48;2;255;231;148m▄\e[38;2;254;230;147;48;2;255;231;148m▄\e[48;2;255;231;148m  \e[m\n";
   cout<<"\e[48;2;255;231;148m       \e[38;2;235;212;136;48;2;255;231;148m▄\e[38;2;52;8;7;48;2;95;83;53m▄\e[38;2;127;0;0;48;2;162;0;0m▄\e[38;2;49;29;7;48;2;204;0;0m▄\e[38;2;41;90;20;48;2;150;7;3m▄\e[38;2;49;119;28;48;2;89;22;7m▄\e[38;2;56;135;33;48;2;44;33;9m▄\e[38;2;59;145;35;48;2;39;39;10m▄\e[38;2;61;149;36;48;2;35;42;11m▄\e[38;2;62;152;37;48;2;34;43;11m▄\e[38;2;15;37;9;48;2;21;11;3m▄\e[38;2;63;154;38;48;2;31;50;11m▄\e[38;2;53;108;24;48;2;88;77;9m▄\e[38;2;108;65;3;48;2;32;21;1m▄\e[38;2;92;51;2;48;2;7;2;0m▄\e[38;2;47;117;27;48;2;148;38;14m▄\e[38;2;53;129;31;48;2;145;28;12m▄\e[38;2;47;115;27;48;2;141;19;9m▄\e[38;2;41;104;21;48;2;130;12;5m▄\e[38;2;40;95;21;48;2;61;6;0m▄\e[38;2;35;86;21;48;2;100;0;1m▄\e[38;2;34;83;20;48;2;45;0;1m▄\e[38;2;32;78;18;48;2;123;0;4m▄\e[38;2;37;77;18;48;2;171;0;0m▄\e[38;2;42;80;18;48;2;117;59;0m▄\e[38;2;92;107;7;48;2;225;203;0m▄\e[38;2;165;154;0;48;2;157;142;0m▄\e[38;2;190;170;0;48;2;217;195;0m▄\e[38;2;184;167;0;48;2;181;164;0m▄\e[38;2;214;194;0;48;2;188;169;0m▄\e[38;2;230;206;0;48;2;249;224;0m▄\e[38;2;94;85;0;48;2;131;122;0m▄\e[38;2;45;2;3;48;2;136;9;6m▄\e[38;2;56;14;4;48;2;119;29;2m▄\e[38;2;107;20;1;48;2;215;4;2m▄\e[38;2;126;33;7;48;2;215;0;4m▄\e[38;2;139;43;10;48;2;161;0;3m▄\e[38;2;100;40;8;48;2;125;1;3m▄\e[38;2;109;60;3;48;2;75;33;0m▄\e[38;2;191;98;0;48;2;146;25;1m▄\e[38;2;98;43;2;48;2;196;1;1m▄\e[38;2;127;60;13;48;2;248;0;0m▄\e[38;2;95;73;17;48;2;222;0;1m▄\e[38;2;71;73;17;48;2;206;0;0m▄\e[38;2;53;74;17;48;2;194;0;0m▄\e[38;2;50;71;17;48;2;193;0;0m▄\e[38;2;41;61;14;48;2;190;0;0m▄\e[38;2;46;72;17;48;2;190;0;0m▄\e[38;2;47;61;15;48;2;173;0;0m▄\e[38;2;48;5;5;48;2;86;33;22m▄\e[38;2;73;66;42;48;2;209;189;120m▄\e[38;2;205;186;119;48;2;255;231;148m▄\e[38;2;254;230;147;48;2;255;231;148m▄\e[m\n";
   cout<<"\e[48;2;255;231;148m       \e[38;2;154;139;89;48;2;192;166;113m▄\e[38;2;75;2;2;48;2;110;2;3m▄\e[38;2;21;53;12;48;2;17;16;4m▄\e[38;2;63;155;38;48;2;44;107;25m▄\e[38;2;62;151;37;48;2;63;154;37m▄\e[38;2;63;154;38;48;2;64;156;38m▄\e[48;2;64;156;38m    \e[38;2;15;37;9;48;2;16;38;9m▄\e[48;2;64;156;38m \e[38;2;50;103;23;48;2;39;79;18m▄\e[38;2;191;105;2;48;2;75;40;0m▄\e[38;2;240;121;0;48;2;96;50;1m▄\e[38;2;25;52;12;48;2;25;64;15m▄\e[38;2;65;150;36;48;2;64;156;38m▄\e[38;2;62;149;36;48;2;64;156;38m▄\e[38;2;62;151;36;48;2;64;156;38m▄\e[48;2;64;156;38m  \e[38;2;64;156;38;48;2;64;155;38m▄\e[38;2;64;156;38;48;2;63;153;37m▄▄\e[38;2;64;156;38;48;2;62;152;37m▄\e[38;2;64;156;38;48;2;58;140;33m▄\e[38;2;64;156;38;48;2;66;129;25m▄\e[38;2;64;156;38;48;2;62;85;11m▄\e[38;2;64;156;38;48;2;60;53;0m▄\e[38;2;64;156;38;48;2;66;68;0m▄\e[38;2;64;156;38;48;2;75;112;16m▄\e[38;2;64;156;38;48;2;52;119;28m▄\e[38;2;64;156;38;48;2;57;126;31m▄\e[38;2;64;156;38;48;2;57;130;32m▄\e[38;2;64;156;38;48;2;59;133;31m▄\e[38;2;64;156;38;48;2;59;137;33m▄\e[38;2;64;156;38;48;2;60;141;34m▄\e[38;2;64;154;37;48;2;55;134;32m▄\e[38;2;52;80;16;48;2;104;90;11m▄\e[38;2;71;39;1;48;2;228;116;0m▄\e[38;2;117;64;0;48;2;192;101;1m▄\e[38;2;165;86;0;48;2;22;54;12m▄\e[38;2;44;26;0;48;2;49;119;28m▄\e[38;2;39;93;21;48;2;61;149;36m▄\e[38;2;64;156;38;48;2;62;152;37m▄\e[38;2;62;151;37;48;2;58;140;34m▄\e[38;2;54;130;31;48;2;48;117;28m▄\e[38;2;60;146;36;48;2;60;146;35m▄\e[38;2;64;156;38;48;2;60;147;36m▄\e[38;2;60;146;35;48;2;25;60;14m▄\e[38;2;32;79;19;48;2;2;2;2m▄\e[38;2;19;28;11;48;2;109;98;65m▄\e[38;2;135;120;69;48;2;244;221;142m▄\e[m\n";
   cout<<"\e[48;2;255;231;148m     \e[38;2;254;230;147;48;2;255;231;148m▄\e[38;2;210;189;121;48;2;255;231;148m▄\e[38;2;79;76;45;48;2;149;133;86m▄\e[38;2;24;58;14;48;2;2;2;1m▄\e[38;2;49;119;27;48;2;40;100;21m▄\e[38;2;63;152;37;48;2;57;139;34m▄\e[38;2;63;154;37;48;2;39;93;22m▄\e[38;2;52;126;30;48;2;47;114;27m▄\e[38;2;39;94;22;48;2;51;123;29m▄\e[38;2;49;119;27;48;2;64;156;38m▄\e[48;2;64;156;38m \e[38;2;57;139;33;48;2;62;152;37m▄\e[38;2;49;117;28;48;2;23;55;13m▄\e[38;2;57;93;19;48;2;59;143;35m▄\e[38;2;160;86;1;48;2;66;95;18m▄\e[38;2;248;125;0;48;2;100;60;2m▄\e[38;2;255;128;0;48;2;148;78;0m▄\e[38;2;255;128;0;48;2;139;74;0m▄\e[38;2;255;128;0;48;2;137;75;0m▄\e[38;2;249;125;0;48;2;97;62;3m▄\e[38;2;185;94;0;48;2;52;70;7m▄\e[38;2;11;29;6;48;2;50;123;30m▄\e[38;2;51;128;28;48;2;64;156;38m▄\e[48;2;64;156;38m  \e[38;2;49;118;29;48;2;64;156;38m▄\e[38;2;37;86;20;48;2;64;156;38m▄\e[38;2;42;97;23;48;2;63;155;38m▄\e[38;2;49;118;28;48;2;60;146;36m▄\e[38;2;45;112;24;48;2;53;130;29m▄\e[38;2;42;102;24;48;2;50;121;29m▄\e[38;2;38;90;19;48;2;55;132;31m▄\e[38;2;32;81;18;48;2;58;142;34m▄\e[38;2;45;110;25;48;2;60;147;35m▄\e[38;2;47;116;26;48;2;60;147;35m▄\e[38;2;49;119;27;48;2;59;144;34m▄\e[38;2;54;129;31;48;2;61;147;35m▄\e[38;2;56;137;33;48;2;61;149;36m▄\e[38;2;47;90;18;48;2;60;145;36m▄\e[38;2;143;81;2;48;2;71;99;21m▄\e[38;2;250;126;0;48;2;183;97;3m▄\e[38;2;255;128;0;48;2;251;126;0m▄\e[48;2;255;128;0m \e[38;2;255;128;0;48;2;253;127;0m▄\e[38;2;243;123;0;48;2;158;81;0m▄\e[38;2;91;47;0;48;2;23;20;2m▄\e[38;2;39;94;22;48;2;61;149;36m▄\e[48;2;64;156;38m \e[38;2;57;141;33;48;2;51;126;29m▄\e[38;2;56;136;32;48;2;57;139;33m▄\e[38;2;39;95;22;48;2;60;147;35m▄\e[38;2;11;30;7;48;2;34;86;18m▄\e[38;2;60;147;36;48;2;60;146;36m▄\e[38;2;48;114;27;48;2;58;138;33m▄\e[38;2;58;139;34;48;2;29;69;16m▄\e[m\n";
   cout<<"\e[48;2;255;231;148m  \e[38;2;252;229;146;48;2;255;231;148m▄\e[38;2;228;207;132;48;2;255;231;148m▄\e[38;2;66;61;33;48;2;250;227;145m▄\e[38;2;37;88;19;48;2;129;120;73m▄\e[38;2;61;149;36;48;2;55;81;31m▄\e[38;2;55;134;32;48;2;54;124;31m▄\e[38;2;33;82;14;48;2;60;146;35m▄\e[38;2;33;83;13;48;2;64;155;38m▄\e[38;2;33;78;13;48;2;64;156;38m▄\e[38;2;33;76;14;48;2;64;156;38m▄\e[38;2;43;104;24;48;2;64;156;38m▄\e[38;2;53;130;31;48;2;61;150;36m▄\e[38;2;60;147;35;48;2;57;139;33m▄\e[38;2;60;146;35;48;2;64;155;38m▄\e[38;2;16;40;7;48;2;46;113;26m▄\e[38;2;50;120;28;48;2;59;144;35m▄\e[38;2;98;59;1;48;2;77;57;2m▄\e[38;2;253;127;0;48;2;240;121;0m▄\e[48;2;255;128;0m    \e[38;2;236;121;0;48;2;255;128;0m▄\e[48;2;254;127;0m \e[38;2;203;103;0;48;2;133;70;1m▄\e[38;2;28;67;16;48;2;31;76;18m▄\e[38;2;64;156;38;48;2;62;152;37m▄\e[38;2;64;156;38;48;2;54;126;30m▄\e[38;2;64;156;38;48;2;52;127;30m▄\e[38;2;64;156;38;48;2;53;129;31m▄\e[38;2;63;153;37;48;2;57;138;33m▄\e[38;2;60;145;36;48;2;52;125;30m▄\e[38;2;64;156;38;48;2;32;81;13m▄\e[38;2;64;156;38;48;2;47;117;27m▄\e[38;2;64;156;38;48;2;52;127;29m▄\e[38;2;64;156;38;48;2;55;136;33m▄\e[38;2;64;156;38;48;2;60;145;35m▄\e[38;2;64;156;38;48;2;60;147;35m▄\e[38;2;64;156;38;48;2;62;149;36m▄\e[38;2;64;156;38;48;2;62;151;37m▄\e[38;2;35;88;16;48;2;61;149;36m▄\e[38;2;151;77;0;48;2;45;37;2m▄\e[38;2;195;98;0;48;2;223;114;0m▄\e[38;2;195;98;0;48;2;241;122;0m▄\e[48;2;255;128;0m    \e[38;2;255;128;0;48;2;220;111;0m▄\e[38;2;49;27;0;48;2;6;16;3m▄\e[38;2;51;126;29;48;2;62;149;36m▄\e[48;2;64;156;38m \e[38;2;64;156;38;48;2;55;134;32m▄\e[38;2;20;47;10;48;2;11;30;6m▄\e[38;2;23;57;13;48;2;25;62;14m▄\e[38;2;28;69;15;48;2;46;114;26m▄\e[38;2;55;135;32;48;2;42;104;24m▄\e[38;2;64;156;38;48;2;63;153;37m▄\e[m\n";
   cout<<"\e[38;2;251;227;145;48;2;255;231;148m▄\e[38;2;172;156;98;48;2;247;223;143m▄\e[38;2;63;92;33;48;2;177;162;100m▄\e[38;2;57;137;33;48;2;54;73;31m▄\e[38;2;64;156;38;48;2;51;123;30m▄\e[48;2;64;156;38m \e[38;2;64;156;38;48;2;61;148;36m▄\e[38;2;64;156;38;48;2;55;132;31m▄\e[38;2;64;156;38;48;2;58;142;34m▄\e[38;2;64;156;38;48;2;57;140;33m▄▄\e[38;2;57;139;34;48;2;58;141;34m▄\e[38;2;37;84;18;48;2;54;132;32m▄\e[38;2;37;89;19;48;2;33;79;17m▄\e[38;2;26;62;13;48;2;23;57;11m▄\e[38;2;30;46;16;48;2;25;57;10m▄\e[38;2;108;96;57;48;2;21;51;10m▄\e[38;2;150;123;80;48;2;13;31;4m▄\e[38;2;180;138;90;48;2;71;51;12m▄\e[38;2;185;140;93;48;2;182;100;14m▄\e[38;2;167;125;82;48;2;185;96;5m▄\e[38;2;151;110;61;48;2;198;100;0m▄\e[38;2;115;65;13;48;2;245;124;0m▄\e[38;2;197;99;0;48;2;255;128;0m▄\e[38;2;250;126;0;48;2;203;105;2m▄\e[38;2;30;15;0;48;2;184;93;0m▄\e[38;2;178;101;3;48;2;204;103;0m▄\e[38;2;43;106;24;48;2;28;67;16m▄\e[48;2;64;156;38m              \e[38;2;52;125;29;48;2;23;60;13m▄\e[38;2;142;77;1;48;2;202;102;0m▄\e[38;2;161;81;0;48;2;174;89;0m▄\e[38;2;152;77;0;48;2;174;89;0m▄\e[38;2;253;127;0;48;2;255;128;0m▄\e[38;2;243;122;0;48;2;253;127;0m▄\e[38;2;156;80;0;48;2;237;119;0m▄\e[48;2;255;128;0m  \e[38;2;131;66;0;48;2;104;53;0m▄\e[38;2;46;111;27;48;2;43;106;25m▄\e[48;2;64;156;38m  \e[38;2;64;156;38;48;2;55;133;32m▄\e[38;2;60;146;35;48;2;25;61;13m▄\e[38;2;44;105;23;48;2;35;85;19m▄\e[38;2;35;88;19;48;2;62;151;36m▄\e[38;2;58;143;34;48;2;64;156;38m▄\e[m\n";
   cout<<"\e[38;2;29;45;13;48;2;189;171;109m▄\e[38;2;60;147;36;48;2;49;80;25m▄\e[38;2;64;156;38;48;2;60;144;35m▄\e[48;2;64;156;38m      \e[38;2;60;146;35;48;2;64;156;38m▄\e[38;2;32;81;18;48;2;56;136;33m▄\e[38;2;45;111;27;48;2;34;76;18m▄\e[38;2;33;51;17;48;2;52;123;29m▄\e[38;2;216;164;112;48;2;45;68;28m▄\e[38;2;255;192;128;48;2;138;115;82m▄\e[38;2;255;192;128;48;2;190;148;103m▄\e[38;2;255;192;128;48;2;240;182;122m▄\e[38;2;255;192;128;48;2;253;191;128m▄\e[48;2;255;192;128m   \e[38;2;255;192;128;48;2;246;186;124m▄\e[38;2;255;192;128;48;2;195;150;103m▄\e[38;2;243;185;126;48;2;101;72;39m▄\e[38;2;94;75;49;48;2;161;82;0m▄\e[38;2;58;35;2;48;2;72;39;0m▄\e[38;2;36;91;19;48;2;70;69;11m▄\e[38;2;64;156;38;48;2;62;152;37m▄\e[48;2;64;156;38m               \e[38;2;63;153;37;48;2;67;103;20m▄\e[38;2;47;115;28;48;2;84;42;1m▄\e[38;2;32;45;8;48;2;124;62;0m▄\e[38;2;202;101;0;48;2;247;124;0m▄\e[38;2;239;120;0;48;2;237;119;0m▄\e[38;2;98;50;0;48;2;102;53;0m▄\e[38;2;178;92;0;48;2;235;119;0m▄\e[38;2;224;114;0;48;2;250;126;0m▄\e[38;2;33;61;12;48;2;82;44;0m▄\e[38;2;64;156;38;48;2;55;132;32m▄\e[48;2;64;156;38m    \e[38;2;64;156;38;48;2;61;150;36m▄\e[38;2;64;156;38;48;2;41;100;24m▄\e[38;2;46;110;27;48;2;17;41;9m▄\e[m\n";
   cout<<"\e[38;2;64;156;38;48;2;51;124;29m▄\e[48;2;64;156;38m       \e[38;2;48;117;28;48;2;64;155;38m▄\e[38;2;30;64;13;48;2;38;92;21m▄\e[38;2;44;111;25;48;2;38;93;21m▄\e[38;2;91;89;52;48;2;41;87;23m▄\e[38;2;244;184;123;48;2;158;122;81m▄\e[48;2;255;192;128m          \e[38;2;250;189;126;48;2;255;192;128m▄\e[38;2;182;139;95;48;2;241;184;126m▄\e[38;2;228;173;117;48;2;83;66;47m▄\e[38;2;67;69;41;48;2;47;114;27m▄\e[38;2;47;115;27;48;2;64;156;38m▄\e[48;2;64;156;38m                \e[38;2;64;156;38;48;2;64;155;38m▄\e[38;2;64;156;38;48;2;46;113;25m▄\e[38;2;56;139;33;48;2;71;73;10m▄\e[38;2;63;104;21;48;2;180;94;0m▄\e[38;2;20;27;4;48;2;63;33;0m▄\e[38;2;36;60;12;48;2;89;53;2m▄\e[38;2;61;139;33;48;2;84;79;9m▄\e[38;2;64;156;38;48;2;59;143;34m▄\e[48;2;64;156;38m        \e[m\n";
   cout<<"\e[48;2;64;156;38m       \e[38;2;51;125;29;48;2;63;155;37m▄\e[38;2;18;44;10;48;2;10;26;4m▄\e[38;2;55;134;33;48;2;53;127;30m▄\e[38;2;66;72;37;48;2;39;79;23m▄\e[38;2;212;161;108;48;2;175;134;90m▄\e[38;2;255;192;128;48;2;254;191;128m▄\e[48;2;255;192;128m         \e[38;2;233;180;125;48;2;255;192;128m▄\e[38;2;255;192;128;48;2;252;190;127m▄\e[38;2;166;129;91;48;2;86;67;47m▄\e[38;2;89;79;52;48;2;221;172;121m▄\e[38;2;225;171;114;48;2;188;144;97m▄\e[38;2;66;69;35;48;2;28;69;17m▄\e[38;2;51;122;28;48;2;64;156;38m▄\e[48;2;64;156;38m                  \e[38;2;63;153;37;48;2;61;149;36m▄\e[38;2;54;133;30;48;2;29;68;12m▄\e[38;2;45;108;25;48;2;32;76;16m▄\e[38;2;49;119;28;48;2;58;141;34m▄\e[48;2;64;156;38m         \e[m\n";
   cout<<"\e[48;2;64;156;38m      \e[38;2;63;153;37;48;2;64;156;38m▄\e[38;2;19;50;10;48;2;37;91;21m▄\e[38;2;61;149;36;48;2;46;108;25m▄\e[38;2;36;79;21;48;2;45;107;26m▄\e[38;2;128;98;67;48;2;111;86;58m▄\e[38;2;255;192;128;48;2;250;189;127m▄\e[48;2;255;192;128m \e[38;2;239;181;122;48;2;255;192;128m▄\e[38;2;234;179;122;48;2;255;192;128m▄\e[48;2;255;192;128m    \e[38;2;252;190;128;48;2;254;191;128m▄\e[38;2;178;136;94;48;2;255;192;128m▄\e[38;2;255;192;128;48;2;254;192;128m▄\e[38;2;193;147;99;48;2;206;158;108m▄\e[38;2;148;114;77;48;2;224;170;114m▄\e[38;2;254;192;128;48;2;247;186;124m▄\e[38;2;20;17;13;48;2;23;18;13m▄\e[38;2;236;180;122;48;2;237;180;121m▄\e[38;2;94;86;59;48;2;92;74;53m▄\e[38;2;47;113;27;48;2;42;103;25m▄\e[48;2;64;156;38m                    \e[38;2;64;156;38;48;2;62;151;37m▄\e[38;2;30;74;16;48;2;40;94;21m▄\e[38;2;63;152;37;48;2;64;156;38m▄\e[48;2;64;156;38m        \e[m\n";
   cout<<"\e[49;38;2;64;156;38m▀▀▀▀▀▀\e[49;38;2;62;149;36m▀\e[49;38;2;19;46;9m▀\e[49;38;2;63;154;37m▀\e[49;38;2;34;79;21m▀\e[49;38;2;122;94;64m▀\e[49;38;2;254;192;128m▀\e[49;38;2;236;181;126m▀\e[49;38;2;142;113;82m▀\e[49;38;2;110;85;58m▀\e[49;38;2;243;183;124m▀\e[49;38;2;254;192;128m▀▀▀▀\e[49;38;2;140;107;72m▀\e[49;38;2;238;183;126m▀\e[49;38;2;189;143;97m▀\e[49;38;2;131;100;66m▀\e[49;38;2;254;192;128m▀\e[49;38;2;51;47;31m▀\e[49;38;2;215;162;107m▀\e[49;38;2;49;78;21m▀\e[49;38;2;62;151;36m▀\e[49;38;2;64;156;38m▀▀▀\e[49;38;2;57;139;34m▀\e[49;38;2;53;127;30m▀\e[49;38;2;52;128;29m▀\e[49;38;2;46;115;24m▀\e[49;38;2;46;111;24m▀▀\e[49;38;2;46;115;24m▀\e[49;38;2;46;111;24m▀▀\e[49;38;2;46;113;24m▀\e[49;38;2;46;111;24m▀\e[49;38;2;47;118;25m▀\e[49;38;2;53;129;30m▀\e[49;38;2;53;127;30m▀\e[49;38;2;57;143;34m▀\e[49;38;2;60;146;35m▀\e[49;38;2;63;154;37m▀\e[49;38;2;64;156;38m▀\e[49;38;2;44;108;26m▀\e[49;38;2;57;139;34m▀\e[49;38;2;64;156;38m▀▀▀▀▀▀▀▀\e[m\n";
}

int main(int argc, char* argv[]) {
    // Parse arguments more efficiently
    unordered_map<string, string> args;
    
    for (int i = 1; i < argc; ++i) {
        string arg = argv[i];
        
        if (arg == "-ms" || arg == "-mp" || arg == "-gp" || arg == "-ep" || 
            arg == "-irf" || arg == "-ird" || arg == "-m" || arg == "-t") {
            if (i + 1 < argc) {
                args[arg] = argv[++i];
            } else {
                std::cerr << "Error: Argument " << arg << " requires a value.\n";
                displayUsage();
                return 1;
            }
        } else if (arg == "-vb" || arg == "-mascot") {
            args[arg] = "true";
        } else {
            std::cerr << "Error: Unknown argument " << arg << ".\n";
            displayUsage();
            return 1;
        }
    }
    
    // Check mandatory arguments
    vector<string> required = {"-irf", "-ird", "-m", "-t"};
    for (const string& req : required) {
        if (args.find(req) == args.end()) {
            std::cerr << "Error: Missing mandatory argument " << req << ".\n";
            displayUsage();
            return 1;
        }
    }
    
    // Create parameters with bounds checking
    try {
        Parameters params(std::stoi(args["-m"]), std::stoi(args["-t"]));
        
        if (args.count("-ms")) params.match_score = std::stoi(args["-ms"]);
        if (args.count("-mp")) params.mismatch_penalty = std::stoi(args["-mp"]);
        if (args.count("-gp")) params.gap_penalty = std::stoi(args["-gp"]);
        if (args.count("-ep")) params.extend_penalty = std::stoi(args["-ep"]);
        if (args.count("-vb")) params.verbose = true;
        
        // Load input data
        ApproxPatternMatching assembler;
        vector<string> reads = loadStrings(args["-ird"]);
        string reference = buildReference(loadStrings(args["-irf"]));
        
        // Run the assembly process
        assembler.run(reference, reads, params);
        
        // Print mascot if requested
        if (args.count("-mascot")) {
            PrintMascot();
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}