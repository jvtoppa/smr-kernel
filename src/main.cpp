#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <algorithm>
#include <set>
#include "libsais.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <divsufsort.h>
#include <unistd.h>
#include <sdsl/lcp.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp.hpp>
#include <sdsl/rmq_support.hpp>
#include <climits>
#include <cmath>
#include <cstdlib>
using namespace std;
using namespace sdsl;
std::string K (std::string&);

struct rlz
{
    int startInR;
    int size;
    char last;
};

std::vector<std::pair<uint32_t, uint32_t>> periods(
	const std::string& T,
	const std::vector<int32_t>& SA,
	const std::vector<int32_t>& LCP 
);

std::vector<std::pair<uint32_t, uint32_t>> supermaximal_repeats(
	const std::string& T,
	const std::vector<int32_t>& SA,
	const std::vector<int32_t>& LCP 
);

std::string merge(
	const std::string& T,
	const std::vector<std::pair<uint32_t, uint32_t>>& p
);

std::string merge_simple(
	const std::string& T,
	const std::vector<std::pair<uint32_t, uint32_t>>& p
);


int recursion_depth(std::string T){

	uint32_t it = 0;
	uint32_t n = T.length();

	while(T.length()>0){

		T = K(T);
		++it;

	}

	return it;

}

vector<rlz> buildRLZ(const string& R, const string& S, int& bestBits) {
    int n = R.size();
    int m = S.size();
    
    if (n == 0 || m == 0) return {};
    
    cout << "Building suffix array (n=" << n << ")..." << endl;
    
    vector<int> sa(n);
    if (divsufsort((const unsigned char*)R.data(), (saidx_t*)sa.data(), n) != 0)
    {
        cerr << "Failed to build suffix array" << endl;
        return {};
    }
    

    
    cout << "Parsing S (m=" << m << ")..." << endl;
    
    vector<rlz> factors;
    factors.reserve(m / 10);
    
    int pos = 0;
    
    while (pos < m) {

        
        int l = 0, r = n-1;
        int len = 0;
        
        while (l <= r && pos + len < m && len < n) {
            char c = S[pos + len];
            
            int lo = l, hi = r;
            int left = n; 
            while (lo <= hi) {
                int mid = lo + (hi - lo) / 2;
                if (mid < 0 || mid >= n) break;
                
                int sa_pos = sa[mid];
                char rc = '\0';
                if (sa_pos + len < n) {
                    rc = R[sa_pos + len];
                }
                
                if (rc >= c) {
                    if (rc == c) left = min(left, mid);
                    hi = mid - 1;
                } else {
                    lo = mid + 1;
                }
            }
            
            lo = l, hi = r;
            int right = -1;
            while (lo <= hi) {
                int mid = lo + (hi - lo) / 2;
                if (mid < 0 || mid >= n) break;
                
                int sa_pos = sa[mid];
                char rc = '\0';
                if (sa_pos + len < n) {
                    rc = R[sa_pos + len];
                }
                
                if (rc <= c) {
                    if (rc == c) right = max(right, mid);
                    lo = mid + 1;
                } else {
                    hi = mid - 1;
                }
            }
            
            if (left > right || left == n || right == -1) {
                break;
            }
            
            l = left;
            r = right;
            len++;
        }
        
        if (l > r || l >= n || r < 0) {
            factors.push_back({-1, 0, S[pos]});
            pos++;
            continue;
        }
        
        int max_len = 0;
        int best_pos = -1;
        
        l = max(0, l);
        r = min(n-1, r);
        
        int max_to_check = 10000;
        if (r - l + 1 > max_to_check) {
            int step = (r - l) / max_to_check;
            for (int i = l; i <= r && i < n; i += max(1, step)) {
                int sa_pos = sa[i];
                int match_len = min(len, n - sa_pos);
                
                bool valid = true;
                for (int j = 0; j < match_len; j++) {
                    if (pos + j >= m || sa_pos + j >= n || 
                        S[pos + j] != R[sa_pos + j]) {
                        valid = false;
                        break;
                    }
                }
                
                if (!valid) continue;
                
                while (pos + match_len < m && 
                       sa_pos + match_len < n &&
                       S[pos + match_len] == R[sa_pos + match_len]) {
                    match_len++;
                }
                
                if (match_len > max_len) {
                    max_len = match_len;
                    best_pos = sa_pos;
                }
            }
        }
        else
        {
            for (int i = l; i <= r && i < n; i++) {
                int sa_pos = sa[i];
                int match_len = min(len, n - sa_pos);
                
                bool valid = true;
                for (int j = 0; j < match_len; j++) {
                    if (pos + j >= m || sa_pos + j >= n || 
                        S[pos + j] != R[sa_pos + j]) {
                        valid = false;
                        break;
                    }
                }
                
                if (!valid) continue;
                
                // Extend match
                while (pos + match_len < m && 
                       sa_pos + match_len < n &&
                       S[pos + match_len] == R[sa_pos + match_len]) {
                    match_len++;
                }
                
                if (match_len > max_len) {
                    max_len = match_len;
                    best_pos = sa_pos;
                }
            }
        }
        
        if (max_len == 0) {
            factors.push_back({-1, 0, S[pos]});
            pos++;
        } else {
            factors.push_back({
                best_pos,
                max_len,
                (pos + max_len < m) ? S[pos + max_len] : '\0'
            });
            pos += max_len + 1;
        }
    }
    
 


	size_t factor_count = factors.size();
	size_t ref_len = R.size();
    size_t input_len = S.length();
    
    if (factor_count > 0) {
        const double DNA_BITS = 8.0;
        double offset_bits = std::max(std::ceil(std::log2(ref_len)), 1.0);
        double avg_len = (double)input_len / factor_count;
        double length_bits = std::max(std::ceil(std::log2(avg_len)), 1.0);
        
        double bits_per_factor = offset_bits + length_bits + 2.0 + DNA_BITS;
        double rlz_bits = factor_count * bits_per_factor;
        double total_bits = ref_len * DNA_BITS + rlz_bits;
        
        std::cout << "Ref bits: " << ref_len * DNA_BITS << "\n";
		std::cout << "RLZ bits: " << rlz_bits << "\n";
        std::cout << "Total bits: " << total_bits << "\n";
		std::cout << "Size of RLZ: " << factor_count << "\n";
		std::cout << "Size of reference: " << ref_len << "\n";
		std::cout << "Total size:" << ref_len + factor_count << "\n\n";

        bestBits = total_bits;
    }
	/*
    for (size_t i = 0; i < factors.size(); i++)
	{
		cout << factors[i].startInR << " " << factors[i].size << " " << factors[i].last << "\n";
	}
	*/
	
    return factors;
}

void process(std::string T) {
    //String<char> seqT = T;
   // std::cout << "iteration\tkernel length" << std::endl;
	std::string originalT = T;
	std::cout << "Iteration 0"  << std::endl;

    T = K(T);
    //String<char> kT = T;
    int bestBits;
	std::cout << "Starting to build RLZ for Iteration 0...\n";
    buildRLZ(T, originalT, bestBits);


    uint32_t it = 1;
    while (T.size() > 0) {
        T = K(T);
        //  std::cout <<"Iteration " <<it << std::endl;
		
        
        std::cout << "Starting to build RLZ for Iteration " << it << "...\n";
        int prevBest = bestBits;
        buildRLZ(T, originalT, bestBits);
        if(bestBits > prevBest) break;
        ++it;
    }
}


void brute_force_max_ratio(){

	std::cout << "n\tratio\tdifference\tstring\tkernel\n";

	std::string max_str;
	std::string max_ker;

	for(int n=1;n<32;++n){

		int max_k=0;
		std::string s(n,'a');
		for(uint32_t i=0;i < (uint32_t(1)<<n) ;++i){
			uint32_t x=i;
			for(int k = 0;k < n; ++k){
				s[k] = x%2 ? 'a' : 'b';
				x = x>>1;
			}
			std::string k = K(s);

			if(k.length() > max_k){
				max_str = s;
				max_ker = k;
				max_k = k.length();
			}

		}

		std::cout << n << "\t" << double(max_k)/n << "\t" << (n-max_k) << "\t" << max_str << "\t" << max_ker << std::endl;
	
	}

}

int main(int argc, char* argv[]) {
    
	std::string T;
	
	if (argc < 2)
	{
		std::ostringstream buffer;
		buffer << std::cin.rdbuf();
		T = buffer.str();
		cout << "Read " << T.size() << "charaters. \n";
	}
    else
	{
		unsigned long max_chars = std::stoul(argv[1]);
		
		
		T.resize(max_chars);
	
		std::cin.read(&T[0], max_chars);
		
		T.resize(std::cin.gcount());

	}
    
    process(T);

    return 0;
}



std::string K(std::string& T){
    uint32_t n = T.size();
	std::cout << "Starting K...\n";
    
	// step 1 compute SA, LCP
    
	std::vector<int32_t> SA(n);
	std::vector<int32_t> LCP(n);
	
	{
	std::vector<int32_t> PLCP(n);
	std::cout << "Building SA and LCP of T (" << n << " bytes)..." << std::endl;
	libsais((uint8_t*)T.data(), SA.data(), n, 0, nullptr);
	libsais_plcp((uint8_t*)T.data(), SA.data(), PLCP.data(), n);
	libsais_lcp(PLCP.data(), SA.data(), LCP.data(), n);
	std::cout << "Success!" << std::endl;
	}

	// compute period of each SMR of T

	std::cout << "Computing periods of SMR of T..." << std::endl;
	auto R = periods(T,SA,LCP);
	std::cout << "Success! found " << R.size() << " periods." << std::endl;
	//for(auto p:R) std::cout << p.first << "," << p.second << std::endl;

	std::string res;

	// step 3 extract characters avoiding overlaps

	return merge(T, R); 

}

std::string merge_simple(
	const std::string& T,
	const std::vector<std::pair<uint32_t, uint32_t>>& p
	) {

	std::string res;
	uint32_t i=0; //last char we appended is T[i-1]

	for(auto x : p){
		for(int j=std::max(x.first,i); j<=x.second;++j ) res += T[j];		
		i = x.second+1;
	}

	return res;

}

//assumes character 0 does not appear in T
std::string merge(
	const std::string& T,
	const std::vector<std::pair<uint32_t, uint32_t>>& p
	) {

	std::string res;
	if (p.empty()) return res;

	// Pre-allocate pi vector to reuse memory and avoid reallocations
	std::vector<int> pi;

	for (const auto& range : p) {
		uint32_t start = range.first;
		uint32_t end = range.second;
		if (start > end || start >= T.length()) continue;

		size_t m = end - start + 1;
		if (res.empty()) {
			res.append(T, start, m);
			continue;
		}

		size_t n = res.length();
		size_t max_possible_overlap = std::min(n, m);

		// Construct: [New Segment] + 0 + [Suffix of res]
		std::string combine;
		combine.reserve(m + 1 + max_possible_overlap);
		combine.append(T, start, m);
		combine.push_back(0);
		combine.append(res, n - max_possible_overlap, max_possible_overlap);

		// Compute KMP prefix function (linear time)
		pi.assign(combine.length(), 0);
		for (int i = 1, j = 0; i < (int)combine.length(); i++) {
			while (j > 0 && combine[i] != combine[j])
				j = pi[j - 1];
			if (combine[i] == combine[j])
				j++;
			pi[i] = j;
		}

		// The last element in pi gives the longest prefix of segment that is a suffix of res
		size_t overlap = pi.back();
		
		// Append only the non-overlapping suffix
		res.append(T, start + overlap, m - overlap);
	}
	return res;
}


/**
 * Computes the period of every supermaximal repeat.
 * Results are returned as (start_position, end_position) where end is inclusive.
 */
std::vector<std::pair<uint32_t, uint32_t>> periods(
	const std::string& T,
	const std::vector<int32_t>& SA,
	const std::vector<int32_t>& LCP
) {
	uint32_t n = static_cast<uint32_t>(T.size());
	std::vector<std::pair<uint32_t, uint32_t>> R;

	for (uint32_t i = 0; i < n; ++i) {
		// Identify an LCP-interval [i, j]
		if (i + 1 < n && LCP[i + 1] > 0) {
			uint32_t j = i + 1;
			uint32_t current_lcp = LCP[j];
			
			while (j + 1 < n && LCP[j + 1] == current_lcp) {
				j++;
			}

			// Local Maximum Check (Right-maximality)
			bool left_boundary = (i == 0 || LCP[i] < current_lcp);
			bool right_boundary = (j + 1 == n || LCP[j + 1] < current_lcp);

			if (left_boundary && right_boundary) {
				std::set<uint8_t> left_chars;
				bool all_distinct = true;
				uint32_t leftmost_pos = n;

				// Verify Left-Maximality and find leftmost position
				for (uint32_t k = i; k <= j; ++k) {
					uint32_t pos = SA[k];
					if (pos < leftmost_pos) leftmost_pos = pos;

					uint8_t c = (pos == 0) ? 0x00 : T[pos - 1];
					
					if (left_chars.count(c)) {
						all_distinct = false;
						break;
					}
					left_chars.insert(c);
				}

				if (all_distinct && left_chars.size() > 1) {
					// --- KMP-based period detection ---
					uint32_t len = current_lcp;
					uint32_t period_len = len;
					
					if (len > 1) {
						std::vector<uint32_t> pi(len, 0);
						for (uint32_t k = 1; k < len; k++) {
							uint32_t m = pi[k - 1];
							while (m > 0 && T[leftmost_pos + k] != T[leftmost_pos + m]) {
								m = pi[m - 1];
							}
							if (T[leftmost_pos + k] == T[leftmost_pos + m]) {
								m++;
							}
							pi[k] = m;
						}
						
						uint32_t longest_border = pi[len - 1];
						period_len = len - longest_border;
							
					}

					// Store as (start, end) where end is inclusive
					R.push_back({leftmost_pos, leftmost_pos + period_len - 1});
				}
			}
			i = j - 1; 
		}
	}

	// Sort results by start position
	std::sort(R.begin(), R.end());

	return R;
}

std::vector<std::pair<uint32_t, uint32_t>> supermaximal_repeats(
    const std::string& T,
    const std::vector<int32_t>& SA,
    const std::vector<int32_t>& LCP
) {
    uint32_t n = static_cast<uint32_t>(T.size());
    std::vector<std::pair<uint32_t, uint32_t>> R;

    for (uint32_t i = 0; i < n; ++i) {
        if (i + 1 < n && LCP[i + 1] > 0) {
            uint32_t j = i + 1;
            uint32_t current_lcp = LCP[j];

            while (j + 1 < n && LCP[j + 1] == current_lcp) {
                j++;
            }

            bool left_boundary = (i == 0 || LCP[i] < current_lcp);
            bool right_boundary = (j + 1 == n || LCP[j + 1] < current_lcp);

            if (left_boundary && right_boundary) {
                std::set<uint8_t> left_chars;
                bool all_distinct = true;
                uint32_t leftmost_pos = n;

                for (uint32_t k = i; k <= j; ++k) {
                    uint32_t pos = SA[k];
                    leftmost_pos = std::min(leftmost_pos, pos);

                    uint8_t c = (pos == 0) ? 0x00 : T[pos - 1];
                    if (!left_chars.insert(c).second) {
                        all_distinct = false;
                        break;
                    }
                }

                if (all_distinct && left_chars.size() > 1) {
                    R.push_back({leftmost_pos,
                                 leftmost_pos + current_lcp - 1});
                }
            }

            i = j - 1;
        }
    }

    std::sort(R.begin(), R.end());
    return R;
}
