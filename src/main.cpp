#include <iostream>
#include <vector>
#include <fstream>
#include <cstdint>
#include <algorithm>
#include "libsais.h"

void K (std::vector<uint8_t>&,std::vector<uint8_t>&);


int main(int argc, char* argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
		return 1;
	}

	std::ifstream file(argv[1], std::ios::binary | std::ios::ate);
	std::streamsize n = file.tellg();
	file.seekg(0, std::ios::beg);
	std::vector<uint8_t> T(n);
	std::vector<uint8_t> ker;
	
	file.read(reinterpret_cast<char*>(T.data()), n);

	uint32_t depth = 1;

	K(T,ker);

	return 0;
}



void K(std::vector<uint8_t>& T, std::vector<uint8_t>& ker){
	
	uint32_t n = T.size();

	std::vector<int32_t> SA(n);
	std::vector<int32_t> LCP(n);
	
	{
	std::vector<int32_t> PLCP(n);

	std::cout << "Building SA and LCP for " << n << " bytes..." << std::endl;

	libsais(T.data(), SA.data(), n, 0, nullptr);

	libsais_plcp(T.data(), SA.data(), PLCP.data(), n);

	libsais_lcp(PLCP.data(), SA.data(), LCP.data(), n);

	std::cout << "Success!" << std::endl;
	}

	std::vector<uint8_t> res;
	

}
