#ifndef INCLUDE_H
#define INCLUDE_H

//intel hexl library 
#include "hexl/hexl.hpp"

//Random Number Generator
#include "RandomNumberGenerator/lac_param.h"
#include "RandomNumberGenerator/rand.hpp"


//LWE
#include "LWEscheme/LWE_32.hpp"
#include "LWEscheme/LWE_64.hpp"

//polynomial ring
#include "PolynomialRing/Poly_ring_64.hpp"
typedef vector<int64_t> polynomial;

//RLWE
#include "RLWEscheme/RLWE_64.hpp"
#include "RLWEscheme/LUT_64.hpp"

//test
#include "test/testRNG.hpp"
#include "test/testLWE.hpp"
#include "test/testPoly.hpp"
#include "test/testRLWE.hpp"


//C++
#include <iostream>
#include <vector>
#include <ctime>
#include<sys/time.h>
#include<time.h>
#include<unistd.h>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <condition_variable>
#include <chrono>
#include <thread>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <immintrin.h>
#include <string>
#include <memory>

#endif