# P2N3 (Privacy-Preserving Neural Network of NTU) - hexl accerlation
## Introduction
P2N3 (**P**rivacy-**P**reserving **N**eural **N**etwork of **N**TU) is a C++ open-source library which implements an optimized fully homomorphic encryption (FHE) scheme for privacy-preserving neural networks.  
Our optimized FHE scheme enables us the ability to perform weighted sums and convolutions on the approximate LWE-based additive homomorphic encryption scheme, and to evaluate non-polynomial activations on FHEW ciphertexts. This FHE scheme has the following properties: 
- It can be applied to neural networks of arbitrary depth.
- It supports many kinds of widely used activations, such as the most popular ReLU. 
- When applied to inference of privacy-preserving neural networks, it is fast and accurate.

For its application in other fields such as pandemic modelling, please refer to our publication [[2]](https://ieeexplore.ieee.org/document/10376224). 

[[1]](https://ieeexplore.ieee.org/document/10398424) K.Y. Lam, X. Lu, L. Zhang, X. Wang, H. Wang and S.Q. Goh. “Efficient FHE-based Privacy-Enhanced Neural Network for Trustworthy AI-as-a-Service”, IEEE Transactions on Dependable and Secure Computing, IEEE.   
[[2]](https://ieeexplore.ieee.org/document/10376224) L. Zhang, X. Wang, J. Wang, R. Pung, H. Wang and K.Y. Lam. “An Efficient FHE-Enabled Secure Cloud-Edge Computing Architecture for IoMTs Data Protection With Its Application to Pandemic Modelling”, IEEE Internet of Things Journal, IEEE.


It is the implementation of publication [[1]](https://ieeexplore.ieee.org/document/10398424) using [hexl](https://github.com/intel/hexl).  
If you prefer not using hexl, please refer to repo [P2N3](https://github.com/HintSight-Technology/P2N3).  

## Prerequisites
- [OpenSSL](https://www.openssl.org/)
- [hexl](https://github.com/intel/hexl)
## Installation
Supported platforms: Linux and MacOS.  

```
mkdir build && cd build
cmake ..
make
```
 

## Document
Please refer to [P2N3_Programmer_Guide.pdf](https://github.com/HintSight-Technology/P2N3/blob/main/P2N3_Programmer_Guide.pdf) from repo [P2N3](https://github.com/HintSight-Technology/P2N3). 

## Examples
We offer some applications of our FHE scheme. It includes inference on MNIST dataset and facial recognition.   
REMARK: All the test data are stored in the form of extracted feature vector (in plaintext). Please refer to the detailed experiments in our publication [[1]](https://ieeexplore.ieee.org/document/10398424).  

[[1]](https://ieeexplore.ieee.org/document/10398424) K.Y. Lam, X. Lu, L. Zhang, X. Wang, H. Wang and S.Q. Goh. “Efficient FHE-based Privacy-Enhanced Neural Network for Trustworthy AI-as-a-Service”, IEEE Transactions on Dependable and Secure Computing, IEEE.   

## License
This software is distributed under the BSD-3-Clause-Clear license. 
