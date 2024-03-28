//code for testing ntt, intt, pmul in FasterNTT
#include <iostream>
using namespace std;

/*
extern void RNG(){
  unsigned int len_out = 10;
  unsigned char seed[SEED_LEN];
  uint64_t out[len_out];
  random_bytes(seed,SEED_LEN);
  uint64_t p = 1099511627776 + 1;
  int a = gen_uniform_int64(out, len_out, p, 41, seed);
   for (int i=0;i<len_out;++i){
    printf("%lld \n",out[i]);
   }
}
*/

inline int64_t modq_64(int64_t number, int64_t q) {
  while (number < 0) {
    number += q;
  }
  /*
  if (q - number >> 63) {
    number %= q;
  }
  else;
  if (number >= (q>>1)) {
    number -= q;
  }
  */
  while (number >= q){
    number -= q;
  }
  if(number >= q/2){
    number -= q;
  }
  return number;
}

// key generation algorithm
// INPUT: dimension n
// OUTPUT: random binary vector x
extern int* LWE64_KeyGen(int n){
  unsigned int len_out = n;
  unsigned char seed[SEED_LEN];
  int *x = new int [n];
  random_bytes(seed,SEED_LEN);
  int r = gen_bernoulli(x,len_out,seed);
  if (r == 1){
    cout <<"Error in LWE Key Generation: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in LWE Key Generation: Length" <<endl;
  }
  return x;
}

//encryption algorithm
// INPUT: modulus q, dimension n, variance 2^(-k), k >= 1, plaintext m, LWE key x
// OUTPUT：LWE ciphertext (vec(a),b)
extern vector<int64_t> LWE64_Enc(int64_t q, int n, int k, int64_t m, const int* x){
  vector<int64_t> ct;

  //generate random array a
  double logq=log(q)/log(2.0);
  int int_logq = (int) logq;
  if(logq > (double)int_logq){
    int_logq++;
  }

  unsigned int len_out = n;
  unsigned char seed[SEED_LEN];
  int64_t *a = new int64_t [len_out];
  random_bytes(seed,SEED_LEN);
  int r = gen_uniform_int64(a, len_out, q, int_logq, seed);
  if (r == 1){
    cout <<"Error in generation random array a of LWE Encryption: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in generation random array a of LWE encryption: modulus" <<endl;
  }

  //generate error array e
  int e[8];
  random_bytes(seed,SEED_LEN);
  r = gen_ternary_var(e, 8, k, seed);
  if (r == 1){
    cout <<"Error in generation a random error of LWE Encryption: NULL "<<endl;
  }
  if(r == 2){
    cout <<"Error in generation a random error of LWE Encryption: modulus" <<endl;
  }

  // b = ( - <a,x> + m + e) mod q
  int64_t b=0;
  for(int i=0;i<n;++i){
    int64_t ax = (a[i]-q/2)*(int64_t)x[i];
    b += ax; 
    b = modq_64(b,q);
    ct.push_back(a[i]-q/2);
  }
  b *= (-1);
  b = modq_64(b,q);
  b += (m+(int64_t)e[0]);
  b = modq_64(b,q);
  ct.push_back(b);
  //return ciphertext
  return ct;

}


// LWE Decryption
// INPUT: modulus q, dimension n, LWE ciphertext c = (vec(a),b), LWE key x
// OUTPUT: Decryption (b + <a,x>) mod q
extern int64_t LWE64_Dec(int64_t q, int n, const vector<int64_t>& c, const int* x) {
  int64_t ip2 = c[n];
  for (int i = 0; i < n; ++i) {
    if (x[i] == -1) {
      ip2 -= c[i];
    }
    else if (x[i] == 1) {
      ip2 += c[i];
    }
    else;
    ip2 = modq_64(ip2, q);
  }
  int64_t m = modq_64(ip2 , q);
  return m;
}

extern int64_t LWE64_Dec_q(int64_t q, int n, const vector<uint64_t>& c, const int* x) {
  int64_t ip2 = (int64_t)c[n];
  for (int i = 0; i < n; ++i) {
    if (x[i] == -1) {
      ip2 -= (int64_t)c[i];
    }
    else if (x[i] == 1) {
      ip2 += (int64_t)c[i];
    }
    else;
    ip2 = modq_64(ip2, q);
  }
  int64_t m = modq_64(ip2 , q);
  return m;
}


// Add two LWE ciphertexts  
// INPUT: modulus q, dimension n, LWE ciphertext ct1 and ct2
// OUTPUT: LWE ciphertext (ct1 + ct2) mod q
vector<int64_t> LWE64_Add_ct(int64_t q, int n, const vector<int64_t>& ct1, const vector<int64_t>& ct2) {
  vector<int64_t> ct(n+1,0);
  for (register int i = 0; i <= n; ++i) {
    ct[i]=modq_64(ct1[i] + ct2[i],q); 
  }
  return ct;
}

vector<uint64_t> LWE64_Add_ct_q(int64_t q, int n, const vector<uint64_t>& ct1, const vector<uint64_t>& ct2) {
  uint64_t modulus = (uint64_t)q;
  vector<uint64_t> output(n+1);
  //output.reserve(n+1);
  intel::hexl::EltwiseAddMod(output.data(), ct1.data(), ct2.data(), n+1, modulus);
  return output;
}

// Add a LWE ciphertext and a plaintext number k  
// INPUT: modulus q, dimension n, LWE ciphertext ct1, plaintext number k
// OUTPUT: LWE ciphertext of dec(ct1) + k
vector<int64_t> LWE64_Plain_Add_ct(int64_t q, int n, const vector<int64_t>& ct1, int64_t k) {
  vector<int64_t> ct=ct1;
  ct[n] = modq_64(ct[n]+k, q);
  return ct;
}

vector<uint64_t> LWE64_Plain_Add_ct_q(int64_t q, int n, const vector<uint64_t>& ct1, int64_t k) {
  vector<uint64_t> ct(ct1);
  uint64_t temp;
  if(k >= 0) {
    temp = ct[n]+(uint64_t)k;
  }
  else{
    temp = ct[n]+(uint64_t)(k+q);
  }
  if(temp >= (uint64_t)q){
    ct[n] = temp-(uint64_t)q;
  }
  else{
    ct[n] = temp;
  }
  return ct;
}


// Multiply a ciphertext and a plaintext number k
// INPUT: modulus q, dimension n, ciphertext ct1, plaintext number k
// OUTPUT: LWE ciphertext of k*dec(ct1)
vector<int64_t> LWE64_Plain_Multi_ct(int64_t q, int n, const vector<int64_t>& ct1, int64_t k) {
  vector<int64_t> ct;
  int size = ct1.size();
  for (register int i = 0; i <size; ++i) {
    ct.emplace_back(modq_64(ct1[i] * k,q));
  }
  return ct;
}

vector<uint64_t> LWE64_Plain_Multi_ct_q(int64_t q, int n, const vector<uint64_t>& ct1, int64_t k) {
  uint64_t modulus = (uint64_t)q;
  if(k >= 0){
    vector<uint64_t> op2(n+1,(uint64_t)k);
    vector<uint64_t> ct(n+1);
    intel::hexl::EltwiseMultMod(ct.data(), ct1.data(), op2.data(), n+1, modulus,1);
    return ct;
  }
  if(k < 0){
    uint64_t tempk = (uint64_t)(-1*k);
    vector<uint64_t> op2(n+1,tempk);
    vector<uint64_t> ct(n+1);
    intel::hexl::EltwiseMultMod(ct.data(), ct1.data(), op2.data(), n+1, modulus,1);
    for(int i = 0 ; i < n+1 ; ++i){
      ct[i] = modulus-ct[i];
    }
    return ct;
  }
}


// modulus Switching, mod q1 -> mod q2  
// INPUT: first modulus q1, second modulus q2, dimension n, LWE ciphertext ct with modulus q1
// OUTPUT: LWE ciphertext with modulus q2
vector<int64_t> LWE64_Rounding(int64_t q1, int64_t q2, int n, const vector<int64_t>& ct) {
  vector<int64_t> ans=ct;
  for (register int i = 0; i <n+1; ++i) {
    ans[i] = q2 * ans[i] / q1;
  }
  return ans;
}


