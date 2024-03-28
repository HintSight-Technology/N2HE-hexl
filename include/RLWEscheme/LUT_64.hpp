#include <iostream>
#include <ctime>
#include<sys/time.h>
#include<time.h>
using namespace std;

typedef vector<int64_t> polynomial;
extern intel::hexl::NTT* g_ntt;

 

// ====================================================== Evaluation key based on RGSW ====================================
// LWE secret key is encrypted under rlwe secret key


// RGSW encryption of x1 * x2
extern vector<vector<polynomial>> key_encrypt_1(int n, int64_t q, int k,int var, int x1,int x2, const polynomial & s, int64_t b, int logb) {
  vector<vector<polynomial>> ans;
  //zero polynomial
  polynomial zero(2,0);
  //zero.exp = 0;
  zero[0] = 0;
  zero[1] = 0;
  //1x^0 polynomial
  polynomial one(2,0);
  //one.exp = 0;
  one[0] = 0;
  one[1] = 1;

  if(n == 2048){
    if (x1*x2 == 0) {
     ans = RGSWct_2048(n, q, k, var, zero, s,b, logb);
    }
   else {
     ans = RGSWct_2048(n, q, k, var, one, s,b, logb);
   }
  }
  else{
    cout <<"undefined input n in key_encrypt_1. "<<endl;
  }
  
  return ans;
}

// RGSW encryption of x1 * (1 - x2)
extern vector<vector<polynomial>> key_encrypt_2(int n, int64_t q, int k, int var, int x1,int x2, const polynomial & s, int64_t b, int logb) {
  vector<vector<polynomial>> ans;
  //zero polynomial
  polynomial zero(2,0);
  //zero.exp = 0;
  zero[0] = 0;
  zero[1] = 0;
  //1x^0 polynomial
  polynomial one(2,0);
  //one.exp = 0;
  one[0] = 0;
  one[1] = 1;

  if(n == 2048){
    if (x1*(1-x2) == 1) {
     ans = RGSWct_2048(n, q, k, var, one, s,b,logb);
    }
    else {
      ans = RGSWct_2048(n, q, k, var, zero, s,b,logb);
    }
  }
    else{
    cout <<"undefined input n in key_encrypt_2. "<<endl;
  }
  
  
  return ans;
}

// RGSW encryption of (1 - x1) * x2
extern vector<vector<polynomial>> key_encrypt_3(int n, int64_t q, int k, int var, int x1, int x2, const polynomial& s, int64_t b, int logb) {
  vector<vector<polynomial>> ans;
  //zero polynomial
  polynomial zero(2, 0);
  //zero.exp = 0;
  zero[0] = 0;
  zero[1] = 0;
  //1x^0 polynomial
  polynomial one(2, 0);
  //one.exp = 0;
  one[0] = 0;
  one[1] = 1;

  if(n==2048){
    if (x2 * (1 - x1) == 1) {
      ans = RGSWct_2048(n, q, k, var, one, s, b, logb);
    }
    else {
      ans = RGSWct_2048(n, q, k, var, zero, s, b, logb);
    }
  }
  else{
    cout <<"undefined input n in key_encrypt_3. "<<endl;
  }
  
  return ans;
}


// RGSW encryption of x1 * x2
extern vector<vector<uint64_t>> key_encrypt_1_ntt(int n, int64_t q, int k,int var, int x1,int x2, const polynomial & s, int64_t b, int logb) {
  vector<vector<uint64_t>> ans;
  //zero polynomial
  polynomial zero(2,0);
  //zero.exp = 0;
  zero[0] = 0;
  zero[1] = 0;
  //1x^0 polynomial
  polynomial one(2,0);
  //one.exp = 0;
  one[0] = 0;
  one[1] = 1;

  if(n == 2048){
    if (x1*x2 == 0) {
     ans = RGSWct_2048_ntt(n, q, k, var, zero, s,b, logb);
    }
   else {
     ans = RGSWct_2048_ntt(n, q, k, var, one, s,b, logb);
   }
  }
  else{
    cout <<"undefined input n in key_encrypt_1_ntt. "<<endl;
  }
  
  return ans;
}

// RGSW encryption of x1 * (1 - x2)
extern vector<vector<uint64_t>> key_encrypt_2_ntt(int n, int64_t q, int k, int var, int x1,int x2, const polynomial & s, int64_t b, int logb) {
  vector<vector<uint64_t>> ans;
  //zero polynomial
  polynomial zero(2,0);
  //zero.exp = 0;
  zero[0] = 0;
  zero[1] = 0;
  //1x^0 polynomial
  polynomial one(2,0);
  //one.exp = 0;
  one[0] = 0;
  one[1] = 1;

  if(n == 2048){
    if (x1*(1-x2) == 1) {
     ans = RGSWct_2048_ntt(n, q, k, var, one, s,b,logb);
    }
    else {
      ans = RGSWct_2048_ntt(n, q, k, var, zero, s,b,logb);
    }
  }
    else{
    cout <<"undefined input n in key_encrypt_2_ntt. "<<endl;
  }
  
  return ans;
}

// RGSW encryption of (1 - x1) * x2
extern vector<vector<uint64_t>> key_encrypt_3_ntt(int n, int64_t q, int k, int var, int x1, int x2, const polynomial& s, int64_t b, int logb) {
  vector<vector<uint64_t>> ans;
  //zero polynomial
  polynomial zero(2, 0);
  //zero.exp = 0;
  zero[0] = 0;
  zero[1] = 0;
  //1x^0 polynomial
  polynomial one(2, 0);
  //one.exp = 0;
  one[0] = 0;
  one[1] = 1;

  if(n==2048){
    if (x2 * (1 - x1) == 1) {
      ans = RGSWct_2048_ntt(n, q, k, var, one, s, b, logb);
    }
    else {
      ans = RGSWct_2048_ntt(n, q, k, var, zero, s, b, logb);
    }
  }
  else{
    cout <<"undefined input n in key_encrypt_3_ntt. "<<endl;
  }
  
  return ans;
}



// ======================================================  Look-up-Table  ==========================================================================

// Activation Function 
// ReLU
extern int64_t Evaluate_T(int64_t x) {
  return (x > 0 ? x : 0);
}

//construct f based on Evaluate_T
// INPUT: dimension n, scaler delta, modulus q
// OUTPUT: Blindrotate polynomial f
extern polynomial p_f(int n, int64_t delta, int64_t q) {

  // construct LUT polynomial f in R_q
  // f.exp = n - 1;
  polynomial f(n + 1,0);
  
  f[0] = (int64_t)n - 1;
  int64_t temp_coef;

  //f: divide q to n parts. 
  for (int j = 0; j < n; j++) {
    if (j == 0) {
      temp_coef = Evaluate_T(0);
    }
    else if (j < (n / 2 )) {
      temp_coef = delta*Evaluate_T(-1 * (j * q / (4 * (int64_t)n)));

    }
    else {
      temp_coef = (-1)* delta * Evaluate_T(((int64_t)n - j) * q / (4 * (int64_t)n));
    }
    f[j + 1] = temp_coef;
  }

  return f;
}

extern polynomial p_f0(int n, int64_t delta, int64_t q) {

  // construct LUT polynomial f in R_q
  // f.exp = n - 1;
  polynomial f(n + 1,0);
  
  f[0] = (int64_t)n - 1;
  int64_t temp_coef;

  //f: divide q to n parts. 
  for (int j = 0; j < n; j++) {
    if (j == 0) {
      temp_coef = Evaluate_T(0);
    }
    else if (j < (n / 2 )) {
      temp_coef = delta*Evaluate_T(-1 * (j * q / (2 * (int64_t)n)));

    }
    else {
      temp_coef = (-1)* delta * Evaluate_T(((int64_t)n - j) * q / (2 * (int64_t)n));
    }
    f[j + 1] = temp_coef;
  }

  return f;
}



// LUT
// Input LWE ciphertext ct_LWE of m
// OUTPUT RLWE ciphertext of X^m * f
//template<int NN, int qNum, enum NTTSelector NTT, enum ImplSelector Impl>
extern vector<int64_t> LUT_2048(int index, int q1, int64_t q, int n, int64_t delta, int k, int N, vector<int> const &ct_LWE, 
  const vector<vector<vector<uint64_t>>> &_ek0, const vector<vector<vector<uint64_t>>> &_ek1, const vector<vector<vector<uint64_t>>> &_ek2, 
  polynomial const& f, int64_t b, int logb) {
  // INPUT: id of the current node index
  // INPUT: modulus q1 and q, dimension n and N, scaler delta, k = log(q)
  // INPUT: Look-up LWE ciphertext ct_LWE
  // INPUT: RGSW Evaluation Key ek0,ek1,ek2, RLWE key s
  // INPUT: Blindrotate polynomial f
  // INPUT: decompostion base b, log(b)
  // OUTPUT: RLWE ciphertext X^(dec(ct_LWE)) * f
  struct timeval tstart1, tend1;
//  gettimeofday(&tstart1,NULL);
  //insert index
  //using R=Rq<NN, qNum, NTT, Impl>;

  uint64_t modulus = (uint64_t)q;
  //intel::hexl::NTT ntta(n, modulus);

  
  //polynomial p_index(2, 0);
  //p_index[1] = index;

  //change the ciphertext modulus of LWE ct from q to 2n
  
  vector<int64_t> a_n(N);
  
  int64_t b_n = 2 * (int64_t)n * ct_LWE[N] / q1;
  int64_t b_n2 = 2 * (int64_t)n * ct_LWE[N] / (q1/2);
  if (b_n2 - 2*b_n >= 1 && b_n2>0) {
    b_n++;
  }
  else if (b_n2 - 2*b_n <= -1 && b_n2 < 0) {
    b_n--;
  }
  
  //int64_t b_n = ct_LWE[N];


  for (int i = 0; i < N; ++i) {
    a_n[i]= 2 * (int64_t)n * ct_LWE[i] / (q1);
    int64_t tempani= 2 * (int64_t)n * ct_LWE[i] / (q1/2);
    if (tempani - 2*a_n[i] >= 1 && tempani>0) {
      a_n[i]++;
    }
    else if (tempani - 2*a_n[i] <= -1 && tempani < 0) {
      a_n[i]--;
    }
    
  }

  //initialize AC0=(f*X^b,0)
  //vector<polynomial> AC0(2);
  vector<vector<uint64_t>> AC0(2, vector<uint64_t>(n,0));
  //compute X^{b_n mod n}
  
  int64_t temp = -1;
  if (b_n >= 0) {
    temp = b_n;
  }
  else {
    while (temp < 0) {
      temp = b_n + n;
    }
  }
  polynomial xb((int)temp+2,0);
  xb[0] = temp;
  if (b_n >= 0) {
    //xb.exp = b_n;
    xb[(int)temp+1] = 1;
  }
  else if (b_n < 0) {
    //xb.exp = n + b_n;
    xb[(int)temp+1] = -1;
  }
  else;
  //compute AC[0]=0X^0, AC0[1]=f*X^b

//  gettimeofday(&tend1,NULL);
//  double time_init = (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);
//  cout <<"LUT init time: "<<time_init<<" us."<<endl;

//  gettimeofday(&tstart1,NULL);

  AC0[1]=(NTTMul_u(f, xb, n, q));

//  gettimeofday(&tend1,NULL);
//  double time_init2 = (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);
//  cout <<"LUT init2 time: "<<time_init2<<" us."<<endl;

  double time_xa = 0.0;
  double time_xaek = 0.0;
  double time_nttac = 0.0;
  double time_odot = 0.0;
  double time_update = 0.0;

  for(int i = 0; i < N/2 ; ++i){

 //   gettimeofday(&tstart1,NULL);

    //generate ntt(x^a-1)
    //R xa1,xa2,xa3;
    vector<uint64_t> xa1;
    xa1.reserve(n);
    vector<uint64_t> xa2;
    xa2.reserve(n);
    vector<uint64_t> xa3;
    xa3.reserve(n);

    XA_minus_1_2048(n,modulus,a_n[2*i]+a_n[2*i+1],xa1);
    XA_minus_1_2048(n,modulus,a_n[2*i],xa2);
    XA_minus_1_2048(n,modulus,a_n[2*i+1],xa3);

    //gettimeofday(&tstart1,NULL);

    g_ntt->ComputeForward(xa1.data(), xa1.data(), 1, 1);
    g_ntt->ComputeForward(xa2.data(), xa2.data(), 1, 1);
    g_ntt->ComputeForward(xa3.data(), xa3.data(), 1, 1); 

  //  gettimeofday(&tend1,NULL);
  //  time_xa += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);

     //compute ntt(X^a-1)*ntt(ek), store in Xa_times_ek

  //  gettimeofday(&tstart1,NULL);

    vector<vector<vector<uint64_t>>> Xa_times_ek(2*(k/logb),vector<vector<uint64_t>>(2, vector<uint64_t>(n)));

    vector<uint64_t> tempxaek1;
        tempxaek1.reserve(n);
        vector<uint64_t> tempxaek2;
        tempxaek2.reserve(n);

    for(int j1 = 0; j1 < 2*(k/logb) ; ++j1){
      for(int j2 = 0; j2 < 2 ; ++j2){


        intel::hexl::EltwiseMultMod(tempxaek1.data(), xa1.data(), _ek0[i][j1*2+j2].data(), n, modulus, 1);
        intel::hexl::EltwiseMultMod(tempxaek2.data(), xa2.data(), _ek1[i][j1*2+j2].data(), n, modulus, 1);

        intel::hexl::EltwiseAddMod(tempxaek1.data(), tempxaek1.data(), tempxaek2.data(), n, modulus);
        intel::hexl::EltwiseMultMod(tempxaek2.data(), xa3.data(), _ek2[i][j1*2+j2].data(), n, modulus, 1);
        intel::hexl::EltwiseAddMod(Xa_times_ek[j1][j2].data(), tempxaek1.data(), tempxaek2.data(), Xa_times_ek[j1][j2].size(), modulus);

      }
    }

    vector<uint64_t>().swap(xa1);
    vector<uint64_t>().swap(xa2);
    vector<uint64_t>().swap(xa3);
    vector<uint64_t>().swap(tempxaek1);
    vector<uint64_t>().swap(tempxaek2);

  //  gettimeofday(&tend1,NULL);
  //  time_xaek += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);

    //decompose AC0 to 2*k/logb

  //  gettimeofday(&tstart1,NULL);
    vector<vector<uint64_t>> de_AC00 = bit_poly_u(k,AC0[0],n,q,b,logb);
    vector<vector<uint64_t>> de_AC01 = bit_poly_u(k,AC0[1],n,q,b,logb);

 //   gettimeofday(&tend1,NULL);
 //   time_nttac += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);

    vector<vector<uint64_t>> ntt_de_AC00(k/logb, vector<uint64_t>(n));
    vector<vector<uint64_t>> ntt_de_AC01(k/logb, vector<uint64_t>(n));


    for(int j1 = 0; j1 < k/logb ; ++j1){
      
      g_ntt->ComputeForward(de_AC00[j1].data(), de_AC00[j1].data(), 1, 1);
      g_ntt->ComputeForward(de_AC01[j1].data(), de_AC01[j1].data(), 1, 1);
      ntt_de_AC01[j1] = move(de_AC01[j1]);
      ntt_de_AC00[j1] = move(de_AC00[j1]);
    }

   // vector<vector<uint64_t>>().swap(de_AC00);
  //  vector<vector<uint64_t>>().swap(de_AC01);
    

   // gettimeofday(&tend1,NULL);
   // time_nttac += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);


  //  gettimeofday(&tstart1,NULL);
    //b <> beta = ntt_de_AC01 <> Xa_times_ek[first k/logb]
    //vector<R> b_beta;
    vector<vector<uint64_t>> b_beta(2,vector<uint64_t>(n));

    vector<uint64_t> tempac01ek0;
    tempac01ek0.reserve(n);
    vector<uint64_t> tempac01ek1;
    tempac01ek1.reserve(n);

    intel::hexl::EltwiseMultMod(tempac01ek0.data(), ntt_de_AC01[0].data(), Xa_times_ek[0][0].data(), n, modulus, 1);
    intel::hexl::EltwiseMultMod(tempac01ek1.data(), ntt_de_AC01[0].data(), Xa_times_ek[0][1].data(), n, modulus, 1);

    for(int j1 = 1 ; j1 < k/logb ; ++j1){
      vector<uint64_t> ac01ek0;
      ac01ek0.reserve(n);
    //  vector<uint64_t> sum0;
    //  sum0.reserve(n);

      intel::hexl::EltwiseMultMod(ac01ek0.data(), ntt_de_AC01[j1].data(), Xa_times_ek[j1][0].data(), n, modulus, 1);
      intel::hexl::EltwiseAddMod(tempac01ek0.data(), tempac01ek0.data(), ac01ek0.data(), n, modulus);

    //  tempac01ek0 = move(sum0);

    //  vector<uint64_t> ac01ek1;
    //  ac01ek1.reserve(n);
   //   vector<uint64_t> sum1;
    //  sum1.reserve(n);

      intel::hexl::EltwiseMultMod(ac01ek0.data(), ntt_de_AC01[j1].data(), Xa_times_ek[j1][1].data(), n, modulus, 1);
      intel::hexl::EltwiseAddMod(tempac01ek1.data(), tempac01ek1.data(), ac01ek0.data(), n, modulus);

     // tempac01ek1 = move(sum1);

    }

    //b_beta[0]=tempac01ek0;
    b_beta[0] = move(tempac01ek0);
   // b_beta[1]=tempac01ek1;
    b_beta[1] = move(tempac01ek1);

    vector<vector<uint64_t>>().swap(ntt_de_AC01);


    vector<vector<uint64_t>> a_alpha(2,vector<uint64_t>(n));

    vector<uint64_t> tempac00ek0;
    tempac00ek0.reserve(n);
    vector<uint64_t> tempac00ek1;
    tempac00ek1.reserve(n);

    intel::hexl::EltwiseMultMod(tempac00ek0.data(), ntt_de_AC00[0].data(), Xa_times_ek[k/logb][0].data(), n, modulus, 1);
    intel::hexl::EltwiseMultMod(tempac00ek1.data(), ntt_de_AC00[0].data(), Xa_times_ek[k/logb][1].data(), n, modulus, 1);

    for(int j1 = 1 ; j1 < k/logb ; ++j1){

      vector<uint64_t> ac00ek0;
      ac00ek0.reserve(n);
      //vector<uint64_t> sum0;
     // sum0.reserve(n);

      intel::hexl::EltwiseMultMod(ac00ek0.data(), ntt_de_AC00[j1].data(), Xa_times_ek[j1+k/logb][0].data(), n, modulus, 1);
      intel::hexl::EltwiseAddMod(tempac00ek0.data(), tempac00ek0.data(), ac00ek0.data(), n, modulus);

      //tempac00ek0 = sum0;
     // tempac00ek0 = move(sum0);

     // vector<uint64_t> ac00ek1;
    //  ac00ek1.reserve(n);
   //   vector<uint64_t> sum1;
    //  sum1.reserve(n);

      intel::hexl::EltwiseMultMod(ac00ek0.data(), ntt_de_AC00[j1].data(), Xa_times_ek[j1+k/logb][1].data(), n, modulus, 1);
      intel::hexl::EltwiseAddMod(tempac00ek1.data(), tempac00ek1.data(), ac00ek0.data(), n, modulus);

     // tempac00ek1 = sum1;
     // tempac00ek1 = move(sum1);


    }

    //a_alpha[0]=tempac00ek0;
    a_alpha[0] = move(tempac00ek0);
   // a_alpha[1]=tempac00ek1;
    a_alpha[1] = move(tempac00ek1);

    vector<vector<uint64_t>>().swap(ntt_de_AC00);
    vector<vector<vector<uint64_t>>>().swap(Xa_times_ek);

    //a_alpha + b_beta

   // vector<vector<uint64_t>> new_AC0(2,vector<uint64_t>(n));
   // intel::hexl::EltwiseAddMod(new_AC0[0].data(), a_alpha[0].data(), b_beta[0].data(), new_AC0[0].size(), modulus);
   // intel::hexl::EltwiseAddMod(new_AC0[1].data(), a_alpha[1].data(), b_beta[1].data(), new_AC0[1].size(), modulus);

    intel::hexl::EltwiseAddMod(a_alpha[0].data(), a_alpha[0].data(), b_beta[0].data(), n, modulus);
    intel::hexl::EltwiseAddMod(a_alpha[1].data(), a_alpha[1].data(), b_beta[1].data(), n, modulus);

    //vector<vector<uint64_t>>().swap(a_alpha);
 //   vector<vector<uint64_t>>().swap(b_beta);
    
 //   gettimeofday(&tend1,NULL);
 //   time_odot += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);


    //intt(new_AC0)

  //  g_ntt->ComputeInverse(new_AC0[0].data(), new_AC0[0].data(), 1, 1);
  //  g_ntt->ComputeInverse(new_AC0[1].data(), new_AC0[1].data(), 1, 1);

    g_ntt->ComputeInverse(a_alpha[0].data(), a_alpha[0].data(), 1, 1);
    g_ntt->ComputeInverse(a_alpha[1].data(), a_alpha[1].data(), 1, 1);

  //  intel::hexl::EltwiseAddMod(AC0[0].data(), AC0[0].data(), new_AC0[0].data(), n, modulus);
  //  intel::hexl::EltwiseAddMod(AC0[1].data(), AC0[1].data(), new_AC0[1].data(), n, modulus);

    intel::hexl::EltwiseAddMod(AC0[0].data(), AC0[0].data(), a_alpha[0].data(), n, modulus);
    intel::hexl::EltwiseAddMod(AC0[1].data(), AC0[1].data(), a_alpha[1].data(), n, modulus);

//    gettimeofday(&tstart1,NULL);
    
  //  vector<vector<uint64_t>>().swap(new_AC0);
//    gettimeofday(&tend1,NULL);
//    time_update += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);

  }

  vector<polynomial> final_AC0(2, polynomial(n+1));

  int AC00_index = 0;
    int AC01_index = 0;
    int64_t halfq = q/2;
    for(int j1 = 0; j1 < n ; ++j1){
      final_AC0[0][j1+1] += (int64_t)AC0[0][j1];
      while(final_AC0[0][j1+1] >= halfq){
        final_AC0[0][j1+1] -= q;
      }
      if(final_AC0[0][j1+1] != 0){
        AC00_index = j1;
      }

      final_AC0[1][j1+1] += (int64_t)AC0[1][j1];
      while(final_AC0[1][j1+1] >= halfq){
        final_AC0[1][j1+1] -= q;
      }
      if(final_AC0[1][j1+1] != 0){
        AC01_index = j1;
      }

    }
    final_AC0[0][0] = AC00_index;
    final_AC0[1][0] = AC01_index;
   // modq_poly(AC0[0],q);
   // modq_poly(AC0[1],q);


  //return AC0
  //final_AC0.push_back(p_index);
  vector<int64_t> ans = Extract0(q, delta, n, final_AC0);
  ans.push_back((int64_t)index);
/*
  cout <<"LUT ntt(x^a-1) time: "<<time_xa<<" us."<<endl;
  cout <<"LUT ntt(x^a-1)ek time: "<<time_xaek<<" us."<<endl;
  cout <<"LUT ntt(de(AC)) time: "<<time_nttac<<" us. "<<endl;
  cout <<"LUT ntt(x^a-1)ekAC time: "<<time_odot<<" us."<<endl;
  cout <<"LUT update AC time: "<<time_update<<" us. "<<endl;
  double lut_total_time = time_init+time_init2+time_xa+time_xaek+time_nttac+time_odot+time_update;
  cout <<"LUT total time: "<<lut_total_time<<" us."<<endl;
  */
  //return final_AC0;
  return ans;

}

// LUT
// Input LWE ciphertext ct_LWE of m
// OUTPUT RLWE ciphertext of X^m * f
//template<int NN, int qNum, enum NTTSelector NTT, enum ImplSelector Impl>
extern vector<uint64_t> LUT_2048_q(int index, int q1, int64_t q, int n, int64_t delta, int k, int N, vector<int> const &ct_LWE, 
  const vector<vector<vector<uint64_t>>> &_ek0, const vector<vector<vector<uint64_t>>> &_ek1, const vector<vector<vector<uint64_t>>> &_ek2, 
  polynomial const& f, int64_t b, int logb) {
  // INPUT: id of the current node index
  // INPUT: modulus q1 and q, dimension n and N, scaler delta, k = log(q)
  // INPUT: Look-up LWE ciphertext ct_LWE
  // INPUT: RGSW Evaluation Key ek0,ek1,ek2, RLWE key s
  // INPUT: Blindrotate polynomial f
  // INPUT: decompostion base b, log(b)
  // OUTPUT: RLWE ciphertext X^(dec(ct_LWE)) * f
  struct timeval tstart1, tend1;
//  gettimeofday(&tstart1,NULL);
  //insert index
  //using R=Rq<NN, qNum, NTT, Impl>;

  uint64_t modulus = (uint64_t)q;
  //intel::hexl::NTT ntta(n, modulus);

  
 // vector<uint64_t> p_index(2, 0);
 // p_index[1] = (uint64_t)index;

  //change the ciphertext modulus of LWE ct from q to 2n
  
  vector<int64_t> a_n(N);
  
  int64_t b_n = 2 * (int64_t)n * ct_LWE[N] / q1;
  int64_t b_n2 = 2 * (int64_t)n * ct_LWE[N] / (q1/2);
  if (b_n2 - 2*b_n >= 1 && b_n2>0) {
    b_n++;
  }
  else if (b_n2 - 2*b_n <= -1 && b_n2 < 0) {
    b_n--;
  }
  
  //int64_t b_n = ct_LWE[N];


  for (int i = 0; i < N; ++i) {
    a_n[i]= 2 * (int64_t)n * ct_LWE[i] / (q1);
    int64_t tempani= 2 * (int64_t)n * ct_LWE[i] / (q1/2);
    if (tempani - 2*a_n[i] >= 1 && tempani>0) {
      a_n[i]++;
    }
    else if (tempani - 2*a_n[i] <= -1 && tempani < 0) {
      a_n[i]--;
    }
    
  }

  //initialize AC0=(f*X^b,0)
  //vector<polynomial> AC0(2);
  vector<vector<uint64_t>> AC0(2, vector<uint64_t>(n,0));
  //compute X^{b_n mod n}
  
  int64_t temp = -1;
  if (b_n >= 0) {
    temp = b_n;
  }
  else {
    while (temp < 0) {
      temp = b_n + n;
    }
  }
  polynomial xb((int)temp+2,0);
  xb[0] = temp;
  if (b_n >= 0) {
    //xb.exp = b_n;
    xb[(int)temp+1] = 1;
  }
  else if (b_n < 0) {
    //xb.exp = n + b_n;
    xb[(int)temp+1] = -1;
  }
  else;
  //compute AC[0]=0X^0, AC0[1]=f*X^b

//  gettimeofday(&tend1,NULL);
//  double time_init = (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);
//  cout <<"LUT init time: "<<time_init<<" us."<<endl;

//  gettimeofday(&tstart1,NULL);

  AC0[1]=(NTTMul_u(f, xb, n, q));

//  gettimeofday(&tend1,NULL);
//  double time_init2 = (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);
//  cout <<"LUT init2 time: "<<time_init2<<" us."<<endl;

  double time_xa = 0.0;
  double time_xaek = 0.0;
  double time_nttac = 0.0;
  double time_odot = 0.0;
  double time_update = 0.0;

  for(int i = 0; i < N/2 ; ++i){

 //   gettimeofday(&tstart1,NULL);

    //generate ntt(x^a-1)
    //R xa1,xa2,xa3;
    vector<uint64_t> xa1;
    xa1.reserve(n);
    vector<uint64_t> xa2;
    xa2.reserve(n);
    vector<uint64_t> xa3;
    xa3.reserve(n);

    XA_minus_1_2048(n,modulus,a_n[2*i]+a_n[2*i+1],xa1);
    XA_minus_1_2048(n,modulus,a_n[2*i],xa2);
    XA_minus_1_2048(n,modulus,a_n[2*i+1],xa3);

    //gettimeofday(&tstart1,NULL);

    g_ntt->ComputeForward(xa1.data(), xa1.data(), 1, 1);
    g_ntt->ComputeForward(xa2.data(), xa2.data(), 1, 1);
    g_ntt->ComputeForward(xa3.data(), xa3.data(), 1, 1); 

  //  gettimeofday(&tend1,NULL);
  //  time_xa += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);

     //compute ntt(X^a-1)*ntt(ek), store in Xa_times_ek

  //  gettimeofday(&tstart1,NULL);

    vector<vector<vector<uint64_t>>> Xa_times_ek(2*(k/logb),vector<vector<uint64_t>>(2, vector<uint64_t>(n)));

    vector<uint64_t> tempxaek1;
        tempxaek1.reserve(n);
        vector<uint64_t> tempxaek2;
        tempxaek2.reserve(n);

    for(int j1 = 0; j1 < 2*(k/logb) ; ++j1){
      for(int j2 = 0; j2 < 2 ; ++j2){


        intel::hexl::EltwiseMultMod(tempxaek1.data(), xa1.data(), _ek0[i][j1*2+j2].data(), n, modulus, 1);
        intel::hexl::EltwiseMultMod(tempxaek2.data(), xa2.data(), _ek1[i][j1*2+j2].data(), n, modulus, 1);

        intel::hexl::EltwiseAddMod(tempxaek1.data(), tempxaek1.data(), tempxaek2.data(), n, modulus);
        intel::hexl::EltwiseMultMod(tempxaek2.data(), xa3.data(), _ek2[i][j1*2+j2].data(), n, modulus, 1);
        intel::hexl::EltwiseAddMod(Xa_times_ek[j1][j2].data(), tempxaek1.data(), tempxaek2.data(), Xa_times_ek[j1][j2].size(), modulus);

      }
    }

    vector<uint64_t>().swap(xa1);
    vector<uint64_t>().swap(xa2);
    vector<uint64_t>().swap(xa3);
    vector<uint64_t>().swap(tempxaek1);
    vector<uint64_t>().swap(tempxaek2);

  //  gettimeofday(&tend1,NULL);
  //  time_xaek += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);

    //decompose AC0 to 2*k/logb

  //  gettimeofday(&tstart1,NULL);
    vector<vector<uint64_t>> de_AC00 = bit_poly_u(k,AC0[0],n,q,b,logb);
    vector<vector<uint64_t>> de_AC01 = bit_poly_u(k,AC0[1],n,q,b,logb);

 //   gettimeofday(&tend1,NULL);
 //   time_nttac += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);

    vector<vector<uint64_t>> ntt_de_AC00(k/logb, vector<uint64_t>(n));
    vector<vector<uint64_t>> ntt_de_AC01(k/logb, vector<uint64_t>(n));


    for(int j1 = 0; j1 < k/logb ; ++j1){
      
      g_ntt->ComputeForward(de_AC00[j1].data(), de_AC00[j1].data(), 1, 1);
      g_ntt->ComputeForward(de_AC01[j1].data(), de_AC01[j1].data(), 1, 1);
      ntt_de_AC01[j1] = move(de_AC01[j1]);
      ntt_de_AC00[j1] = move(de_AC00[j1]);
    }

   // vector<vector<uint64_t>>().swap(de_AC00);
  //  vector<vector<uint64_t>>().swap(de_AC01);
    

   // gettimeofday(&tend1,NULL);
   // time_nttac += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);


  //  gettimeofday(&tstart1,NULL);
    //b <> beta = ntt_de_AC01 <> Xa_times_ek[first k/logb]
    //vector<R> b_beta;
    vector<vector<uint64_t>> b_beta(2,vector<uint64_t>(n));

    vector<uint64_t> tempac01ek0;
    tempac01ek0.reserve(n);
    vector<uint64_t> tempac01ek1;
    tempac01ek1.reserve(n);

    intel::hexl::EltwiseMultMod(tempac01ek0.data(), ntt_de_AC01[0].data(), Xa_times_ek[0][0].data(), n, modulus, 1);
    intel::hexl::EltwiseMultMod(tempac01ek1.data(), ntt_de_AC01[0].data(), Xa_times_ek[0][1].data(), n, modulus, 1);

    for(int j1 = 1 ; j1 < k/logb ; ++j1){
      vector<uint64_t> ac01ek0;
      ac01ek0.reserve(n);
    //  vector<uint64_t> sum0;
    //  sum0.reserve(n);

      intel::hexl::EltwiseMultMod(ac01ek0.data(), ntt_de_AC01[j1].data(), Xa_times_ek[j1][0].data(), n, modulus, 1);
      intel::hexl::EltwiseAddMod(tempac01ek0.data(), tempac01ek0.data(), ac01ek0.data(), n, modulus);

    //  tempac01ek0 = move(sum0);

    //  vector<uint64_t> ac01ek1;
    //  ac01ek1.reserve(n);
   //   vector<uint64_t> sum1;
    //  sum1.reserve(n);

      intel::hexl::EltwiseMultMod(ac01ek0.data(), ntt_de_AC01[j1].data(), Xa_times_ek[j1][1].data(), n, modulus, 1);
      intel::hexl::EltwiseAddMod(tempac01ek1.data(), tempac01ek1.data(), ac01ek0.data(), n, modulus);

     // tempac01ek1 = move(sum1);

    }

    //b_beta[0]=tempac01ek0;
    b_beta[0] = move(tempac01ek0);
   // b_beta[1]=tempac01ek1;
    b_beta[1] = move(tempac01ek1);

    vector<vector<uint64_t>>().swap(ntt_de_AC01);


    vector<vector<uint64_t>> a_alpha(2,vector<uint64_t>(n));

    vector<uint64_t> tempac00ek0;
    tempac00ek0.reserve(n);
    vector<uint64_t> tempac00ek1;
    tempac00ek1.reserve(n);

    intel::hexl::EltwiseMultMod(tempac00ek0.data(), ntt_de_AC00[0].data(), Xa_times_ek[k/logb][0].data(), n, modulus, 1);
    intel::hexl::EltwiseMultMod(tempac00ek1.data(), ntt_de_AC00[0].data(), Xa_times_ek[k/logb][1].data(), n, modulus, 1);

    for(int j1 = 1 ; j1 < k/logb ; ++j1){

      vector<uint64_t> ac00ek0;
      ac00ek0.reserve(n);
      //vector<uint64_t> sum0;
     // sum0.reserve(n);

      intel::hexl::EltwiseMultMod(ac00ek0.data(), ntt_de_AC00[j1].data(), Xa_times_ek[j1+k/logb][0].data(), n, modulus, 1);
      intel::hexl::EltwiseAddMod(tempac00ek0.data(), tempac00ek0.data(), ac00ek0.data(), n, modulus);

      //tempac00ek0 = sum0;
     // tempac00ek0 = move(sum0);

     // vector<uint64_t> ac00ek1;
    //  ac00ek1.reserve(n);
   //   vector<uint64_t> sum1;
    //  sum1.reserve(n);

      intel::hexl::EltwiseMultMod(ac00ek0.data(), ntt_de_AC00[j1].data(), Xa_times_ek[j1+k/logb][1].data(), n, modulus, 1);
      intel::hexl::EltwiseAddMod(tempac00ek1.data(), tempac00ek1.data(), ac00ek0.data(), n, modulus);

     // tempac00ek1 = sum1;
     // tempac00ek1 = move(sum1);


    }

    //a_alpha[0]=tempac00ek0;
    a_alpha[0] = move(tempac00ek0);
   // a_alpha[1]=tempac00ek1;
    a_alpha[1] = move(tempac00ek1);

    vector<vector<uint64_t>>().swap(ntt_de_AC00);
    vector<vector<vector<uint64_t>>>().swap(Xa_times_ek);

    //a_alpha + b_beta

   // vector<vector<uint64_t>> new_AC0(2,vector<uint64_t>(n));
   // intel::hexl::EltwiseAddMod(new_AC0[0].data(), a_alpha[0].data(), b_beta[0].data(), new_AC0[0].size(), modulus);
   // intel::hexl::EltwiseAddMod(new_AC0[1].data(), a_alpha[1].data(), b_beta[1].data(), new_AC0[1].size(), modulus);

    intel::hexl::EltwiseAddMod(a_alpha[0].data(), a_alpha[0].data(), b_beta[0].data(), n, modulus);
    intel::hexl::EltwiseAddMod(a_alpha[1].data(), a_alpha[1].data(), b_beta[1].data(), n, modulus);

    //vector<vector<uint64_t>>().swap(a_alpha);
 //   vector<vector<uint64_t>>().swap(b_beta);
    
 //   gettimeofday(&tend1,NULL);
 //   time_odot += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);


    //intt(new_AC0)

  //  g_ntt->ComputeInverse(new_AC0[0].data(), new_AC0[0].data(), 1, 1);
  //  g_ntt->ComputeInverse(new_AC0[1].data(), new_AC0[1].data(), 1, 1);

    g_ntt->ComputeInverse(a_alpha[0].data(), a_alpha[0].data(), 1, 1);
    g_ntt->ComputeInverse(a_alpha[1].data(), a_alpha[1].data(), 1, 1);

  //  intel::hexl::EltwiseAddMod(AC0[0].data(), AC0[0].data(), new_AC0[0].data(), n, modulus);
  //  intel::hexl::EltwiseAddMod(AC0[1].data(), AC0[1].data(), new_AC0[1].data(), n, modulus);

    intel::hexl::EltwiseAddMod(AC0[0].data(), AC0[0].data(), a_alpha[0].data(), n, modulus);
    intel::hexl::EltwiseAddMod(AC0[1].data(), AC0[1].data(), a_alpha[1].data(), n, modulus);

//    gettimeofday(&tstart1,NULL);
    
  //  vector<vector<uint64_t>>().swap(new_AC0);
//    gettimeofday(&tend1,NULL);
//    time_update += (tend1.tv_sec-tstart1.tv_sec)*1000000+(tend1.tv_usec-tstart1.tv_usec);

  }

  vector<vector<uint64_t>> final_AC0(2, vector<uint64_t>(n+1));

  int AC00_index = 0;
  int AC01_index = 0;
  for(int j1 = 0; j1 < n ; ++j1){
      final_AC0[0][j1+1] += AC0[0][j1];
      if(final_AC0[0][j1+1] != 0){
        AC00_index = j1;
      }

      final_AC0[1][j1+1] += AC0[1][j1];
      if(final_AC0[1][j1+1] != 0){
        AC01_index = j1;
      }

    }
    final_AC0[0][0] = AC00_index;
    final_AC0[1][0] = AC01_index;
   // modq_poly(AC0[0],q);
   // modq_poly(AC0[1],q);


  //return AC0
  //final_AC0.push_back(p_index);
  vector<uint64_t> ans = Extract0_q(q, delta, n, final_AC0);
  ans.push_back((uint64_t)index);
/*
  cout <<"LUT ntt(x^a-1) time: "<<time_xa<<" us."<<endl;
  cout <<"LUT ntt(x^a-1)ek time: "<<time_xaek<<" us."<<endl;
  cout <<"LUT ntt(de(AC)) time: "<<time_nttac<<" us. "<<endl;
  cout <<"LUT ntt(x^a-1)ekAC time: "<<time_odot<<" us."<<endl;
  cout <<"LUT update AC time: "<<time_update<<" us. "<<endl;
  double lut_total_time = time_init+time_init2+time_xa+time_xaek+time_nttac+time_odot+time_update;
  cout <<"LUT total time: "<<lut_total_time<<" us."<<endl;
  */
  //return final_AC0;
  return ans;

}


