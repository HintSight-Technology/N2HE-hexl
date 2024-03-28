
#include <iostream>
#include "include.hpp"

using namespace std;

// Switches
const bool ENABLE_MULTITHREADING = true;

intel::hexl::NTT* g_ntt = nullptr;


// ================================================

// parameters for Neural Network
                     
const int nh2 = 30;                     // #Nodes in Hidden Layer
const int nh = 30;
const int classes = 30;

int scaler0 = 16;

int scaler1 = 16;                       // Scaler for Weights between Input Layer and Hidden Layer, and Bias in Hidden Layer

int scaler2 = 4;                        // Scaler for Weights between Input Layer and Hidden Layer, and Bias in Hidden Layer

int ww1[512][nh2];                        // (int) Weights between Input Layer and Hidden Layer

int64_t ww2[nh2][classes];                        // (int) Weights between Output Layer and Hidden Layer

int ibias1[nh2];                          // (int) Bias in Hidden Layer

int64_t ibias2[classes];                          // (int) Bias in Output Layer

int64_t hat = -3000;                      // threshold for unknown people (max < hat)
int64_t threshold = 4000;               // threshold for unknown people (max - 2nd max)


// Test set 

vector< vector<int> > test_vectors;

vector<int> test_labels;

//parameter for NTT

const int64_t mod = 576460752213245953;
const int64_t delta = 1099511627776;       //delta * ip1 < q/4 && delta*ip2 < q/2

// ================================================

// Define Type polynomial: 
// polynomial p is an array of size (deg(p) + 2): [deg(p), coef_0(p), coef_1(p), coef_2(p),..., coef_{deg(p)}(p)]
typedef vector<int64_t> polynomial;

vector<vector<vector<uint64_t>>> _ek0, _ek1, _ek2;
int * lwe_s;
// ================================================
struct Params_ip1 {
  int index;
  int q1;
  int NN;
  vector<vector<int>> const &ct_LWE; 
};


struct Params_lut {
  int n;
  int k;
  int NN;
  int q1;
  int64_t b;
  int logb;
  int64_t delta;
  int64_t q2;
  vector<int> const &ct_LWE;
  polynomial const &f;
};

struct Params_ip2 {
  uint64_t index;
  int64_t q2;
  int64_t delta;
  int n;
  vector<vector<uint64_t>> const & ct_lut; 
};

// Define to sort
bool cmp(const vector<int> &a, const vector<int> &b) {
   return a[513]< b[513];
 }

bool cmp3(const vector<uint64_t> &a, const vector<uint64_t> &b) {
   return a[2049]< b[2049];
 }


// ======================================================  Neural Network  ==========================================================================

//read test vectors and labels
void read_vector() {
  ifstream infile;
  infile.open("../FacicalRecognition_data/testImgAll.txt");
  if(!infile.is_open()){
    cout <<"Cannot open file testImgAll.txt"<<endl;
  }
   int num = 0;
   infile >> num;
  //int num = 1;
  for (int i = 0; i < num; ++i) {
    vector<int> ivector;
    for (int j = 0; j < 512; ++j) {
      double t;
      infile >> t;
      t = t * scaler0;
      int tt=(int)t;
      if (t > 0 && t - tt > 0.5) {
        ivector.push_back( tt + 1);
      }
      else if (t < 0 && tt - t>0.5) {
        ivector.push_back( tt - 1);
      }
      else {
        ivector.push_back( tt);
      }
    }
    test_vectors.push_back(ivector);
  }
  infile.close();
}


// Import trained Neural Network
// Turn it into int Neural Network

void read_w() {

  ifstream infile;
  infile.open("../FacicalRecognition_data/cnn30.txt");
  if(!infile.is_open()){
    cout <<"Cannot open file cnn30.txt"<<endl;
  }
  double round_w = 1.0;
  double tempw;
  for (int i = 0; i < nh2; i++)
  {
    for (int j = 0; j < 512; j++)
    {
      infile >> tempw;
      double t = tempw * scaler1;
      ww1[j][i] = ((int)(tempw * scaler1));
      if (t - ww1[j][i] > 0.5 && t > 0) {
        ww1[j][i]++;
      }
      else if (t < 0 && ww1[j][i] - t>0.5) {
        ww1[j][i]--;
      }
    }
  }
  for (int j = 0; j < classes; j++)
  {
    for (int i = 0; i < nh2; i++)
    {
      infile >> tempw;
      double t = tempw * scaler2;
      int64_t tt = (int64_t)(tempw * scaler2);
      if (t > 0 && t - tt > 0.5) {
        ww2[i][j] = tt + 1;
      }
      else if (t < 0 && tt - t>0.5) {
        ww2[i][j] = tt - 1;
      }
      else {
        ww2[i][j] = tt;
      }
    }
  }
  for (int i = 0; i < nh2; i++)
  {
    infile >> tempw;
    double t = tempw * scaler1;
    int tt = (int)(tempw * scaler1);
    if (t > 0 && t - tt > 0.5) {
      ibias1[i] = tt + 1;
    }
    else if (t < 0 && tt - t>0.5) {
      ibias1[i] = tt - 1;
    }
    else {
      ibias1[i] = tt;
    }
  }
  for (int i = 0; i < classes; i++)
  {
    infile >> tempw;
    double t = tempw * scaler2;
    int64_t tt = (int64_t)(tempw * scaler2);
    if (t > 0 && t - tt > 0.5) {
      ibias2[i] = tt + 1;
    }
    else if (t < 0 && tt - t>0.5) {
      ibias2[i] = tt - 1;
    }
    else {
      ibias2[i] = tt;
    }
  }
  infile.close();
}

// label index
string print_name(int index) {
  return "name"+to_string(index);
}


// END OF Import trained Neural Network

vector<int> IP1(int index, int q1, int N,  vector<vector<int>> const &ct_LWE){
  if(index % 4 == 0){
    vector<int> ct_ip1 = LWE32_Plain_Multi_ct(q1, N, ct_LWE[0], ww1[0][index/4]);
    for (register int i = 1; i < 128; ++i) {
      vector<int> ct_iip0 = LWE32_Plain_Multi_ct(q1, N, ct_LWE[i], ww1[i][index/4]);
      ct_ip1 = LWE32_Add_ct(q1, N, ct_iip0, ct_ip1);
    }
    //ct_ip1 = LWE32_Plain_Add_ct(q1, N, ct_ip1, 2*ibias1[index]);
    ct_ip1.emplace_back(index);
    return ct_ip1;
  }
  if(index % 4 == 1){
    vector<int> ct_ip1 = LWE32_Plain_Multi_ct(q1, N, ct_LWE[0+128], ww1[0+128][index/4]);
    for (register int i = 1; i < 128; ++i) {
      vector<int> ct_iip0 = LWE32_Plain_Multi_ct(q1, N, ct_LWE[i+128], ww1[i+128][index/4]);
      ct_ip1 = LWE32_Add_ct(q1, N, ct_iip0, ct_ip1);
    }
    //ct_ip1 = LWE32_Plain_Add_ct(q1, N, ct_ip1, 2*ibias1[index]);
    ct_ip1.emplace_back(index);
    return ct_ip1;
  }
  if(index % 4 == 2){
    vector<int> ct_ip1 = LWE32_Plain_Multi_ct(q1, N, ct_LWE[0+256], ww1[0+256][index/4]);
    for (register int i = 1; i < 128; ++i) {
      vector<int> ct_iip0 = LWE32_Plain_Multi_ct(q1, N, ct_LWE[i+256], ww1[i+256][index/4]);
      ct_ip1 = LWE32_Add_ct(q1, N, ct_iip0, ct_ip1);
    }
    //ct_ip1 = LWE32_Plain_Add_ct(q1, N, ct_ip1, 2*ibias1[index]);
    ct_ip1.emplace_back(index);
    return ct_ip1;
  }
  else{
    vector<int> ct_ip1 = LWE32_Plain_Multi_ct(q1, N, ct_LWE[0+384], ww1[0+384][index/4]);
    for (register int i = 1; i < 128; ++i) {
      vector<int> ct_iip0 = LWE32_Plain_Multi_ct(q1, N, ct_LWE[i+384], ww1[i+384][index/4]);
      ct_ip1 = LWE32_Add_ct(q1, N, ct_iip0, ct_ip1);
    }
    ct_ip1 = LWE32_Plain_Add_ct(q1, N, ct_ip1, 2*scaler0*ibias1[index/4]);
    ct_ip1.emplace_back(index);
    return ct_ip1;
  }
}

vector<uint64_t> LUT_q(int q1, int64_t q2, int n, int64_t delta, int k, int N, vector<int> const &ct_ip1, 
    polynomial const& f, int64_t b, int logb) {
    int index = ct_ip1[N+1];
    vector<int> ct_ip1_t (ct_ip1.begin(),ct_ip1.end()-1);
    return LUT_2048_q(index,q1,q2, n, delta, k, N,ct_ip1_t, _ek0, _ek1,_ek2, f,b,logb);
}

vector<uint64_t> IP2_q(uint64_t index, int64_t q2, int64_t delta, int n2, vector<vector<uint64_t>> const & ct_lut){
   vector<uint64_t> u_temp_ip2 = LWE64_Plain_Multi_ct_q(q2, n2, ct_lut[0], ww2[0][index]);
      for (int j = 1; j < nh; ++j) {
        vector<uint64_t> u_ct_iipj = LWE64_Plain_Multi_ct_q(q2, n2, ct_lut[j], ww2[j][index]);
        u_temp_ip2 = LWE64_Add_ct_q(q2, n2, u_temp_ip2, u_ct_iipj);
      }
      u_temp_ip2 = LWE64_Plain_Add_ct_q(q2, n2, u_temp_ip2, 2*ibias2[index]*delta*scaler0 );
      u_temp_ip2.emplace_back(index);
      return u_temp_ip2;
}



int main(){
  
  cout <<"Task: Privacy-preserving facial recognition. "<<endl;
  int start_index = 0;
  int test_count = 37;
  
  // ======================== Initialization and Information Output ==========================
  read_vector();
  read_w();

  //set n,q,t for LWE scheme
  cout << "Parameters for LWE scheme:" << endl;
  const int n1 = 512;                        //vector length
  const int q1 = 1<<13;                      //ciphertext modulus
  const int var1 = 3;                          //variance of error
  
  cout << "vector length = " << n1 <<", ciphertext modulus = " << q1 <<", variance of error = 2^ -"<<var1 << endl;
  cout << "----------" << endl;

  
  //set n,q,k for RLWE scheme
  cout << "Parameters for RLWE scheme:" << endl;
  const int n2 = 1<<11; //vector length=exp of polynomial=2048
  const int n3 = n2;
  const int var2=3;                      //variance of error
  const int k2 = 60; //k=log(q)
  const int64_t b = 1<<30;
  const int logb = 30; //logb=log2(b)-1
  const int64_t q2 = mod; //ciphertext modulus
  
  cout << "polynomial degree = " << n2 << ", ciphertext modulus: " << q2 <<", variance of error = 2^ -"<<var2  <<", external product based on " << logb << " bits extension" << endl;
  cout << "----------" << endl;
   struct timeval tstart, tend;

  g_ntt = new intel::hexl::NTT(2048,576460752213245953);
  _ek0 = vector<vector<vector<uint64_t>>>(n1/2, vector<vector<uint64_t>>(4*(k2/logb), vector<uint64_t>(n2)));
  _ek1 = vector<vector<vector<uint64_t>>>(n1/2, vector<vector<uint64_t>>(4*(k2/logb), vector<uint64_t>(n2)));
  _ek2 = vector<vector<vector<uint64_t>>>(n1/2, vector<vector<uint64_t>>(4*(k2/logb), vector<uint64_t>(n2)));

  int *x = new int [n1];
  x=LWE32_KeyGen(n1);
  //generate LWE sk and read test images



  // ======================== Evaluation Key Generation ======================================

  //encrypt the RGSW key
  polynomial s = RLWE64_KeyGen(n2);
  int *lwe_s = new int [n2];
  for (int i = (int)s[0] + 1; i >= 1; --i) {
    lwe_s[(int)s[0] + 1 - i]=s[i];
  }
 

  uint64_t modulus = (uint64_t)q2;
    //ntt ek0,ek1,ek2
   gettimeofday(&tstart,NULL);
  
  //intel::hexl::NTT ntt(n2, modulus);


  for (int i = 0; i < n1/2; ++i) {
    _ek0[i] = key_encrypt_1_ntt(n2, q2, k2,var2, x[2*i],x[2*i+1], s, b, logb);
  }

  for(int i = 0 ; i < n1/2 ; ++i){
    _ek1[i] = key_encrypt_2_ntt(n2, q2, k2,var2, x[2*i],x[2*i+1], s, b, logb);
  }

  for(int i = 0 ; i < n1/2 ; ++i){
    _ek2[i] = key_encrypt_3_ntt(n2, q2, k2,var2, x[2*i],x[2*i+1], s, b, logb);
  }

   gettimeofday(&tend,NULL);
  double  time_ek = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;
  cout <<"Time for evaluation key generation: "<<time_ek<<"s. "<<endl;

  double ave_time = 0;
  double ave_ip1_time = 0;
  double ave_lut_time = 0;
  double ave_ip2_time = 0;
 // struct timeval tstart, tend;
  int correct_num=0;
  int correct_num_test = 0;

  // =========================== Start Testing ===============================================
  cout <<"Test start from photo No. "<<start_index <<", to No. "<< start_index+test_count<<endl;
  for (int tt = start_index; tt < start_index + test_count; ++tt) {

    cout <<"Test no: "<<tt <<": ";
    //encrypt image i, scale from {0,1} to {0,2}

    vector<vector<int>> ct_i;
    for (int j = 0; j < 512; ++j) {
      vector<int>ctj = LWE32_Enc(q1, n1,var1, 2*test_vectors[tt][j], x);
      ct_i.push_back(ctj);
    }
    int true_label = 0;

    //compute the first inner-product
    //lut, construct polynomial f
    polynomial f = p_f(n2, delta, q1);

    
    vector<vector<uint64_t>>u_ct_lut(nh, vector<uint64_t>(n2+1));
    //vector<vector<uint64_t>> u_ct_ip2(classes, vector<uint64_t>(n2+1));
    vector<vector<uint64_t>> ip2_results(classes);
    double time_ip1 = 0;
    double time_lut = 0;
    double time_ip2 = 0;


    if (ENABLE_MULTITHREADING == 1){
      //multithreading

      //construct inputs
      gettimeofday(&tstart,NULL);

      vector<vector<int>> ip1_results_t(4*nh);
      vector<vector<int>> ip1_results(nh);
      vector<vector<uint64_t>> lut_results(nh);
      

      vector<Params_ip1> inputs_ip1;
      for(int i = 0 ; i < 4*nh ; ++i){
        Params_ip1 a = {i, q1, n1, ct_i};
        inputs_ip1.emplace_back(a);
      }

       vector<thread> w;

      for (int ii = 0 ; ii < 4*nh ; ++ii){
        w.emplace_back([&,ii]() {

          auto r = IP1(inputs_ip1[ii].index, inputs_ip1[ii].q1, inputs_ip1[ii].NN, inputs_ip1[ii].ct_LWE);

          ip1_results_t[ii] = move(r);

          });
      }

      for (auto& th : w) {
        th.join();
      }

      sort(ip1_results_t.begin(), ip1_results_t.end(), cmp);

      for(int i = 0 ; i < nh ; ++i){
        ip1_results[i] = LWE32_Add_ct(q1, n1+1, ip1_results_t[4*i], ip1_results_t[4*i+1]);
        ip1_results[i] = LWE32_Add_ct(q1, n1+1, ip1_results[i], ip1_results_t[4*i+2]);
        ip1_results[i] = LWE32_Add_ct(q1, n1+1, ip1_results[i], ip1_results_t[4*i+3]);
/*
        ip1_results[i].pop_back();
        cout<<LWE32_Dec(q1,n1,ip1_results[i],x)<<" ";
        ip1_results[i].push_back(i);
*/
        ip1_results[i][n1+1]=i;
      }
 //     cout <<"---------"<<endl;

       gettimeofday(&tend,NULL);
      time_ip1 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;

      gettimeofday(&tstart,NULL);

      vector<Params_lut> inputs_lut;
      for (int i = 0; i < nh; ++i) {
        Params_lut a = { n2,  k2, n1, q1, b,logb, delta, q2, ip1_results[i], f };
        inputs_lut.emplace_back(a);
      }

      vector<thread> w2;

      for (int ii = 0 ; ii < nh ; ++ii){
        w2.emplace_back([&,ii]() {

          auto r = LUT_q(inputs_lut[ii].q1, inputs_lut[ii].q2, inputs_lut[ii].n, inputs_lut[ii].delta, inputs_lut[ii].k, inputs_lut[ii].NN, 
            inputs_lut[ii].ct_LWE, inputs_lut[ii].f, inputs_lut[ii].b, inputs_lut[ii].logb);

          lut_results[ii] = move(r);

          });
      }

      for (auto& th : w2) {
        th.join();
      }

      sort(lut_results.begin(), lut_results.end(), cmp3);

      for (int i = 0; i < nh; ++i) {
        lut_results[i].pop_back();
        u_ct_lut[i] = lut_results[i];
      //  cout<<LWE64_Dec_q(q2,n2,u_ct_lut[i],lwe_s)/delta<<endl;
      }

      gettimeofday(&tend,NULL);
      time_lut = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;      // time cost of LUT (multithreading)

      gettimeofday(&tstart,NULL);

      vector<Params_ip2> inputs_ip2;
      for (int i = 0; i < classes; ++i) {
        Params_ip2 a = { (uint64_t)i, q2, delta, n2, u_ct_lut};
        inputs_ip2.emplace_back(a);
      }

      vector<thread> w3;

      for (int ii = 0 ; ii < classes ; ++ii){
        w3.emplace_back([&,ii]() {

          auto r = IP2_q(inputs_ip2[ii].index, inputs_ip2[ii].q2, inputs_ip2[ii].delta, inputs_ip2[ii].n, inputs_ip2[ii].ct_lut);

          ip2_results[ii] = move(r);

          });
      }

      for (auto& th : w3) {
        th.join();
      }

      gettimeofday(&tend,NULL);
      time_ip2 = tend.tv_sec-tstart.tv_sec+(tend.tv_usec-tstart.tv_usec)/1000000.0;

      
    }

    

    //decrypt
    int64_t max = -1*mod;
    int64_t max2 = max;
    int cindex = 0;
    int cindex2 = 0;
    double likelihood[classes];

    //check vector
    // ofstream disc_out;

    // disc_out.open("result.txt",ios::app);
    for (int i = 0; i < classes; ++i) {
      uint64_t index = ip2_results[i][2049];
      ip2_results[i].pop_back();
      int64_t u_tempr = LWE64_Dec_q(q2, n2, ip2_results[i], lwe_s);
      likelihood[i] = 1.0 * u_tempr / delta;
      if (likelihood[i] > max) {
        max = (int64_t) likelihood[i];
        cindex = (int)index;
      }
    }
    
    for (int i = 0; i < classes; ++i) {
      if (likelihood[i] > max2 && i != cindex) {
        max2 = (int64_t) likelihood[i];
        cindex2 = i;
      }
    }

  //  cout <<"max: "<<max<<" , max2: "<<max2<<endl;

    double total_time = time_ip1+time_lut+time_ip2;
    ave_time += total_time;
    ave_ip1_time += time_ip1;
    ave_lut_time += time_lut;
    ave_ip2_time += time_ip2;


    // output 
     if (max < hat && max-max2<threshold) {
      cout <<  "Prediction: This people is not in dataset, time: " << total_time << "s. " << endl;
     } 
     else {
      cout <<  "Prediction: " << print_name(cindex) << ", time: " << total_time << "s. " << endl;
     }

  }
  
  // =========================== End of Testing ==============================================
  cout << "----------" << endl;
  cout <<"Average time: "<<ave_time/test_count<<"s. "<<endl;
  cout <<"Average time for ip1: "<<ave_ip1_time/test_count<<"s. "<<endl;
  cout <<"Average time for lut: "<<ave_lut_time/test_count<<"s. "<<endl;
  cout <<"Average time for ip2: "<<ave_ip2_time/test_count<<"s. "<<endl;
  return 0;

}



