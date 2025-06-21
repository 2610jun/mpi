#include <iostream>
#include "polar.h"
// #include "pac.h"
#include <random>
#include <cstdarg>
#include <iomanip>

using namespace std;

// function declaration

// console outputput package
string format(const char* fmt, ...);

//  for EbN0 (x axis)
vector<double> linspace(double start, double ed, int num);

// decoder 
void run_SC_parallel(int world_rank, int world_size, int N, int K, int EBN0_BEGIN, int EBN0_END, double EBN0_INTERVAL, int N_MAX_BLK, int TARGET_FRAME_ERR, string RM_METHOD);
void run_CA_SCL_parallel(int world_rank, int world_size, int N, int K, int NUM_LIST, string CRC_TYPE, int EBN0_BEGIN, int EBN0_END, double EBN0_INTERVAL, int N_MAX_BLK, int TARGET_FRAME_ERR, string RM_METHOD);
//void run_CONV();
//void run_PAC();

int main(int argc, char** argv){
    //////////////////////////// MPI INITIALIZE ////////////////////////////////
   
    // MPI_initialize
    MPI_Init(&argc, &argv);

    // Get the number of processes
    int world_size; // number of process?
    MPI_Comm_size(MPI_COMM_WORLD, &world_size); 
    
    // Get the rank of the process
    int world_rank;   
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    

    // Running parameter
    string ENCODER_TYPE = argv[1];
    string DECODER_TYPE = argv[2];
    
    int N = stoi(argv[3]);
    int K = stoi(argv[4]);
    int NUM_LIST = stoi(argv[5]);
    string CRC_TYPE = argv[6];

    int EBN0_BEGIN = stoi(argv[7]);
    int EBN0_END = stoi(argv[8]);
    double EBN0_INTERVAL = stod(argv[9]);

    int N_MAX_BLK = stoi(argv[10]);
    int TARGET_FRAME_ERR = stoi(argv[11]);

    string RM_METHOD = argv[12];
    //////////////////////////// MPI INITIALIZE ////////////////////////////////
   
    // Test code to setup the MPI process.
    
     if(ENCODER_TYPE == "POLAR" && DECODER_TYPE == "SC"){
               if(world_rank == 0) cout<<"MPI Polar SC simulation started with " << world_size << " processes."<<endl; 
               run_SC_parallel(world_rank, world_size, N, K, EBN0_BEGIN, EBN0_END, EBN0_INTERVAL, N_MAX_BLK, TARGET_FRAME_ERR, RM_METHOD); // CA_SCL에 MPI code 나중
     }
     else if (ENCODER_TYPE == "POLAR" && DECODER_TYPE == "CA_SCL"){
               
               if(world_rank == 0) cout<<"MPI Polar CA-SCL simulation started with " << world_size << " processes."<<endl; 
               run_CA_SCL_parallel(world_rank, world_size, N, K, NUM_LIST, CRC_TYPE, EBN0_BEGIN, EBN0_END, EBN0_INTERVAL, N_MAX_BLK, TARGET_FRAME_ERR, RM_METHOD); 
     }

     MPI_Finalize();
}


void run_SC_parallel(int world_rank, int world_size, int N, int K, int EBN0_BEGIN, int EBN0_END, double EBN0_INTERVAL, int N_MAX_BLK, int TARGET_FRAME_ERR, string RM_METHOD){
    
    double R = (double)K/N;

    // instantiate a POLAR object
    POLAR polar = POLAR(N,K, 0, RM_METHOD);
     
    default_random_engine random_engine(world_rank + time(0)); // random vectror를 각 rank별로 차이를 줘야함
    bernoulli_distribution  bern_dist;

    //cout<<"random ok?"<<endl;

    // Running parameters
    /**
    * Es: 하나의 심볼을 전송하는데 필요한 에너지
    * Eb: 하나의 비트를 전송하는데 필요한 에너지
    * N0: Noise power density
    * Es=k*Eb
    * k=심볼당 비트수. (BPSK->1, QPSK->2, 16-QAM->4)
    * BER을 비교할때는 대부분 EbN0 사용
    * Modulation symbol과 관련된 경우에는 EsN0를 사용한다
    * 
    * 이번 데모에서는 EbN0 사용하면 될거같음. 그리고, 어차피 BPSK 사용할거라 두개 같은 값임
    * Polar code의 경우 BPSK가 보편적인 modulation 방법이다.
    */


    // parameters about simulation 
    double linspace_num_point = (EBN0_END-EBN0_BEGIN)/EBN0_INTERVAL + 1;
    //cout<<"EBN0 Interval: "<<EBN0_INTERVAL<<endl;   
    //cout<<"num_point: "<<linspace_num_point<<endl;
    vector<double> EbN0_dB = linspace(EBN0_BEGIN, EBN0_END, linspace_num_point); // 7 point in 0dB~3dB
    string decoder_name = "SC";

    vector<double> EbN0;
    vector<double> sigma;
    //vector<double> N0;
    for (auto e : EbN0_dB){
        double x = pow(10.0, e/10);
        EbN0.push_back(x);
        //cout<<x<<endl;
        sigma.push_back(sqrt(1/(2*R*x)));
        //cout<<sqrt(1/(2*R*x))<<endl;
    }

    vector<double> ber(EbN0.size(), 0), bler(EbN0.size(), 0);
    vector<int> n_bit_errs_local(EbN0.size(), 0), n_blk_errs_local(EbN0.size(), 0);
    vector<int> n_bit_errs_total(EbN0.size(), 0), n_blk_errs_total(EbN0.size(), 0);


    // 헤더 출력
    if(world_rank == 0){
    cout << "============================================================\n";
    cout << "Eb/N0             : " 
                                    << EBN0_BEGIN << "[dB] ~ " 
                                    << EBN0_END << "[dB], " 
                                    << EBN0_INTERVAL << "[dB]\n";

    cout << "Target Frame Error: " << TARGET_FRAME_ERR<<"\n";
    cout << "Max num frame     : " << N_MAX_BLK<<"\n";
    cout << "N                 : " << N << "\n";
    cout << "K                 : " << K << "\n";
    cout << "Decoder           : " << decoder_name << "\n";
    cout << "============================================================\n\n";

    cout << "|| Eb/N0 ||     BER     ||     BLER     || n_bit_errs || n_blk_errs || n_blk_done || elapsed time ||" << endl;
    cout << "====================================================================================================" << endl;
    }
    // Loop for each EbN0 (SNR)
    for(size_t i = 0; i<EbN0.size(); i++){
        //cout<<str<<endl;
        string str;// 모르겠다
        int print_len = 0; // 추후에 출력할때 지우고 다시 쓰는 동작을 구현하기 위해 필요하다

        int n_blks_done_local = 0;
        int n_blks_done_total = 0;
        clock_t ts = clock();
        while((n_blks_done_total<N_MAX_BLK)&&(n_blk_errs_total[i]<TARGET_FRAME_ERR)){
            // generate random bit stream
            //cout << "[Rank " << world_rank << "] 현재 출력 메시지입니다." << std::endl;
            vector<bool> msg;

            for (int j = 0; j < K; j++)
                msg.push_back(bern_dist(random_engine));
            
            // polar encoding
            vector<bool> enc = polar.encoder(&msg);
            // rate matching (Puncturing)
            vector<bool> rm_enc = enc;
            // Punctucing 구현한부분. Puncturing 구현한다음 나중에 추가하자
            /*
                vector<bool> rm_enc = poalr.rate_matching(&enc);
            */

            // BPSK+ AWGN
            normal_distribution<double> norm_dist(0, sigma[i]);
            vector<double> r; r.reserve(N);
            for (auto e : rm_enc)
                r.push_back(1 - 2.0 * e + norm_dist(random_engine));
            
            // Compute soft bits as LLR(Log Likelihood Ratio)
            vector <double> llr; llr.reserve(N);
            for (auto e : r){
                llr.push_back(2.0 * e / (sigma[i]*sigma[i])); // 여기는 한번 더 체크해보자 
            }

            // Rate recovery (Puncturing)
            // Not yet
            //vector<double> rr_llr = polar.rate_recovery(&llr); 나중에 puncturing 공부하게되면 적용해보자.
            vector<double> rr_llr = llr; 
            
            // SC decoder
            vector<bool> msg_cap = polar.sc_decoder(&rr_llr);

            // SCL decoder 
            //vector<bool> msg_cap = polar.scl_decoder(&rr_llr);

            // count errors
            int n_errs = 0;
            for(int j=0; j<K; j++){
                if(msg[j]!=msg_cap[j]){
                    n_errs++;
                }
            }

            // n_blk_errs는 하나라도 error가 생겼다면 카운트해준다
            if (n_errs){
                n_bit_errs_local[i] += n_errs;
                n_blk_errs_local[i]++;
            }

            n_blks_done_local += 1;

            // ALLReduce
            MPI_Allreduce(&n_bit_errs_local[i], &n_bit_errs_total[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&n_blk_errs_local[i], &n_blk_errs_total[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&n_blks_done_local, &n_blks_done_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            
            // BER and BLER
     
            ber[i] = n_bit_errs_total[i] * 1.0 / (K * n_blks_done_total);
            bler[i] = n_blk_errs_total[i] * 1.0 / n_blks_done_total;
            //cout<< "rank 0"<<endl;
          
            // 실시간 진행상태 출력
            if (world_rank == 0 && n_blks_done_total % 10 == 0) {
                //double elapsed_time = (clock() - ts) / (double)CLOCKS_PER_SEC;
                str = format("|| %.1f  || %12.5e || %12.5e || %10d || %10d || %10d || %9.1f s  ||",
                             EbN0_dB[i], ber[i], bler[i], n_bit_errs_total[i], n_blk_errs_total[i], n_blks_done_total, (clock()-ts)/(double)CLOCKS_PER_SEC);
                cout << string(print_len, '\b');
                cout << str << flush;
                print_len = str.length();
            }
         }
         // print progress when one SNR is finished
         if(world_rank == 0){
         str = format("|| %.1f  || %12.5e || %12.5e || %10d || %10d || %10d || %9.1f s  ||",
         EbN0_dB[i], ber[i], bler[i], n_bit_errs_total[i], n_blk_errs_total[i], n_blks_done_total, (clock()-ts)/(double)CLOCKS_PER_SEC);
         cout << string(print_len, '\b') << str << flush << endl;
         }
    }

    // 파일 경로 및 이름 지정
     if(world_rank == 0){     
     string folder = "result/SC/";
     string filename = decoder_name + "_" + to_string(N) + "_" + to_string(K) + ".txt";
     string filepath = folder + filename;

     ofstream fout(filepath);
     if (!fout.is_open()) {
          cerr << "Failed to open the file: " << filepath << endl;
     }

     // 헤더 출력
     fout << "==================================================\n";
     fout << "Eb/N0             : " 
                                    << EBN0_BEGIN << "[dB] ~ " 
                                    << EBN0_END << "[dB], " 
                                    << EBN0_INTERVAL << "[dB]\n";

     fout << "Target Frame Error: " << TARGET_FRAME_ERR<<"\n";
     fout << "Max num frame     : " << N_MAX_BLK <<"\n";
     fout << "N                 : " << N << "\n";
     fout << "K                 : " << K << "\n";
     fout << "Decoder           : " << decoder_name << "\n";
     fout << "==================================================\n\n";
     // BER 출력
     fout << "BER: \n";
     fout << fixed << setprecision(8);
     for (double b : ber)
          fout << b << "\n";
     fout << "\n";

     // FER 출력
     fout << "FER: \n";
     for (double f : bler)
          fout << f << "\n";

     fout.close();
     cout << "File is uploaded!: " << filepath << endl;
     }    
}



void run_CA_SCL_parallel(int world_rank, int world_size, int N, int K, int NUM_LIST, string CRC_TYPE, int EBN0_BEGIN, int EBN0_END, double EBN0_INTERVAL, int N_MAX_BLK, int TARGET_FRAME_ERR, string RM_METHOD){
    
    
    // crc polynomial
    vector<bool> crc_g = crc_polynomial(CRC_TYPE); // crc_generator 뭐쓸지는 여기서 바꾸기. 나중에 generator matric 종류에 따른 차이도 비교?
    int crc_len = crc_g.size()-1; // 왜냐하면, crc_g로 나누는 나머지가 crc 값이기에 crc_g의 승수보다 1 작아야한다.
    int K_pad = K + crc_len; // crc + info
    double R = (double)K/N;

    // instantiate a POLAR object
    POLAR polar = POLAR(N, K, crc_len, RM_METHOD);

    // 랜덤한 info bit들 생성하기
    default_random_engine random_engine(world_rank + time(0));
    bernoulli_distribution  bern_dist;

    // parameters about simulation 
    int linspace_num_point = (EBN0_END - EBN0_BEGIN)/EBN0_INTERVAL + 1;
    
    
    vector<double> EbN0_dB = linspace(EBN0_BEGIN,EBN0_END,linspace_num_point); // 7 point in 0dB~3dB
    string decoder_name = "CA_SCL";

    vector<double> EbN0;
    vector<double> sigma;
    //vector<double> N0;
    for (auto e : EbN0_dB){
        double x = pow(10.0, e/10);
        EbN0.push_back(x);
        sigma.push_back(sqrt(1/(2*R*x)));
    }

    vector<double> ber(EbN0.size(), 0), bler(EbN0.size(), 0);
    vector<int> n_bit_errs_local(EbN0.size(), 0), n_blk_errs_local(EbN0.size(), 0);
    vector<int> n_bit_errs_total(EbN0.size(), 0), n_blk_errs_total(EbN0.size(), 0);

    if(world_rank == 0){
    // 헤더 출력
    cout << "============================================================\n";
    cout << "Eb/N0             : " 
                                    << EBN0_BEGIN << "[dB] ~ " 
                                    << EBN0_END << "[dB], " 
                                    << EBN0_INTERVAL << "[dB]\n";

    cout << "Target Frame Error: " << TARGET_FRAME_ERR<<"\n";
    cout << "Max num frame     : " << N_MAX_BLK<<"\n";
    cout << "N                 : " << N << "\n";
    cout << "K                 : " << K << "\n";
    cout << "CRC               : " << CRC_TYPE << "\n";
    cout << "Decoder           : " << decoder_name << "\n";
    cout << "L                 : " << NUM_LIST <<"\n";
    cout << "============================================================\n\n";

    cout << "|| Eb/N0 ||     BER     ||     BLER    || n_bit_errs || n_blk_errs || n_blk_done || elapsed time ||" << endl;
    cout << "===================================================================================================" << endl;
    }

    // Loop for each EbN0 (SNR)
    for(size_t i = 0; i<EbN0.size(); i++){
        // print progress
        string str = format("\nNow running EbN0: %.2f dB [%d of %d]", EbN0_dB[i], i+1, EbN0.size()); // format for console display.  str은 콘솔 출력 위한 버퍼로 계속 사용할것.
        //cout<<str<<endl;

        int print_len = 0; // 추후에 출력할때 지우고 다시 쓰는 동작을 구현하기 위해 필요하다

        int n_blks_done_local = 0;
        int n_blks_done_total = 0;
        clock_t ts = clock();
        // 1. max_block까지 돌아가거나 block error가 target
        while((n_blks_done_total<N_MAX_BLK)&&(n_blk_errs_total[i])<TARGET_FRAME_ERR){
            //cout << "n_blks_done: " << n_blks_done << ", n_blk_errs[i]: " << n_blk_errs[i] << endl;
            // generate random bit stream
            vector<bool> msg;
            msg.reserve(K);
            for (int j = 0; j < K; j++)
                msg.push_back(bern_dist(random_engine));

            // generate CRC and attach it to the message
            vector<bool> crc = POLAR::crc_gen(&msg, crc_g); // crc 계산
            vector<bool> crc_msg = msg; crc_msg.reserve(K_pad); // 메모리 미리 확보 (size K_pad)
            crc_msg.insert(crc_msg.end(), crc.begin(), crc.end()); // msg 뒤에 crc 붙여준다

            // polar encoding
            vector<bool> enc = polar.encoder(&crc_msg);
            
            // rate matching (Puncturing).. Not yet
            vector<bool> rm_enc = enc;
            // Punctucing 구현한부분. Puncturing 구현한다음 나중에 추가하자
            /*
                vector<bool> rm_enc = poalr.rate_matching(&enc);
            */

            // BPSK+ AWGN
            normal_distribution<double> norm_dist(0, sigma[i]);
            vector<double> r; r.reserve(N);
            for (auto e : rm_enc)
                r.push_back(1 - 2.0 * e + norm_dist(random_engine));
            
            // Compute soft bits as LLR(Log Likelihood Ratio)
            vector <double> llr; llr.reserve(N);
            for (auto e : r){
                llr.push_back(2.0 * e / (sigma[i]*sigma[i])); // 여기는 한번 더 체크해보자 
            }

            // Rate recovery (Puncturing)
            // Not yet
            //vector<double> rr_llr = polar.rate_recovery(&llr); 나중에 puncturing 공부하게되면 적용해보자.
            vector<double> rr_llr = llr; 
            
            // SCL decoder
            // cout<<"OK?"<<endl;
            // option 1. 아직 최적화 안된 decoder
            //vector<bool> msg_cap = polar.scl_decoder(&rr_llr, crc_g, num_List);
            // option 2. 최적화 진행되었음. 
            vector<bool> msg_cap = polar.scl_decoder_opt(rr_llr, crc_g, NUM_LIST);
            // vector<bool> msg_cap = polar.scs_decoder(&rr_llr, crc_g, num_List);
            //display_vec(msg_cap);
            

            // SCL decoder 
            //vector<bool> msg_cap = polar.scl_decoder(&rr_llr);

            // count errors
            int n_errs = 0;
            for(int j=0; j<K; j++){
                if(msg[j]!=msg_cap[j]){
                    n_errs++;
                }
            }

            // n_blk_errs는 하나라도 error가 생겼다면 카운트해준다
            if (n_errs){
                n_bit_errs_local[i] += n_errs;
                n_blk_errs_local[i]++;
            }

            n_blks_done_local += 1;

            MPI_Allreduce(&n_bit_errs_local[i], &n_bit_errs_total[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);            
            MPI_Allreduce(&n_blk_errs_local[i], &n_blk_errs_total[i], 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);            
            MPI_Allreduce(&n_blks_done_local, &n_blks_done_total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);            

            ber[i] = n_bit_errs_total[i] * 1.0 / (K * n_blks_done_total);
            bler[i] = n_blk_errs_total[i] * 1.0 / n_blks_done_total;
                  
            // print progress for every 10 blocks
            // flush 사용하여 실시간 업데이트
            //if(n_blks_done % 10 == 0){
            //    str = format("Elapsed time: %.1f seconds, # tx blocks: %d, # error blocks:%d, BER: %.10f, BLER: %.10f", (clock()-ts)/(double)CLOCKS_PER_SEC, n_blks_done, n_blk_errs[i], ber[i], bler[i]);
            //    cout<<string(print_len, '\b'); //print_len 길이만큼 백스페이스 문자 출력하여 이전에 출력된 문자열 지우기
            //    cout<<str<<flush; //새로운 진행상태 출력. flush는 버퍼를 즉시 비워서 출력이 실시간으로 반영되도록 보장해준다.
            //    print_len = str.length();// 현재 출력된 문자열의 길이를 print_len에 저장하여, 다음 출력 시 이전 문자열을 정확히 덮어쓰도록 보장해준다.
            //}



            // 실시간 진행상태 출력
            if (world_rank == 0 && n_blks_done_total % 10 == 0) {
                //double elapsed_time = (clock() - ts) / (double)CLOCKS_PER_SEC;
                str = format("|| %.1f  || %12.5e || %12.5e || %10d || %10d || %10d || %9.1f s  ||",
                             EbN0_dB[i], ber[i], bler[i], n_bit_errs_total[i], n_blk_errs_total[i], n_blks_done_total, (clock()-ts)/(double)CLOCKS_PER_SEC);
                cout << string(print_len, '\b');
                cout << str << flush;
                print_len = str.length();
            }
            

        }
        if(world_rank == 0){
        str = format("|| %.1f  || %12.5e || %12.5e || %10d || %10d || %10d || %9.1f s  ||",
            EbN0_dB[i], ber[i], bler[i], n_bit_errs_total[i], n_blk_errs_total[i], n_blks_done_total, (clock()-ts)/(double)CLOCKS_PER_SEC);
        cout << string(print_len, '\b') << str << flush << endl;
        }
    }


    // 파일 경로 및 이름 지정
    if(world_rank == 0){
    string folder = "result/CA_SCL/";
    string filename = decoder_name + "_" + to_string(N) + "_" + to_string(K) + "_"+ to_string(NUM_LIST) + "_" + to_string(crc_len) + ".txt";
    string filepath = folder + filename;

    ofstream fout(filepath);
    if (!fout.is_open()) {
        cerr << "Failed to open the file: " << filepath << endl;
    }

    // 헤더 출력
    fout << "==================================================\n";
    fout << "Eb/N0             : " 
                                    << EBN0_BEGIN << "[dB] ~ " 
                                    << EBN0_END << "[dB], " 
                                    << EBN0_INTERVAL << "[dB]\n";

    fout << "Target Frame Error: " << TARGET_FRAME_ERR<<"\n";
    fout << "Max num frame     : " << N_MAX_BLK<<"\n";
    fout << "N                 : " << N << "\n";
    fout << "K                 : " << K << "\n";
    fout << "CRC               : " << CRC_TYPE << "\n";
    fout << "Decoder           : " << decoder_name << "\n";
    fout << "L                 : " << NUM_LIST <<"\n";
    fout << "==================================================\n\n";

    // BER 출력
    fout << "BER: \n";
    fout << fixed << setprecision(8);
    for (double b : ber)
        fout << b << "\n";
    fout << "\n";

    // FER 출력
    fout << "FER: \n";
    for (double f : bler)
        fout << f << "\n";

    fout.close();
    cout << "File is uploaded!: " << filepath << endl;

    }
}


string format(const char* fmt, ...) {
    int size = 512;
    // Chatgpt 추천 수정코드
    //char* buffer = 0;
    //buffer = new char[size];
    char* buffer = new char[size];
    va_list vl;
    va_start(vl, fmt);
    int nsize = vsnprintf(buffer, size, fmt, vl);
    va_end(vl); // 
    if (size <= nsize) { //fail delete buffer and try again
        delete[] buffer;
        buffer = 0;
        buffer = new char[nsize + 1]; //+1 for /0
        nsize = vsnprintf(buffer, size, fmt, vl);
    }
    string ret(buffer);
    //va_end(vl);
    delete[] buffer;
    return ret;
}

vector<double> linspace(double start, double ed, int num) {
    // catch rarely, throw often
    if (num < 2) {
        // 쳇지피티가 고치라고 한다
        // throw new exception();
        throw exception();
    }
    int partitions = num - 1;
    vector<double> pts;
    // length of each segment
    double length = (ed - start) / partitions;
    // first, not to change
    pts.push_back(start);
    for (int i = 1; i < num - 1; i++) {
        pts.push_back(start + i * length);
    }
    // last, not to change
    pts.push_back(ed);
    return pts;
}



