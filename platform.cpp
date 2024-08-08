#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <mpi.h>

#include <calculategk.hpp>
#include <instruction.hpp>

#define NUM_THREADS 1
/*
class instruction {
public:
    char flag = '@';
    char inst_name[4];
    int qubit0 = -1;
    int qubit1 = -1;
    float angle = 100;
};
*/

//double calculate_gk() {
//}

using namespace std;

int main(int argc, char *argv[]) {

    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//if(rank == 0){

    int epoch = 100;

    instruction inst_list[1000];
    int inst_num = 0;
    int traning_parameter_num = 0;

    srand((unsigned)time(NULL));

    //reading network file
    FILE *fp=fopen(argv[1],"r");
    //FILE *fp=fopen("network_26.dat","r");
    FILE *fp_data = fopen(argv[2],"r");
    //FILE *fp_data = fopen("26_feature","r");
//    std::cout << fp << std::endl;
    const int buf_size = 30;
    char buf[buf_size];
    if(fp == NULL) {
        printf("can't read network file.\n");
        exit(0);
    }
    int count = 0;

    while(fgets(buf, buf_size, fp) != NULL) {
	
	instruction tmp_inst;
	int last_blank_postion = 0;
	int count_blank = 0;

	if(buf[0] == '#') traning_parameter_num++;

        for(int i = 0; i < buf_size; i++) {
	    if(buf[0] != '*' && buf[0] != '#' && buf[0] != '&' && buf[0] != '$') break;
            if(buf[i] == ' ' || buf[i] == '\n') {
                char tmp_char[10];
	        for(int t = 0; t < 10; t++) tmp_char[t] = 0;
	        int tmp_char_long = i - last_blank_postion;
                for(int j = 0; j < tmp_char_long; j++) {
                    tmp_char[j] = buf[last_blank_postion + j];
                }
		if(count_blank == 0) {
		    inst_list[count].flag = tmp_char[0];
		}
		else if(count_blank == 1) {
		    for(int j = 0; j < 10; j++) inst_list[count].inst_name[j] = tmp_char[j];
		}
		else if(count_blank == 2) {
		    inst_list[count].qubit0 = atoi(tmp_char);//tmp_char;//atoi(tmp_char);
		}
		else if(count_blank == 3) {
		    if(strcmp(inst_list[count].inst_name, "CNOT") == 0 || strcmp(inst_list[count].inst_name, "CZ") == 0 || strcmp(inst_list[count].inst_name, "SWAP") == 0) {
		        inst_list[count].qubit1 = atoi(tmp_char);//tmp_char;//atoi(tmp_char);
		    }
		    else if(strcmp(inst_list[count].inst_name, "RX") == 0 || strcmp(inst_list[count].inst_name, "RY") == 0 || strcmp(inst_list[count].inst_name, "RZ") == 0) {
			if(inst_list[count].flag == '*') {
			    inst_list[count].feature_index = atoi(tmp_char);
			}
			else {
			    inst_list[count].angle = atof(tmp_char);
			}
		    }
		    else {
	                inst_list[count].angle = atof(tmp_char);//tmp_char;//atof(tmp_char);
		    }
		}

	        last_blank_postion = i + 1;
		count_blank++;

		if(buf[i] == '\n') {
		    count++;
		    break;
		}

	    }
	}
    }
    fclose(fp);
    //end reading network file
if(rank == 0) {
    for(int i = 0; i < count; i++){
	std::cout << inst_list[i].flag << ' ' << inst_list[i].inst_name << ' ' << inst_list[i].qubit0  << ' ' << inst_list[i].qubit1 << ' ' << inst_list[i].angle << std::endl;
    }
    std::cout << std::endl;
}
//printf("debug 3\n");
    //read feature
    //FILE *fp_data = fopen(argv[2],"r");
    //pair<vector<vector<double>>,int> data_set;
    vector<pair<vector<double>,int>> data_set;
    vector<pair<vector<double>,int>> training_data_set;
    vector<pair<vector<double>,int>> test_data_set;

    /*const int feature_num = 2;
    const int sample_num = 50;
    float feature_array[sample_num][feature_num];
    int lable_array[sample_num];*/

    const int data_buf_size = 10 * 100;//(feature_num + 1);
    char data_buf[data_buf_size];
    if(fp == NULL) {
        printf("can't read dataset file.\n");
        exit(0);
    }
    int sample_count = 0;

    //printf("debug 4\n");
    /*
printf("debug 4\n");
std::cout << argv[2] << std::endl;
std::cout << fp_data << std::endl;
char* test = fgets(data_buf, data_buf_size, fp_data);
printf("debug 3\n");*/


    while(fgets(data_buf, data_buf_size, fp_data) != NULL) {
	int last_blank_postion = 0;
	int blank_count = 0;
//printf("debug 2\n");	
	vector<double> feature;
	int lable = 0;
	for(int i = 0; i < data_buf_size; i++) {
	    if(data_buf[i] == ' ') {
	        char tmp_char[10];
		for(int t = 0; t < 10; t++) tmp_char[t] = 0;
		int tmp_char_long = i - last_blank_postion;
		for(int t = 0; t < tmp_char_long; t++) {
		    tmp_char[t] = data_buf[last_blank_postion + t];
		}
		float tmp_feature = atof(tmp_char);// * (1 / 1.57) * 90.0;
		feature.push_back(tmp_feature);
		//feature_array[sample_count][blank_count] = feature;
		//blank_count++;
		last_blank_postion = i;
	    }
            else if(data_buf[i] == '\r' || data_buf[i] == '\n') {
	        char tmp_char[10];
		for(int t = 0; t < 10; t++) tmp_char[t] = 0;
		tmp_char[0] = data_buf[i - 1];
		lable = atoi(tmp_char);
                //int lable = atoi(tmp_char);
                //lable_array[sample_count] = lable;
		break;
	    }
	    //push_back(make_pair(feature,lable));
	}
	data_set.push_back(make_pair(feature,lable));
	sample_count++;
    }

    fclose(fp_data);
    for(int i = 0; i < data_set.size(); i++) {
	if(i < 0.8 * data_set.size()) {//划分训练集、测试集
	    training_data_set.push_back(data_set[i]);
	}
	else {
	    test_data_set.push_back(data_set[i]);
	}
    }
    //end read feature
//printf("debug 1\n");
/*    for(int i = 0; i < training_data_set.size(); i++) {
	for(int j = 0; j < training_data_set[i].first.size(); j++) {
	    std::cout << training_data_set[i].first[j] << ' ';
	}
	std::cout << training_data_set[i].second << std::endl;
    }
    std::cout << std::endl;
*/

//    for(int i = 0; i < sample_count; i++) {
//	for(int j = 0; j < feature_num; j++) {
//	    std::cout << feature_array[i][j] << ' ';
//	}
//	std::cout << lable_array[i] << std::endl;
//    }
//    std::cout << std::endl;

    //intial traning parameter 

//int s_num = 100;
//while(s_num--) {
//if(rank == 0){
    //float *traning_parameter_array = NULL;//new float[traning_parameter_num]();
    //int group = training_data_set.size() / size + 1;
    int group = traning_parameter_num / size;
    if(float(traning_parameter_num)/size - traning_parameter_num / size > 0.0) group += 1;
    //float* g_array = NULL;
    //int param_num;
    //int stop = 0;

if(rank == 0){
    float *traning_parameter_array = new float[traning_parameter_num]();
    for(int i = 0; i < traning_parameter_num; i++) {
	traning_parameter_array[i] = rand()/(RAND_MAX + 1.0) * 180.0 - 90.0;//初始化参数
	std::cout << traning_parameter_array[i] << ' ';
    }

    float* g_array = new float[training_data_set.size()]();
//    MPI_Bcast(traning_parameter_array, traning_parameter_num, MPI_FLOAT, 0, MPI_COMM_WORLD);
    std::cout << "traning_parameter_num: " << traning_parameter_num << std::endl << std::endl;

    //traning
    //rmdir("training");
    //mkdir("training",S_IRWXU);

    //double learning_rate = 0.9;
    double decay1 = 0.9;
    double decay2 = 0.999;
    double stepsize = 0.01;
    double eps = 1.0e-8;

    double m(0), v(0), mk(0), vk(0);

    double* last_m = new double[traning_parameter_num]();
    double* last_v = new double[traning_parameter_num]();

    //FILE *result = fopen("result_loss.data","w+");

//    int group = training_data_set.size() / size + 1;

    for(int i = 0; i < epoch; i++) {
	//float* g_array = new float[training_data_set.size()]();
	float LCE = 0;
	int training_correct_count = 0;
	int test_correct_count = 0;

	//char dump_name[30] = {0};
        //sprintf(dump_name,"training/epoch%d_dump.data",i);
        //FILE *dump = fopen(dump_name,"w+");

	//double rhoC = 0.0;

	for(int sample = 0; sample < training_data_set.size(); sample++) {
	    float *feature_array = new float[training_data_set[sample].first.size()]();//[feature_num];
	    for(int j = 0; j < training_data_set[sample].first.size(); j++) {
		feature_array[j] = training_data_set[sample].first[j];
	    }
	    //cout<<"debug1"<<endl;
	    //int lable = lable_array[sample];
	    int lable = training_data_set[sample].second;
	    //float N_prob = calculate_prob(inst_list,feature_array[sample],traning_parameter_array,sample,i);
	    float N_prob = calculate_prob(inst_list,feature_array,traning_parameter_array,-2,sample,i,0);
	    //N_prob = resu.real;
	    //rhoC += resu.imag;
	    g_array[sample] = N_prob;

	    //fprintf(dump,"epoch%d no_shift sample%d: %f\n",i,sample,N_prob);

	    delete feature_array;

            if(training_data_set[sample].second/*lable_array[sample]*/ == 0) {
                LCE -= log(N_prob);
	    }
	    else {
		LCE -= log(1 - N_prob);
	    }
	    int pred_lable = -1;
	    if(N_prob >= 0.5) {
		pred_lable = 0;
	    }
	    else {
		pred_lable = 1;
	    }
	    if(pred_lable == training_data_set[sample].second/*lable_array[sample]*/) training_correct_count++;
	}

	float training_accuracy = (float)training_correct_count / training_data_set.size();
	//if(i == 0 && (training_accuracy < 0.1 || training_accuracy > 0.9)) assert(0);

/*	FILE *fp_sc = fopen("S_C.csv","a+");
        fprintf(fp_sc,"%f\t", rhoC / 90);
        if(i == epoch -1) fprintf(fp_sc,"\n");
        fclose(fp_sc);*/

//	close(result);

        for(int sample = 0; sample < test_data_set.size(); sample++) {
            float *feature_array = new float[test_data_set[sample].first.size()]();//[feature_num];
            for(int j = 0; j < test_data_set[sample].first.size(); j++) {
                feature_array[j] = test_data_set[sample].first[j];
            }
            //int lable = lable_array[sample];
            int lable = test_data_set[sample].second;
            //float N_prob = calculate_prob(inst_list,feature_array[sample],traning_parameter_array,sample,i);
            float N_prob = calculate_prob(inst_list,feature_array,traning_parameter_array,-1,sample,i,0);
	    //N_prob = resu.real;;
            //g_array[sample] = N_prob;

            delete feature_array;

            int pred_lable = -1;
            if(N_prob >= 0.5) {
                pred_lable = 0;
            }
            else {
                pred_lable = 1;
            }
            if(pred_lable == test_data_set[sample].second/*lable_array[sample]*/) test_correct_count++;
        }

	//float training_accuracy = (float)training_correct_count / training_data_set.size();
	float test_accuracy = (float)test_correct_count / test_data_set.size();
	printf("epoch:%d loss:%f training accuracy:%f test accuracy:%f\n",i,LCE,training_accuracy,test_accuracy);

        time_t nowtime;
	time(&nowtime); //获取1970年1月1日0点0分0秒到现在经过的秒数
	tm* p = localtime(&nowtime); //将秒数转换为本地时间,年从1900算起,需要+1900,月为0-11,所以要+1
	printf("%04d:%02d:%02d %02d:%02d:%02d\n", p->tm_year + 1900, p->tm_mon + 1, p->tm_mday,p->tm_hour,p->tm_min,p->tm_sec);
                        
	//double parameter_dev[32] = {0.0};

	    int total_send_size = training_data_set.size() + traning_parameter_num + 1;
            float* send_array = new float[total_send_size]();
            for(int a = 0; a < total_send_size - 2; a++) {
                if(a < training_data_set.size()) {
                    send_array[a] = g_array[a];
                }
                else if(a >= training_data_set.size() && a < training_data_set.size() + traning_parameter_num) {
                    send_array[a] = traning_parameter_array[a - training_data_set.size()];
                }
            }
            //send_array[total_send_size - 2] = param;

	    int stop = 0;
            if(i == epoch - 1) stop = 1;
            send_array[total_send_size - 1] = stop;
            //MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
            for(int c = 1; c < rank; c++) {
                MPI_Send(send_array,total_send_size,MPI_FLOAT,c,98,MPI_COMM_WORLD);
            }

	    //for(int sample = rank * group; sample < min<int>(rank * group + group, training_data_set.size()); sample++) {
//float parameter_dev[32] = {0.0};
float *parameter_dev = new float[traning_parameter_num]();
	    //traning_parameter_num
        for(int param = rank * group; param < min<int>(rank * group + group, traning_parameter_num); param++) {
	//for(int param = 0; param < traning_parameter_num; param++) {

	    float pcost = 0;
	    //float LCE = 0;
	    /*int total_send_size = training_data_set.size() + traning_parameter_num + 2;
	    float* send_array = new float[total_send_size]();
	    for(int a = 0; a < total_send_size - 2; a++) {
		if(a < training_data_set.size()) {
	            send_array[a] = g_array[a];
		}
		else if(a >= training_data_set.size() && a < training_data_set.size() + traning_parameter_num) {
		    send_array[a] = traning_parameter_array[a - training_data_set.size()];
		}
	    }
	    send_array[total_send_size - 2] = param;

	    int stop = 0;
	    if(i == epoch - 1 && param == traning_parameter_num - 1) stop = 1;
	    send_array[total_send_size - 1] = stop;
	    //MPI_Bcast(&stop, 1, MPI_INT, 0, MPI_COMM_WORLD);
            for(int c = 1; c < rank; c++) {
		MPI_Send(send_array,total_send_size,MPI_FLOAT,c,98,MPI_COMM_WORLD);
            }*/

	    float *pcost_array = new float[training_data_set.size()]();

	    #pragma omp parallel for num_threads(NUM_THREADS)
	    for(int sample = 0; sample < training_data_set.size(); sample++) {
	    //for(int sample = rank * group; sample < min<int>(rank * group + group, training_data_set.size()); sample++) {
		
		float *feature_array = new float[training_data_set[sample].first.size()]();//[feature_num];
                for(int j = 0; j < training_data_set[sample].first.size(); j++) {
                    feature_array[j] = training_data_set[sample].first[j];
                }

		float gkp = calculate_prob(inst_list,feature_array,traning_parameter_array, param,sample,i,3.141492/2);
		//gkp = gkp_resu.real;
		//fprintf(dump,"epoch%d P_shift param%d sample%d: %f\n",i,param,sample,gkp);
		float gkm = calculate_prob(inst_list,feature_array,traning_parameter_array, param,sample,i,-3.141492/2);
		//gkm = gkm_resu.real;
		//fprintf(result,"epoch%d P_shift param%d sample%d: %f",i,param,sample,gkp);
		//fprintf(dump,"epoch%d M_shift param%d sample%d: %f\n",i,param,sample,gkm);
		float gk = (gkp - gkm)/2;

		if(training_data_set[sample].second == 0) {
		    //pcost = pcost - (1 / g_array[sample]) * gk;
		    pcost_array[sample] = (1 / g_array[sample]) * gk;
                }
                else {
		    //pcost = pcost - (1 / (1 - g_array[sample])) * ( -gk);
		    pcost_array[sample] = (1 / (1 - g_array[sample])) * ( -gk);
                }
		delete feature_array;
	    }

	    for(int sample = 0; sample < training_data_set.size(); sample++) {
                pcost -= pcost_array[sample];
            }

            delete pcost_array;

	    /*for(int c = 1; c < rank; c++) {
		float pcost_recv;
		MPI_Status status;
		MPI_Recv(&pcost_recv,1,MPI_FLOAT,c,99,MPI_COMM_WORLD,&status);
		pcost -= pcost_recv;
            }*/
	    //delete feature_array;
	    pcost = 0.25 * pcost / sample_count;
	    
	    double grad = pcost;

	    m = decay1 * last_m[param] + (1.0 - decay1)*grad;
            v = decay2 * last_v[param] + (1.0 - decay2)*grad*grad;
	    int t = epoch;
            mk = m / (1.0 - pow(decay1, t));
            vk = v / (1.0 - pow(decay2, t));

	    //adam
	//    traning_parameter_array[param] -= stepsize * (mk / (sqrt(vk) + eps));
	    //sgd
            //traning_parameter_array[param] -= 0.02 * pcost;

	//double tmp_dump = stepsize * (mk / (sqrt(vk) + eps));
	parameter_dev[param] = stepsize * (mk / (sqrt(vk) + eps));
	//printf("param_dump:%f ",tmp_dump);

	    last_m[param] = m;
	    last_v[param] = v;

	}
        //delete g_array;

	    for(int c = 1; c < rank; c++) {
                float *traning_parameter_array_recv = new float[traning_parameter_num]();

                //float pcost_recv;
                MPI_Status status;
                MPI_Recv(traning_parameter_array_recv,traning_parameter_num,MPI_FLOAT,c,99,MPI_COMM_WORLD,&status);
                //pcost -= pcost_recv;

		for(int param_recv = c * group; param_recv < min<int>(c * group + group, traning_parameter_num); param_recv++) {
		    traning_parameter_array[param_recv] -= traning_parameter_array_recv[param_recv];
		}

		delete traning_parameter_array_recv;
            }

	for(int i = 0; i < group; i++) {
            traning_parameter_array[i] -= parameter_dev[i];
        }

	delete parameter_dev;

  }
    delete g_array;

    delete last_m;
    delete last_v;

    delete traning_parameter_array;
}//end rank == 0
else {
    double decay1 = 0.9;
    double decay2 = 0.999;
    double stepsize = 0.01;
    double eps = 1.0e-8;

    double m(0), v(0), mk(0), vk(0);

    double* last_m = new double[traning_parameter_num]();
    double* last_v = new double[traning_parameter_num]();

    //cout << "debug101" << endl;
    while(1) {
	float* g_array = new float[training_data_set.size()]();
	float *traning_parameter_array = new float[traning_parameter_num]();

	int total_send_size = training_data_set.size() + traning_parameter_num + 1;
        float* send_array = new float[total_send_size]();
//cout << "debug102" << endl;
        MPI_Status status;
	MPI_Recv(send_array,total_send_size,MPI_FLOAT,0,98,MPI_COMM_WORLD,&status);
//cout << "debug103" << endl;
        //float parameter_dev[32] = {0.0};
        float *parameter_dev = new float[traning_parameter_num]();

	for(int a = 0; a < total_send_size - 1; a++) {
                if(a < training_data_set.size()) {
                    g_array[a] = send_array[a];
                }
                else if(a >= training_data_set.size() && a < training_data_set.size() + traning_parameter_num) {
                    traning_parameter_array[a - training_data_set.size()] = send_array[a];
                }
        }
//	int param = send_array[total_send_size - 2];
        int stop = send_array[total_send_size - 1];

	
    for(int param = rank * group; param < min<int>(rank * group + group, traning_parameter_num); param++) {
        //for(int param = 0; param < traning_parameter_num; param++) {

            //float pcost = 0;
        float pcost = 0;
        float *pcost_array = new float[training_data_set.size()]();	
	//int param = param_num;
	#pragma omp parallel for num_threads(NUM_THREADS)
	for(int sample = 0; sample < training_data_set.size(); sample++) {
	//for(int sample = rank * group; sample < min<int>(rank * group + group, training_data_set.size()); sample++) {


                float *feature_array = new float[training_data_set[sample].first.size()]();//[feature_num];
                for(int j = 0; j < training_data_set[sample].first.size(); j++) {
                    feature_array[j] = training_data_set[sample].first[j];
                }

                float gkp = calculate_prob(inst_list,feature_array,traning_parameter_array, param,sample,0,3.141492/2);
                //gkp = gkp_resu.real;
                //fprintf(dump,"epoch%d P_shift param%d sample%d: %f\n",i,param,sample,gkp);
                float gkm = calculate_prob(inst_list,feature_array,traning_parameter_array, param,sample,0,-3.141492/2);
                //gkm = gkm_resu.real;
                //fprintf(result,"epoch%d P_shift param%d sample%d: %f",i,param,sample,gkp);
                //fprintf(dump,"epoch%d M_shift param%d sample%d: %f\n",i,param,sample,gkm);
                float gk = (gkp - gkm)/2;

                if(training_data_set[sample].second == 0) {
                    //pcost = pcost - (1 / g_array[sample]) * gk;
		    pcost_array[sample] = (1 / g_array[sample]) * gk;
                }
                else {
                    //pcost = pcost - (1 / (1 - g_array[sample])) * ( -gk);
		    pcost_array[sample] = (1 / (1 - g_array[sample])) * ( -gk);
                }
                delete feature_array;
        }

	    for(int sample = 0; sample < training_data_set.size(); sample++) {
                pcost -= pcost_array[sample];
            }

            delete pcost_array;

	pcost = 0.25 * pcost / sample_count;

            double grad = pcost;

            m = decay1 * last_m[param] + (1.0 - decay1)*grad;
            v = decay2 * last_v[param] + (1.0 - decay2)*grad*grad;
            int t = epoch;
            mk = m / (1.0 - pow(decay1, t));
            vk = v / (1.0 - pow(decay2, t));

            //adam
            //traning_parameter_array[param] -= stepsize * (mk / (sqrt(vk) + eps));
            //sgd
            //traning_parameter_array[param] -= 0.02 * pcost;

        //double tmp_dump = stepsize * (mk / (sqrt(vk) + eps));
        //parameter_dev[param] = tmp_dump;
        //printf("param_dump:%f ",tmp_dump);

	    parameter_dev[param] = stepsize * (mk / (sqrt(vk) + eps));

            last_m[param] = m;
            last_v[param] = v;


    }

	MPI_Send(parameter_dev,traning_parameter_num,MPI_FLOAT,0,99,MPI_COMM_WORLD);
	delete g_array;
	delete traning_parameter_array;
	delete send_array;
	delete parameter_dev;

	if(stop) break;
    }

    delete last_m;
    delete last_v;
}
    return 0;
}
