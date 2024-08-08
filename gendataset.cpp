#include <iostream>
#include <time.h>
#include <stdlib.h>

int main(int argc, char const *argv[]){


    int feature_num = atoi(argv[1]);//2;
    int sample_num = atoi(argv[2]);//50;

    FILE *fp=fopen(argv[3],"w");

    srand((unsigned)time(NULL));

    for(int j = 0; j < sample_num; j++){
	int lable = 0;
	float x = rand()/(RAND_MAX + 1.0);
	float y = rand()/(RAND_MAX + 1.0);
	if((x > 0.5 && y > 0.5) || (x > 0.5 && y < 0.5 && ((x - 1)*(x - 1) + (y - 0.5)*(y - 0.5) < 0.25)) || (x < 0.5 && y > 0.5 && ((x - 0.5)*(x - 0.5) + (y -1)*(y -1) > 0.25))) {
	    lable = 0;
	}
	else {
	    lable = 1;
	}

	for(int k = 0; k < feature_num; k++) {
                float tmp_rand = rand()/(RAND_MAX + 1.0);
                float sample = 0;
                if(k == 0) {
		    sample = x;
		}
		else if(k == 1) {
		    sample = y;
		}
                fprintf(fp,"%f ",sample);
            }
            fprintf(fp,"%d\n",lable);

	//fprintf(fp,"%f %f %d\n",x,y,lable);

      /*  int lable = rand()%2;
	//float tmp_rand = rand()/(RAND_MAX + 1.0);
	if(lable == 0) {
	    
	    //float sample = 0.2 + tmp_rand/5;
	    //fprintf(fp,"%f %d\n",sample,lable);
	    for(int k = 0; k < feature_num; k++) {
		float tmp_rand = rand()/(RAND_MAX + 1.0);
		float sample = 0;
		if(k < 2) sample = tmp_rand/2;
		fprintf(fp,"%f ",sample);
	    }
	    fprintf(fp,"%d\n",lable);
	}
	else {
	    //float sample = 0.6 + tmp_rand/5;
	    for(int k = 0; k < feature_num; k++) {
		float tmp_rand = rand()/(RAND_MAX + 1.0);
		float sample = 0;
		if(k < 2) sample = 0.5 + tmp_rand/2;
                fprintf(fp,"%f ",sample);
            }
            fprintf(fp,"%d\n",lable);

	    //fprintf(fp,"%f %d\n",sample,lable);
	    //std::cout << sample << ' ' << lable << std::endl;
	}
    */
    }
    fclose(fp);

    return 0;
}
