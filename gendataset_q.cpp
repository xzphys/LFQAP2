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
//	float y = rand()/(RAND_MAX + 1.0);
	if(x > 0.5) {
	    fprintf(fp,"%f ",0.0);
	    fprintf(fp,"%f ",90.0);
	    fprintf(fp,"%f ",90.0);
	    fprintf(fp,"%f ",0.0);
	    fprintf(fp,"%d\n",0);
	}
	else {
	    fprintf(fp,"%f ",0.0);//rand()/(RAND_MAX + 1.0));
	    fprintf(fp,"%f ",0.0);//rand()/(RAND_MAX + 1.0));
	    fprintf(fp,"%f ",0.0);//rand()/(RAND_MAX + 1.0));
	    fprintf(fp,"%f ",0.0);//rand()/(RAND_MAX + 1.0));
	    fprintf(fp,"%d\n",1);
	}

    }
    fclose(fp);

    return 0;
}
