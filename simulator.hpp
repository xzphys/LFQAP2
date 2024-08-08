#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <vector>
#include <algorithm>

using namespace std;

typedef struct complextmp {
    double real, imag;
    double assist;
}Complex;

typedef struct qr {
    unsigned short size;
    Complex *matrix;
}quReg;

void Qprint(quReg *qr);
quReg* simulator(vector <vector<char>> fp_name);
quReg* newQuReg(size_t n);
//quReg* simulator(FILE *fp);

#endif
