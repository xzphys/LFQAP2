#include <instruction.hpp>
#include <simulator.hpp>

float calculate_pg(instruction inst_list[], float feature_array[], float *traning_parameter_array, int param, int sample, int epoch, float shift);
float calculate_prob(instruction inst_list[], float feature_array[], float *traning_parameter_array, int param, int sample, int epoch, float shift);
