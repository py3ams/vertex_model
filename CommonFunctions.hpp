#ifndef COMMONFUNCTIONS_H
#define COMMONFUNCTIONS_H

#include<cmath>
#include<vector>
#define pi 3.1415926535897932384626433832795

using namespace std;

// find straight line distance (i.e. curvature not accounted for) between two points
double findStraightLineDistanceBetweenTwoNodes(double* node_1_position, double* node_2_position) {
	double sum_squared_differences = 0;
	
	for(unsigned dim=0;dim<2;dim++) {
		sum_squared_differences += pow(node_2_position[dim]-node_1_position[dim], 2);
	}
	
	return sqrt(sum_squared_differences);	
}

// find the unit vector between two positions
void findUnitVector(double* node_1_position, double* node_2_position, double* unit_vector) {
	
	double distance_between_nodes = findStraightLineDistanceBetweenTwoNodes(node_1_position, node_2_position);
	
	for(unsigned dim=0;dim<2;dim++) {
		unit_vector[dim] = node_2_position[dim] - node_1_position[dim];
		unit_vector[dim] /= distance_between_nodes;
	}
}

// find the norm of a vector
double findNorm(double* vector) {
	double sum_squares = 0;
	
	for(unsigned dim=0;dim<2;dim++) {
		sum_squares += pow(vector[dim], 2);
	}
	
	return sqrt(sum_squares);
}

	
#endif
