/*
 * matrixLib.h
 *
 *  Created on: 14.10.2020
 *      Author: Lorenz Bornholdt
 */
#include <stdint.h>
#include <stdio.h>
#ifndef SIMPLEMATRIX_H_
#define SIMPLEMATRIX_H_

typedef struct {
	double ** cells;
	size_t nmbrOfRows;
	size_t nmbrOfColums;
} Matrix;


Matrix * initMatrix(size_t numberOfRows, size_t numberOfColums);
void fillMatrixWithRandomValues(Matrix * m);
void freeMatrix(Matrix * m);
void printMatrix(Matrix * m);
Matrix * matrixInverse(Matrix *mA) ;
#endif /* MATRIXLIB_H_ */
