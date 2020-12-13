/*
 ============================================================================
 Name        : matrixLib.c
 Author      : Lorenz

 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <limits.h>
#include "simple-matrix.h"
#include "mathHelper.h"
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#define CHECK_IF_DOUBLE_ZERO(nmbr) ((nmbr<0.00000000001) ? TRUE : FALSE)

static inline int getRandomWithRange(int min, int max) {
	return rand() % (max + 1 - min) + min;

}

//ZEILEN,SPALTEN
Matrix * initMatrix(size_t numberOfRows, size_t numberOfColums) {

	Matrix * newMatrix = NULL;
	newMatrix = (Matrix*) malloc(sizeof(Matrix));
	newMatrix->cells = NULL;
	newMatrix->nmbrOfColums = numberOfColums;
	newMatrix->nmbrOfRows = numberOfRows;
	newMatrix->cells = (double **) malloc(numberOfRows * sizeof(double *));
	for (int i = 0; i < numberOfRows; i++) {
		newMatrix->cells[i] = (double *) malloc(
				numberOfColums * sizeof(double));
	}
	return newMatrix;

}


void fillMatrixWithRandomValues(Matrix * m) {
	const int max = 100;
	const int min = 1;

	for (size_t i = 0; i < m->nmbrOfRows; i++) {
		for (size_t j = 0; j < m->nmbrOfColums; j++) {
			m->cells[i][j] = getRandomWithRange(min, max);

		}
	}

}

void freeMatrix(Matrix * m) {
	free(m);
}
void printMatrix(Matrix * m) {
	if (m == NULL) {
		printf("\n");
		return;
	}
	for (size_t i = 0; i < m->nmbrOfRows; i++) {
		for (size_t j = 0; j < m->nmbrOfColums; j++) {
			printf("%f\t", m->cells[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

Matrix * matrixProduct(Matrix *mA, Matrix *mB) {

	if (mA->nmbrOfColums != mB->nmbrOfRows) {
		return NULL;
	}
	Matrix * mC = initMatrix(mA->nmbrOfRows, mB->nmbrOfColums);
	for (size_t i = 0; i < mC->nmbrOfRows; i++) {
		for (size_t j = 0; j < mC->nmbrOfColums; j++) {
			for (size_t k = 0; k < mA->nmbrOfColums; k++) {
				mC->cells[i][j] = mC->cells[i][j]
						+ mA->cells[i][k] * mB->cells[k][j];
			}
		}

	}

	return mC;

}

bool copyMatrix(Matrix * dest, Matrix * src) {
	if (dest->nmbrOfColums != src->nmbrOfColums) {
		return FALSE;
	}
	if (dest->nmbrOfRows != src->nmbrOfRows) {
		return FALSE;
	}

	for (size_t i = 0; i < src->nmbrOfRows; i++) {
		for (size_t j = 0; j < src->nmbrOfColums; j++) {
			dest->cells[i][j] = src->cells[i][j];
		}
	}
	return TRUE;
}

Matrix * matrixInverse(Matrix *mA) {

	if (mA->nmbrOfColums != mA->nmbrOfRows) {
		return NULL;
	}

	size_t n = mA->nmbrOfColums;
	double ratio = 0.0;
	Matrix * mInv = NULL;
	Matrix * m = NULL;
	mInv = initMatrix(n, n);
	m = initMatrix(n, n);
	copyMatrix(m, mA);
	/*Identity matrix*/
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			mInv->cells[i][j] = (i == j) ? 1 : 0;
		}
	}
	/*Gauﬂ Jordan Elimination*/
	for (size_t i = 0; i < n; i++) {
		if (CHECK_IF_DOUBLE_ZERO(m->cells[i][i])) {
			printf("Error\n");
			return NULL;
		}
		for (size_t j = 0; j < n; j++) {
			if (i != j) {
				ratio = m->cells[j][i] / m->cells[i][i];

				for (size_t k = 0; k < n; k++) {

					m->cells[j][k] = m->cells[j][k] - ratio * m->cells[i][k];

				}
				for (size_t k = 0; k < n; k++) {

					mInv->cells[j][k] = mInv->cells[j][k]
							- ratio * mInv->cells[i][k];
				}
			}
		}
	}
	/*Normalize*/
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			mInv->cells[i][j] = mInv->cells[i][j] / m->cells[i][i];
		}

	}

	freeMatrix(m);
	return mInv;

}

Matrix * transpose(Matrix * m) {
	if (m == NULL) {
		return NULL;
	}
	if (m->nmbrOfColums != m->nmbrOfRows) {
		return NULL;
	}
	size_t n = m->nmbrOfColums;
	Matrix * transposeMatrix = initMatrix(n, n);

	for (size_t i = 0; i < m->nmbrOfRows; i++) {
		for (size_t j = 0; j < m->nmbrOfColums; j++) {
			if (i != j) {
				transposeMatrix->cells[j][i] = m->cells[i][j];
			} else {
				transposeMatrix->cells[i][j] = m->cells[i][j];
			}
		}
	}
	return transposeMatrix;
}
static double getCofactor(Matrix * m, Matrix * temp, int p, int q, int n) {

	size_t i = 0;
	size_t j = 0;

	for (size_t r = 0; r < n; r++) {
		for (size_t c = 0; c < n; c++) {

			if (r != p && q != c) {
				temp->cells[i][j++] = m->cells[r][c];
			}
			if (j == n - 1) {
				j = 0;
				i++;
			}
		}
	}
}

double calculateDeterminante(Matrix ** m, int recCounter) {
	Matrix * pntr = *m;

	if (pntr->nmbrOfColums != pntr->nmbrOfRows) {
		return 0;
	}

	int n = pntr->nmbrOfColums - recCounter;
	recCounter++;

	/*Just one Element*/
	if (n == 1) {
		return pntr->cells[0][0];
	}
	/*2x2*/
	if (n == 2) {
		return pntr->cells[0][0] * pntr->cells[1][1]
				- pntr->cells[0][1] * pntr->cells[1][0];
	}
	/*3x3: Sarrus*/
	if (n == 3) {
		return pntr->cells[0][0] * pntr->cells[1][1] * pntr->cells[2][2]
				+ pntr->cells[0][1] * pntr->cells[1][2] * pntr->cells[2][0]
				+ pntr->cells[0][2] * pntr->cells[1][0] * pntr->cells[2][1]
				- pntr->cells[0][2] * pntr->cells[1][1] * pntr->cells[2][0]
				- pntr->cells[0][0] * pntr->cells[1][2] * pntr->cells[2][1]
				- pntr->cells[0][1] * pntr->cells[1][0] * pntr->cells[2][2];
	}
	double det = 0;
	int sign = 1;
	Matrix * temp = NULL;
	temp = initMatrix(n, n);

	for (size_t i = 0; i < n; i++) {
		getCofactor(pntr, temp, 0, i, n);
		det += sign * pntr->cells[0][i]
				* calculateDeterminante(&temp, recCounter);
		sign = -sign;
	}

	freeMatrix(temp);

	return det;
}


static bool notInArray(double * array, double nmbr, int index, int len) {

	for (int i = 0; i < len; i++) {
		if (i != index && CHECK_IF_DOUBLE_ZERO(fabs(array[i] - nmbr))) {
			return FALSE;
		}
	}
	return TRUE;
}



//int main(void) {
//	srand(time(NULL));
//	Matrix * matrixA = NULL;
//	matrixA = initMatrix(3, 3);
//	fillMatrixWithRandomValues(matrixA);
//	printMatrix(matrixA);
//	printf("Det: %f\n", calculateDeterminante(&matrixA, 0));
//	getZeroOfFunction(3, 3.0, 6.0, -4.0);
//
////	Matrix * matrixB = NULL;
////	matrixB = initMatrix(4, 1);
////	fillMatrixWithRandomValues(matrixB);
////
////
////	Matrix * matrixC = NULL;
////	matrixC = matrixProduct(matrixA, matrixB);
////
////
////	Matrix * matrixD = matrixInverse(matrixA);
////	printMatrix(matrixD);
////	Matrix * matrixE = transpose(matrixD);
////	printMatrix(matrixE);
////	freeMatrix(matrixA);
////	freeMatrix(matrixB);
////	freeMatrix(matrixC);
////	freeMatrix(matrixD);
////	freeMatrix(matrixE);
//	return EXIT_SUCCESS;
//}
