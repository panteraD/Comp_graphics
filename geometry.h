//
// Created by andrey on 14.05.15.
//

#ifndef COMP_GRAPHICS_GEOMETRY_H
#define COMP_GRAPHICS_GEOMETRY_H

#include <math.h>
#include "model.h"
#include "stdio.h"

typedef double Vec4[4];
typedef double Matrix[16];

Vec3 * vectorAdd3D(Vec3 *t0, Vec3 *t1, Vec3 *res) ;

Vec3 * vectorSub3D(Vec3 *t0, Vec3 *t1, Vec3 *res);

Vec3 * vectorMult3D (Vec3 *t0, Vec3 *t1, Vec3 *res);

Vec3 * vectorMultSimple3D(Vec3 *t0, double a, Vec3 *res) ;
double  vectorMultScalar (Vec3 *t0, Vec3 *t1);
Vec3 * normalize(Vec3 * t0);

void Round(Vec3 *t0);

void swapVec3(Vec3 *t0, Vec3 *t1);
void printVec3(Vec3 *t0);

double  norm(Vec3 *t0);

Vec3 * initVec3(double x, double y, double z, Vec3 * result);




void swap(int *a, int *b);
void swapD(double *a, double *b);
int abs(int a);


Vec4  * transform3D4D(Vec3 *v3, Vec4 *v4);
Vec3 * transform4D3D( Vec4 *v4,Vec3 *v3 );
Vec4 *matrixMultiplyMatrixVec4(Matrix *matrix, Vec4 *vec4, Vec4 *result);
Matrix  *matrixMultuplyMatrixMatrix(Matrix *A, Matrix *B, Matrix * result);
Matrix * lookAt(Vec3 *eye, Vec3 *center, Vec3 *up, Matrix *ModelView);
Matrix * viewport(int x, int y, int width, int height, int depth, Matrix * viewPort);
Matrix * transponse(Matrix *x);
Matrix * inversion(Matrix *x);
Matrix * identity(Matrix *x);
Matrix * clone(Matrix *x, Matrix *y);


#endif //COMP_GRAPHICS_GEOMETRY_H