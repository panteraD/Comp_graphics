//
// Created by andrey on 14.05.15.
//

#include "geometry.h"
#include <math.h>

Vec3 * vectorAdd3D(Vec3 *t0, Vec3 *t1, Vec3 *res) {

    (*res)[0] = (*t0)[0] + (*t1)[0];
    (*res)[1] = (*t0)[1] + (*t1)[1];
    (*res)[2] = (*t0)[2] + (*t1)[2];
    return res;
}



Vec3 * vectorSub3D(Vec3 *t0, Vec3 *t1, Vec3 *res) {
    (*res)[0] = (*t0)[0] - (*t1)[0];
    (*res)[1] = (*t0)[1] - (*t1)[1];
    (*res)[2] = (*t0)[2] - (*t1)[2];
    return res;
}




Vec3 * vectorMultSimple3D(Vec3 *t0, double a, Vec3 *res) {

    (*res)[0] = (*t0)[0] * a;
    (*res)[1] = (*t0)[1] * a;
    (*res)[2] = (*t0)[2] * a;
    return res;
}



Vec3 * vectorMult3D (Vec3 *t0, Vec3 *t1, Vec3 *res) {

    (*res)[0] = (*t0)[1] * (*t1)[2] - (*t0)[2] * (*t1)[1];   //aybz-azby
    (*res)[1] = (*t0)[2] * (*t1)[0] - (*t0)[0] * (*t1)[2];   //azbx-axbz
    (*res)[2] = (*t0)[0] * (*t1)[1] - (*t0)[1] * (*t1)[0];   //axby-aybx
    return res;
}

double vectorMultScalar (Vec3 *t0, Vec3 *t1) {

    return ((*t0)[0] * (*t1)[0] + (*t0)[1] * (*t1)[1] + (*t0)[2] * (*t1)[2]);
}



double  norm(Vec3 *t0) {
    return (double) sqrt((*t0)[0] * (*t0)[0] + (*t0)[1] * (*t0)[1] + (*t0)[2] * (*t0)[2]);

}


Vec3 * normalize(Vec3 * t0) {

    double n = (double) (1.0f / norm(t0));
    vectorMultSimple3D(t0, n, t0);
    return t0;
}


void Round(Vec3 *t0) {
    (*t0)[0] = (int) ((*t0)[0] + 0.5);
    (*t0)[1] = (int) ((*t0)[1] + 0.5);
    (*t0)[2] = (int) ((*t0)[2] + 0.5);
}



//TODO delete?
void line (tgaImage *image,
           int x0, int y0,
           int x1, int y1,
           tgaColor color)
{
    int steep = 0;
    if (abs(y1 - y0) > abs(x1 - x0)) {
        steep = 1;
        swap(&x0, &y0);
        swap(&x1, &y1);
    }

    if (x0 > x1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
    }

    int x;
    double y;
    double k = ((double)(y1 - y0))/(x1 - x0);
    for (x = x0, y = y0; x <= x1; ++x, y += k) {
        if (steep) {
            tgaSetPixel(image, y, x, color);
        } else {
            tgaSetPixel(image, x, y, color);
        }
    }
}

void swap(int *a, int *b) {
    int t = *a;
    *a = *b;
    *b = t;
}

void swapVec3(Vec3 *t0, Vec3 *t1) {
    double temp;
    for (int i = 0; i < 3; i++) {
        temp = (*t0)[i];
        (*t0)[i] = (*t1)[i];
        (*t1)[i] = temp;
    }
}

int abs(int a) {
    return (a >= 0) ? a : -a;
}

void printVec3(Vec3 *t0){
    printf("%f %f %f \n" ,(*t0)[0],(*t0)[1],(*t0)[2]);

}

Vec4  * transform3D4D(Vec3 *v3, Vec4 *v4){
    for(int i=0;i<3; i++){
        (*v4)[i]=(*v3)[i];
    }
    (*v4)[3]=1.;
    return v4;
}

Vec3 * transform4D3D( Vec4 *v4,Vec3 *v3 ){
    for(int i=0;i<3; i++){
        (*v3)[i]=(*v4)[i]/(*v4)[3];
    }

    return v3;
}

Vec4 * MatrixMuliply( Matrix *matrix, Vec4 *vec4, Vec4 *result){
    for(int i=0; i< 4;i++){
        double sum=0;
        for(int j=0;j<4;j++){
            sum+=(*matrix)[j+i*4]*(*vec4)[j];
        }
        (*result)[i]=sum;
    }
    return result;

}

void lookAt(Vec3 *eye, Vec3 *center, Vec3 *up) {
    Vec3 z;
    vectorSub3D(eye, center, &z);
    normalize(&z);
    
}
/*
 * Matrix Matrix::transpose() {
+    Matrix result(cols, rows);
+    for(int i=0; i<rows; i++)
+        for(int j=0; j<cols; j++)
+            result[j][i] = m[i][j];
+    return result;
+}
+
+Matrix Matrix::inverse() {
+    assert(rows==cols);
+    // augmenting the square matrix with the identity matrix of the same dimensions a => [ai]
+    Matrix result(rows, cols*2);
+    for(int i=0; i<rows; i++)
+        for(int j=0; j<cols; j++)
+            result[i][j] = m[i][j];
+    for(int i=0; i<rows; i++)
+        result[i][i+cols] = 1;
+    // first pass
+    for (int i=0; i<rows-1; i++) {
+        // normalize the first row
+        for(int j=result.cols-1; j>=0; j--)
+            result[i][j] /= result[i][i];
+        for (int k=i+1; k<rows; k++) {
+            float coeff = result[k][i];
+            for (int j=0; j<result.cols; j++) {
+                result[k][j] -= result[i][j]*coeff;
+            }
+        }
+    }
+    // normalize the last row
+    for(int j=result.cols-1; j>=rows-1; j--)
+        result[rows-1][j] /= result[rows-1][rows-1];
+    // second pass
+    for (int i=rows-1; i>0; i--) {
+        for (int k=i-1; k>=0; k--) {
+            float coeff = result[k][i];
+            for (int j=0; j<result.cols; j++) {
+                result[k][j] -= result[i][j]*coeff;
+            }
+        }
+    }
+    // cut the identity matrix back
+    Matrix truncate(rows, cols);
+    for(int i=0; i<rows; i++)
+        for(int j=0; j<cols; j++)
+            truncate[i][j] = result[i][j+cols];
+    return truncate;
+}
 * */


