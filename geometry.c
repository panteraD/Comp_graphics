//
// Created by andrey on 14.05.15.
//

#include "geometry.h"

Vec3 *vectorAdd3D(Vec3 *t0, Vec3 *t1, Vec3 *res) {

    (*res)[0] = (*t0)[0] + (*t1)[0];
    (*res)[1] = (*t0)[1] + (*t1)[1];
    (*res)[2] = (*t0)[2] + (*t1)[2];
    return res;
}


Vec3 *vectorSub3D(Vec3 *t0, Vec3 *t1, Vec3 *res) {
    (*res)[0] = (*t0)[0] - (*t1)[0];
    (*res)[1] = (*t0)[1] - (*t1)[1];
    (*res)[2] = (*t0)[2] - (*t1)[2];
    return res;
}


Vec3 *vectorMultSimple3D(Vec3 *t0, double a, Vec3 *res) {

    (*res)[0] = (*t0)[0] * a;
    (*res)[1] = (*t0)[1] * a;
    (*res)[2] = (*t0)[2] * a;
    return res;
}


Vec3 *vectorMult3D(Vec3 *t0, Vec3 *t1, Vec3 *res) {

    (*res)[0] = (*t0)[1] * (*t1)[2] - (*t0)[2] * (*t1)[1];   //aybz-azby
    (*res)[1] = (*t0)[2] * (*t1)[0] - (*t0)[0] * (*t1)[2];   //azbx-axbz
    (*res)[2] = (*t0)[0] * (*t1)[1] - (*t0)[1] * (*t1)[0];   //axby-aybx
    return res;
}

double vectorMultScalar(Vec3 *t0, Vec3 *t1) {

    return ((*t0)[0] * (*t1)[0] + (*t0)[1] * (*t1)[1] + (*t0)[2] * (*t1)[2]);
}


double norm(Vec3 *t0) {
    return (double) sqrt((*t0)[0] * (*t0)[0] + (*t0)[1] * (*t0)[1] + (*t0)[2] * (*t0)[2]);

}


Vec3 *normalize(Vec3 *t0) {

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
void line(tgaImage *image,
          int x0, int y0,
          int x1, int y1,
          tgaColor color) {
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
    double k = ((double) (y1 - y0)) / (x1 - x0);
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

void printVec3(Vec3 *t0) {
    printf("%f %f %f \n", (*t0)[0], (*t0)[1], (*t0)[2]);

}

Vec4 *transform3D4D(Vec3 *v3, Vec4 *v4) {
    for (int i = 0; i < 3; i++) {
        (*v4)[i] = (*v3)[i];
    }
    (*v4)[3] = 1.;
    return v4;
}

Vec3 *transform4D3D(Vec4 *v4, Vec3 *v3) {
    for (int i = 0; i < 3; i++) {
        (*v3)[i] = (*v4)[i] / (*v4)[3];
    }

    return v3;
}

Vec4 *MatrixMuliply(Matrix *matrix, Vec4 *vec4, Vec4 *result) {
    for (int i = 0; i < 4; i++) {
        double sum = 0;
        for (int j = 0; j < 4; j++) {
            sum += (*matrix)[j + i * 4] * (*vec4)[j];
        }
        (*result)[i] = sum;
    }
    return result;

}

// ???
Matrix * lookAt(Vec3 *eye, Vec3 *center, Vec3 *up, Matrix *ModelView) {
    Vec3 z;
    vectorSub3D(eye, center, &z);
    normalize(&z);
    Vec3 x;
    vectorMult3D(up, &z, &x);
    normalize(&x);
    Vec3 y;
    vectorMult3D(&z, &x, &y);
    normalize(&y);
    Matrix Minv;
    identity(&Minv);
    Matrix Tr;
    identity(&Tr);
    for (int i = 0; i < 3; ++i) {
        Minv[i+0*4]= x[i];
        Minv[i+1*4]= y[i];
        Minv[i+2*4]= z[i];
        Tr[i+3*4]=-(*center)[i];

    }
    MatrixMuliply(&Minv,&Tr,ModelView);
    return ModelView;

}



Matrix *transponse(Matrix *x) {
    Matrix copy;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            copy[i + j * 4] = (*x)[j + i * 4];
        }
    }
    //copy matrix  values into ref
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            (*x)[j + i * 4] = copy[j + i * 4];
        }
    }
    return x;

}

//not working
Matrix *inverse(Matrix *x) {
    Matrix result, E;
    //copy + init E
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[j + i * 4] = (*x)[j + i * 4];
            if (i == j) {
                E[j + i * 4] = 1.;
            }
            else {
                E[j + i * 4] = 0.;
            }
        }
    }
    //fp
    for (int i = 0; i < 3; i++) {

        for (int j = 0; j < 4; j++) {
            result[j + i * 4] /= result[i + i * 4];
            E[j + i * 4] /= result[i + i * 4];
        }

        for (int k = i + 1; k < 4; k++) {
            double coeff = result[i + k * 4];
            for (int j = 0; j < 4; j++) {
                result[j + k * 4] -= result[j + i * 4] * coeff;
                E[j + k * 4] -= E[j + i * 4] * coeff;
            }
        }
        //normalize
        for (int l = 0; l < 4; ++l) {
            E[l + 3 * 4] /= result[3 + 3 * 4];
        }

        //second pass
        for (int i = 3; i > 0; i--) {
            for (int k = i - 1; k >= 0; k--) {
                double c = result[i + k * 4];
                for (int j = 0; j < 4; j++) {
                    result[j + k * 4] -= result[j + i * 4] * c;
                    E[j + k * 4] -= E[j + i * 4] * c;
                }
            }
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                (*x)[j + i * 4] = E[j + i * 4];
            }
        }
        return x;


    }


}

Matrix * inversion(Matrix *x, int N) {
    double temp;
    //TODO add proffs
    Matrix E, result;


    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[j + i * 4] = (*x)[j + i * 4];
            if (i == j) {
                E[j + i * 4] = 1.;
            }
            else {
                E[j + i * 4] = 0.;
            }
        }
    }


    for (int k = 0; k < N; k++)
    {
        temp = result[k+k*4];

        for (int j = 0; j < 4; j++)
        {
            result[j+k*4] /= temp;
            E[j+k*4] /= temp;
        }

        for (int i = k + 1; i < N; i++)
        {
            temp = result[k+i*4];

            for (int j = 0; j < N; j++)
            {
                result[j+i*4] -= result[j+k*4] * temp;
                E[j+i*4] -= E[j+k*4] * temp;
            }
        }
    }

    for (int k = N - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = result[k+i*4];

            for (int j = 0; j < N; j++)
            {
                result[j+i*4] -= result[j+k*4] * temp;
                E[j+i*4] -= E[j+k*4] * temp;
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            (*x)[j + i * 4] = E[j + i * 4];
        }
    }
    return x;

}

Matrix * identity(Matrix *x){

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {

            if (i == j) {
                (*x)[j + i * 4] = 1.;
            }
            else {
                (*x)[j + i * 4] = 0.;
            }
        }
    }
    return x;

}

Matrix * viewport(int x, int y, int width, int height, int depth, Matrix * viewPort){
    identity(viewPort);
    (*viewPort[3+0*4])=x+width/2.;
    (*viewPort[3+1*4])=y+height/2.;
    (*viewPort[3+2*4])=depth/2.;

    (*viewPort[0+0*4])=width/2.;
    (*viewPort[1+1*4])=height/2.;
    (*viewPort[2+2*4])=depth/2.f;

    return viewPort;
}






