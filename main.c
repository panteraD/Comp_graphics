#include <stdio.h>
#include <stdlib.h>
#include "model.h"
#include "geometry.h"
#include <omp.h>


void line(tgaImage *image,
          int x0, int y0,
          int x1, int y1,
          tgaColor color);


void triangle2(tgaImage *image, double *zbuffer, Vec3 *t0, Vec3 *t1, Vec3 *t2, Vec3 *uv0, Vec3 *uv1, Vec3 *uv2,
               double intensity, Model *cat);


const unsigned width = 1000;
const unsigned int height = 1000;
const unsigned int depth = 255;// base 255

//закраска фонга
int main(int argc, char **argv) {


    int rv = 0;


    if (argc < 2) {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        return -1;
    }
    tgaImage *image = tgaNewImage(width, height, RGB);


    Model *model = loadFromObj(argv[2]);
    loadDiffuseMap(model, argv[3]);

    Vec3 screen_coords[3];
    Vec3 uv_coords[3];
    Vec3 *uv[3];
    Vec3 *memory_coords[3];
    Vec3 *world_coords[3];
    //Vec3 *light_dir = getVertex(cat,rand()%1500,rand()%2);

    Vec3 left;
    Vec3 right;
    Vec3 n;

    /*
//    Matrix test = {3. , 2., 8., 5.,
//                    11., 135., 56., 84.,
//                    55., 9., 2., 43.,
//                    2., 66., 8., 5., };

//    Matrix test = {3., 2., 8., 5.,
//                   11., 15., 5., 8.,
//                   5., 9., 2., 3.,
//                   2., 6., 8., 5.,};
//    inversion(&test, 4);
     */


    double *zbuffer = malloc(sizeof(double) * width * height); //изменить тип на short

    for (int i = 0; i < width * height; i++) {
        zbuffer[i] = -3000000;
    }

    double c = 10;//20-30



    /*
     * PROGRAM VARIABLES:
             Vec3f light_dir = Vec3f(1,-1,1).normalize();
             Vec3f eye(1,1,3);
             Vec3f center(0,0,0);
     */

    Vec3 light_dir, eye, center;
    //initVec3(0., 0., -1., &light_dir);
    initVec3(1., -1., 1., &light_dir);
    normalize(&light_dir);

    initVec3(1., 1., 3., &eye);

    initVec3(0., 0., 0., &center);

    /*
    Matrix A = {1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.};
    Matrix B = {10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160.};
    Matrix C;
    matrixMultuplyMatrixMatrix(&A, &B, &C);
     */

    /*
     *      Matrix ModelView  = lookat(eye, center, Vec3f(0,1,0));
            Matrix Projection = Matrix::identity(4);
            Matrix ViewPort   = viewport(width/8, height/8, width*3/4, height*3/4);
            Projection[3][2] = -1.f/(eye-center).norm();

            std::cerr << ModelView << std::endl;
            std::cerr << Projection << std::endl;
            std::cerr << ViewPort << std::endl;
            Matrix z = (ViewPort*Projection*ModelView);
            std::cerr << z << std::endl;
     */


    Matrix ModelView;
    Vec3 up;
    initVec3(0., 1., 0., &up);
    lookAt(&eye, &center, &up, &ModelView);


    Matrix projection;
    identity(&projection);

    Vec3 temp2;
    vectorSub3D(&eye, &center, &temp2);

    double dd = norm(&temp2);
    projection[2 + 3 * 4] = -1. /dd;



    Matrix ViewPort;
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4, depth, &ViewPort);

    Vec4 tttt,rrrr;

    omp_set_num_threads(4);

    //#pragma omp parallel for
    for (int i = 0; i < model->nface; i++) {

        printf("face %d of %d\n", i + 1, model->nface);

        for (int j = 0; j < 3; j++) {
            memory_coords[j] = getVertex(model, i, j); //нужно ли 2 пременных? //TODO:remove this
            world_coords[j] = getVertex(model, i, j);

            // версия перспективного искажения

            transform3D4D(memory_coords[j], &tttt);
            matrixMultiplyMatrixVec4(&projection, &tttt, &rrrr);
            transform4D3D(&rrrr, &screen_coords[j]);




            //IN CIRCLE:screen_coords[j] =  Vec3f(ViewPort*Projection*ModelView*Matrix(v));

//            Vec4 mv, mv2;
//            transform3D4D(memory_coords[j], &mv); //[x,y,z] => [x,y,z,1];
//            Matrix matrix, matrix2;
//            matrixMultuplyMatrixMatrix(&ViewPort, &projection, &matrix);
//            matrixMultuplyMatrixMatrix(&matrix, &ModelView, &matrix2);
//            matrixMultiplyMatrixVec4(&matrix2, &mv, &mv2);
//            transform4D3D(&mv2, &screen_coords[j]);



            (screen_coords[j])[0] = ((screen_coords[j])[0] + 1.0) * (double) (width) / 2.;
            (screen_coords[j])[1] = ((screen_coords[j])[1] + 1.0) * (double) (height) / 2.;
            (screen_coords[j])[2] = ((screen_coords[j])[2] + 1.0) * (double) (depth) / 2.;

            printVec3(&screen_coords[j]);

            uv[j] = getDiffuseUV(model, i, j);
            uv_coords[j][0] = (*uv[j])[0];
            uv_coords[j][1] = (*uv[j])[1];
            uv_coords[j][2] = (*uv[j])[2];
        }
        vectorSub3D(world_coords[2], world_coords[0], &left);
        vectorSub3D(world_coords[1], world_coords[0], &right);
        vectorMult3D(&left, &right, &n);

        normalize(&n);

        double intensity = vectorMultScalar(&n, &light_dir);


        // if (intensity > 0.0f)    // убрать
        triangle2(image, zbuffer, &screen_coords[0], &screen_coords[1], &screen_coords[2], &uv_coords[0],
                  &uv_coords[1], &uv_coords[2],
                  intensity, model);


        char str[10];
        if (i < 0) {
            sprintf(str, "head%d.tga", i);

            tgaSaveToFile(image, str);
        }


    }


    tgaFlipVertically(image);


    if (-1 == tgaSaveToFile(image, argv[1])) {
        perror("tgaSateToFile");
        rv = -1;
    }

    tgaImage *zimage = tgaNewImage(width, height, RGB);
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            tgaSetPixel(zimage, i, j, tgaRGB(zbuffer[i + j * width], zbuffer[i + j * width], zbuffer[i + j * width]));
        }
    }
    tgaFlipVertically(zimage);
    tgaSaveToFile(zimage, "zim.tga");


    tgaFreeImage(image);

    free(zbuffer);
    printf("\ndone!\n");
    return rv;
}


//IT WORKS!
void triangle2(tgaImage *image, double *zbuffer, Vec3 *t0, Vec3 *t1, Vec3 *t2, Vec3 *uv0, Vec3 *uv1, Vec3 *uv2,
               double intensity, Model *cat) {
    if ((*t0)[1] == (*t1)[1] && (*t0)[1] == (*t2)[1]) return;
    if ((*t0)[1] > (*t1)[1]) {
        swapVec3(t0, t1);
        swapVec3(uv0, uv1);
    }
    if ((*t0)[1] > (*t2)[1]) {
        swapVec3(t0, t2);
        swapVec3(uv0, uv2);
    }
    if ((*t1)[1] > (*t2)[1]) {
        swapVec3(t1, t2);
        swapVec3(uv1, uv2);
    }


    double total_height = ((*t2)[1] - (*t0)[1]);

    Vec3 A;
    Vec3 B;
    Vec3 uvA;
    Vec3 uvB;

    int p0, p1, p2;
    for (double i = 0; i < total_height; i += 0.5) {
        int second_half = 0; //flag
        if (i > (*t1)[1] - (*t0)[1]) {
            second_half = 1;
        }
        else if ((*t1)[1] == (*t0)[1]) {
            second_half = 1;
        }
        double segment_height = second_half == 1 ? (*t2)[1] - (*t1)[1] : (*t1)[1] - (*t0)[1];
        double alpha = i / total_height;

        double beta = (i - (second_half == 1 ? (*t1)[1] - (*t0)[1] : 0)) / segment_height;


        (vectorAdd3D(t0, vectorMultSimple3D(vectorSub3D(t2, t0, &A), alpha, &A), &A));

        (second_half == 1) ? vectorAdd3D(t1, vectorMultSimple3D(vectorSub3D(t2, t1, &B), beta, &B), &B) :
        vectorAdd3D(t0, vectorMultSimple3D(vectorSub3D(t1, t0, &B), beta, &B), &B);

        vectorAdd3D(uv0, vectorMultSimple3D(vectorSub3D(uv2, uv0, &uvA), alpha, &uvA), &uvA);

        (second_half == 1) ? vectorAdd3D(uv1, vectorMultSimple3D(vectorSub3D(uv2, uv1, &uvB), beta, &uvB), &uvB) :
        vectorAdd3D(uv0, vectorMultSimple3D(vectorSub3D(uv1, uv0, &uvB), beta, &uvB), &uvB);


        if ((A)[0] > (B)[0]) {
            swapVec3(&A, &B);
            swapVec3(&uvA, &uvB);
        }


        for (double j = (A)[0]; j <= (B)[0]; j++) {

            double phi = ((B)[0] == (A)[0]) ? 1.0 : (j - A[0]) / (B[0] - A[0]);
            Vec3 P;
            Vec3 uP;
            Vec3 T1, T2, uT1, uT2;


            vectorSub3D(&B, &A, &T1);
            //Round(&T1);
            vectorMultSimple3D(&T1, phi, &T2);
            //Round(&T2);
            vectorAdd3D(&A, &T2, &P);
            //Round(&P);

            // vectorAdd3D(&uvA, vectorMultSimple3D(vectorSub3D(&uvB, &uvA, &uP), phi, &uP), &uP);

            vectorSub3D(&uvB, &uvA, &uT1);
            vectorMultSimple3D(&uT1, phi, &uT2);
            vectorAdd3D(&uvA, &uT2, &uP);


            p0 = P[0];
            p1 = P[1];
            p2 = P[2];


            int idx = (int) (p0 + 0.5) + (int) ((p1) * width + 0.5); //TODO: deal with it

            if (idx > (3800 * 3800 - 1) || idx < 0) {
                continue;
            }

            if (P[2] > zbuffer[idx]) {
                zbuffer[idx] = P[2];

                tgaColor color = getDiffuseColor(cat, &uP);


                tgaSetPixel(image, (int) (P[0]), (int) (P[1]),
                            tgaRGB(Red(color) * intensity, Green(color) * intensity, Blue(color) * intensity)); //2
                // tgaSetPixel(image, p0, p1, tgaRGB(255 * intensity, 255 * intensity, 255 * intensity)); //2

            }

            //сделлать единый тип double вместо float
            //проверить алгоритм  трегольника
            //попробовать изменить параметр depth
            //попробовать сделать swap для указателей

        }


    }

}






