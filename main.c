#include <stdio.h>
#include <stdlib.h>
#include "model.h"
#include "geometry.h"
#include <omp.h>


void triangle(tgaImage *image, double *zbuffer, Vec3 *t0, Vec3 *t1, Vec3 *t2, Vec3 *uv0, Vec3 *uv1, Vec3 *uv2,
              Vec3 *inten, Vec3 *light_dir, Model *cat);


const unsigned width = 3800;
const unsigned int height = 3800;
const unsigned int depth = 255;


int main(int argc, char **argv) {


    int rv = 0;


    if (argc < 2) {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        return -1;
    }
    tgaImage *image = tgaNewImage(width, height, RGB);


    Model *model = loadFromObj(argv[2]);
    loadDiffuseMap(model, argv[3]);
    loadNormalMap(model, argv[4]);


    Vec3 screen_coords[3];
    Vec3 uv_coords[3];
    Vec3 *uv[3];
    Vec3 *memory_coords[3];
    Vec3 *intensive[3];

    /*
    for(unsigned int i=0; i< width;i++){
        for(unsigned int j = 0; j< height; j++){
            tgaSetPixel(image, i, j, tgaRGB(128 , 128 , 128 ));
        }
    }
     */


    double *zbuffer = malloc(sizeof(double) * width * height); //изменить тип на short

    for (int i = 0; i < width * height; i++) {
        zbuffer[i] = -3000000;
    }


    Vec3 light_dir, eye, center;

    Vec3 light1, light2;
    initVec3(0., 0., 1, &light1);
    initVec3(0., -0., 0., &light2);
    vectorAdd3D(&light1, &light2, &light_dir);


    // initVec3(0., 0., 1., &light_dir);
    //  initVec3(1., -1., 1., &light_dir);
    normalize(&light_dir);

    initVec3(1., 1., 3., &eye);
    initVec3(0., 0., 0., &center);


    Matrix ModelView;
    Vec3 up;
    initVec3(0., 1., 0., &up);
    lookAt(&eye, &center, &up, &ModelView);


    Matrix projection;
    identity(&projection);

    Vec3 C;
    vectorSub3D(&eye, &center, &C);

    double dd = norm(&C);
    projection[2 + 3 * 4] = -1. / dd;


    Matrix ViewPort;
    viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4, depth, &ViewPort);

    Matrix M, temp; // camera transformation
    matrixMultuplyMatrixMatrix(&ViewPort, &projection, &temp);
    matrixMultuplyMatrixMatrix(&temp, &ModelView, &M);






    // #pragma omp parallel for
    // #pragma omp for
#pragma omp  for
    for (int i = 0; i < model->nface; i++) {

        printf("face %d of %d\n", i + 1, model->nface);
        Vec3 inten, sv;


        for (int j = 0; j < 3; j++) {
            memory_coords[j] = getVertex(model, i, j);


            /*
             * Vec3f n = (world_coords[1] - world_coords[0]) ^ (world_coords[2] - world_coords[0]); // поменять местами множители
                n.normalize();
                освященность?
             */

            /*
            (screen_coords[j])[0] = ((*memory_coords[j])[0] + 1.0) * width / 2;
            (screen_coords[j])[1] = ((*memory_coords[j])[1] + 1.0) * height / 2;
            (screen_coords[j])[2] = ((*memory_coords[j])[2] + 1.0) * depth / 2;
             */


            //преобразования для поворта камеры
            Vec4 mv, mv2;
            transform3D4D(memory_coords[j], &mv); //[x,y,z] => [x,y,z,1];
            matrixMultiplyMatrixVec4(&M, &mv, &mv2);
            transform4D3D(&mv2, &screen_coords[j]);

            //printVec3(&screen_coords[j]);

            uv[j] = getDiffuseUV(model, i, j);
            uv_coords[j][0] = (*uv[j])[0];
            uv_coords[j][1] = (*uv[j])[1];
            uv_coords[j][2] = (*uv[j])[2];


            /*
            //преобразовния нормалей для тонировки Гуро
            //normals transforamtions for Gourand
            intensive[j]=getNorm(model, i,j);
            sv[0]=-(*intensive)[j][0];
            sv[1]=-(*intensive)[j][1];
            sv[2]=-(*intensive)[j][2];


            Vec4 C, temp3;
            transform3D4D(&sv, &temp);
            temp[3]=0;
            Matrix M1;
            clone(&M,&M1);
            transponse(&M1);
            inversion(&M1);
            matrixMultiplyMatrixVec4(&M1,&temp,&temp3);

            transform4D3D(&temp3, &sv);
            normalize(&sv);


            inten[j]=vectorMultScalar(&sv,&light_dir);
             */

        }


        triangle(image, zbuffer, &screen_coords[0], &screen_coords[1], &screen_coords[2],
                 &uv_coords[0], &uv_coords[1], &uv_coords[2],
                 &inten, &light_dir, model);


        char str[10];
        if (i < 0) {
            sprintf(str, "head%d.tga", i);

            tgaSaveToFile(image, str);
        }


    }
#pragma omp barrier


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


void triangle(tgaImage *image, double *zbuffer, Vec3 *t0, Vec3 *t1, Vec3 *t2,
              Vec3 *uv0, Vec3 *uv1, Vec3 *uv2,
              Vec3 *intnen, Vec3 *light_dir, Model *cat) {
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

    int p0, p1;
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

        double ityA = (*intnen)[0] + ((*intnen)[2] - (*intnen)[0]) * alpha;
        double ityB = second_half ?
                      (*intnen)[1] + ((*intnen)[2] - (*intnen)[1]) * beta :
                      (*intnen)[1] + ((*intnen)[1] - (*intnen)[2]) * beta;


        if ((A)[0] > (B)[0]) {
            swapVec3(&A, &B);
            swapVec3(&uvA, &uvB);
            swapD(&ityA, &ityB);
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

            //интреполяция текстурных коррдинтат
            vectorSub3D(&uvB, &uvA, &uT1);
            vectorMultSimple3D(&uT1, phi, &uT2);
            vectorAdd3D(&uvA, &uT2, &uP);

            //гуро
            double ityP = ityA + (ityB - ityA) * phi;


            p0 = P[0];
            p1 = P[1];


            int idx = (int) (p0 + 0.5) + (int) ((p1) * width + 0.5); //TODO: deal with it

            if (idx > (width * height - 1) || idx < 0) {
                continue;
            }

            if (P[2] > zbuffer[idx]) {
                zbuffer[idx] = P[2];

                tgaColor color;
#pragma omp critical
                {
                    color = getDiffuseColor(cat, &uP);
                }

                Vec3 n;

#pragma omp critical
                {
                    getNormal(cat, &n, &uP);
                }
                normalize(&n);
                n[0] *= -1;
                n[1] *= -1;
                n[2] *= -1;
                double ityFhong = vectorMultScalar(&n, light_dir);  //?
                if (ityFhong < 0.)
                    ityFhong = 0.;


                //tgaSetPixel(image, p0, p1, tgaRGB(255 * ityP, 255 * ityP, 255 * ityP)); //Gourad

                // tgaSetPixel(image, p0, p1, tgaRGB(255 * ityFhong, 255 * ityFhong, 255 * ityFhong));//Phong

                double lightHack = 5.;

                tgaSetPixel(image, p0, p1,
                            tgaRGB(Red(color) * ityFhong + lightHack, Green(color) * ityFhong + lightHack,
                                   Blue(color) * ityFhong + lightHack)); //Phong text , +0.5 more light hack

                /*tgaSetPixel(image, p0, p1, tgaRGB(Red(color) * ityFhong, Green(color) * ityFhong,
                                                  Blue(color) * ityFhong)); //Phong text , les light
                                                  */


            }

        }

    }

}






