#include <stdio.h>
#include <stdlib.h>
#include "model.h"
#include "geometry.h"
#include "tga.h"


void line(tgaImage *image,
          int x0, int y0,
          int x1, int y1,
          tgaColor color);

void triangle(tgaImage *image, double *zbuffer, Vec3 *t0, Vec3 *t1, Vec3 *t2, Vec3 *uv0, Vec3 *uv1, Vec3 *uv2,
              double intensity, Model *cat);

void triangle2(tgaImage *image, double *zbuffer, Vec3 *t0, Vec3 *t1, Vec3 *t2, Vec3 *uv0, Vec3 *uv1, Vec3 *uv2,
               double intensity, Model *cat);


const unsigned width = 3800;
const unsigned int height = 3800;
const unsigned int depth = 255;// base 255

int main(int argc, char **argv) {


    int rv = 0;


    if (argc < 2) {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        return -1;
    }
    tgaImage *image = tgaNewImage(width, height, RGB);
    //читаем файл с модели

    Model *cat = loadFromObj(argv[2]);
    loadDiffuseMap(cat, argv[3]);

    Vec3 screen_coords[3];
    Vec3 uv_coords[3];
    Vec3 *uv[3];
    Vec3 *memory_coords[3];
    Vec3 *world_coords[3];
    //Vec3 *light_dir = getVertex(cat,rand()%1500,rand()%2);
    Vec3 light_dir;
    light_dir[0] = 0.0;
    light_dir[1] = 0.0;
    light_dir[2] = -1.0;
    Vec3 left;
    Vec3 right;
    Vec3 n;
    printf("nfaces %d", cat->nface);

    double *zbuffer = malloc(sizeof(double) * width * height); //изменить тип на short

    for (int i = 0; i < width * height; i++) {
        zbuffer[i] = -3000000;
    }

    double c = 10;//20-30

    Matrix projection = {1., 0., 0., 0.,
                         0., 1., 0., 0.,
                         0., 0., 1., 0.,
                         0., 0., -1. / c, 1.};

    Vec4 t,r;
   // for(int e=0;e<100; e+=10) {


    for (int i = 0; i < cat->nface; i++) {
        if(i==1179 || i==1188 || i==1189 || i== 2396 || i==2386 ||i==2397){
            continue;
        }

        printf("face %d of %d\n", i + 1, cat->nface);

        for (int j = 0; j < 3; j++) {
            memory_coords[j] = getVertex(cat, i, j); //нужно ли 2 пременных? //TODO:remove this
            world_coords[j] = getVertex(cat, i, j);


            transform3D4D(memory_coords[j],&t);
            MatrixMuliply(&projection,&t,&r);
            transform4D3D(&r,&screen_coords[j]);

            (screen_coords[j])[0] = ((screen_coords[j])[0] + 1.0) * (double) (width) / 2.;
            (screen_coords[j])[1] = ((screen_coords[j])[1] + 1.0) * (double) (height) / 2.;
            (screen_coords[j])[2] = ((screen_coords[j])[2] + 1.0) * (double) (depth) / 2.;






            uv[j] = getDiffuseUV(cat, i, j);
            uv_coords[j][0] = (*uv[j])[0];
            uv_coords[j][1] = (*uv[j])[1];
            uv_coords[j][2] = (*uv[j])[2];
        }
        vectorSub3D(world_coords[2], world_coords[0], &left);
        vectorSub3D(world_coords[1], world_coords[0], &right);
        vectorMult3D(&left, &right, &n);

        normalize(&n);

        double intensity = vectorMultScalar(&n, &light_dir);


        if (intensity > 0.0f)    // убрать
            triangle2(image, zbuffer, &screen_coords[0], &screen_coords[1], &screen_coords[2], &uv_coords[0],
                      &uv_coords[1], &uv_coords[2],
                      intensity, cat);


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

            int idx = (int) (p0 + 0.5) + (int) ((p1) * width + 0.5);

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


//use valgrind
//draws triangle by coordinates
void triangle(tgaImage *image, double *zbuffer, Vec3 *t0, Vec3 *t1, Vec3 *t2, Vec3 *uv0, Vec3 *uv1, Vec3 *uv2,
              double intensity, Model *cat) {
    if ((*t0)[1] == (*t1)[1] && (*t0)[1] == (*t2)[1]) return;
    // sort the vertices, t0, t1, t2 lower-to-upper
    if ((*t0)[1] > (*t1)[1]) { swapVec3(t0, t1); }
    if ((*t0)[1] > (*t2)[1]) { swapVec3(t0, t2); }
    if ((*t1)[1] > (*t2)[1]) { swapVec3(t1, t2); }
    int total_height = (int) ((*t2)[1] - (*t0)[1]); ///???? привести?


    Vec3 A;
    Vec3 B;
    Vec3 uvA;
    Vec3 uvB;

    int p0, p1, p2;
    for (int i = 0; i < total_height; i++) {
        int second_half = 0;
        if (i > (*t1)[1] - (*t0)[1]) {
            second_half = 1;
        }
        else if ((*t1)[1] == (*t0)[1]) {
            second_half = 1;
        }
        int segment_height = second_half == 1 ? (*t2)[1] - (*t1)[1] : (*t1)[1] - (*t0)[1];
        double alpha = (double) i / (double) total_height;

        double beta = (i - (second_half == 1 ? (*t1)[1] - (*t0)[1] : 0)) / (double) segment_height;

        // Vec3i A = t0 + Vec3f(t2-t0)*alpha;
        //Vec3i B = second_half ? t1 + Vec3f(t2-t1)*beta : t0 + Vec3f(t1-t0)*beta;
        // if (A.x>B.x) std::swap(A, B);

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

        (A)[0] = (int) ((A)[0] + 0.5);
        (A)[1] = (int) ((A)[1] + 0.5);
        (A)[2] = (int) ((A)[2] + 0.5);
        (B)[0] = (int) ((B)[0] + 0.5);
        (B)[1] = (int) ((B)[1] + 0.5);
        (B)[2] = (int) ((B)[2] + 0.5);


        for (int j = (A)[0]; j <= (B)[0]; j++) {

            double phi = ((B)[0] == (A)[0]) ? 1.0 : (j - A[0]) / (B[0] - A[0]);
            Vec3 P;
            Vec3 uP;
            Vec3 T1, T2;


            // vectorAdd(&A,vectorMultSimple3D(vectorSub3D(&B,&A,&P),phi,&P),&P);

            // Vec3i P = Vec3f(A) + Vec3f(B-A)*phi;


            vectorSub3D(&B, &A, &T1);
            //Round(&T1);
            vectorMultSimple3D(&T1, phi, &T2);
            //Round(&T2);
            vectorAdd3D(&A, &T2, &P);
            //Round(&P);



            vectorAdd3D(&uvA, vectorMultSimple3D(vectorSub3D(&uvB, &uvA, &uP), phi, &uP), &uP);


            p0 = P[0];
            p1 = P[1];
            p2 = P[2];




            /*

           //int idx = j + (p1)*width;  //3

        //   int idx = j + (*t0)[1]*width; //spiderman

            int idx = j + ((*t0)[1])*width; //1

            if( p2 > zbuffer[idx] ) {
                zbuffer[idx] = p2;
                tgaSetPixel(image,j,(*t0)[1]+i, color); //1 хуже

            }
            */




            /////////////////////////////
            //2
            //int idx = p0 + (p1)*width;  //2
            int idx = (int) (p0 + 0.5) + (int) ((p1) * width + 0.5);
            if(idx>3800*3800-1){
                continue;
            }

            if (P[2] > zbuffer[idx]) {
                zbuffer[idx] = P[2];
                tgaColor color = getDiffuseColor(cat, &uP);

                 tgaSetPixel(image,p0,p1, tgaRGB( Red(color)*intensity, Green(color)*intensity, Blue(color)*intensity)); //2
                //tgaSetPixel(image, p0, p1, tgaRGB(255 * intensity, 255 * intensity, 255 * intensity)); //2

            }



            //сделлать единый тип double вместо float
            //проверить алгоритм разторизапции трегольника
            //попробовать изменить параметр depth


            /*

           tgaColor color = getDiffuseColor(cat,&uP);

            tgaSetPixel(image,j,(*t0)[1]+i, tgaRGB( Red(color)*intensity, Green(color)*intensity, Blue(color)*intensity));
           //tgaSetPixel(image,j,(*t0)[1]+i, color); //simple gives pic with no artifacts without z buff

            */
        }


    }
    printf("  %d", p2);
}





