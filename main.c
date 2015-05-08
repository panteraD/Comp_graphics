#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "model.h"
#include "tga.h"


void swap(int *a, int *b);
int abs(int a);
void line (tgaImage *image, 
           int x0, int y0,
           int x1, int y1,
           tgaColor color);
void triangle(tgaImage *image, double *zbuffer, Vec3 *t0, Vec3 *t1, Vec3 *t2,Vec3 *uv0,Vec3 *uv1,Vec3 *uv2, double intensity, Model *cat);

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
 
const unsigned width = 1600;
const unsigned int height = 1600;
const unsigned int depth = 800;// base 255

int main(int argc, char **argv)
{
    int rv = 0;
    
    
    if (argc < 2) {
        fprintf(stderr, "Usage: %s outfile\n", argv[0]);
        return -1;
    }
    tgaImage * image = tgaNewImage(width, height, RGB);
    //читаем файл с модели
   
	Model *cat = loadFromObj(argv[2]);
    loadDiffuseMap(cat,argv[3]);
  
    Vec3 screen_coords[3];
    Vec3 *memory_coords[3];
    Vec3 *world_coords[3];
    //Vec3 *light_dir = getVertex(cat,rand()%1500,rand()%2);
    Vec3 light_dir; 
    light_dir[0]=0.0;
    light_dir[1]=0.0;
    light_dir[2]=-1.0;
    Vec3 left;
    Vec3 right;
    Vec3 n;
    printf("nfaces %d",cat->nface);
    
    double *zbuffer = malloc (sizeof (double) * width*height); //изменить тип на short
   
    for(int i = 0; i < width*height; i++) {
        zbuffer[i]  = -3000000 ; 
    }
   
	for (int i = 0; i < cat->nface; i++) {
        
        printf("face %d of %d\n",i+1, cat->nface);	
		
		for (int j=0; j<3; j++) {
		  memory_coords[j] = getVertex(cat,i,j); //нужно ли 2 пременных?
          world_coords[j] = getVertex(cat,i,j);
		  (screen_coords[j])[0]=((*memory_coords[j])[0]+1.0)*width/2;
		  (screen_coords[j])[1]=((*memory_coords[j])[1]+1.0)*height/2;
          (screen_coords[j])[2]=((*memory_coords[j])[2]+1.0)*depth/2;
		}
		vectorSub3D(world_coords[2], world_coords[0], &left);
        vectorSub3D(world_coords[1], world_coords[0], &right);
        vectorMult3D(&left,&right,&n);
      
        normalize(&n);
      
        double intensity = vectorMultScalar(&n,&light_dir);
        
        //printf("%f\n",intensity); //laggy comment
        
       
        
        Vec3 *uv[3];
        for(int k = 0; k< 3;k++){
            uv[k]=getDiffuseUV(cat,i,k);
            
        }
        
        
		//if (intensity > 0.0f)    // убрать
            triangle(image, zbuffer ,&screen_coords[0],&screen_coords[1],&screen_coords[2], uv[0], uv[1], uv[2],intensity,cat);
                     
        
        
        
       
        char str[10];
        if(i<0){
                sprintf(str,"head%d.tga",i);
                
                tgaSaveToFile(image, str);
        }
       
			
    }
    
   
    
    
		
	tgaFlipVertically(image);
    

	if (-1 == tgaSaveToFile(image, argv[1])) {
	  perror("tgaSateToFile");
	  rv = -1;
	}
	
	tgaImage *zimage = tgaNewImage(width, height, RGB);
    for(int i = 0; i< width; i++) {
        for( int j = 0; j < height ; j++) {
             tgaSetPixel(zimage, i,j, tgaRGB(zbuffer[i+j*width], zbuffer[i+j*width] ,zbuffer[i+j*width] ));
        }
    }
    tgaFlipVertically(zimage);
	tgaSaveToFile(zimage, "zim.tga");
	
      
    tgaFreeImage(image);
  
    free(zbuffer);
	printf("\ndone!\n");
    return rv;
}

//use valgrind
//draws triangle by coordinates
void triangle(tgaImage *image,  double  *zbuffer, Vec3 *t0, Vec3 *t1, Vec3 *t2,Vec3 *uv0,Vec3 *uv1,Vec3 *uv2, double intensity, Model *cat) {
    if ((*t0)[1]==(*t1)[1] && (*t0)[1]==(*t2)[1]) return;
    // sort the vertices, t0, t1, t2 lower-to-upper 
    if ((*t0)[1]>(*t1)[1]) { swapVec3(t0, t1); }
    if ((*t0)[1]>(*t2)[1]) { swapVec3(t0, t2); }
    if ((*t1)[1]>(*t2)[1]) { swapVec3(t1, t2); }
    int total_height = (int)((*t2)[1]-(*t0)[1]); ///???? привести?
	

    Vec3 A;
    Vec3 B;
    Vec3 uvA;
    Vec3 uvB;
	
    int p0,p1,p2;
    for (int i=0; i < total_height; i++) {		
        int second_half = 0;
		if (i>(*t1)[1]-(*t0)[1]) {
		  second_half=1;
		}
		else if ( (*t1)[1]==(*t0)[1]) {
		  second_half =1;
		}
        int segment_height = second_half==1 ? (*t2)[1]-(*t1)[1] : (*t1)[1]-(*t0)[1];
        double alpha = (double)i/(double)total_height;
       
        double beta  = (double)(i-(second_half == 1 ? (*t1)[1]-(*t0)[1] : 0))/(double)segment_height; 
        
       // Vec3i A = t0 + Vec3f(t2-t0)*alpha;
        //Vec3i B = second_half ? t1 + Vec3f(t2-t1)*beta : t0 + Vec3f(t1-t0)*beta;
       // if (A.x>B.x) std::swap(A, B);
      
		(vectorAdd3D(t0, vectorMultSimple3D(vectorSub3D(t2, t0,&A),alpha,&A), &A)); 
		
		( second_half == 1)? vectorAdd3D(t1,vectorMultSimple3D(vectorSub3D(t2,t1,&B),beta,&B),&B):
							 vectorAdd3D(t0,vectorMultSimple3D(vectorSub3D(t1,t0,&B),beta,&B),&B); 
		
       vectorAdd3D(uv0,vectorMultSimple3D(vectorSub3D(uv2,uv0,&uvA),alpha,&uvA),&uvA); 
       ( second_half == 1)? vectorAdd3D(uv1,vectorMultSimple3D(vectorSub3D(uv2,uv1,&uvB),beta,&uvB),&uvB):
                            vectorAdd3D(uv0,vectorMultSimple3D(vectorSub3D(uv1,uv0,&uvB),beta,&uvB),&uvB); 

        
        
       
       
        if ((A)[0]>(B)[0]) { 
            swapVec3(&A, &B) ;
            swapVec3(&uvA,&uvB);
        }
        
        (A)[0]=(int)((A)[0]+0.5);
        (A)[1]=(int)((A)[1]+0.5);
        (A)[2]=(int)((A)[2]+0.5);
        (B)[0]=(int)((B)[0]+0.5);
        (B)[1]=(int)((B)[1]+0.5);
        (B)[2]=(int)((B)[2]+0.5);
        
        
        
        
          
        for (int j=(A)[0]; j<=(B)[0]; j++) {
           
            double phi = ((B)[0]==(A)[0]) ? 1.0 : (double)(j-A[0])/(double)(B[0]-A[0]);
            Vec3 P; 
            Vec3 uP;
            Vec3 T1,T2;
      
            
           // vectorAdd(&A,vectorMultSimple3D(vectorSub3D(&B,&A,&P),phi,&P),&P);
            
           // Vec3i P = Vec3f(A) + Vec3f(B-A)*phi;
                        
            
            vectorSub3D(&B,&A,&T1);
            //Round(&T1);
            vectorMultSimple3D(&T1,phi,&T2);
            //Round(&T2);
            vectorAdd3D(&A,&T2,&P);
            //Round(&P);
            
           
             
           vectorAdd3D(&uvA,vectorMultSimple3D(vectorSub3D(&uvB,&uvA,&uP),phi,&uP),&uP);  
            
           
            
            p0=P[0];
            p1=P[1];         
            p2=P[2]; 
            
           
            
            
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
              int idx = (int)(p0+0.5) + (int)( (p1)*width +0.5);
            
            if( P[2] > zbuffer[idx] ) {
                zbuffer[idx] = P[2];
                tgaColor color = getDiffuseColor(cat,&uP);
                
               // tgaSetPixel(image,p0,p1, tgaRGB( Red(color)*intensity, Green(color)*intensity, Blue(color)*intensity)); //2
                tgaSetPixel(image,p0,p1, tgaRGB( 255*intensity, 255*intensity, 255*intensity)); //2
                           
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
    printf("  %d",p2);
}



Vec3 * vectorAdd3D(Vec3 *t0, Vec3 *t1, Vec3 *res) {
  
  (*res)[0]=(*t0)[0]+(*t1)[0];
  (*res)[1]=(*t0)[1]+(*t1)[1];
  (*res)[2]=(*t0)[2]+(*t1)[2];
  return res;
}



Vec3 * vectorSub3D(Vec3 *t0, Vec3 *t1, Vec3 *res) {
 // printf("sub %f \n",(*t0)[0]);//,(*t1)[0],(*res)[0]); //
  (*res)[0]=(*t0)[0]-(*t1)[0];
  (*res)[1]=(*t0)[1]-(*t1)[1];
  (*res)[2]=(*t0)[2]-(*t1)[2];
  return res;
}




Vec3 * vectorMultSimple3D(Vec3 *t0, double a, Vec3 *res) {

  (*res)[0]=(*t0)[0]*a;
  (*res)[1]=(*t0)[1]*a;
  (*res)[2]=(*t0)[2]*a;
  return res;
}



Vec3 * vectorMult3D (Vec3 *t0, Vec3 *t1, Vec3 *res) {
    
    (*res)[0] = (*t0)[1]*(*t1)[2]-(*t0)[2]*(*t1)[1];   //aybz-azby
    (*res)[1] = (*t0)[2]*(*t1)[0]-(*t0)[0]*(*t1)[2];   //azbx-axbz
    (*res)[2] = (*t0)[0]*(*t1)[1]-(*t0)[1]*(*t1)[0];   //axby-aybx
    return res;
}

double vectorMultScalar (Vec3 *t0, Vec3 *t1) {
  
    return ((*t0)[0]*(*t1)[0]+(*t0)[1]*(*t1)[1]+(*t0)[2]*(*t1)[2]);
}



double  norm(Vec3 *t0){
    return (double)sqrt((*t0)[0]*(*t0)[0]+(*t0)[1]*(*t0)[1]+(*t0)[2]*(*t0)[2]);
    
}


Vec3 * normalize(Vec3 * t0) {
       
        double n = (double)(1.0f/norm(t0));
        vectorMultSimple3D(t0,  n,t0 );
        return t0;
}


void Round(Vec3 *t0) {
     (*t0)[0]=(int)((*t0)[0]+0.5); 
     (*t0)[1]=(int)((*t0)[1]+0.5);     
     (*t0)[2]=(int)((*t0)[2]+0.5); 
}




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

void swapVec3(Vec3 *t0, Vec3 *t1){
	double temp;
	for(int i=0;i<3;i++){
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

