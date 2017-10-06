#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>

#define TRUE (0==0)
#define FALSE (0==1)

#define MAX_BYTE 255
#define MAX_DEPTH 3
#define N_SPHERES 2
#define N_LIGHTS 2

#define AREA 500
#define SPHERE_RAY 10

//floating point init data
#define AMBIENT 0.35f
#define LIGHT_INT 3.5f

//scale texture
#define TEXTURE_SCALE 0.5f

#define WID 1280
#define HEI 800

#define PI ((float) (355/113))
#define FOV         ((float) (PI/4))        /* field of view in rads (pi/4) */
#define HALF_FOV        (FOV * 0.5)
#define NRAN    1024
#define MASK    (NRAN - 1)
#define DIVX 4
#define DIVY 2
#define RAY_MAG 1000



//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

//STRUCTURES


struct c_color
{
    float r;
    float g;
    float b;
    float alpha;

}typedef color;

struct p_point
{
    float x;
    float y;
    float z;
}typedef point;

struct r_ray
{
    point o;
    point d;
}typedef ray;

struct v_viewport
{
    int width;
    int height;
}typedef viewport;

struct c_camera
{
    point eye;
    point lookat;
    point up,u,v,w;
    float distance;
    viewport view;
}typedef camera;

struct s_sphere
{
    point center;
    point maxBB;
    point minBB;
    color c;
    float cr;
    float kd; //[0,1] 
    float kf; //[0,1]
    float r;

    int e;
    float ks;
}typedef sphere;

struct o_open_cylinder
{
    float bottom;
    float top;
    float r;

    point maxBB;
    point minBB;
    color c;
    float cr;
    float kd; //[0,1] 
    float kf; //[0,1]

    int e;
    float ks;
}typedef open_cylinder;

struct r_rectangle
{
    point p0;
    point a;
    point b;

    color c;
    float cr;
    float kd; //[0,1] 
    float kf; //[0,1]

    int e;
    float ks;
}typedef rectangle;

struct p_plane
{
    point normal;
    color c;
    float d;

    float cr;
    float kd; //[0,1] 
    float kf; //[0,1]
    int e;
    float ks;

    point m_uAxis;
    point m_vAxis;
}typedef plane;

struct p_light
{
    float ls; //light intensity
    point l; //light coordinates
    color c; //light color
}typedef ponctual_light;

//struct g_grid
//{
//  point p0;
//  point p1;
//  int nx, ny, nz;
//  int total;
//}typedef grid;

enum object {
    NONE,
    SPHERE,
    PLANE,
    OPEN_CYLINDER,
    RECTANGLE,
    CONE,
    BOX
};

typedef unsigned char uchar;

struct varFromLoop{
    int samples;
    int s;
    float rcp_samples;
    uchar* image;
    camera c;
} typedef vFLoop;

struct retangulo {
    int inicio_x;
    int inicio_y;
    int fim_x;
    int fim_y;
} typedef retangulo;

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

// FUNCTIONS 
ray get_primary_ray(camera *cam, int x, int y, int sample);
point get_sample_pos(int x, int y, int sample);
point jitter(int x, int y, int s);
void generateScene();

void addVec(point *a, point *b);
void subVec(point *a, point *b);
void scalarMulVec(point *a,float sc);

//int inside_grid(grid *gr, point * p);
void initImage(camera *c,uchar * frame);
//void setupGrid(grid *g,point max, point min, int nx, int ny, int nz);
ray * generatePrimaryRays(camera *c);
float dist(point *o , point *d);
void setupCamera(camera *c);
void initCamera(camera *c, point eye, point lookat, int width, int height);
void printPrimaryRays(ray *r, int size);
void normalize(point *a);
void crossProduct(point *a, point *b, point *r);
float dotProduct(point *a, point *b);
float distancia(point a, point b);
float clamp(float x, float min, float max);
float m_max(float a,float b);
float m_min(float a,float b);
//void intersectBoundingBox(grid * gr, ray * raio, point *t_min, point *t_max, float *t0, float *t1);
int save_bmp(char * file_name, camera *c, uchar *frame);

int intersectSphere(sphere * s, ray * r, float *t);
int intersectPlane(plane *p, ray *r, float *t);
int intersectOpenCylinder(open_cylinder *cylinder, ray *r, float *t_res);
int intersectRectangle(rectangle *rec, ray *r, float *t);

enum object intersect3Dscene(int *idx,ray *r, float *t);

void gamut(float *r, float *g, float *b);
float checkerTexture(float u, float v);

//generate random spheres into the global 'data' array and lights into 'lights' array
void generateRandomSpheres();
void generateRandomLightSources();

//conversion from Float to Uchar
uchar floatToIntColor(float c);

//actual ray tracing CORE functions

color trace(camera cam, ray *raio, int iter);
color shade(camera cam, point *incid , enum object obj, int index, point *p, int iter);



// CONSTANTS AND GLOBAL VARIABLES

const float epsilon = 0.111f;

sphere data[N_SPHERES]; //data array with spheres
ponctual_light lights[N_LIGHTS]; //data array with lights

plane ground = {{0.0f,1.0f,0.0f},{4.4f}};

open_cylinder cylinder;
rectangle rect;

color black = {0.0f,0.0f,0.0f,0.0f}; //default black color
color white = {1.0f,1.0f,1.0f,1.0f}; //default white color
color red = {1.0f,0.0f,0.0f,0.0f}; //default white color

//double aspect = 1.333333f;
double aspect = WID/HEI;

point urand[NRAN];
int irand[NRAN];


#define DIV 20

void raytracerLoop(vFLoop *vl, int inicio_x, int inicio_y, int fim_x, int fim_y){
    int i,j,s;
    for(i = inicio_x ; i < fim_x ; i++)
    {
        for(j = inicio_y ; j < fim_y ; j++)
        {
            float r, g, b;
            r = g = b = 0.0;
            
            for(s = 0; s < vl->samples; s++) {
                ray rr = get_primary_ray(&(vl->c), i, j, s);    
                color col = trace(vl->c,&rr,0);
                r += col.r;
                g += col.g;
                b += col.b;
            }

            r = r * vl->rcp_samples;
            g = g * vl->rcp_samples;
            b = b * vl->rcp_samples;

            //ray rr = get_primary_ray(&c, i, j, samples); 
            //color clr = trace(c,&rr,0);

            //red green blue color components
            vl->image[ 3* (i * vl->c.view.height + j) + 0] = floatToIntColor(r);
            vl->image[ 3* (i * vl->c.view.height + j) + 1] = floatToIntColor(g);
            vl->image[ 3* (i * vl->c.view.height + j) + 2] = floatToIntColor(b);            
        }
    }
    printf("PASSO AQ DENTRO\n");
}

void comecaraytracerloopcoordenadas(vFLoop *vl, int num_divisoes_x, int num_divisoes_y){
    int width = vl->c.view.width;
    int height = vl->c.view.height;

    int divisao_horizontal = width/num_divisoes_x;
    int divisao_vertical = height/num_divisoes_y;

    int num_retangulos = num_divisoes_x * num_divisoes_y;
    retangulo *ret = (retangulo *) malloc(num_retangulos * sizeof(retangulo));

    for(int i = 0; i < num_divisoes_x; i++){
        for(int j = 0; j < num_divisoes_y; j++){
            ret[i * num_divisoes_y + j].inicio_x = i * divisao_horizontal;
            ret[i * num_divisoes_y + j].inicio_y = j * divisao_vertical;

            if(i == num_divisoes_x - 1)
                ret[i * num_divisoes_y + j].fim_x = width;
            else
                ret[i * num_divisoes_y + j].fim_x = ((i + 1) * divisao_horizontal);
            if(j == num_divisoes_y - 1)
                ret[i * num_divisoes_y + j].fim_y = height;
            else
                ret[i * num_divisoes_y + j].fim_y = ((j + 1) * divisao_vertical);
        } 
    }

    for(int i = 0; i < num_retangulos; i++)
        raytracerLoop(vl, ret[i].inicio_x, ret[i].inicio_y, ret[i].fim_x, ret[i].fim_y);
}

void comecaraytracerloop(vFLoop *vl, int num_divisoes){
    int diferenca = num_divisoes + 1;
    int num_divisoes_x = num_divisoes_x;
    int num_divisoes_y = 1;
    for(int i = 1; i <= num_divisoes; i++){
        if(num_divisoes % i == 0){
            int diferenca_atual = i - num_divisoes/i;
            if(diferenca_atual < diferenca){
                num_divisoes_x = i;
                num_divisoes_y = num_divisoes/i;
            }
        }
    }
    comecaraytracerloopcoordenadas(vl, num_divisoes_x, num_divisoes_y);
}

int main(int argc, char ** argv)
{
    int i,j;
    uchar *image;
    camera c;
    point eye;
    point lookat;
    int samples; 
    int s;    
    float rcp_samples;// = 1.0 / (float)samples;
    //char fname[20];
    //ray * rays;
    //color cor;

    //srand ( time(NULL) );

    //---init virtual camera---
    //point eye = {10.0f,400.0f,1000.0f};
    //point eye = {0.0f,2.0f,-20.0f};
    eye.x = 0.0f;
    eye.y = 2.0f;
    eye.z = -20.0f;

    //point lookat = {0.5f,0.0f,0.0f};
    lookat.x = 0.5f;
    lookat.y = 0.0f;
    lookat.z = 0.0f;

    initCamera(&c,eye,lookat,WID,HEI);
    setupCamera(&c);

    //---malloc the image frame---
    image = (uchar *) malloc(c.view.width * c.view.height * 3 * sizeof(uchar));
    if(image == NULL)
    {
        fprintf(stderr,"Error. Cannot malloc image frame.\n");
        return 0;
    }

    //---just init the image frame with some data---
    initImage(&c,image);

    //---insert random N_SPHERES into the 'data' array
    //generateRandomSpheres();
    generateScene();

    //---insert random N_LIGHTS into the 'lights' array
    generateRandomLightSources();

    //---create a 1D array with primary rays coordinates
    //rays = generatePrimaryRays(&c);

    for(i=0; i<NRAN; i++) urand[i].x = (double)rand() / RAND_MAX - 0.5;
    for(i=0; i<NRAN; i++) urand[i].y = (double)rand() / RAND_MAX - 0.5;
    for(i=0; i<NRAN; i++) irand[i] = (int)(NRAN * ((double)rand() / RAND_MAX));

    //-------------------------------------------ray tracing loop-----------------------------------------------------------
    vFLoop vl;

    vl.samples = 8;
    vl.s = 0;
    vl.rcp_samples = 1.0 / (float)vl.samples;

    vl.image = image;
    vl.c = c;

    //Loop aqui
    if(argc == 2){
        comecaraytracerloop(&vl, atoi(argv[1]));
    } else {
        comecaraytracerloop(&vl, DIV);
    }

    // int tamPieceLinear = c.view.width*c.view.height/DIV;
    // uchar* imagePiece = (uchar *) malloc(tamPieceLinear*sizeof(uchar));
    // uchar * aux = imagePiece;

    // int fimPiece = 0;
    // int inicioPiece;
    // int cont;

    // for(int y = 0; y < DIV; y++)
    // {
    //     int inicioPiece = fimPiece;
    //     int fimPiece = inicioPiece + tamPieceLinear;
        
    //     *(imagePiece) = raytracerLoop(vl,y*c.view.height/DIV);
    //     cont = 0;
    //     imagePiece = aux;
    //     for(int i=inicioPiece; i<fimPiece;i++,cont++){
    //         *(image + i) = *(imagePiece + cont);
    //     }

    //     imagePiece = aux;
    // }
    

    //*image = raytracerLoop(vl,0,c.view.height);
   
    //printPrimaryRays(rays,c.view.width*c.view.height); //for testing only

    if(save_bmp("output_rt.bmp",&c,vl.image) != 0)
    {
        fprintf(stderr,"Cannot write image 'output.bmp'.\n");
        return 0;
    }

    //---freeing data---
    //free(rays);
    free(image);

    //---exit---
    return 0;
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

ray get_primary_ray(camera *cam, int x, int y, int sample) {
    ray ray;
    float m[3][3];
    point i, j = {0, 1, 0}, k, dir, orig, foo;

    k.x = cam->lookat.x - cam->eye.x;
    k.y = cam->lookat.y - cam->eye.y;
    k.z = cam->lookat.z - cam->eye.z;

    normalize(&k);  
    //NORMALIZE(k);

    crossProduct(&j,&k,&i);
    crossProduct(&k,&i,&j);

    //i = cross_product(j, k);
    //j = cross_product(k, i);
    m[0][0] = i.x; m[0][1] = j.x; m[0][2] = k.x;
    m[1][0] = i.y; m[1][1] = j.y; m[1][2] = k.y;
    m[2][0] = i.z; m[2][1] = j.z; m[2][2] = k.z;

    ray.o.x = ray.o.y = ray.o.z = 0.0;
    ray.d = get_sample_pos(x, y, sample);
    ray.d.z = 1.0 / HALF_FOV;
    ray.d.x *= RAY_MAG;
    ray.d.y *= RAY_MAG;
    ray.d.z *= RAY_MAG;

    dir.x = ray.d.x + ray.o.x;
    dir.y = ray.d.y + ray.o.y;
    dir.z = ray.d.z + ray.o.z;
    foo.x = dir.x * m[0][0] + dir.y * m[0][1] + dir.z * m[0][2];
    foo.y = dir.x * m[1][0] + dir.y * m[1][1] + dir.z * m[1][2];
    foo.z = dir.x * m[2][0] + dir.y * m[2][1] + dir.z * m[2][2];

    orig.x = ray.o.x * m[0][0] + ray.o.y * m[0][1] + ray.o.z * m[0][2] + cam->eye.x;
    orig.y = ray.o.x * m[1][0] + ray.o.y * m[1][1] + ray.o.z * m[1][2] + cam->eye.y;
    orig.z = ray.o.x * m[2][0] + ray.o.y * m[2][1] + ray.o.z * m[2][2] + cam->eye.z;

    ray.o = orig;
    ray.d.x = foo.x + orig.x;
    ray.d.y = foo.y + orig.y;
    ray.d.z = foo.z + orig.z;

    return ray;
}


point get_sample_pos(int x, int y, int sample) {
    point pt;
    //double xsz = 2.0, ysz = WID / aspect;
    static double sf = 0.0;

    if(sf == 0.0) {
        sf = 2.0 / (double)WID;
    }

    pt.x = ((double)x / (double)WID) - 0.5;
    pt.y = -(((double)y / (double)HEI) - 0.65) / aspect;

    if(sample) {
        point jt = jitter(x, y, sample);
        pt.x += jt.x * sf;
        pt.y += jt.y * sf / aspect;
    }
    return pt;
}

/* jitter function taken from Graphics Gems I. */
point jitter(int x, int y, int s) {
    point pt;
    pt.x = urand[(x + (y << 2) + irand[(x + s) & MASK]) & MASK].x;
    pt.y = urand[(y + (x << 2) + irand[(y + s) & MASK]) & MASK].y;
    return pt;
}

void gamut(float *r, float *g, float *b)
{
    float max = m_max(*r,m_max(*g,*b));
    if(max > 1.0f)
    {
        *r = *r/max;
        *g = *g/max;
        *b = *b/max;
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void addVec(point *a, point *b)
{
    a->x = a->x + b->x;
    a->y = a->y + b->y;
    a->z = a->z + b->z;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void subVec(point *a, point *b)
{
    a->x = a->x - b->x;
    a->y = a->y - b->y;
    a->z = a->z - b->z;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void scalarMulVec(point *a,float sc)
{
    a->x = a->x*sc;
    a->y = a->y*sc;
    a->z = a->z*sc;
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


//int inside_grid(grid *gr, point * p) {
//  if ((p->x > gr->p0.x) && (p->x < gr->p1.x)) {
//      if ((p->y > gr->p0.y) && (p->y < gr->p1.y)) {
//          if ((p->z > gr->p0.z) && (p->z < gr->p1.z)) {
//              return TRUE;
//          }
//      }
//  }
//  return FALSE;
//}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

float dist(point *o , point *d)
{
    return sqrt((d->x - o->x)*(d->x - o->x) + (d->y - o->y)*(d->y - o->y) + (d->z - o->z)*(d->z - o->z));
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void normalize(point *a)
{
    float m = a->x*a->x + a->y*a->y + a->z*a->z;
    m = sqrt(m);
    a->x = a->x/m;
    a->y = a->y/m;
    a->z = a->z/m;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void crossProduct(point *a, point *b, point *r)
{
    r->x = a->y*b->z - a->z*b->y;
    r->y = a->z*b->x - a->x*b->z;
    r->z = a->x*b->y - a->y*b->x;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

float dotProduct(point *a, point *b)
{
    return a->x*b->x + a->y*b->y + a->z*b->z;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

float clamp(float x, float min, float max) {
    return (x < min ? min : (x > max ? max : x));
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void setupCamera(camera *c)
{
    //don't use anymore

    c->w.x = c->eye.x - c->lookat.x;
    c->w.y = c->eye.y - c->lookat.y;
    c->w.z = c->eye.z - c->lookat.z;

    normalize(&c->w);

    crossProduct(&c->w,&c->up,&c->u);
    normalize(&c->u);

    crossProduct(&c->u,&c->w,&c->v);
    normalize(&c->v);   

    //printf("u: %f %f %f\n",c->u.x,c->u.y,c->u.z);
    //printf("v: %f %f %f\n",c->v.x,c->v.y,c->v.z);
    //printf("w: %f %f %f\n",c->w.x,c->w.y,c->w.z);

}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


void initCamera(camera *c, point eye, point lookat, int width, int height)
{
    //c->eye.x = -200;
    //c->eye.y = 200;
    //c->eye.z = 600;

    c->eye = eye;
    c->lookat = lookat;
    //c->lookat.y = c->lookat.y + 100;

    c->up.x = 0;
    c->up.y = 1;
    c->up.z = 0;

    c->u.x = 0;
    c->u.y = 0;
    c->u.z = 0;

    c->v.x = 0;
    c->v.y = 0;
    c->v.z = 0;

    c->w.x = 0;
    c->w.y = 0;
    c->w.z = 0;

    c->view.width = width;
    c->view.height = height;

    //c->distance = dist(&c->eye,&c->lookat); 
    c->distance = 1.0f/HALF_FOV;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


ray * generatePrimaryRays(camera *c)
{
    int i,j = 0;
    float xv = 0;
    float yv = 0;

    ray * res;

    res = (ray *) malloc(c->view.width * c->view.height * sizeof(ray));
    for(i = 0 ; i < c->view.width ; i ++)
    {
        for(j = 0 ; j < c->view.height ; j ++)
        {
            xv = i-(c->view.width/2);
            yv = (c->view.height/2)-j;

            res[j+i*c->view.height].d.x = c->u.x*xv + c->v.x*yv - c->w.x*c->distance;
            res[j+i*c->view.height].d.y = c->u.y*xv + c->v.y*yv - c->w.y*c->distance;
            res[j+i*c->view.height].d.z = c->u.z*xv + c->v.z*yv - c->w.z*c->distance;

            res[j+i*c->view.height].o.x = c->eye.x;
            res[j+i*c->view.height].o.y = c->eye.y;
            res[j+i*c->view.height].o.z = c->eye.z;
        }
    }

    return res;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void printPrimaryRays(ray *r, int size)
{
    int i;

    for(i = 0 ; i < size ; i ++)
    {
        printf(" %f %f %f -> %f %f %f\n",r[i].o.x,r[i].o.y,r[i].o.z,r[i].d.x,r[i].d.y,r[i].d.z);
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------


//void setupGrid(grid *g,point max, point min, int nx, int ny, int nz)
//{
//    g->p0 = min;
//    g->p1 = max;
//    g->nx = nx;
//    g->ny = ny;
//    g->nz = nz;
//    g->total = nx*ny*nz;
//}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

float m_max(float a, float b)
{
    return a>b?a:b;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

float m_min(float a, float b)
{
    return a<b?a:b;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

//void intersectBoundingBox(grid * gr, ray * r, point *t_min, point *t_max, float *t0, float *t1)
//{
//        float a;
//        float b;
//        float c;
//
//  point rd;
//
//  rd.x = r->d.x - r->o.x;
//  rd.y = r->d.y - r->o.y;
//  rd.z = r->d.z - r->o.z;
//  normalize(&rd);
//
//  a = 1.0f/rd.x;
//  b = 1.0f/rd.y;
//  c = 1.0f/rd.z;
//
//  if (a >= 0) {
//      t_min->x = (gr->p0.x - r->o.x)*a;
//      t_max->x = (gr->p1.x - r->o.x)*a;
//  } else {
//      t_min->x = (gr->p1.x - r->o.x)*a;
//      t_max->x = (gr->p0.x - r->o.x)*a;
//  }
//
//  if (b >= 0) {
//      t_min->y = (gr->p0.y - r->o.y)*b;
//      t_max->y = (gr->p1.y - r->o.y)*b;
//  } else {
//      t_min->y = (gr->p1.y - r->o.y)*b;
//      t_max->y = (gr->p0.y - r->o.y)*b;
//  }
//
//  if (c >= 0) {
//      t_min->z = (gr->p0.z - r->o.z)*c;
//      t_max->z = (gr->p1.z - r->o.z)*c;
//  } else {
//      t_min->z = (gr->p1.z - r->o.z)*c;
//      t_max->z = (gr->p0.z - r->o.z)*c;
//  }
//
//  *t0 = m_max(m_max(t_min->x, t_min->y), t_min->z);
//  *t1 = m_min(m_min(t_max->x, t_max->y), t_max->z);
//}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void initImage(camera *c,uchar * frame)
{
    int i;
    for(i = 0 ; i < c->view.height*c->view.width ; i ++)
    {
        if(i % 10 == 0)
        {
            frame[3*i + 0] = 250;
            frame[3*i + 1] = 250;
            frame[3*i + 2] = 250;
        }
        else
        {
            frame[3*i + 0] = 0;
            frame[3*i + 1] = 0;
            frame[3*i + 2] = 0;
        }
    }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

int save_bmp(char * file_name, camera *c, uchar *frame)
{
    FILE *fpBMP;

    int i, j;

    // Header and 3 bytes per pixel
    unsigned long ulBitmapSize;
    char ucaBitmapSize[4];
    ulBitmapSize = (c->view.height * c->view.width *3 )+54;



    //printf("x_s,y_s = %d , %d\n",x_size,y_size);
    ucaBitmapSize[3]= (ulBitmapSize & 0xFF000000) >> 24;
    ucaBitmapSize[2]= (ulBitmapSize & 0x00FF0000) >> 16;
    ucaBitmapSize[1]= (ulBitmapSize & 0x0000FF00) >> 8;
    ucaBitmapSize[0]= (ulBitmapSize & 0x000000FF);

    //remove(file_name);

    //Create bitmap file
    fpBMP=fopen(file_name,"wb");
    if(fpBMP == NULL)
        return -1;

    //Write header
    //All values are in big endian order (LSB first)
    // BMP signature + filesize
    fprintf(fpBMP,"%c%c%c%c%c%c%c%c%c%c", 66, 77, ucaBitmapSize[0],
            ucaBitmapSize[1], ucaBitmapSize[2], ucaBitmapSize[3], 0, 0,
            0, 0);


    // Image offset, infoheader size, image width
    fprintf(fpBMP,"%c%c%c%c%c%c%c%c%c%c", 54, 0, 0, 0, 40, 0 , 0, 0,
            (c->view.width & 0x00FF), (c->view.width & 0xFF00)>>8);

    // Image height, number of panels, num bits per pixel
    fprintf(fpBMP,"%c%c%c%c%c%c%c%c%c%c", 0, 0, (c->view.height & 0x00FF),
            (c->view.height & 0xFF00) >> 8, 0, 0, 1, 0, 24, 0);

    // Compression type 0, Size of image in bytes 0 because uncompressed
    fprintf(fpBMP,"%c%c%c%c%c%c%c%c%c%c", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    fprintf(fpBMP,"%c%c%c%c%c%c%c%c%c%c", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    fprintf(fpBMP,"%c%c%c%c", 0, 0 ,0, 0);

    //ok, mas a imagem esta invertida
    for(i = c->view.height -1 ; i >= 0 ; i--)
    {
        for(j = 0 ; j< c->view.width ; j++)
        {
            putc(frame[ 3* (j * c->view.height + i) + 0], fpBMP);
            putc(frame[ 3* (j * c->view.height + i) + 1], fpBMP);
            putc(frame[ 3* (j * c->view.height + i) + 2], fpBMP);

        }
        for(j = 0 ; j < c->view.width % 4 ; j++)
            putc(0,fpBMP);
    }

    /*for(i = c->view.height - 1; i >= 0; i--) {
    //in bitmaps the bottom line of the image is at the beginning of the file
    for(j = 0; j < c->view.width; j++){
    putc(frame[ 3* (i * c->view.width + j) + 0], fpBMP); //b
    putc(frame[ 3* (i * c->view.width + j) + 1], fpBMP); //g
    putc(frame[ 3* (i * c->view.width + j) + 2], fpBMP); //r
    }
    for(j = 0; j < c->view.width % 4; j++)
    putc(0, fpBMP);
    }*/

    fclose(fpBMP);
    return 0;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

enum object intersect3Dscene(int *idx, ray *r, float *t)
{
    enum object obj = NONE;
    int k = 0;

    float dist1 = 0;
    float max_dist1 = FLT_MAX;
    float lowest_dist1 = FLT_MAX;

    float dist2 = 0;
    float max_dist2 = FLT_MAX;
    float lowest_dist2 = FLT_MAX;

    float dist3 = 0;
    float max_dist3 = FLT_MAX;
    float lowest_dist3 = FLT_MAX;

    float dist4 = 0;
    float max_dist4 = FLT_MAX;
    float lowest_dist4 = FLT_MAX;

    //first check for RAY-SPHERE intersection
    for(k = 0 ; k < N_SPHERES ; k++)
    {
        if(intersectSphere(&data[k],r,&dist1) == TRUE)
        {
            dist1 = dist1 - epsilon;
            if(dist1 < max_dist1 && dist1 > 0)
            {
                max_dist1 = dist1; //stores the max distance
                lowest_dist1 = dist1; //stores the lowest distance
                *idx = k; //stores the index of the sphere
                obj = SPHERE;
            }
        }
    }

    //check for RAY-PLANE intersection
    if(intersectPlane(&ground,r,&dist2) == TRUE)
    {
        dist2 = dist2 - epsilon;
        if(dist2 < max_dist2 && dist2 > 0)
        {
            max_dist2 = dist2;
            lowest_dist2 = dist2;
            obj = PLANE;
        }
    }

    //check for RAY-OPENCYLINDER intersection
    if(intersectOpenCylinder(&cylinder,r,&dist3) == TRUE)
    {
        dist3 = dist3 - epsilon;
        if(dist3 < max_dist3 && dist3 > 0)
        {
            max_dist3 = dist3;
            lowest_dist3 = dist3;
            obj = OPEN_CYLINDER;
        }
    }

    //check for RAY-RECTANGLE intersection
    if(intersectRectangle(&rect,r,&dist4) == TRUE)
    {
        dist4 = dist4 - epsilon;
        if(dist4 < max_dist4 && dist4 > 0)
        {
            max_dist4 = dist4;
            lowest_dist4 = dist4;
            obj = RECTANGLE;
        }
    }

    //get the lowest distance
    if(lowest_dist1 < lowest_dist2 && lowest_dist1 < lowest_dist3 && lowest_dist1 < lowest_dist4)
    {
        obj = SPHERE;
        *t = lowest_dist1;
    }
    else if(lowest_dist2 < lowest_dist1 && lowest_dist2 < lowest_dist3 && lowest_dist2 < lowest_dist4)
    {
        obj = PLANE;
        *t = lowest_dist2;
    }
    else if(lowest_dist3 < lowest_dist1 && lowest_dist3 < lowest_dist2 && lowest_dist3 < lowest_dist4)
    {
        obj = OPEN_CYLINDER;
        *t = lowest_dist3;
    }
    else if(lowest_dist4 < lowest_dist1 && lowest_dist4 < lowest_dist2 && lowest_dist4 < lowest_dist1)
    {
        obj = RECTANGLE;
        *t = lowest_dist4;
    }

    return obj;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------- 

color trace(camera cam,ray *r, int iter)
{
    float max_dist = FLT_MAX;
    float lowest_dist;
    int index = -1;

    enum object obj = NONE;

    if(iter > MAX_DEPTH)
    {
        return black;
    }
    else
    {
        index = -1;
        max_dist = FLT_MAX;

        //FALTA O INDEX (tem que ser passado para o SHADE)
        obj = intersect3Dscene(&index,r,&lowest_dist);

        if(obj != NONE && lowest_dist > epsilon)
        {
            //lowest_dist = lowest_dist - epsilon;
            point intersection;
            point dir;

            dir.x = r->d.x - r->o.x;
            dir.y = r->d.y - r->o.y;
            dir.z = r->d.z - r->o.z;
            normalize(&dir);

            intersection.x = r->o.x + lowest_dist*(dir.x);
            intersection.y = r->o.y + lowest_dist*(dir.y);
            intersection.z = r->o.z + lowest_dist*(dir.z);

            //call the shade function and, after that, returns the result of the function
            return shade(cam,&dir,obj,index,&intersection,iter);

            //return white; //just for testing and learning what the shade function actually does. 
        }
        else
        {
            //background color
            return black;
        }
    }
    return black;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

color shade(camera cam,point *incid ,enum object obj, int index, point *p, int iter)
{
    enum object ob;
    point normal;
    point ldir;
    point ldir2;
    float dot2;
    ray r_sombra;
    ray r_sec;
    point w0;
    int il;
    color res2;
    int idx=-1;
    int sombra = FALSE;
    float dist;
    float ka = AMBIENT;
    float dot;
    float sample_f;
    float u,v; //texture coordinates
    float r,g,b;
    color res;
    point reflexao;
    point normal_aux;// = normal;
    point normal_aux2;// = normal;
    point light_reflection; 

    if(obj == SPHERE)
    {
        //ambient
        r = data[index].c.r*ka;
        g = data[index].c.g*ka;
        b = data[index].c.b*ka;

        //normal
        normal.x = p->x - data[index].center.x;
        normal.y = p->y - data[index].center.y;
        normal.z = p->z - data[index].center.z;
        normalize(&normal);
    }
    else if(obj == PLANE)
    {
        //ambient
        r = ground.c.r*ka;
        g = ground.c.g*ka;
        b = ground.c.b*ka;

        //normal
        normal = ground.normal;

        //texture
        u = dotProduct(p,&ground.m_uAxis);
        v = dotProduct(p,&ground.m_vAxis);

    }
    else if(obj == OPEN_CYLINDER)
    {
        //ambient
        r = cylinder.c.r*ka;
        g = cylinder.c.g*ka;
        b = cylinder.c.b*ka;

        //normal
        //sr.normal = Normal((ox + t * dx) * inv_radius, 0.0, (oz + t * dz) * inv_radius);

        normal.x = p->x * (1.0f/cylinder.r);
        normal.y = 0.0f;
        normal.z = p->z * (1.0f/cylinder.r);
        normalize(&normal);
    }
    else if(obj == RECTANGLE)
    {
        //ambient
        r = rect.c.r*ka;
        g = rect.c.g*ka;
        b = rect.c.b*ka;

        //normal
        crossProduct(&rect.a,&rect.b,&normal);
        normalize(&normal);
    }


    //auxiliary vectors
    reflexao = *incid;
    normal_aux = normal;
    normal_aux2 = normal;

    //reflection vector, to be used later
    scalarMulVec(&normal_aux,2*dotProduct(incid,&normal_aux));
    subVec(&reflexao,&normal_aux);
    normalize(&reflexao);

    for(il = 0 ; il < N_LIGHTS ; il++)
    {
        //light direction
        ldir.x = lights[il].l.x - p->x;
        ldir.y = lights[il].l.y - p->y;
        ldir.z = lights[il].l.z - p->z;
        normalize(&ldir);

        //opposite light direction
        ldir2.x = p->x - lights[il].l.x;
        ldir2.y = p->y - lights[il].l.y;
        ldir2.z = p->z - lights[il].l.z;
        normalize(&ldir2);

        //light reflection
        normal_aux2 = normal;
        scalarMulVec(&normal_aux2,2*dotProduct(&normal,&ldir2));
        light_reflection = ldir2;
        subVec(&light_reflection,&normal_aux2);
        normalize(&light_reflection);

        //shadow ray
        r_sombra.o = *p;
        r_sombra.d.x = lights[il].l.x;
        r_sombra.d.y = lights[il].l.y;
        r_sombra.d.z = lights[il].l.z;

        idx=-1;

        ob = intersect3Dscene(&idx,&r_sombra,&dist);
        if(ob != NONE)
            sombra = TRUE;

        //if no shadow exists, compute ambient + diffuse + specular light
        if(!sombra)
        {
            dot = dotProduct(&ldir,&normal);

            if(dot > 0)
            {
                //difuse
                if(obj == SPHERE)
                {
                    r = r + ((lights[il].c.r * lights[il].ls)*(data[index].c.r * data[index].kd/PI))*dot;
                    g = g + ((lights[il].c.g * lights[il].ls)*(data[index].c.g * data[index].kd/PI))*dot;
                    b = b + ((lights[il].c.b * lights[il].ls)*(data[index].c.b * data[index].kd/PI))*dot;
                }
                else if(obj == PLANE)
                {
                    r = r + ((lights[il].c.r * lights[il].ls)*(ground.c.r * ground.kd/PI))*dot;
                    g = g + ((lights[il].c.g * lights[il].ls)*(ground.c.g * ground.kd/PI))*dot;
                    b = b + ((lights[il].c.b * lights[il].ls)*(ground.c.b * ground.kd/PI))*dot;
                }
                else if(obj == OPEN_CYLINDER)
                {
                    r = r + ((lights[il].c.r * lights[il].ls)*(cylinder.c.r * cylinder.kd/PI))*dot;
                    g = g + ((lights[il].c.g * lights[il].ls)*(cylinder.c.g * cylinder.kd/PI))*dot;
                    b = b + ((lights[il].c.b * lights[il].ls)*(cylinder.c.b * cylinder.kd/PI))*dot;
                }
                else if(obj == RECTANGLE)
                {
                    r = r + ((lights[il].c.r * lights[il].ls)*(rect.c.r * rect.kd/PI))*dot;
                    g = g + ((lights[il].c.g * lights[il].ls)*(rect.c.g * rect.kd/PI))*dot;
                    b = b + ((lights[il].c.b * lights[il].ls)*(rect.c.b * rect.kd/PI))*dot;
                }

                //specular

                w0.x = cam.eye.x - p->x;
                w0.y = cam.eye.y - p->y;
                w0.z = cam.eye.z - p->z;
                normalize(&w0);

                dot2 = dotProduct(&light_reflection,&w0);
                if(dot2 > 0)
                {
                    if(obj == SPHERE)
                    {
                        r = r + data[index].ks * pow(dot2,data[index].e) * (lights[il].c.r * lights[il].ls) *dot;
                        g = g + data[index].ks * pow(dot2,data[index].e) * (lights[il].c.r * lights[il].ls) *dot;
                        b = b + data[index].ks * pow(dot2,data[index].e) * (lights[il].c.r * lights[il].ls) *dot;                   
                    }
                    else if(obj == PLANE)
                    {
                        r = r + ground.ks * pow(dot2,ground.e) * (lights[il].c.r * lights[il].ls) *dot;
                        g = g + ground.ks * pow(dot2,ground.e) * (lights[il].c.r * lights[il].ls) *dot;
                        b = b + ground.ks * pow(dot2,ground.e) * (lights[il].c.r * lights[il].ls) *dot;
                    }
                    else if(obj == OPEN_CYLINDER)
                    {
                        r = r + cylinder.ks * pow(dot2,cylinder.e) * (lights[il].c.r * lights[il].ls) *dot;
                        g = g + cylinder.ks * pow(dot2,cylinder.e) * (lights[il].c.r * lights[il].ls) *dot;
                        b = b + cylinder.ks * pow(dot2,cylinder.e) * (lights[il].c.r * lights[il].ls) *dot;
                    }
                    else if(obj == RECTANGLE)
                    {
                        r = r + rect.ks * pow(dot2,rect.e) * (lights[il].c.r * lights[il].ls) *dot;
                        g = g + rect.ks * pow(dot2,rect.e) * (lights[il].c.r * lights[il].ls) *dot;
                        b = b + rect.ks * pow(dot2,rect.e) * (lights[il].c.r * lights[il].ls) *dot;
                    }
                }
            }
        }
    }

    //if object is reflexive (like a mirror)

    //for the ground plane
    if(obj == PLANE)
    {
        if(ground.kf > 0)
        {

            r_sec.o.x = p->x;
            r_sec.o.y = p->y;
            r_sec.o.z = p->z;

            r_sec.d.x = reflexao.x + p->x;
            r_sec.d.y = reflexao.y + p->y;
            r_sec.d.z = reflexao.z + p->z;

            //Recursive call to TRACE a new ray
            res2 = trace(cam, &r_sec,iter+1);
            if(res2.r != 0.0f || res2.g != 0.0f || res2.b != 0.0f  ) //if not black
            {
                sample_f = ((ground.cr * ground.kf)/(dotProduct(&normal,&reflexao)));
                r += ground.c.r * sample_f * res2.r * dotProduct(&normal,&reflexao);
                g += ground.c.g * sample_f * res2.g * dotProduct(&normal,&reflexao);
                b += ground.c.b * sample_f * res2.b * dotProduct(&normal,&reflexao);

            }
        }

        r = r * checkerTexture(u,v);
        g = g * checkerTexture(u,v);
        b = b * checkerTexture(u,v);

    }
    //for spheres...
    else if(obj == SPHERE)
    {
        if(data[index].kf > 0)
        {

            r_sec.o.x = p->x;
            r_sec.o.y = p->y;
            r_sec.o.z = p->z;

            r_sec.d.x = reflexao.x + p->x;
            r_sec.d.y = reflexao.y + p->y;
            r_sec.d.z = reflexao.z + p->z;

            //Recursive call to TRACE a new ray
            res2 = trace(cam, &r_sec,iter+1);
            if(res2.r != 0.0f || res2.g != 0.0f || res2.b != 0.0f  ) //if not black
            {
                sample_f = ((data[index].cr * data[index].kf)/(dotProduct(&normal,&reflexao)));
                r += sample_f * res2.r * dotProduct(&normal,&reflexao);
                g += sample_f * res2.g * dotProduct(&normal,&reflexao);
                b += sample_f * res2.b * dotProduct(&normal,&reflexao);
            }
        }
    }
    else if(obj == OPEN_CYLINDER)
    {
        if(cylinder.kf > 0)
        {

            r_sec.o.x = p->x;
            r_sec.o.y = p->y;
            r_sec.o.z = p->z;

            r_sec.d.x = reflexao.x + p->x;
            r_sec.d.y = reflexao.y + p->y;
            r_sec.d.z = reflexao.z + p->z;

            //Recursive call to TRACE a new ray
            res2 = trace(cam, &r_sec,iter+1);
            if(res2.r != 0.0f || res2.g != 0.0f || res2.b != 0.0f  ) //if not black
            {
                sample_f = ((cylinder.cr * cylinder.kf)/(dotProduct(&normal,&reflexao)));
                r += sample_f * res2.r * dotProduct(&normal,&reflexao);
                g += sample_f * res2.g * dotProduct(&normal,&reflexao);
                b += sample_f * res2.b * dotProduct(&normal,&reflexao);
            }
        }
    }
    else if(obj == RECTANGLE)
    {
        if(rect.kf > 0)
        {

            r_sec.o.x = p->x;
            r_sec.o.y = p->y;
            r_sec.o.z = p->z;

            r_sec.d.x = reflexao.x + p->x;
            r_sec.d.y = reflexao.y + p->y;
            r_sec.d.z = reflexao.z + p->z;

            //Recursive call to TRACE a new ray
            res2 = trace(cam, &r_sec,iter+1);
            if(res2.r != 0.0f || res2.g != 0.0f || res2.b != 0.0f  ) //if not black
            {
                float sample_f = ((rect.cr * rect.kf)/(dotProduct(&normal,&reflexao)));
                r += sample_f * res2.r * dotProduct(&normal,&reflexao);
                g += sample_f * res2.g * dotProduct(&normal,&reflexao);
                b += sample_f * res2.b * dotProduct(&normal,&reflexao);
            }
        }
    }

    //gamut ok...
    gamut(&r,&g,&b);

    res.r = r;
    res.g = g;
    res.b = b;

    return res;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

int intersectOpenCylinder(open_cylinder *cylinder, ray *r, float *t_res)
{
    float t,a,b,c,disc,e,yhit,denom;
    float ox = r->o.x;
    float oy = r->o.y;
    float oz = r->o.z;

    float dx = r->d.x - r->o.x;
    float dy = r->d.y - r->o.y;
    float dz = r->d.z - r->o.z;

    float temp = dx*dx + dy*dy + dz*dz;
    temp = sqrt(temp);

    dx = dx/temp;
    dy = dy/temp;
    dz = dz/temp;

    a = dx * dx + dz * dz;      
    b = 2.0f * (ox * dx + oz * dz);                 
    c = ox * ox + oz * oz - cylinder->r * cylinder->r;
    disc = b * b - 4.0f * a * c ;

    if (disc < 0.0f)
        return FALSE;
    else {  
        e = sqrt(disc);
        denom = 2.0f * a;
        t = (-b - e) / denom;    // smaller root

        if (t > 0) {
            yhit = oy + t * dy;

            if (yhit > cylinder->bottom && yhit < cylinder->top) {
                //tmin = t;
                *t_res = t;

                // test for hitting from inside
                //if (-ray.d * sr.normal < 0.0)
                //  sr.normal = -sr.normal;
                return TRUE;
            }
        } 

        t = (-b + e) / denom;    // larger root
        *t_res = t; 

        if (t > 0) {
            yhit = oy + t * dy;

            if (yhit > cylinder->bottom && yhit < cylinder->top ) {
                //tmin = t;
                *t_res = t;
                // test for hitting inside surface
                //if (-ray.d * sr.normal < 0.0)
                //  sr.normal = -sr.normal;

                return TRUE;
            }
        } 
    }
    return FALSE;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

//Requires Suqare Root

int intersectSphere(sphere * s, ray * r, float *t)
{
    float A,B,C,disc,t0,t1,l,m,n;
    float xd = r->d.x - r->o.x;
    float yd = r->d.y - r->o.y;
    float zd = r->d.z - r->o.z;

    float temp = xd*xd + yd*yd + zd*zd;
    temp = sqrt(temp);

    xd = xd/temp;
    yd = yd/temp;
    zd = zd/temp;

    l = (r->o.x - s->center.x);
    m = (r->o.y - s->center.y);
    n = (r->o.z - s->center.z);

    A = 1;//xd*xd + yd*yd + zd*zd;
    B = 2*(xd*l + yd*m + zd*n);
    C = l*l + m*m + n*n - s->r*s->r;

    //printf("A = %f\n",A);
    //printf("B = %f\n",B);
    //printf("C = %f\n",C);

    disc = B*B - 4*A*C;
    //printf("disc = %f\n",disc);
    if(disc < 0)
    {
        return FALSE;
    }

    t0 = (float) ((-B - sqrt(disc))/2*A);

    if(t0 > epsilon)
    {
        *t = t0;
        return TRUE;
    }

    t1 = (float) ((-B + sqrt(disc))/2*A);

    if(t1 > epsilon)
    {
        *t = t1;
        return TRUE;
    }

    return FALSE;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void generateScene()
{
    int i = 0;

    data[0].center.x = -1.5f;
    data[0].center.y = 1.3f;
    data[0].center.z = -3.0f;
    data[0].r = 0.7f;
    data[0].c.r = 1.0f;
    data[0].c.g = 0.2f;
    data[0].c.b = 0.05f;
    data[0].e = 50.0f;
    data[0].cr = 0.3f;
    data[0].kf = 0.8f; 


    data[1].center.x = 2.5f;
    data[1].center.y = 1.4f;
    data[1].center.z = 0.0f;
    data[1].r = 0.6f;
    data[1].c.r = 0.1f;
    data[1].c.g = 0.85f;
    data[1].c.b = 1.0f;
    data[1].e = 50.0f;
    data[1].cr = 0.4f;
    data[1].kf = 0.4f; 

    data[2].center.x = 0;
    data[2].center.y = 1;
    data[2].center.z = 2;
    data[2].r = 1.0f;
    data[2].c.r = 1.0f;
    data[2].c.g = 0.5f;
    data[2].c.b = 0.1f;
    data[2].e = 60.0f;
    data[2].cr = 0.7f;
    data[2].kf = 0.7f; 

    for(i = 0 ; i < N_SPHERES ; i++)
    {
        data[i].maxBB.x = data[i].center.x + data[i].r;
        data[i].maxBB.y = data[i].center.y + data[i].r;
        data[i].maxBB.z = data[i].center.z + data[i].r;

        data[i].minBB.x = data[i].center.x - data[i].r;
        data[i].minBB.y = data[i].center.y - data[i].r;
        data[i].minBB.z = data[i].center.z - data[i].r;

        data[i].kd = 0.5f;  
        data[i].ks = 0.25f;
    }

    //cylinder
    cylinder.bottom = 0.1f;
    cylinder.top = 2.5f;
    cylinder.r = 1.0f;

    cylinder.c.r = 0.5f;
    cylinder.c.g = 0.5f;
    cylinder.c.b = 1.0f;
    cylinder.e = 60.0f;
    cylinder.cr = 0.7f;
    cylinder.kf = 0.5f; 

    //rectangle
    //data[1].center.x = 2.5f;
    //data[1].center.y = 1.4f;
    //data[1].center.z = 0.0f;

    rect.p0.x = 1.5f;
    rect.p0.y = 0.5f;
    rect.p0.z = 3.0f;

    rect.a.x = 2.5f;
    rect.a.y = 0.0f;
    rect.a.z = -2.0f;

    rect.b.x = 0.0f;
    rect.b.y = 2.0f;
    rect.b.z = 0.0f;

    rect.c.r = 0.8f;
    rect.c.g = 0.8f;
    rect.c.b = 0.8f;
    rect.e = 20.0f;
    rect.cr = 0.7f;
    rect.kf = 0.3f; 

    //also generate a flat ground
    ground.normal.x = 0.0f;
    ground.normal.y = 1.0f;
    ground.normal.z = 0.0f;
    ground.d = 0.0f; //distance

    ground.c.r = 0.3f;
    ground.c.g = 0.3f;
    ground.c.b = 0.3f;

    ground.kf = 1.0f;
    ground.kd = 0.8f;
    ground.cr = 1.0f; //what does it? don't remember!
    ground.e = 1.0f;
    ground.ks = 0.25f;

    //u axis coordinates for texture
    ground.m_uAxis.x = ground.normal.y;
    ground.m_uAxis.y = ground.normal.z;
    ground.m_uAxis.z = -ground.normal.x;

    //v axis coordinates for texture
    crossProduct(&ground.m_uAxis,&ground.normal,&ground.m_vAxis); 
}


//will genereate N_SPHERES into the global 'data' array
void generateRandomSpheres()
{
    int i;

    // uncoment for fixed sphere locations **** 

    //float x = (rand() % AREA) - (rand() % AREA);
    //float y = SPHERE_RAY + 5; // (rand() % AREA) - (rand() % AREA);
    //float z = (rand() % AREA) - (rand() % AREA);

    for(i = 0; i < N_SPHERES ; i++)
    {
        //random spehere locations
        data[i].center.x = (rand() % AREA) - (rand() % AREA);
        data[i].center.y = SPHERE_RAY + 5;//(rand() % AREA) - (rand() % AREA);
        data[i].center.z = (rand() % AREA) - (rand() % AREA);

        // uncoment for fixed sphere locations ****
        /*
           if(i == 0)
           {
           data[i].center.x = x;
           data[i].center.y = y;
           data[i].center.z = z;
           }
           else
           {
           data[i].center.x = data[i-1].center.x + 3*SPHERE_RAY;
           data[i].center.y = data[i-1].center.y; 
           data[i].center.z = data[i-1].center.z;
           }
           */

        //ray is fixed, but can also be random
        data[i].r = SPHERE_RAY;

        //creating bounding boxes
        data[i].maxBB.x = data[i].center.x + data[i].r;
        data[i].maxBB.y = data[i].center.y + data[i].r;
        data[i].maxBB.z = data[i].center.z + data[i].r;

        data[i].minBB.x = data[i].center.x - data[i].r;
        data[i].minBB.y = data[i].center.y - data[i].r;
        data[i].minBB.z = data[i].center.z - data[i].r;


        //RAND_MAX is a constant defined in <cstdlib>. 
        //Its default value may vary between implementations
        //but it is granted to be at least 32767.

        //random colors
        data[i].c.r = (float) rand()/RAND_MAX;
        data[i].c.g = (float) rand()/RAND_MAX;        
        data[i].c.b = (float) rand()/RAND_MAX;

        //fixed colors (from 0 to 1)
        //data[i].c.r = 0.7f;
        //data[i].c.g = 0.7f;        
        //data[i].c.b = 0.7f;

        //printf("data[%d].kf = %f\n",h,data[h].kf);
        data[i].kf = 1.0f;
        data[i].kd = 0.5f;  
        data[i].cr = 1.0f;
        data[i].e = 100;
        data[i].ks = 0.25f;

    }

    //also generate a flat ground
    ground.normal.x = 0.0f;
    ground.normal.y = 1.0f;
    ground.normal.z = 0.0f;
    ground.d = 10.0f; //distance

    ground.c.r = 0.9f;
    ground.c.g = 0.9f;
    ground.c.b = 0.9f;

    ground.kf = 0.85f;
    ground.kd = 0.8f;
    ground.cr = 1.0f; //what does it? don't remember!
    ground.e = 100;
    ground.ks = 0.25f;

    //u axis coordinates for texture
    ground.m_uAxis.x = ground.normal.y;
    ground.m_uAxis.y = ground.normal.z;
    ground.m_uAxis.z = -ground.normal.x;

    //v axis coordinates for texture
    crossProduct(&ground.m_uAxis,&ground.normal,&ground.m_vAxis); 
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

uchar floatToIntColor(float c)
{
    return (uchar) (c * MAX_BYTE);
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

void generateRandomLightSources()
{

    /*int i;

      for(i = 0; i < N_LIGHTS ; i++)
      {
      lights[i].l.x = (rand() % AREA) - (rand() % AREA);
      lights[i].l.y = 1000.0f; //(rand() % AREA) - (rand() % AREA);
      lights[i].l.z = (rand() % AREA) - (rand() % AREA);

    //white color
    lights[i].c = white;

    //light intensity
    lights[i].ls = LIGHT_INT;
    }*/

    lights[0].l.x = -50;
    lights[0].l.y = 200;
    lights[0].l.z = -50;

    //white color
    lights[0].c = white;

    //light intensity
    lights[0].ls = LIGHT_INT;

    lights[1].l.x = 40;
    lights[1].l.y = 40;
    lights[1].l.z = 200;

    //white color
    lights[1].c = white;

    //light intensity
    lights[1].ls = LIGHT_INT;


}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

int intersectRectangle(rectangle *rec, ray *r, float *t_res)
{
    float dx,dy,dz,temp,t1,t2,t,ddota,ddotb,a_len_squared,b_len_squared; 
    point normal,aux1,aux2,p,d;
    crossProduct(&rec->a,&rec->b, &normal);

    dx = r->d.x - r->o.x;
    dy = r->d.y - r->o.y;
    dz = r->d.z - r->o.z;

    temp = dx*dx + dy*dy + dz*dz;
    temp = sqrt(temp);

    dx = dx/temp;
    dy = dy/temp;
    dz = dz/temp;

    aux1.x = rec->p0.x - r->o.x;
    aux1.y = rec->p0.y - r->o.y;
    aux1.z = rec->p0.z - r->o.z;

    t1 = dotProduct(&aux1,&normal);

    aux2.x = dx;
    aux2.y = dy;
    aux2.z = dz;

    t2 = dotProduct(&aux2,&normal);
    t = t1/t2;
    /*float t = (p0 - ray.o) * normal / (ray.d * normal); */

    if (t <= epsilon)
        return FALSE;

    p.x = r->o.x + t*dx;
    p.y = r->o.y + t*dy;
    p.z = r->o.z + t*dz;

    d.x = p.x - rec->p0.x;
    d.y = p.y - rec->p0.y;
    d.z = p.z - rec->p0.z;

    //Point3D p = ray.o + t * ray.d;
    //Vector3D d = p - p0;

    ddota = dotProduct(&d,&rec->a);

    a_len_squared = rec->a.x * rec->a.x + rec->a.y * rec->a.y + rec->a.z * rec->a.z;
    b_len_squared = rec->b.x * rec->b.x + rec->b.y * rec->b.y + rec->b.z * rec->b.z;

    if (ddota < 0.0 || ddota > a_len_squared)
        return FALSE;

    ddotb = dotProduct(&d,&rec->b);

    if (ddotb < 0.0 || ddotb > b_len_squared)
        return FALSE;

    *t_res = t;

    return TRUE;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

int intersectPlane(plane *p, ray *r, float *t)
{
    point rd;
    float d,dist;

    rd.x = r->d.x - r->o.x;
    rd.y = r->d.y - r->o.y;
    rd.z = r->d.z - r->o.z;
    normalize(&rd);

    d = dotProduct(&p->normal, &rd );
    if (d != 0)
    {
        dist = -(dotProduct(&p->normal, &r->o) + p->d) / d;
        if (dist > epsilon && dist < FLT_MAX)
        {
            *t = dist;
            return TRUE;
        }
    }
    return FALSE;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

float checkerTexture(float u, float v)
{
    float ff = fmod(floor(u*TEXTURE_SCALE) + floor(v*TEXTURE_SCALE),2.0f);
    if(ff < 1.0f && ff > -1.0f)
    {
        return 0;
    }
    else
    {
        return 1;
    }

    //return (fmod(floor(u*TEXTURE_SCALE) + floor(v*TEXTURE_SCALE),2.0f) < 1.0f ? 0:1);

    //return 1.0f;
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------



//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------