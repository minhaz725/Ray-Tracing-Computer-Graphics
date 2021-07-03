#include<bits/stdc++.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/glut.h>
#include "bitmap_image.hpp"
using namespace std;
#define pi (2*acos(0.0))

double cameraHeight;
double cameraAngle;
int window_width = 500;
int window_height = 500;
int drawgrid;
int drawaxes;
double angle;
double theta = pi/180;
double fovY;
double move_dist = 2;
double bighemiradius=30;
double smallhemi_and_cylinder_radius=10;
double cylinder_height=100;

double totalpart;
double totalpart_minus_first_hemi;
double totalpart_minus_two_hemis;
double onlybarrelpart;
int bulletcount;


unsigned int image_width = 768;

struct point
{
    double x,y,z;

    point(){
    }
    point(double a,double b,double c)
    {
        x=a;
        y=b;
        z=c;
    }

    point scalermul(double value)
    {
        return point(x*value, y*value ,z*value);
    }

    double dot(point p)
    {
        return x*p.x+y*p.y+z*p.z;
    }
    point cross(point p)
    {
        return point(
                y*p.z-z*p.y,
                z*p.x-x*p.z,
                x*p.y-y*p.x
        );
    }

    double magnitude(point p)
    {
        return sqrt(x*x+y*y+z*z);
    }

    point add(point p)
    {
        return point(x+p.x,y+p.y,z+p.z);
    }
    point sub(point p)
    {
        return point(x-p.x,y-p.y,z-p.z);
    }

    point rotation(point p, double rotationangle)
    {
        point rotated;
        point self(x,y,z);
        rotated = self.scalermul(cos(rotationangle));
        rotated = rotated.add( p.cross(self).scalermul(sin(rotationangle)) );
        return rotated;

    }

};

struct Ray
{
    point start, dir;
    double color[3];
    Ray()
    {
        start = point(0, 0, 0);
        dir = point(0, 0, 0);
        color[0] = 0;
        color[1] = 1; //green
        color[2] = 0;
    }

    Ray(point start, point dir)
    {
        this->start = start;
        this->dir = dir;

       // this->dir.normalize(); magnitude
    }
    void setColor( double r,double g, double b)
    {
        color[0] = r; color[1] =g; color[2]=b;
    }
};


point bulletlocations[1000];
point pos(100,100,100), u(0,0,1), r(-1/sqrt(2),1/sqrt(2),0),l(-1/sqrt(2),-1/sqrt(2),0);

class Object
{
public:
    double color[3];
    double height, width, length;
    double ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff; //Double co_efficients[4]
    double specular_exponent; // shine


    void setAmbientCoeff(double ambientCoeff) {
        ambient_coeff = ambientCoeff;
    }

    void setDiffuseCoeff(double diffuseCoeff) {
        diffuse_coeff = diffuseCoeff;
    }

    void setSpecularCoeff(double specularCoeff) {
        specular_coeff = specularCoeff;
    }

    void setReflectionCoeff(double reflectionCoeff) {
        reflection_coeff = reflectionCoeff;
    }

    void setSpecularExponent(double specularExponent) {
        specular_exponent = specularExponent;
    }

    void setColor( double r,double g, double b)
    {
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }
    point center; //ref point
    double radius;

    point leftCorner;
    double base;

    point A, B, C;

    Object() {}
    virtual void draw() {}
    virtual double intersect(Ray ray, double *color, int level)
    {
        return -1;
    }

    virtual double intersecting_point(Ray ray)
    {
        return -1;
    }
};



class Sphere : public Object
{
public:

    Sphere(point Center, double Radius){
        center=Center;
        radius=Radius;
    }

    void draw()
    {
        glPushMatrix();
        glTranslated(center.x, center.y, center.z);
        glColor3f(color[0], color[1], color[2]);
        glutSolidSphere(radius, 90, 90);
        glPopMatrix();

    }

};


class Floor : public Object {
public:
    double floorWidth, tileWidth;
    point origin;
    int tile_quantity;

    Floor(double floorWidth, double tileWidth) {
        //    this->ambient_coeff = 0.4;
        //    this->diffuse_coeff = this->specular_coeff = this->reflection_coeff = 0.2;
        //    this->specular_exponent = 1.0;

        this->floorWidth = floorWidth;
        this->tileWidth = tileWidth;

        //left corner
        this->origin = point(-floorWidth / 2, -floorWidth / 2, 0);
        this->tile_quantity = this->floorWidth / tileWidth;
    }

    void draw() {
        glBegin(GL_QUADS);
        {
            for (int i = 0; i < tile_quantity; i++) {
                for (int j = 0; j < tile_quantity; j++) {
                    glColor3f((i + j) % 2, (i + j) % 2, (i + j) % 2);

                    glVertex3f(origin.x + i * tileWidth, origin.y + j * tileWidth, origin.z);
                    glVertex3f(origin.x + i * tileWidth, origin.y + (j + 1) * tileWidth, origin.z);
                    glVertex3f(origin.x + (i + 1) * tileWidth, origin.y + (j + 1) * tileWidth, origin.z);
                    glVertex3f(origin.x + (i + 1) * tileWidth, origin.y + j * tileWidth, origin.z);
                }
            }
        }
        glEnd();
    }

};
vector <Object*> objects;
vector <point> lights;

Object *board;
point L,R,U;

void loadTestData()
{
    board = new Floor(3000, 30);
    objects.push_back(board);
    Object *temp;
    point Center(0,0,10);
    double Radius = 10;
    temp=new Sphere(Center, Radius); // Center(0,0,10), Radius 10
    temp->setColor(1,0,0);
    //temp->setCoEfficients(0.4,0.2,0.2,0.2)
    //temp->setShine(1)
    objects.push_back(temp);
    point light1(-50,50,50);
    lights.push_back(light1);

}


void capture()
{
    bitmap_image image(image_width, image_width);

    double plane_dist = (window_height / 2) / tan(theta*(fovY / 2)); //CHK
    point topLeft, l, r, u;

    l = L.scalermul(plane_dist);
    r = R.scalermul(window_width / 2);
    u = U.scalermul( window_height / 2);

    topLeft = u.add(l.add(pos)).sub(r);

    double du = window_width/image_width;
    double dv = window_height/image_width;
// Choose middle of the grid cell

   // r = R.scalermul(0.5*du);
   // u = U.scalermul( 0.5*dv);
   // topLeft = topLeft.add(r).sub(u);

    int nearest;
    double t, tMin;

    point corner;
    double *dummyColor = new double[3];

    for(int i = 0; i < image_width; i++)
    {
        for(int j = 0; j < image_width; j++)
        {
            corner = topLeft.add(R.scalermul( i * du).sub( U.scalermul( j * dv))) ;
            Ray ray(pos, corner.sub(pos));

            nearest = -1;
            tMin = 1e4 * 1.0;
            for(int k = 0; k < objects.size(); k++)
            {
                //by giving level 0 we denote that we  only want to know the nearest object
                t = objects[k]->intersect(ray, dummyColor, 0);

                if(t > 0 && t < tMin)
                    tMin = t, nearest = k;
            }

            if(nearest != -1)
            {
                t = objects[nearest]->intersect(ray, dummyColor, 1);

                for(int c = 0; c < 3; c++)
                {
                    if(dummyColor[c] < 0.0)
                        dummyColor[c] = 0.0;

                    else if(dummyColor[c] > 1.0)
                        dummyColor[c] = 1.0;
                }
            }

            else
                dummyColor[0] = dummyColor[1] = dummyColor[2] = 0.0;

            image.set_pixel(i, j, 255 * dummyColor[0], 255 * dummyColor[1], 255 * dummyColor[2]);
        }
    }

    image.save_image("1605093_rayTracing.bmp");
    image.clear();

    cout << "image captured\n";

}


void drawAxes()
{
    if(drawaxes==1)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);{
            glColor3f(1.0, 0.0, 0.0);
            glVertex3f( 100,0,0);
            glVertex3f(-100,0,0);

            glColor3f(0.0, 1.0, 0.0);
            glVertex3f(0,-100,0);
            glVertex3f(0, 100,0);

            glColor3f(0.0, 0.0, 1.0);
            glVertex3f(0,0, 100);
            glVertex3f(0,0,-100);
        }glEnd();
    }
}


void drawGrid()
{
    int i;
    if(drawgrid==1)
    {
        glColor3f(0.6, 0.6, 0.6);	//grey
        glBegin(GL_LINES);{
            for(i=-8;i<=8;i++){

                if(i==0)
                    continue;	//SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i*10, -90, 0);
                glVertex3f(i*10,  90, 0);
                //lines parallel to X-axis
                glVertex3f(-90, i*10, 0);
                glVertex3f( 90, i*10, 0);



                glVertex3f(i*10, 0, -90);
                glVertex3f(i*10, 0,  90);

                glVertex3f(-90, 0, i*10);
                glVertex3f(90, 0, i*10);


                //glVertex3f(0, -90, i*10);
                //glVertex3f(0,  90, i*10);

                //glVertex3f(0, i*10,-90);
                //glVertex3f(0, i*10,90);


            }
        }glEnd();
    }
}

void drawSquare(double a)
{
    //glColor3f(1.0,0.0,0.0);
    glBegin(GL_QUADS);{
        glVertex3f( a, a,2);
        glVertex3f( a,-a,2);
        glVertex3f(-a,-a,2);
        glVertex3f(-a, a,2);
    }glEnd();
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
            glVertex3f(points[i].x,points[i].y,0);
            glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCylinder(double radius,double height,int segments)
{
    int i;
    double shade;
    int altercolor = 1;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        //glColor3f(shade,shade,shade);

        altercolor = 1 - altercolor;
        glColor3f(altercolor,altercolor,altercolor);
        glBegin(GL_QUADS);
        {
            glVertex3f(points[i].x,points[i].y,0);
            glVertex3f(points[i+1].x,points[i+1].y,0);
            glVertex3f(points[i].x,points[i].y,height);
            glVertex3f(points[i+1].x,points[i+1].y,height);
        }
        glEnd();
    }
}


void drawSphere(double radius,int slices,int stacks)
{
    struct point points[100][100];
    int i,j;
    double h,r;
    //generate points
    for(i=0;i<=stacks;i++)
    {
        h=radius*sin(((double)i/(double)stacks)*(pi/2));
        r=radius*cos(((double)i/(double)stacks)*(pi/2));
        for(j=0;j<=slices;j++)
        {
            points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
            points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
            points[i][j].z=h;
        }
    }
    //draw quads using generated points
    for(i=0;i<stacks;i++)
    {
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
        for(j=0;j<slices;j++)
        {
            glBegin(GL_QUADS);{
                //upper hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
            }glEnd();
        }
    }
}



void drawHemiSphere(double radius,int slices,int stacks)
{
    struct point points[100][100];
    int i,j;
    double h,r;
    int altercolor = 1;
    //generate points
    for(i=0;i<=stacks;i++)
    {
        h=radius*sin(((double)i/(double)stacks)*(pi/2));
        r=radius*cos(((double)i/(double)stacks)*(pi/2));

        for(j=0;j<=slices;j++)
        {
            points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
            points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
            points[i][j].z=h;
        }
    }
    //draw quads using generated points
    for(i=0;i<stacks;i++)
    {
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
        for(j=0;j<slices;j++)
        {
            altercolor = 1 - altercolor;
            glColor3f(altercolor,altercolor,altercolor);
            glBegin(GL_QUADS);{
                //upper hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

            }glEnd();
        }
    }
}



void drawFlower(double radius,int slices,int stacks)
{
    struct point points[100][100];
    int i,j;
    double h,r;
    int altercolor=1;
    //generate points
    for(i=0;i<=stacks;i++)
    {
        h=radius*sin(((double)i/(double)stacks)*(pi/2));
        r=radius*cos(((double)i/(double)stacks)*(pi/2));
        r = 2 * radius - r;
        for(j=0;j<=slices;j++)
        {
            points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
            points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
            points[i][j].z=h;
        }
    }
    //draw quads using generated points
    for(i=0;i<stacks;i++)
    {
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
        for(j=0;j<slices;j++)
        {
            altercolor = 1 - altercolor;
            glColor3f(altercolor,altercolor,altercolor);
            glBegin(GL_QUADS);{
                //upper hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

            }glEnd();
        }
    }
}

void drawLightSources()
{
   // if(!display_light_source)
   //     return;

    glColor3f(1, 1, 1);
    for(int i = 0; i < lights.size(); i++)
    {
        glPushMatrix();
        glTranslated(lights[i].x, lights[i].y, lights[i].z);
        glutSolidSphere(1, 90, 90);
        glPopMatrix();
    }
}

void drawObjects()
{
    drawLightSources();
    for(auto & object : objects)
        object->draw();
}





void drawSS()
{

    glRotatef(90,0,1,0);
    //white board
    glPushMatrix();
    {
        glTranslatef(0.0,0.0,-800.0);
        glColor3f(1.0,1.0,1.0);
        drawSquare(300);
    }
    glPopMatrix();
    glPushMatrix();
    {

        {
            //first hemi
            glRotatef(totalpart,1,0,0);
            drawHemiSphere(bighemiradius,100,20);
        }

        {
            //2nd hemi

            glRotatef(totalpart_minus_first_hemi,0,1,0);

            glRotatef(180,0,1,0);
            drawHemiSphere(bighemiradius,80,20);
        }

        {
            //3rd hemi
            //3rder center shuru 2nd tar center + 3rder rad
            glRotatef(onlybarrelpart,0,0,1);
            glTranslatef(0,0,bighemiradius);
            glRotatef(totalpart_minus_two_hemis,0,1,0);
            glTranslatef(0,0,-bighemiradius);
            glTranslatef(0,0,bighemiradius+smallhemi_and_cylinder_radius);
            glRotatef(180,0,1,0);
            drawHemiSphere(smallhemi_and_cylinder_radius,80,20);
        }


        {
            //cyclinder
            glRotatef(180,0,1,0);
            drawCylinder(smallhemi_and_cylinder_radius,cylinder_height,80);

        }


        {
            //flower
            glTranslatef(0,0,cylinder_height);
            drawFlower(smallhemi_and_cylinder_radius,80,20);
        }


    }
    glPopMatrix();



    for(int i=0;i<bulletcount;i++)
    {
        glPushMatrix();
        {

            glRotatef(bulletlocations[i].x,1,0,0);
            glRotatef(bulletlocations[i].y,0,1,0);
            glRotatef(bulletlocations[i].z,0,0,1);
            glTranslatef(0.0,0.0,-795.0);
            glColor3f(1.0,0.0,0.0);
            drawSquare(10);
        }
        glPopMatrix();
    }




}

void keyboardListener(unsigned char key, int x,int y){

    switch(key){


        case '1':
            r=r.rotation(u,theta);
            l=l.rotation(u,theta);
            break;
        case '2':
            r=r.rotation(u,-theta);
            l=l.rotation(u,-theta);
            break;
        case '3':
            l=l.rotation(r,theta);
            u=u.rotation(r,theta);
            break;
        case '4':
            l=l.rotation(r,-theta);
            u=u.rotation(r,-theta);
            break;
        case '5':
            r=r.rotation(l,theta);
            u=u.rotation(l,theta);
            break;
        case '6':
            r=r.rotation(l,-theta);
            u=u.rotation(l,-theta);

            break;
        case '9':
            drawgrid=1-drawgrid;
            break;
        case 'q':
            if(totalpart<-45) break;
            totalpart=(totalpart-1);

            break;
        case 'w':
            if(totalpart>45) break;
            totalpart=(totalpart+1);

            break;
        case 'e':
            if(totalpart_minus_first_hemi<-45) break;
            totalpart_minus_first_hemi=(totalpart_minus_first_hemi-1);

            break;
        case 'r':
            if(totalpart_minus_first_hemi>45) break;
            totalpart_minus_first_hemi=(totalpart_minus_first_hemi+1);

            break;
        case 'a':
            if(totalpart_minus_two_hemis>45) break;
            totalpart_minus_two_hemis=(totalpart_minus_two_hemis+1);

            break;
        case 's':
            if(totalpart_minus_two_hemis<-45) break;
            totalpart_minus_two_hemis=(totalpart_minus_two_hemis-1);


            break;
        case 'd':
            if(onlybarrelpart>45) break;
            onlybarrelpart=(onlybarrelpart+1);

            break;
        case 'f':
            if(onlybarrelpart<-45) break;
            onlybarrelpart=(onlybarrelpart-1);

            break;
        case '0':
            capture();
            break;


        default:
            break;
    }
}


void specialKeyListener(int key, int x,int y){
    switch(key){
        case GLUT_KEY_DOWN:		//down arrow key
            pos = pos.sub(l.scalermul(move_dist));
            break;
        case GLUT_KEY_UP:		// up arrow key
            //cameraHeight += 3.0;
            pos = pos.add(l.scalermul(move_dist));
            break;

        case GLUT_KEY_RIGHT:
            //cameraAngle += 0.03;
            pos = pos.add(r.scalermul(move_dist));
            break;
        case GLUT_KEY_LEFT:
            //cameraAngle -= 0.03;
            pos = pos.sub(r.scalermul(move_dist));
            break;

        case GLUT_KEY_PAGE_UP:
            pos = pos.add(u.scalermul(move_dist));
            break;
        case GLUT_KEY_PAGE_DOWN:
            pos = pos.sub(u.scalermul(move_dist));
            break;

        case GLUT_KEY_INSERT:
            break;

        case GLUT_KEY_HOME:
            break;
        case GLUT_KEY_END:
            break;

        default:
            break;
    }
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
    switch(button){
        case GLUT_LEFT_BUTTON:
            if(state == GLUT_DOWN){
                double  angleLimit = (atan2 (300, 800))*180/pi;
                // if ( fabs(totalpart) < angleLimit && fabs(totalpart_minus_first_hemi+ totalpart_minus_two_hemis) < angleLimit )
                //    {
                bulletlocations[bulletcount].x=totalpart;
                bulletlocations[bulletcount].y=totalpart_minus_first_hemi;
                bulletlocations[bulletcount].z=totalpart_minus_two_hemis;
                bulletcount++;
                //      }
            }
            break;

        case GLUT_RIGHT_BUTTON:
            if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
                drawaxes=1-drawaxes;
            }
            break;

        case GLUT_MIDDLE_BUTTON:
            //........
            break;

        default:
            break;
    }
}



void display(){

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
    //gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    gluLookAt(pos.x,pos.y,pos.z,	pos.x+l.x, pos.y+l.y, pos.z+l.z,	  u.x,u.y,u.z);


    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawObjects();
    drawGrid();
    //drawSS();
    glutSwapBuffers();
}


void animate(){
    //angle+=0.05;
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init(){
    //codes for initialization
    drawgrid=0;
    drawaxes=1;
    cameraHeight=150.0;
    cameraAngle=1.0;
    angle=0;
    fovY = 80;
    //clear the screen
    glClearColor(0,0,0,0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(fovY,	1,	1,	1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv){
    glutInit(&argc,argv);
    glutInitWindowSize(window_width , window_height );
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("Ray Tracing");

    loadTestData();
    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}
