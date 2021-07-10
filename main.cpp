#include<bits/stdc++.h>
#include<math.h>
#include "H1605093.h"
#include <windows.h>
#include <GL/glut.h>
#include "bitmap_image.hpp"
using namespace std;
#define pi (2*acos(0.0))
#define EPSILON  0.0000001
#define NP 1
#define FP 1000

double cameraHeight;
double cameraAngle;
double window_width = 500;
double window_height = 500;
int drawgrid;
int drawaxes;
double theta = pi/180;
double fovY;
double move_dist = 2;
int levelOfRecursion;
int numberOfLights;

unsigned int image_width;
int numberOfObjects;


point pos(100, 90, 20),u(0, 0, 1), r(-1/sqrt(2.0),  1/sqrt(2.0), 0), l(-1/sqrt(2.0),  -1/sqrt(2.0), 0);

bool inNearAndFar(double t)
{
    return (t >= NP && t <= FP);
}
vector <Object*> objects;
vector <point> lights;
Object *board;


class Sphere : public Object
{
public:

    Sphere(point Center, double Radius){
        center=Center;
        radius=Radius;
    }
    Sphere(){
    }

    void draw()
    {
        glPushMatrix();
        glTranslated(center.x, center.y, center.z);
        glColor3f(color[0], color[1], color[2]);
        glutSolidSphere(radius, 90, 90);
        glPopMatrix();

    }

    double intersecting_point(Ray ray) override
    {
        point r0 = ray.start, rd = ray.dir;
        //point r0 (0,100,0), rd(0,1,0)
        double a = ray.dir.dot(ray.dir);     //as rd is unit
        double b = 2 * rd.dot(r0.sub(center));
        double c = r0.sub(center).dot(r0.sub(center)) - (radius * radius);
        double d = (b * b) - (4 * a * c);

        if(d < 0)
            return -1;

        d = sqrt(d);
        double t1 = (- b - d) / (2 * a);
        double t2 = (- b + d) / (2 * a);

        return min(t1, t2);
    }

    double intersect(Ray ray, double *current_color, int level) override {
        double t = intersecting_point(ray);
        if (t <= 0)
            return -1;

        if (!inNearAndFar(t))
            return -1;

        if (level == 0)
            return t;

        for (int c = 0; c < 3; c++)
            current_color[c] = color[c] * ambient_coeff;



        point intersectionPoint = ray.start.add(ray.dir.scalermul(t));
        point normal = intersectionPoint.getNormal(center);
        point reflection = ray.dir.getReflection(normal);


        for (int i = 0; i < lights.size(); i++) {
            point L = lights[i].sub(intersectionPoint);
            L.setColor(lights[i].color[0],lights[i].color[1],lights[i].color[2]);
            L.normalize();
            //cout << L.color[0] << L.color[1] << L.color[2] <<endl;
            point start = intersectionPoint.add(L.scalermul(EPSILON));
            Ray sunLight(start, L);

            point N = intersectionPoint.getNormal(center);
            point R = L.getReverseReflection(N);

            point V = ray.start.sub(intersectionPoint);
            V.normalize();


            bool obscured = false;
            for (int j = 0; j < objects.size(); j++) {
                double temp = objects[j]->intersecting_point(sunLight);

                if (temp > 0) {
                    obscured = true;
                    break;
                }
            }

            if (!obscured) {
                double cosTheta = max(0.0, L.dot(N));
                double cosPhi = max(0.0, R.dot(V));

                double lambart = diffuse_coeff * cosTheta;
                double phong = pow(cosPhi, shine) * specular_coeff;

                for (int c = 0; c < 3; c++)
                {
                    current_color[c] += L.color[c]*lambart * color[c];
                    current_color[c] += L.color[c]*(phong*1.0) * color[c];
                }

            }
        }
        int nearest, tMin, t2;

        if(level < levelOfRecursion)
        {
            point start = intersectionPoint.add(reflection);
            Ray reflectionRay(start, reflection);

            nearest = -1;
            tMin = 1e4;
            double *reflected_color = new double[3];
            reflected_color[0] = reflected_color[1] = reflected_color[2] = 0.0;

            for(int k = 0; k < objects.size(); k++)
            {
                t2 = objects[k]->intersect(reflectionRay, reflected_color, 0);

                if(t2 > 0 && t2 < tMin)
                    tMin = t2, nearest = k;
            }

            if(nearest != -1)
            {
                t2 = objects[nearest]->intersect(reflectionRay, reflected_color, level + 1);

                for(int c = 0; c < 3; c++)
                    current_color[c] += (reflected_color[c] * reflection_coeff);
            }

            delete[] reflected_color;
        }

        return t;
    }


};


class Floor : public Object {
public:
    double floorWidth, tileWidth;
    point origin;
    int tile_quantity;

    Floor(double floorWidth, double tileWidth) {
        this->ambient_coeff = 0.4;
        this->diffuse_coeff = this->specular_coeff = this->reflection_coeff = 0.2;
        this->shine = 1.0;

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

    point getNormal()
    {
        return point(0, 0, 1);
    }

    bool CheckSurf(point p)
    {
        point dist = p.sub(origin);

        if(dist.x < 0 || dist.x > floorWidth || dist.y < 0 || dist.y > floorWidth)
            return false;

        return true;
    }

    void setColor(point p)
    {
        point dist = p.sub(origin);

        int tileX = (dist.x / tileWidth);
        int tileY = (dist.y / tileWidth);

        for(int c = 0; c < 3; c++)
            color[c] = (tileX + tileY) % 2;
    }

    double intersecting_point(Ray ray)
    {
        if(ray.dir.z == 0)
            return -1;

        double t = (-ray.start.z / ray.dir.z);

        point intersectionPoint = ray.start.add(ray.dir.scalermul(t));

        CheckSurf(intersectionPoint);
        setColor(intersectionPoint);

        return t;
    }

    double intersect(Ray ray, double *current_color, int level) {
        double t = intersecting_point(ray);

        if (t <= 0)
            return -1;

        if (!inNearAndFar(t))
            return -1;

        if (level == 0)
            return t;

        point normal = getNormal();

        point intersectionPoint = ray.start.add(ray.dir.scalermul(t));

        for (int c = 0; c < 3; c++)
            current_color[c] = color[c] * ambient_coeff;

        point reflection = ray.dir.getReflection( normal);


        for (int i = 0; i < lights.size(); i++) {
            point L = lights[i].sub( intersectionPoint);
            L.setColor(lights[i].color[0],lights[i].color[1],lights[i].color[2]);
            L.normalize();

            point start = intersectionPoint.add(L.scalermul( EPSILON));
            Ray sunLight(start, L);

            point N = getNormal();
            point R = L.getReverseReflection( N);

            point V = ray.start.sub(intersectionPoint);
            V.normalize();

            bool obscured = false;
            for (int j = 0; j < objects.size(); j++) {
                double temp = objects[j]->intersecting_point(sunLight);

                if (temp > 0) {
                    obscured = true;
                    break;
                }
            }

            if (!obscured) {
                double cosTheta = max(0.0, L.dot(N));
                double cosPhi = max(0.0, R.dot( V));

                double lambart = diffuse_coeff * cosTheta;
                double phong = pow(cosPhi, shine) * specular_coeff;

                for (int c = 0; c < 3; c++)
                {
                    current_color[c] +=  L.color[c]*lambart * color[c];
                    current_color[c] += L.color[c]*(phong * 1.0) * color[c];
                }

            }

        }


        int nearest, tMin, t2;

        if(level < levelOfRecursion)
        {
            point start = intersectionPoint.add(reflection);
            Ray reflectionRay(start, reflection);

            nearest = -1;
            tMin = 1e4;
            double *reflected_color = new double[3];
            reflected_color[0] = reflected_color[1] = reflected_color[2] = 0.0;

            for(int k = 0; k < objects.size(); k++)
            {
                t2 = objects[k]->intersect(reflectionRay, reflected_color, 0);

                if(t2 > 0 && t2 < tMin)
                    tMin = t2, nearest = k;
            }

            if(nearest != -1)
            {
                t2 = objects[nearest]->intersect(reflectionRay, reflected_color, level + 1);

                for(int c = 0; c < 3; c++)
                    current_color[c] += (reflected_color[c] * reflection_coeff);
            }

            delete[] reflected_color;
        }
        return t;
    }

};

void loadTestData()
{
    board = new Floor(3000, 30);
    objects.push_back(board);

    Object *temp;
    cin >> levelOfRecursion;
    cin >> image_width;
    cin >> numberOfObjects;

    string str;
    for(int i = 0; i < numberOfObjects; i++)
    {
        cin >> str;

        if(str == "sphere")
        {
            temp = new Sphere();
            cin >> temp->center.x >> temp->center.y >> temp->center.z;
            cin >> temp->radius;
        }


        cin >> temp->color[0] >> temp->color[1] >> temp->color[2];
        cin >> temp->ambient_coeff >> temp->diffuse_coeff >> temp->specular_coeff >> temp->reflection_coeff;
        cin >> temp->shine;

        objects.push_back(temp);
    }

    cin >> numberOfLights;;

    point p;
    for(int i = 0; i < numberOfLights; i++)
    {
        cin >> p.x >> p.y >> p.z;
        cin>>p.color[0]>>p.color[1]>>p.color[2];
        lights.push_back(p);
    }

}


void capture()
{
    bitmap_image image(image_width, image_width);

    double plane_dist = (window_height / 2) / tan(theta*(fovY / 2));
    point topLeft, L, R, U;

    L = l.scalermul(plane_dist);
    R = r.scalermul(window_width / 2);
    U = u.scalermul( window_height / 2);

    topLeft = pos.add(L).add(U).sub(R);

    double du = window_width/image_width;
    double dv = window_height/image_width;
// Choose middle of the grid cell

    R = r.scalermul(0.5*du);
    U = u.scalermul( 0.5*dv);
    topLeft = topLeft.add(R).sub(U);

    int nearest;
    double t, tMin;

    point corner;
    double *dummyColor = new double[3];

    for(int i = 0; i < image_width; i++)
    {
        for(int j = 0; j < image_width; j++)
        {

            R = r.scalermul( i * du);  //
            U = u.scalermul( j * dv);
            corner = topLeft.add(R).sub(U);
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

void drawLightSources()
{
    // if(!display_light_source)
    //     return;

    glColor3f(1, 1, 1);
    for(int i = 0; i < lights.size(); i++)
    {
        glPushMatrix();
        glColor3f(lights[i].color[0], lights[i].color[1], lights[i].color[2]);
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
    fovY = 90;
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
    freopen("F:\\RTC\\input.txt", "r", stdin);
    glutInit(&argc,argv);
    glutInitWindowSize(window_width , window_height );
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("Ray Tracing 1605093");

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
