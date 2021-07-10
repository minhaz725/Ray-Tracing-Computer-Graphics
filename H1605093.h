#include<bits/stdc++.h>
#include<stdlib.h>
#include<math.h>

#ifndef RTC_H1605093_H
#define RTC_H1605093_H


struct point
{
    double x,y,z;
    double color[3];

    point(){
    }
    point(double a,double b,double c)
    {
        x=a;
        y=b;
        z=c;
        color[0]=color[1]=color[2]=1;
    }

    void setColor( double r,double g, double b)
    {
        color[0] = r;
        color[1] = g;
        color[2] = b;
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

    void normalize()
    {
        double sqroot = sqrt(x * x + y * y + z * z);

        x /= sqroot;
        x *= 1.0;
        y /= sqroot;
        y *= 1.0;
        z /= sqroot;
        z *= 1.0;
    }

    point getReflection(point normal)
    {
        point original_vec(x,y,z);
        double coeff = original_vec.dot(normal) * 2;
        point reflected_vec = original_vec.sub(normal.scalermul( coeff));

        reflected_vec.normalize();

        return reflected_vec;
    }

    point getNormal(point sec)
    {
        point fir(x,y,z);
        point temp = fir.sub(sec);
        temp.normalize();

        return temp;
    }

    point getReverseReflection(point normal)
    {
        point original_vec(x,y,z);
        double coeff = original_vec.dot( normal) * 2;
        point reflected_vec = normal.scalermul( coeff).sub(original_vec);
        reflected_vec.normalize();

        return reflected_vec;
    }


};

struct Ray
{
    point start, dir;
    Ray()
    {
        start = point(0, 0, 0);
        dir = point(0, 0, 0);
    }

    Ray(point start, point dir)
    {
        this->start = start;
        this->dir = dir;
        this->dir.normalize();
    }
};

class Object
{
public:
    double color[3];
    double ambient_coeff, diffuse_coeff, specular_coeff, reflection_coeff;
    double shine;

    void setCoeffs(double ambientCoeff,double diffuseCoeff,double specularCoeff, double reflectionCoeff)
    {
        ambient_coeff = ambientCoeff;
        diffuse_coeff = diffuseCoeff;
        specular_coeff = specularCoeff;
        reflection_coeff = reflectionCoeff;
    }

    void setShine(double shinevar) {
        shine = shinevar;
    }

    void setColor( double r,double g, double b)
    {
        color[0] = r;
        color[1] = g;
        color[2] = b;
    }
    point center; //ref point
    double radius;

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


#endif //RTC_H1605093_H
