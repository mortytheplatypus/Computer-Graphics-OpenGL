#include<bits/stdc++.h>
#include <GL/glut.h>
#include "bitmap_image.hpp"
using namespace std;

#define PI acos(-1.0)
#define DEGTORAD(degree) (degree * PI / 180.0)
#define INF 100

enum ObjectType {
    NONE, CUBE, PYRAMID, SPHERE 
};

void drawAxes() {
    glLineWidth(2);
    glBegin(GL_LINES);
        // x axis := RED
        glColor3f(1, 0, 0);  
        glVertex3f(-INF, 0, 0);
        glVertex3f(INF, 0, 0);
        
        // y axis := GREEN
        glColor3f(0, 1, 0);   
        glVertex3f(0, -INF, 0);
        glVertex3f(0, INF, 0);

        // z axis := BLUE
        glColor3f(0, 0, 1); 
        glVertex3f(0, 0, -INF);
        glVertex3f(0, 0, INF);
    glEnd();
}

struct Point {
	double x, y, z;

    Point(double x=0, double y=0, double z=0) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Point(const Point &p) {
        this->x = p.x;
        this->y = p.y;
        this->z = p.z;
    }

    double length() {
        return sqrt(x*x + y*y + z*z);
    }

    Point normalize() {
        double len = this->length();
        x /= len;
        y /= len;
        z /= len;

        return Point(x, y, z);
    }

    Point operator+(Point p) {
        return Point(x+p.x, y+p.y, z+p.z);
    }

    Point operator-(Point p) {
        return Point(x-p.x, y-p.y, z-p.z);
    }

    Point operator*(double a) {
        return Point(x*a, y*a, z*a);
    }

    Point operator/(double a) {
        return Point(x/a, y/a, z/a);
    }

    friend istream& operator>>(istream &in, Point &p) {
        in >> p.x >> p.y >> p.z;
        return in;
    }

    friend ostream& operator<<(ostream &out, Point &p) {
        out << fixed << setprecision(7) << p.x << " " << p.y << " " << p.z;
        return out;
    }
};

struct Color {
    double r, g, b;

    Color(double r=0.0, double g=0.0, double b=0.0) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
};

struct Ray {
    Point origin, dir;

    Ray(Point origin=Point(), Point dir=Point(1, 0, 0)) {
        this->origin = origin;
        dir.normalize();
        this->dir = dir;
    }
};

double dotProduct(Point p1, Point p2) {
    return p1.x * p2.x + p1.y * p2.y + p1.z * p2.z;
}

Point crossProduct(Point p1, Point p2) {
    double x, y, z;

    x = p1.y * p2.z - p1.z * p2.y; 
    y = p1.z * p2.x - p1.x * p2.z;
    z = p1.x * p2.y - p1.y * p2.x;

    return Point(x, y, z);
}

Point rotate3D(Point axis, Point rotationAxis, double angle) {
    Point cross = crossProduct(axis, rotationAxis);

    return axis*cos(angle) + cross*sin(angle);
}

struct Component {
    Color color;
    double ambient, diffuse, specular, reflection;
    double shininess;

    Component(Color color=Color(), double ambient=0.0, double diffuse=0.0, double specular=0.0, double reflection=0.0, double shininess=0.0) {
        this->color = color;
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->specular = specular;
        this->reflection = reflection;
        this->shininess = shininess;
    }
};

struct Light {
    Point position;
    double falloff;

    Light(const Point position=Point(), double falloff=0.0) {
        this->position = position;
        this->falloff = falloff;
    }

    void draw() {
        double rgb[3] = {250.0/255.0, 250.0/255.0, 50.0/255.0};
        glColor3f(rgb[0], rgb[1], rgb[2]);

        int slices = 50, stacks = 50;
        glPushMatrix();
            glTranslatef(position.x, position.y, position.z);
            glutSolidSphere(7.5, slices, stacks);
        glPopMatrix();
    }
};

struct Spotlight {
    Point position, look;
    double falloff, cutoff;

    Spotlight(const Point position=Point(), const Point look=Point(), double falloff=0.0, double cutoff=0.0) {
        this->position = position;
        this->look = look;
        this->falloff = falloff;
        this->cutoff = cutoff;
    }

    void draw() {
        double rgb[3] = {1.0, 1.0, 155.0/255.0};
        glColor3f(rgb[0], rgb[1], rgb[2]);

        int slices = 50, stacks = 50;
        glPushMatrix();
            glTranslatef(position.x, position.y, position.z);
            glutSolidSphere(5, slices, stacks);
        glPopMatrix();
    }
    
};

// extern vector <Light> lights;
// extern vector <Spotlight> spotlights;
// extern vector <Sphere> spheres;
// extern vector <Cube> cubes;
// extern vector <Pyramid> pyramids;


struct Cube {
    Point lowestPoint;
    double side;
    Component component;

    Cube(Point lowestPoint=Point(), double side=0.0) {
        this->lowestPoint = lowestPoint;
        this->side = side;
        this->component = Component();
    }

    Cube(const Cube &c) {
        this->lowestPoint = c.lowestPoint;
        this->side = c.side;
        this->component = c.component;
    }

    Ray getNormal(Point point, Ray incidentRay) {
        // xy planes
        if (point.z == lowestPoint.z) {
            return Ray(point, Point(0, 0, -1));
        } else if (point.z == lowestPoint.z + side) {
            return Ray(point, Point(0, 0, 1));
        }

        // yz planes
        if (point.x == lowestPoint.x) {
            return Ray(point, Point(-1, 0, 0));
        } else if (point.x == lowestPoint.x + side) {
            return Ray(point, Point(1, 0, 0));
        }

        // xz planes
        if (point.y == lowestPoint.y) {
            return Ray(point, Point(0, -1, 0));
        } else if (point.y == lowestPoint.y + side) {
            return Ray(point, Point(0, 1, 0));
        }
    }

    void draw() {
        double x = lowestPoint.x, y = lowestPoint.y, z = lowestPoint.z;

        glColor3f(component.color.r, component.color.g, component.color.b);

        for (int i=0; i<2; i++) {
            glBegin(GL_QUADS);
            {
                glVertex3f(x, y, z + i*side);
                glVertex3f(x + side, y, z + i*side);
                glVertex3f(x + side, y + side, z + i*side);
                glVertex3f(x, y + side, z + i*side);
            }
            glEnd();
        }

        for (int i=0; i<2; i++) {
            glBegin(GL_QUADS);
            {
                glVertex3f(x, y + i*side, z);
                glVertex3f(x + side, y + i*side, z);
                glVertex3f(x + side, y + i*side, z + side);
                glVertex3f(x, y + i*side, z + side);
            }
            glEnd();
        }

        for (int i=0; i<2; i++) {
            glBegin(GL_QUADS);
            {
                glVertex3f(x + i*side, y, z);
                glVertex3f(x + i*side, y + side, z);
                glVertex3f(x + i*side, y + side, z + side);
                glVertex3f(x + i*side, y, z + side);
            }
            glEnd();
        }
    }

    double getIntersectingT(Ray ray) {
        /*
            1. Two planes parallel to xy plane
            2. Two planes parallel to yz plane
            3. Two planes parallel to xz plane

            P(x, y, z), n(A, B, C) 
            n.P + D = 0
        */

        Point n, R0 = ray.origin, Rd = ray.dir;
        double tMin = INF, D, t;

        // xy planes := n = (0, 0, 1) and (0, 0, -1)
        n = Point(0, 0, -1);
        for (int i=0; i<2; i++) {
            if (dotProduct(Rd, n) == 0) {
                continue;
            }

            D = lowestPoint.z + i*side;
            t = (D - dotProduct(R0, n)) / dotProduct(Rd, n);

            n = Point(0, 0, 1);

            Point p(R0 + Rd*t);
            if (p.x < lowestPoint.x || p.x > lowestPoint.x + side || p.y < lowestPoint.y || p.y > lowestPoint.y + side) {
                continue;
            }

            if (t>0 && t<tMin) {
                tMin = t;
            }
        }

        // yz planes := n = (1, 0, 0) and (-1, 0, 0)
        n = Point(-1, 0, 0);
        for (int i=0; i<2; i++) {
            if (dotProduct(Rd, n) == 0) {
                continue;
            }

            D = lowestPoint.x + i*side;
            t = (D - dotProduct(R0, n)) / dotProduct(Rd, n);

            n = Point(1, 0, 0);

            Point p(R0 + Rd*t);
            if (p.y < lowestPoint.y || p.y > lowestPoint.y + side || p.z < lowestPoint.z || p.z > lowestPoint.z + side) {
                continue;
            }

            if (t>0 && t<tMin) {
                tMin = t;
            }
        }

        // xz planes := n = (0, 1, 0) and (0, -1, 0)
        n = Point(0, -1, 0);
        for (int i=0; i<2; i++) {
            if (dotProduct(Rd, n) == 0) {
                continue;
            }

            D = lowestPoint.y + i*side;
            t = (D - dotProduct(R0, n)) / dotProduct(Rd, n);
            
            n = Point(0, 1, 0);

            Point p(R0 + Rd*t);
            if (p.x < lowestPoint.x || p.x > lowestPoint.x + side || p.z < lowestPoint.z || p.z > lowestPoint.z + side) {
                continue;
            }

            if (t>0 && t<tMin) {
                tMin = t;
            }
        }

        if (tMin == INF) {
            return -1;
        }

        return tMin;
    }
};

struct Pyramid {
    Point lowest, top; // lowest := lower left corner of base
    vector<Point> corners;
    double width, height;
    Component component;

    Pyramid(const Point lowest=Point(), double width=0.0, double height=0.0, const Component component=Component()) {
        this->lowest = lowest;
        this->width = width;
        this->height = height;
        this->component = component;
    }

    Pyramid(const Pyramid &p) {
        this->lowest = p.lowest;
        this->width = p.width;
        this->height = p.height;
        this->component = p.component;
    }

    void draw() {
        double x = lowest.x, y = lowest.y, z = lowest.z;
        // top = Point(x + width/2.0, y + width/2.0, z + height);
        this->evaluateVariables();

        double rgb[3] = {component.color.r, component.color.g, component.color.b};
        glColor3f(rgb[0], rgb[1], rgb[2]);
        
        // draw base square
        glBegin(GL_QUADS);
        {
            for (int i=0; i<4; i++) {
                glVertex3f(corners[i].x, corners[i].y, corners[i].z);
            }
        }
        glEnd();

        // draw triangles
        for (int i=0; i<4; i++) {
            glBegin(GL_TRIANGLES);
            {
                glVertex3f(corners[i].x, corners[i].y, corners[i].z);
                glVertex3f(top.x, top.y, top.z);
                glVertex3f(corners[(i+1)%4].x, corners[(i+1)%4].y, corners[(i+1)%4].z);
            }
            glEnd();
        }
    }

    void evaluateVariables() {
        this->top = Point(lowest.x + width/2.0, lowest.y + width/2.0, lowest.z + height);

        this->corners.push_back(Point(lowest.x, lowest.y, lowest.z));
        this->corners.push_back(Point(lowest.x + width, lowest.y, lowest.z));
        this->corners.push_back(Point(lowest.x + width, lowest.y + width, lowest.z));
        this->corners.push_back(Point(lowest.x, lowest.y + width, lowest.z));
    }

    double determinant(double matrix[3][3]) {
        double v1 = matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]);
        double v2 = matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]);
        double v3 = matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);
        return v1 - v2 + v3;
    }

    double getTriangleT(Point a, Point b, Point c, Ray ray) {
        Point R0 = ray.origin, Rd = ray.dir;
        
        double AMat[3][3] {
            {a.x - b.x, a.x - c.x, Rd.x},
            {a.y - b.y, a.y - c.y, Rd.y},
            {a.z - b.z, a.z - c.z, Rd.z}
        };
        double betaMat[3][3] = {
            {a.x - R0.x, a.x - c.x, Rd.x},
            {a.y - R0.y, a.y - c.y, Rd.y},
            {a.z - R0.z, a.z - c.z, Rd.z}
        };
        double gammaMat[3][3] = {
            {a.x - b.x, a.x - R0.x, Rd.x},
            {a.y - b.y, a.y - R0.y, Rd.y},
            {a.z - b.z, a.z - R0.z, Rd.z}
        };
        double tMat[3][3] = {
            {a.x - b.x, a.x - c.x, a.x - R0.x},
            {a.y - b.y, a.y - c.y, a.y - R0.y},
            {a.z - b.z, a.z - c.z, a.z - R0.z}
        };

        double A = determinant(AMat);
        if (A == 0) {
            return -1;
        }

        double beta = determinant(betaMat) / A;
        double gamma = determinant(gammaMat) / A;
        double t = determinant(tMat) / A;

        if (beta>0 && gamma>0 && beta+gamma<1 && t>0) {
            return t;
        }

        return -1;
    }

    double getIntersectingT(Ray ray) {
        /*
            1. Base parallal to xy plane
            2. Four triangles
        */
        evaluateVariables();
        Point R0 = ray.origin, Rd = ray.dir;
        double tMin = INF, D, t, beta, gamma, A;

        // base := xy plane
        Point n(0, 0, 1);
        if (dotProduct(Rd, n) != 0) {
            D = - (lowest.z);
            t = - (dotProduct(R0, n) + D) / dotProduct(Rd, n);

            Point p(R0 + Rd*t);
            if (p.x >= lowest.x && p.x <= lowest.x + width && p.y >= lowest.y && p.y <= lowest.y + width) {
                if (t>0 && t<tMin) {
                    tMin = t;
                }
            }
        }

        // four triangles
        for (int i=0; i<4; i++) {
            t = getTriangleT(corners[i], corners[(i+1)%4], top, ray);

            if (t>0 && t<tMin) {
                // std::cout << ":') ";
                tMin = t;
            }
        }

        if (tMin == INF) {
            return -1; 
        }
        return tMin;
    }
};

struct Sphere {
    Point center;
    double radius;
    Component component;

    Sphere(Point center=Point(), double radius=0.0) {
        this->center = center;
        this->radius = radius;
        this->component = Component();
    }

    Sphere(const Sphere &s) {
        this->center = s.center;
        this->radius = s.radius;
        this->component = s.component;
    }

    Point getNormal(Point point) {
        return point-center;
    }

    void draw() {
        double rgb[3] = {component.color.r, component.color.g, component.color.b};
        glColor3f(rgb[0], rgb[1], rgb[2]);
        
        int slices = 100, stacks = 100;
        glPushMatrix();
            glTranslatef(center.x, center.y, center.z);
            glutSolidSphere(radius, slices, stacks);
        glPopMatrix();
    }

    double getIntersectingT(Ray ray) {
        Point R0 = ray.origin - this->center, Rd = ray.dir;

        /*
            Rd.Rd t^2 + 2Rd.R0 t + R0.R0 - r^2 = 0

            a = Rd.Rd, b = 2Rd.R0, c = R0.R0 - r^2
            a = 1 [Rd is normalized]

            D = 4(Rd.R0)^2 - 4(R0.R0 - r^2)
              = 4((Rd.R0)^2 - (R0.R0 - r^2))

            sqrt(D) = 2*sqrt((Rd.R0)^2 - (R0.R0 - r^2))

        */

        double b = dotProduct(R0, Rd);
        double c = dotProduct(R0, R0) - radius*radius;
        double D = b*b - c; // discriminent (a = 1)

        if (D < 0) {
            return -1;
        }

        double t1 = -b+sqrt(D), t2 = -b-sqrt(D);

        if (t2 < 0) {
            if (t1 > 0) {
                return t1;
            } else {
                return -INF;
            }
        }

        return t2;
    }

    double intersect(Ray ray, Color &color, int level, int recursionLevel, vector<Sphere> spheres, vector<Cube> cubes, vector<Pyramid> pyramids, vector<Light> lights, vector<Spotlight> spotlights) {
        double t = getIntersectingT(ray);

        if (t < 0) {
            return -1;
        }

        // if (level >= recursionLevel) {
        //     return t;
        // }

        // std::cout << level << " " << ray.origin << " " << ray.dir << " " << t << "\n";

        if (level == 0) {
            return t;
        }

        Point intersectionPoint = ray.origin + ray.dir*t;
        
        double multiplier = component.ambient;
        double lambert = 0.0, phong = 0.0;

        // Normal lights
        for (Light light : lights) {
            Point lightPos = light.position;
            Point lightDir = (intersectionPoint - lightPos).normalize();
            Point toSource = (lightPos - intersectionPoint).normalize();

            Ray lightRay = Ray(lightPos, lightDir);
            // Ray ray = Ray(intersectionPoint, toSource);
            Point normal = getNormal(intersectionPoint);

            // check if incident ray is obstructed by any object
            double d = (intersectionPoint - lightPos).length();
            if (d < 1e-5) {
                continue;
            }

            bool notObscured = true;
            if (notObscured) {
                for (Sphere sphere : spheres) {
                    double t3 = sphere.getIntersectingT(lightRay);
                    if (t3 > 0 && t3 < d) {
                        notObscured = false;
                        break;
                    }
                }
            }
            if (notObscured) {
                for (Cube cube : cubes) {
                    double t3 = cube.getIntersectingT(lightRay);
                    if (t3 > 0 && t3 < d) {
                        notObscured = false;
                        break;
                    }
                }
            }
            if (notObscured) {
                for (Pyramid pyramid : pyramids) {
                    double t3 = pyramid.getIntersectingT(lightRay);
                    if (t3 > 0 && t3 < d) {
                        notObscured = false;
                        break;
                    }
                }
            }

            if (notObscured) {
                double scale = exp(-d*d*light.falloff);
                lambert += (dotProduct(toSource, normal)*scale);

                Point Rd = (intersectionPoint - lightPos).normalize();
                Point reflection = (Rd - normal*2*dotProduct(Rd, normal)).normalize();
                phong += pow(dotProduct(reflection, toSource), component.shininess)*scale;
            }
        }

        // Spotlights
        for (Spotlight spotlight : spotlights) {
            Point lightPos = spotlight.position;
            Point lightDir = (intersectionPoint - lightPos).normalize();
            Point toSource = (lightPos - intersectionPoint).normalize();

            Ray lightRay = Ray(lightPos, lightDir);
            // Ray ray = Ray(intersectionPoint, toSource);
            Point normal = getNormal(intersectionPoint);
            double angle = acos(dotProduct(lightDir, normal))*180/PI;

            if (fabs(angle) > spotlight.cutoff) {
                continue;
            }

            // check if incident ray is obstructed by any object
            double d = (intersectionPoint - lightPos).length();
            if (d < 1e-5) {
                continue;
            }

            bool notObscured = true;
            if (notObscured) {
                for (Sphere sphere : spheres) {
                    double t3 = sphere.getIntersectingT(lightRay);
                    if (t3 > 0 && t3 < d) {
                        notObscured = false;
                        break;
                    }
                }
            }
            if (notObscured) {
                for (Cube cube : cubes) {
                    double t3 = cube.getIntersectingT(lightRay);
                    if (t3 > 0 && t3 < d) {
                        notObscured = false;
                        break;
                    }
                }
            }
            if (notObscured) {
                for (Pyramid pyramid : pyramids) {
                    double t3 = pyramid.getIntersectingT(lightRay);
                    if (t3 > 0 && t3 < d) {
                        notObscured = false;
                        break;
                    }
                }
            }

            if (notObscured) {
                double scale = exp(-d*d*spotlight.falloff);
                lambert += (dotProduct(toSource, normal)*scale);

                Point Rd = (intersectionPoint - lightPos).normalize();
                Point reflection = (Rd - normal*2*dotProduct(Rd, normal)).normalize();
                phong += pow(dotProduct(reflection, toSource), component.shininess)*scale;
            }
        }

        multiplier += component.diffuse*lambert + component.specular*phong;

        component.color.r *= multiplier;
        component.color.g *= multiplier;
        component.color.b *= multiplier;

        // Reflection
        if (level < recursionLevel) {
            Point normal = getNormal(intersectionPoint);
            Ray reflectedRay = Ray(intersectionPoint, ray.dir - normal*2*dotProduct(ray.dir, normal));

            reflectedRay.origin = reflectedRay.origin + reflectedRay.dir*1e-5;
            double t3 = -1, tMin = 1000007;

            

            for (Sphere sphere : spheres) {
                Color colorTemp;
                double t3 = sphere.intersect(reflectedRay, color, level+1, recursionLevel, spheres, cubes, pyramids, lights, spotlights);
                if (t3 > 0) {
                    component.color.r += colorTemp.r*sphere.component.reflection;
                    component.color.g += colorTemp.g*sphere.component.reflection;
                    component.color.b += colorTemp.b*sphere.component.reflection;
                }
            }
        }

        if (component.color.r > 1.0) {
            component.color.r = 1.0;
        }
        if (component.color.g > 1.0) {
            component.color.g = 1.0;
        } 
        if (component.color.b > 1.0) {
            component.color.b = 1.0;
        }

        return t;
    }
};

struct Flooring {
    double width, ambient, diffuse, reflection;
    bitmap_image texture_w;
    bitmap_image texture_b;

    Flooring(double width=0.5, double ambient=0.0, double diffuse=0.0, double reflection=0.0) {
        this->width = width;
        this->ambient = ambient;
        this->diffuse = diffuse;
        this->reflection = reflection;
        this->texture_b = bitmap_image("texture_b.bmp");
        this->texture_w = bitmap_image("texture_w.bmp");
    }

    Ray getNormal(Ray incidentRay, Point point) {
        if(incidentRay.dir.z > 0) {
            return Ray(point, Point(0, 0, 1));
        }
        else {
            return Ray(point, Point(0, 0, -1));
        }
    }

    void draw() {
        int noOfTiles = INF;
        double startX = -width*INF/2, startY = -width*INF/2;
        
        for (int i = 0; i < noOfTiles; i++) {
            for (int j = 0; j < noOfTiles; j++) {
                int isBlack = (i+j) % 2;
                int rgb[3] = {isBlack, isBlack, isBlack};
                
                glColor3f(rgb[0], rgb[1], rgb[2]);

                glBegin(GL_QUADS);
                {
                    glVertex3f(startX + i * width, startY + j * width, 0);
                    glVertex3f(startX + (i+1) * width, startY + j * width, 0);
                    glVertex3f(startX + (i+1) * width, startY + (j+1) * width, 0);
                    glVertex3f(startX + i * width, startY + (j+1) * width, 0);
                }
                glEnd();
            }
        }
    }

    Color getTextureColor(const Point& p, const bitmap_image& image) { 
        unsigned char r, g, b;
        int x = std::ceil(p.x / width), y = std::ceil(p.y / width); 
        
        double dx = -(x-1) * width + p.x, dy = -p.y + (y) * width;
        x = dx * image.width() / width;
        y = dy * image.height() / width;
        
        image.get_pixel(x, y, r, g, b);

        return Color(r/255.0, g/255.0, b/255.0);
    }

    Color getColorAt(Ray ray, double t) {
        Point p = ray.origin + ray.dir*t;

        double start = -width*INF/2, end = width*INF/2;
        double x, y;
        x = p.x; y = p.y; 
        
        if (p.x < 0) {
            x -= start;
        }
        if (p.y < 0) {
            y -= start;
        }

        if (x < start || x > end || y < start || y > end) {
            return Color(0, 0, 0);
        }

        bool isWhite = !((int) (x/width) % 2) == ((int) (y/width) % 2);
        if (isWhite) {
            return getTextureColor(p, texture_w);
        } else {
            return getTextureColor(p, texture_b);
        }

        int color = isWhite;
        return Color(color*ambient, color*ambient, color*ambient);
    }

    double getIntersectingT(Ray ray) {
        /*
            P := (x, y, z) == (x, y, 0) [xy plane, z=0]
            n := (A, B, C) == (0, 0, 1) [floor is parallel to xy plane]
            A = 0, B = 0, C = 1

            H(P) := n.P + D = 0
            => Ax + By + Cz + D = 0
            => z + D = 0;
            => D = 0 [floor is at z=0]
        */

        Point n(0, 0, 1);
        double t = - dotProduct(n, ray.origin) / dotProduct(n, ray.dir);

        if (t < 0) {
            return -INF;
        }

        return t;
    }
}; 


