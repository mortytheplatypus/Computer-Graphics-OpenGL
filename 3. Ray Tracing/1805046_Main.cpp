#include <bits/stdc++.h>
#include <GL/glut.h>
#include "1805046_Header.h"
using namespace std;

#define WINDOW_WIDTH 700
#define WINDOW_HEIGHT 600
#define WINDOW_X 100
#define WINDOW_Y 50

double fovY, fovX, aspect, zNear, zFar, height, width, dh, dw;
int recursionLevel, screenSize, objectCount, lightCount, spotlightCount;
Point midpoint;

Flooring flooring;
vector<Cube> cubes;
vector<Pyramid> pyramids;
vector<Sphere> spheres;
vector<Light> lights;
vector<Spotlight> spotlights;
vector< vector<Point> > pointBuffer;
vector< vector<Ray> > rays;

Point eye(0, -145, 50), k_(0, 0, 1);
Point look = (Point(0, 0, 0) - eye).normalize();
Point rightV = crossProduct(look, k_).normalize();
Point up = crossProduct(rightV, look).normalize();

void initialize();
void calculatePointbuffer();
void parseInput();
void keyboardListener(unsigned char , int, int);
void specialKeyListener(int, int, int);
void display();
void animate();
void captureBMP(string);

int main(int argc, char **argv){
	glutInit(&argc, argv);
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
    glutInitWindowPosition(WINDOW_X, WINDOW_Y);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutCreateWindow("Ray Tracing");

    parseInput();

    glutDisplayFunc(display);
	glutIdleFunc(animate);

    glutKeyboardFunc(keyboardListener);
	glutSpecialFunc(specialKeyListener);

    initialize();
    glutMainLoop();

	return 0;
}

void initialize() {
	glClearColor(0, 0, 0, 0);

	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();
    glEnable(GL_DEPTH_TEST);

	gluPerspective(fovY, aspect, zNear, zFar);
}

void keyboardListener(unsigned char key, int x, int y) {
    /* 
        Controls:
            0 := ray tracing and output bmp
            1 := rotate/look left
            2 := rotate/look right
            3 := look up
            4 := look down
            5 := tilt counterclockwise
            6 := tilt clockwise
    */

    double dt = PI/180; // d_theta

    if (key == '0') {
        captureBMP("output.bmp");
    } else if (key == '1') {
        rightV = rotate3D(rightV, up, -dt);
        look = rotate3D(look, up, -dt);
    } else if (key == '2') {
        rightV = rotate3D(rightV, up, dt);
        look = rotate3D(look, up, dt);
    } else if (key == '3') {
        up = rotate3D(up, rightV, -dt);
        look = rotate3D(look, rightV, -dt);
    } else if (key == '4') {
        up = rotate3D(up, rightV, dt);
        look = rotate3D(look, rightV, dt);
    } else if (key == '5') {
        rightV = rotate3D(rightV, look, -dt);
        up = rotate3D(up, look, -dt);
    } else if (key == '6') {
        rightV = rotate3D(rightV, look, dt);
        up = rotate3D(up, look, dt);
    }
}

void specialKeyListener(int key, int x, int y) {
    /*
        Contols: 
            up arrow    :=  move forward
            down arrow  :=  move backward
            right arrow :=  move right
            left arrow  :=  move left
            page up     :=  move up
            page down   :=  move down
    */

    double dr = 10; // dr

    if (key == GLUT_KEY_UP) {
        eye = eye + look * dr;
    } else if (key == GLUT_KEY_DOWN) {
        eye = eye - look * dr;
    } else if (key == GLUT_KEY_RIGHT) {
        eye = eye + rightV * dr;
    } else if (key == GLUT_KEY_LEFT) {
        eye = eye - rightV * dr;
    } else if (key == GLUT_KEY_PAGE_UP) {
        eye = eye + up * dr;
    } else if (key == GLUT_KEY_PAGE_DOWN) {
        eye = eye - up * dr;
    } 
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // control viewing (or camera)
    Point center = eye + look;
    gluLookAt(eye.x, eye.y, eye.z, center.x, center.y, center.z, up.x, up.y, up.z);

    drawAxes();
    flooring.draw();

    // for (Cube cube : cubes) {
    //     cube.draw();
    // } 

    // for (Pyramid pyramid : pyramids) {
    //     pyramid.draw();
    // }

    for (Sphere sphere : spheres) {
        sphere.draw();
    }

    for (Light light : lights) {
        light.draw();
    }

    for (Spotlight spotlight : spotlights) {
        spotlight.draw();
    }

    glFlush();
    glutSwapBuffers(); 
}

void animate() {
	glutPostRedisplay();
}

void parseInput() {
    ifstream fin("description.txt");

    fin >> zNear >> zFar >> fovY >> aspect; 
    fovX = aspect*fovY;

    fin >> recursionLevel;
    fin >> screenSize;

    fin >> flooring.width;
    fin >> flooring.ambient >> flooring.diffuse >> flooring.reflection;

    fin >> objectCount;
    for (int i=0; i<objectCount; i++) {
        string objectType; 
        fin >> objectType;
        
        if (objectType == "cube") {
            Cube cube;

            fin >> cube.lowestPoint.x >> cube.lowestPoint.y >> cube.lowestPoint.z;
            fin >> cube.side;
            fin >> cube.component.color.r >> cube.component.color.g >> cube.component.color.b;
            fin >> cube.component.ambient >> cube.component.diffuse >> cube.component.specular >> cube.component.reflection;
            fin >> cube.component.shininess;

            cubes.push_back(cube);
        } else if (objectType == "pyramid") {
            Pyramid pyramid;

            fin >> pyramid.lowest.x >> pyramid.lowest.y >> pyramid.lowest.z;
            fin >> pyramid.width >> pyramid.height;
            fin >> pyramid.component.color.r >> pyramid.component.color.g >> pyramid.component.color.b;
            fin >> pyramid.component.ambient >> pyramid.component.diffuse >> pyramid.component.specular >> pyramid.component.reflection;
            fin >> pyramid.component.shininess;

            pyramids.push_back(pyramid);
        } else if (objectType == "sphere") {
            Sphere sphere;
            fin >> sphere.center.x >> sphere.center.y >> sphere.center.z;
            fin >> sphere.radius;
            fin >> sphere.component.color.r >> sphere.component.color.g >> sphere.component.color.b;
            fin >> sphere.component.ambient >> sphere.component.diffuse >> sphere.component.specular >> sphere.component.reflection;
            fin >> sphere.component.shininess;

            spheres.push_back(sphere);
        } else {
            std::cout << "Unexpected error while taking input\n";
            return;
        }
    }

    fin >> lightCount;
    for (int i=0; i<lightCount; i++) {
        Light light;

        fin >> light.position.x >> light.position.y >> light.position.z;
        fin >> light.falloff;

        lights.push_back(light);
    }

    fin >> spotlightCount;
    for (int i=0; i<spotlightCount; i++) {
        Spotlight spotlight;

        fin >> spotlight.position.x >> spotlight.position.y >> spotlight.position.z;
        fin >> spotlight.falloff;
        fin >> spotlight.look.x >> spotlight.look.y >> spotlight.look.z >> spotlight.cutoff;

        spotlights.push_back(spotlight);
    }
}

void calculatePointbuffer() {
    height = 2*zNear*tan(fovY/2.0); 
    width = 2*zNear*tan(fovX/2.0); 

    dh = height / screenSize; 
    dw = width / screenSize; 

    midpoint = eye + (look*zNear);

    pointBuffer.resize(screenSize);
    rays.resize(screenSize);
    for (int i=0; i<screenSize; i++) {
        pointBuffer[i].resize(screenSize);
        rays[i].resize(screenSize);
    }

    Point bottomRightScreen = midpoint + (rightV*width/2) - (up*height/2);
    double x = bottomRightScreen.x, y = bottomRightScreen.y, z = bottomRightScreen.z;

    for (int i=0; i<screenSize; i++) {
        for (int j=0; j<screenSize; j++) {
            pointBuffer[i][j] = bottomRightScreen - (rightV*j*dw) + (up*i*dh);
            rays[i][j] = Ray(pointBuffer[i][j], (pointBuffer[i][j]-eye));
        }
    }
}

void captureBMP(string fileName) {
    bitmap_image image(screenSize, screenSize);
    std::cout << "Capturing image...\n";

    calculatePointbuffer();

    for (int i=0; i<screenSize; i++) {
        for (int j=0; j<screenSize; j++) {
            Ray ray = Ray(pointBuffer[i][j], (pointBuffer[i][j]-eye).normalize()); 

            //check for intersection
            double tMin = 1000007, t = tMin+1;
            ObjectType objectType = NONE;
            Color color;

            // spheres
            for (Sphere sphere : spheres) {
                // double t = sphere.getIntersectingT(ray);
                double t = sphere.intersect(ray, color, 1, recursionLevel, spheres, cubes, pyramids, lights, spotlights);
                
                // Cube cube = cubes[0];
                if (t > 0 && t < tMin) {
                    tMin = t;
                    objectType = SPHERE;
                    color = sphere.component.color;
                }
            }

            // cubes
            for (Cube cube : cubes) {
                double t = cube.getIntersectingT(ray);
                if (t > 0 && t < tMin) {
                    tMin = t;
                    objectType = CUBE;
                    color = cube.component.color;
                }
            }

            // pyramids
            for (Pyramid pyramid : pyramids) {
                double t = pyramid.getIntersectingT(ray);
                if (t > 0 && t < tMin) {
                    tMin = t;
                    objectType = PYRAMID;
                    color = pyramid.component.color;
                }
            }

            // flooring
            t = flooring.getIntersectingT(ray);
            if (t > 0 && t < tMin) {
                tMin = t;
                color = flooring.getColorAt(ray, t);
            }

            image.set_pixel(j, i, 255*color.r, 255*color.g, 255*color.b);
        }
    }

    image.save_image(fileName);
    std::cout << "Image saved as " << fileName << "\n";
}