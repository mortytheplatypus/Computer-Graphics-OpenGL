#include <math.h>
#include <vector>
#include <GL/glut.h>

void initialize();
void display();
void reshape(int, int);
void keyboard(unsigned char, int, int);
void special(int, int, int);

void drawSpheres();
void drawSphereUtil();
void drawCylinderSegmentUtil();
void drawCylinderSegments();
void drawTriangleUtil();
void drawTriangles();

const float ROOT_TWO = sqrt(2.0), ROOT_THREE = sqrt(3.0), PI = acos(-1.0);
const float ANGLE = acos(-1.0) - acos(-1.0/3.0);
const int subdivision = 5;

float eyex = 4.0, eyey = 4.0, eyez = 4.0;
float centerx = 0, centery = 0, centerz = 0;
float upx = 0, upy = 1, upz = 0;
float rotation_axis = 0.0, rotation_x = 0.0, rotation_y = 0.0, rotation_z = 0.0;
float translate_x = 0.0, translate_y = 0.0, translate_z = 0.0;

float r = 0.5, dr = 0.01, scale = 2*r;

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(680, 680);
    glutCreateWindow("3D Octahedron");
    initialize();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);

    glutMainLoop();

    return 0;
}

void initialize() {
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glEnable(GL_DEPTH_TEST);
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glTranslatef(0.0f, 0.0f, -5.0f);

    gluLookAt(eyex, eyey, eyez,
                centerx, centery, centerz,
                upx, upy, upz);

    glTranslated(translate_x, translate_y, translate_z);

    glRotated(rotation_x, 1, 0, -1);
    glRotated(rotation_y, -1, 1, 0);
    glRotated(rotation_z, 1, 0, 1);

    glRotated(rotation_axis, 0, 1, 0);

    glScaled(2.0, 2.0, 2.0);
    drawSpheres();
    drawCylinderSegments();
    drawTriangles();
    
    glutSwapBuffers();
}

void reshape(int height, int width) {
    if (height == 0)
        height = 1;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (float)width / height, 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
}

void keyboard(unsigned char key, int x, int y) {
    if (key == '.') {
        r -= dr;
        if (r <= 0) {
            r = 0;
        }
    } else if (key == ',') {
        r += dr;
        if (r >= 1/ROOT_TWO) {
            r = 1/ROOT_TWO;
        }
    } else if (key == 'a') {
        rotation_axis -= 2;
    } else if (key == 'd') {
        rotation_axis += 2;
    } else if (key == '1') {
        rotation_y += 1;
    } else if (key == '2') {
        rotation_y -= 1;
    } else if (key == '3') {
        rotation_x += 1;
    } else if (key == '4') {
        rotation_x -= 1;
    } else if (key == '5') {
        rotation_z += 1;
    } else if (key == '6') {
        rotation_z -= 1;
    }
    
    glutPostRedisplay();
}

void special(int key, int x, int y) {
    float PI_BY_FOUR = PI / 4.0, r = 0.1;
    if (key == GLUT_KEY_UP) {
        translate_z += r * cos(PI_BY_FOUR);
        translate_x += r * sin(PI_BY_FOUR);
    } else if (key == GLUT_KEY_DOWN) {
        translate_z -= r * cos(PI_BY_FOUR);
        translate_x -= r * sin(PI_BY_FOUR);
    } else if (key == GLUT_KEY_RIGHT) {
        translate_x += r * cos(PI_BY_FOUR);
        translate_z -= r * sin(PI_BY_FOUR);
    } else if (key == GLUT_KEY_LEFT) {
        translate_x -= r * cos(PI_BY_FOUR);
        translate_z += r * sin(PI_BY_FOUR);
    } else if (key == GLUT_KEY_PAGE_UP) {
        translate_y += r;
    } else if (key == GLUT_KEY_PAGE_DOWN) {
        translate_y -= r;
    }

    glutPostRedisplay();
}

void drawSphereUtil() {
    std::vector<float> vertices;
    float n1[3], n2[3], v[3]; // direction vector intersecting 2 planes, n1 x n2
    float a1, a2; // latitudinal angle along Z-axis

    int verticesPerRow = 1;
    for (int i=0; i<subdivision; i++) {
        verticesPerRow = verticesPerRow << 1;
    }
    verticesPerRow++;

    // rotate latitudinal plane from 45 to -45 degrees along Z-axis (top-to-bottom)
    for (int i=0; i < verticesPerRow; i++) {
        // normal for latitudinal plane: if latitude angle is 0, then normal vector of latitude plane is n2=(0,1,0); therefore, it is rotating (0,1,0) vector by latitude angle a2
        a2 = (45.0 - i * 90.0 / (verticesPerRow - 1));
        a2 *= (PI / 180.0);
        n2[0] = -sin(a2);
        n2[1] = cos(a2);
        n2[2] = 0;

        // rotate longitudinal plane from -45 to 45 along Y-axis (left-to-right)
        for(int j = 0; j < verticesPerRow; j++) {
            // normal for longitudinal plane: if longitude angle is 0, then normal vector of longitude is n1=(0,0,-1); therefore, it is rotating (0,0,-1) vector by longitude angle a1
            a1 = (-45.0 + j * 90.0 / (verticesPerRow - 1));
            a1 *= (PI / 180.0);
            n1[0] = -sin(a1);
            n1[1] = 0;
            n1[2] = -cos(a1);

            // n1 CROSS n2
            v[0] = n1[1] * n2[2] - n1[2] * n2[1];
            v[1] = n1[2] * n2[0] - n1[0] * n2[2];
            v[2] = n1[0] * n2[1] - n1[1] * n2[0];

            float scale_of_vector = 1 / sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
            v[0] *= scale_of_vector;
            v[1] *= scale_of_vector;
            v[2] *= scale_of_vector;

            vertices.push_back(v[0]);
            vertices.push_back(v[1]);
            vertices.push_back(v[2]);
        }
    }

    glBegin(GL_QUADS);
    for (int i=1; i<verticesPerRow; i++) {
        for (int j=1; j<verticesPerRow; j++) {
            glVertex3f(vertices[3*(verticesPerRow*i + j)], vertices[3*(verticesPerRow*(i) + j) + 1], vertices[3*(verticesPerRow*i + j)+2]);
            glVertex3f(vertices[3*(verticesPerRow*(i-1) + j)], vertices[3*(verticesPerRow*(i-1) + j) + 1], vertices[3*(verticesPerRow*(i-1) + j)+2]);
            glVertex3f(vertices[3*(verticesPerRow*(i-1) + j-1)], vertices[3*(verticesPerRow*(i-1) + j-1) + 1], vertices[3*(verticesPerRow*(i-1) + j-1)+2]);
            glVertex3f(vertices[3*(verticesPerRow*i + j-1)], vertices[3*(verticesPerRow*i + j-1) + 1], vertices[3*(verticesPerRow*i + j-1)+2]);
        }
    }
    glEnd();
}

void drawSpheres() {
    glPushMatrix();
        glColor3f(210.0/256.0, 4.0/256.0, 45.0/256.0);
        glRotated(0.0, 0.0, 0.0, 0.0);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glScaled(r, r, r);
        drawSphereUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(255.0/256.0, 172.0/256.0, 28.0/256.0);
        glRotated(180.0, 0.0, 0.0, 1.0);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glScaled(r, r, r);
        drawSphereUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(65.0/256.0, 105.0/256.0, 225.0/256.0);
        glRotated(90.0, 0.0, 1.0, 0.0);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glScaled(r, r, r);
        drawSphereUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(34.0/256.0, 139.0/256.0, 34.0/256.0);
        glRotated(-90.0, 0.0, 1.0, 0.0);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glScaled(r, r, r);
        drawSphereUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(251.0/256.0, 236.0/256.0, 93.0/256.0);
        glRotated(90.0, 0.0, 0.0, 1.0);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glScaled(r, r, r);
        drawSphereUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(250.0/256.0, 249.0/256.0, 246.0/256.0);
        glRotated(-90.0, 0.0, 0.0, 1.0);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glScaled(r, r, r);
        drawSphereUtil();
    glPopMatrix();
}

void drawCylinderSegmentUtil() {
    float h = ROOT_TWO - 2*r;
    double tempx = r, tempz = 0;
    double currx, curry;
    int segments = 200;
    
    glBegin(GL_QUADS);
        for (int i = -segments/2; i <= segments/2; i++) {
            double theta = i * ANGLE / segments;
            currx = r * cos(theta);
            curry = r * sin(theta);

            glVertex3f(currx, h, curry);
            glVertex3f(currx, 0, curry);

            glVertex3f(tempx, 0, tempz);
            glVertex3f(tempx, h, tempz);

            tempx = currx;
            tempz = curry;
        }
    glEnd();
}

void drawCylinderSegments() {
    glColor3f(92.0/256.0, 64.0/256.0, 51.0/256.0);

    // First 4
    glPushMatrix();
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glRotated(45, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    glPushMatrix();
        glRotated(180, 0, 0, 1);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glRotated(45, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    glPushMatrix();
        glTranslated(0, 1 - ROOT_TWO * r, 0);
        glRotated(135, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    glPushMatrix();
        glRotated(180, 0, 0, 1);
        glTranslated(0, 1 - ROOT_TWO * r, 0);
        glRotated(135, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    // Next 4
    glPushMatrix();
        glRotated(90, 0, 1, 0);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glRotated(45, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    glPushMatrix();
        glRotated(90, 0, 1, 0);
        glRotated(180, 0, 0, 1);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glRotated(45, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    glPushMatrix();
        glRotated(90, 0, 1, 0);
        glTranslated(0, 1 - ROOT_TWO * r, 0);
        glRotated(135, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    glPushMatrix();
        glRotated(90, 0, 1, 0);
        glRotated(180, 0, 0, 1);
        glTranslated(0, 1 - ROOT_TWO * r, 0);
        glRotated(135, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    // Last 4
    glPushMatrix();
        glRotated(90, 1, 0, 0);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glRotated(45, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    glPushMatrix();
        glRotated(90, 1, 0, 0);
        glRotated(180, 0, 0, 1);
        glTranslated(1 - ROOT_TWO * r, 0, 0);
        glRotated(45, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    glPushMatrix();
        glRotated(90, 1, 0, 0);
        glTranslated(0, 1 - ROOT_TWO * r, 0);
        glRotated(135, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    glPushMatrix();
        glRotated(90, 1, 0, 0);
        glRotated(180, 0, 0, 1);
        glTranslated(0, 1 - ROOT_TWO * r, 0);
        glRotated(135, 0, 0, 1);
        drawCylinderSegmentUtil();
    glPopMatrix();

    glColor3f(1.0, 1.0, 1.0);
}

void drawTriangleUtil() {    
    glBegin(GL_TRIANGLES);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.0, 1.0);
    glEnd();
}

void drawTriangles() {
    glPushMatrix();
        glColor3f(0.8, 0.8, 0.0);
        glTranslated(r/ROOT_THREE, r/ROOT_THREE, r/ROOT_THREE);
        glScaled(1-ROOT_TWO * r, 1-ROOT_TWO * r, 1-ROOT_TWO * r);
        drawTriangleUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(0.8, 0.8, 0.0);
        glRotated(180.0, 0, 1, 0);
        glTranslated(r/ROOT_THREE, r/ROOT_THREE, r/ROOT_THREE);
        glScaled(1-ROOT_TWO * r, 1-ROOT_TWO * r, 1-ROOT_TWO * r);
        drawTriangleUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(222.0/256.0, 49.0/256.0, 99.0/256.0);
        glRotated(90.0, 0, 1, 0);
        glTranslated(r/ROOT_THREE, r/ROOT_THREE, r/ROOT_THREE);
        glScaled(1-ROOT_TWO * r, 1-ROOT_TWO * r, 1-ROOT_TWO * r);
        drawTriangleUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(222.0/256.0, 49.0/256.0, 99.0/256.0);
        glRotated(-90.0, 0, 1, 0);
        glTranslated(r/ROOT_THREE, r/ROOT_THREE, r/ROOT_THREE);
        glScaled(1-ROOT_TWO * r, 1-ROOT_TWO * r, 1-ROOT_TWO * r);
        drawTriangleUtil();
    glPopMatrix();

    // Next 4
    glPushMatrix();
        glColor3f(0.8, 0.8, 0.0);
        glRotated(180.0, 0, 0, 1);
        glTranslated(r/ROOT_THREE, r/ROOT_THREE, r/ROOT_THREE);
        glScaled(1-ROOT_TWO * r, 1-ROOT_TWO * r, 1-ROOT_TWO * r);
        drawTriangleUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(0.8, 0.8, 0.0);
        glRotated(180.0, 0, 0, 1);
        glRotated(180.0, 0, 1, 0);
        glTranslated(r/ROOT_THREE, r/ROOT_THREE, r/ROOT_THREE);
        glScaled(1-ROOT_TWO * r, 1-ROOT_TWO * r, 1-ROOT_TWO * r);
        drawTriangleUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(222.0/256.0, 49.0/256.0, 99.0/256.0);
        glRotated(180.0, 0, 0, 1);
        glRotated(90.0, 0, 1, 0);
        glTranslated(r/ROOT_THREE, r/ROOT_THREE, r/ROOT_THREE);
        glScaled(1-ROOT_TWO * r, 1-ROOT_TWO * r, 1-ROOT_TWO * r);
        drawTriangleUtil();
    glPopMatrix();

    glPushMatrix();
        glColor3f(222.0/256.0, 49.0/256.0, 99.0/256.0);
        glRotated(180.0, 0, 0, 1);
        glRotated(-90.0, 0, 1, 0);
        glTranslated(r/ROOT_THREE, r/ROOT_THREE, r/ROOT_THREE);
        glScaled(1-ROOT_TWO * r, 1-ROOT_TWO * r, 1-ROOT_TWO * r);
        drawTriangleUtil();
    glPopMatrix();
    
}
