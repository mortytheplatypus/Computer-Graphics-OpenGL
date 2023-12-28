#include <iostream>
#include <time.h>
#include <math.h>
#include <GL/glut.h>
using namespace std;

#define CENTER_X 0
#define CENTER_Y 0.3
#define HOUR_HAND 1
#define MINUTE_HAND 2
#define SECOND_HAND 3

struct Time {
    int hour;
    int minute;
    int second;
} current_time;

void display();
void idle();
void timer(int);

void drawCircle(float);
void drawClockCircle();
void drawHand(int);
void drawClockHands();
void drawOctagon(float, float);
void drawHexadecagon(float, float);
void drawPendulum();
void drawClockBody();
void calculateThetas();
void updateTime();

const float PI = acos(-1.0), TIMER_INTERVAL = 1000.0 / 360.0;
float theta_s, theta_m, theta_h;
float theta_p = 45.0, omega = 0.0, g = 9.8, L = 0.3;

int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitWindowSize(680, 680);
    glutInitWindowPosition(100, 50);
    glutCreateWindow("Clock");
    glutDisplayFunc(display);
    glutTimerFunc(0, timer, 0);
    glutMainLoop();

    return 0;
}

void display() {
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity(); 

    glPushMatrix();
        drawClockBody();
        drawClockCircle();
        drawClockHands();
        drawPendulum();
    glPopMatrix();

    glFlush();
    glutSwapBuffers(); 
}

void timer(int value) {
    float theta = - g / L * sin(theta_p * PI / 180.0);
    omega += theta * 0.01;
    theta_p += omega * 0.01;
    
    glutPostRedisplay();
    glutTimerFunc(TIMER_INTERVAL, timer, 0);
}

void drawClockCircle() {
    drawCircle(0.61);
    drawCircle(0.57);
    drawCircle(0.53);

    glColor3f(245.0/256.0, 222.0/256.0, 179.0/256.0);
    glLineWidth(10);
    for (int i=-1; i < 2; i+=2) {
        glBegin(GL_LINES);
            glVertex2f(CENTER_X + i*0.49, CENTER_Y);
            glVertex2f(CENTER_X + i*0.35, CENTER_Y);
        glEnd();

        glBegin(GL_LINES);
            glVertex2f(CENTER_X, CENTER_Y + i*0.49);
            glVertex2f(CENTER_X, CENTER_Y + i*0.35);
        glEnd();
    }

    float theta;
    for (int i=1; i<12; i++) {
        if (i%3 == 0) {
            continue;
        }

        glBegin(GL_QUADS);
            theta = i * PI / 6.0  + 0.02;
            glVertex2f(CENTER_X + 0.41 * cos(theta), CENTER_Y + 0.41 * sin(theta));
            glVertex2f(CENTER_X + 0.49 * cos(theta), CENTER_Y + 0.49 * sin(theta));
            theta = theta - 0.04;
            glVertex2f(CENTER_X + 0.49 * cos(theta), CENTER_Y + 0.49 * sin(theta));
            glVertex2f(CENTER_X + 0.41 * cos(theta), CENTER_Y + 0.41 * sin(theta));
        glEnd();
    }
    glColor3f(1.0, 1.0, 1.0);
}

void drawCircle(float radius) {
    const int no_of_points = 1200;
    float x, y, theta;

    glColor3f(111.0/256.0, 78.0/256.0, 55.0/256.0);
    glPointSize(5);
    glBegin(GL_POINTS);
        for(int i=0; i < no_of_points; i++) {
            theta = i * 2.0 * PI / no_of_points;

            x = CENTER_X + radius * cos(theta);
            y = CENTER_Y + radius * sin(theta);
            
            glVertex3f(x, y, 0);
        }
    glEnd();
    glColor3f(1.0, 1.0, 1.0);
}

void drawClockHands() {
    updateTime();
    float x, y;

    // printf("%02d:%02d:%02d\n", current_time.hour, current_time.minute, current_time.second);

    calculateThetas();

    // hours hand
    glColor3f(1.0, 0.0, 0.0);
    drawHand(HOUR_HAND);

    // minutes hand
    glColor3f(0.8, 0.8, 0.0);
    drawHand(MINUTE_HAND);

    // seconds hand
    glColor3f(150.0/256.0, 105.0/256.0, 25.0/256.0);
    drawHand(SECOND_HAND);
    
    glColor3f(0.5, 0.5, 0.5);
    drawHexadecagon(CENTER_X, CENTER_Y);

    glColor3f(1.0, 1.0, 1.0);
}

void drawHand(int hand_type) {
    float theta, length, angle;
    if (hand_type == HOUR_HAND) {
        theta = theta_h;
        length = 0.25;
        angle = 0.35;
    } else if (hand_type == MINUTE_HAND) {
        theta = theta_m;
        length = 0.35;
        angle = 0.15;
    } else if (hand_type == SECOND_HAND) {
        theta = theta_s;
        length = 0.45;
        angle = 0.05;
    }

    float l = 2.0 * length / 5.0, x, y;

    glBegin(GL_QUADS);
        glVertex2f(CENTER_X, CENTER_Y);
        
        if (hand_type == HOUR_HAND) {
            glColor3f(0.5, 0.5, 0.0);
        } else if (hand_type == MINUTE_HAND) {
            glColor3f(124.0/256.0, 48.0/256.0, 48.0/256.0);
        } else if (hand_type == SECOND_HAND) {
            glColor3f(240.0/256.0, 220.0/256.0, 100.0/256.0);
        }

        x = CENTER_X + l * cos(theta + angle);
        y = CENTER_Y + l * sin(theta + angle);
        glVertex2f(x, y);

        if (hand_type == HOUR_HAND) {
            glColor3f(1.0, 1.0, 0.0);
        } else if (hand_type == MINUTE_HAND) {
            glColor3f(0.5, 0.0, 0.0);
        } else if (hand_type == SECOND_HAND) {
            glColor3f(229.0/256.0, 170.0/256.0, 112.0/256.0);
        }

        x = CENTER_X + length * cos(theta);
        y = CENTER_Y + length * sin(theta);
        glVertex2f(x, y);

        if (hand_type == HOUR_HAND) {
            glColor3f(0.5, 0.5, 0.0);
        } else if (hand_type == MINUTE_HAND) {
            glColor3f(124.0/256.0, 48.0/256.0, 48.0/256.0);
        } else if (hand_type == SECOND_HAND) {
            glColor3f(194.0/256.0, 42.0/256.0, 42.0/256.0);
        }

        x = CENTER_X + l * cos(theta - angle);
        y = CENTER_Y + l * sin(theta - angle);
        glVertex2f(x, y);
    glEnd();

    glColor3f(1.0, 1.0, 1.0);
}

void drawOctagon(float cx, float cy) {
    float r = 0.075;
    glColor3f(0.7, 0.7, 0.2);
    glBegin(GL_POLYGON);
        for (int i=0; i<8; i++) {
            glVertex2f(cx + r * cos(i * PI / 4.0), cy + r * sin(i * PI / 4.0));
        }
    glEnd();
    glColor3f(1.0, 1.0, 1.0);
}

void drawHexadecagon(float cx, float cy) {
    float r = 0.04;
    glColor3f(82.0/256.0, 64.0/256.0, 51.0/256.0);
    glBegin(GL_POLYGON);
        for (int i=0; i<16; i++) {
            glVertex2f(cx + r * cos(i * PI / 8.0), cy + r * sin(i * PI / 8.0));
        }
    glEnd();
    glColor3f(1.0, 1.0, 1.0);
}

void drawPendulum() {
    const float cx = CENTER_X, cy = CENTER_Y - 0.75;

    glLineWidth(10);
    glPushMatrix();
        glTranslated(cx, cy, 0);
        glRotated(theta_p, 0.0f, 0.0f, 1.0f);
        glBegin(GL_LINES);
            glColor3f(0.9, 0.2, 0.2);
            glVertex2f(0.0, 0.0);
            glColor3f(0.2, 0.8, 0.0);
            glVertex2f(0.0, -L);
            glColor3f(0.3, 0.3, 0.8);
        glEnd();
        drawOctagon(0.0, -L);
    glPopMatrix();
    
    drawHexadecagon(cx, cy);
}

void drawClockBody() {
    glColor3f(245.0/256.0, 222.0/256.0, 179.0/256.0);

    glLineWidth(10);

    glBegin(GL_LINES);
        glVertex2f(-0.5, 0.95);
        glVertex2f(0.5, 0.95);

        glVertex2f(0.5, 0.965);
        glVertex2f(0.8, 0.05);

        glVertex2f(0.8, 0.05);
        glVertex2f(0.5, -0.865);

        glVertex2f(0.5, -0.85);
        glVertex2f(-0.5, -0.85);

        glVertex2f(-0.5, -0.865);
        glVertex2f(-0.8, 0.05);

        glVertex2f(-0.8, 0.05);
        glVertex2f(-0.5, 0.965);
    glEnd();
}

void calculateThetas() {
    theta_s = 2.0 * PI * current_time.second / 60.0;
    theta_m = 2.0 * PI * current_time.minute / 60.0;
    theta_h = 2.0 * PI * current_time.hour / 12.0;
    theta_h += ((current_time.minute * PI) / 180.0); // pi_by_three * (min / 60)

    theta_s -= (PI / 2.0);
    theta_m -= (PI / 2.0);
    theta_h -= (PI / 2.0);

    theta_s *= (-1);
    theta_m *= (-1);
    theta_h *= (-1);
}

void updateTime() {
    time_t t = time(NULL);
    struct tm tm = *localtime(&t);

    current_time.hour = tm.tm_hour;
    current_time.minute = tm.tm_min;
    current_time.second = tm.tm_sec;

    if (current_time.hour > 12) {
        current_time.hour -= 12;
    }
}

