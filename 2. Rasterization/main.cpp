#include <iostream>
#include "utils.h"
#include "bitmap_image.hpp"
using namespace std;

double getMin(double, double, double); 
double getMax(double, double, double); 
Vector calculateXandZ(double, Point, Point, Point);
bitmap_image createBlackBMP();

void ModelingTransformation();
void ViewTransformation();
void ProjectionTransformation();
void ClippingAndScanConversionUsingZBuffer();

int no_of_triangles = 0;
Point eye, look, up;
double fovY, aspect, near, far;
int screenWidth, screenHeight;
double minX, maxX, minY, maxY;
double zMin = -1.0, zMax = 1.0;

int main() {
    ModelingTransformation(); 
    cout << "Stage 1 completed...\n\n";

    ViewTransformation();
    cout << "Stage 2 completed...\n\n";

    ProjectionTransformation();
    cout << "Stage 3 completed...\n\n";

    ClippingAndScanConversionUsingZBuffer();
    cout << "Stage 4 completed...\n\n";

    cout << "Rasterization completed.\n\n";
    return 0;
}

// Stage 1: Modeling Transformation
void ModelingTransformation() {
    ifstream fin("scene.txt");
    ofstream fout("stage1.txt");

    stack<Matrix> matrixStack;
    matrixStack.push(createIdentityMatrix());

    fin >> eye >> look >> up;
    fin >> fovY >> aspect >> near >> far;

    while (true) {
        string cmd;
        fin >> cmd;

        if (cmd == "triangle") {
            Matrix M = matrixStack.top();
            no_of_triangles++;
            
            Point p1, p2, p3;
            fin >> p1 >> p2 >> p3;

            p1 = multiplyMatrixVector(M, p1);
            p2 = multiplyMatrixVector(M, p2);
            p3 = multiplyMatrixVector(M, p3);

            fout << p1 << "\n" << p2 << "\n" << p3 << "\n\n";
        } else if (cmd == "translate") {
            Point p; fin >> p;
            Matrix T = createTranslationMatrix(p), M = matrixStack.top();

            // cout << "\nM: \n"; printMatrix(M);
            // cout << "\nT: \n"; printMatrix(T);

            M = multiplyMatrices(M, T);

            matrixStack.pop();
            matrixStack.push(M);

            // cout << "\nM*T: \n"; printMatrix(M);
        } else if (cmd == "scale") {
            Point p; fin >> p;
            Matrix S = createScaleMatrix(p), M = matrixStack.top();

            // cout << "\nM: \n"; printMatrix(M);
            // cout << "\nS: \n"; printMatrix(S);

            M = multiplyMatrices(M, S);

            matrixStack.pop();
            matrixStack.push(M);

            // cout << "\nM*S: \n"; printMatrix(M);
        } else if (cmd == "rotate") {
            double theta; Point p;
            fin >> theta >> p;

            Matrix R = createRotationMatrix(theta, p), M = matrixStack.top();

            // cout << "\nM: \n"; printMatrix(M);
            // cout << "\nS: \n"; printMatrix(S);

            M = multiplyMatrices(M, R);

            matrixStack.pop();
            matrixStack.push(M);

            // cout << "\nM*R: \n"; printMatrix(M);
        } else if (cmd == "push") {
            matrixStack.push(matrixStack.top());
        } else if (cmd == "pop") {
            if (matrixStack.empty()) {
                cout << "Error: Stack is empty\n";
                return;
            }

            matrixStack.pop();
        } else if (cmd == "end") {
            break;
        } else {
            cout << "Invalid command\n";
            break;
        }
    }

    fin.close();
    fout.close();
}

// Stage 2: View Transformation
void ViewTransformation() {
    ifstream fin("stage1.txt");
    ofstream fout("stage2.txt");

    Matrix V = createViewMatrix(eye, look, up);

    for (int i=0; i<no_of_triangles; i++) {
        Point p1, p2, p3;
        fin >> p1 >> p2 >> p3;

        p1 = multiplyMatrixVector(V, p1);
        p2 = multiplyMatrixVector(V, p2);
        p3 = multiplyMatrixVector(V, p3);

        fout << p1 << "\n" << p2 << "\n" << p3 << "\n\n";
    }

    fin.close();
    fout.close();
}

// Stage 3: Projection Transformation
void ProjectionTransformation() {
    ifstream fin("stage2.txt");
    ofstream fout("stage3.txt");

    Matrix P = createProjectionMatrix(fovY, aspect, near, far);

    for (int i=0; i<no_of_triangles; i++) {
        Point p1, p2, p3;
        fin >> p1 >> p2 >> p3;

        p1 = multiplyMatrixVector(P, p1);
        p2 = multiplyMatrixVector(P, p2);
        p3 = multiplyMatrixVector(P, p3);

        fout << p1 << "\n" << p2 << "\n" << p3 << "\n\n";
    }

    fin.close();
    fout.close();
}

// Stage 4: Clipping and Scan Conversion using Z-Buffer Algorithm
void ClippingAndScanConversionUsingZBuffer() {
    ifstream fin("config.txt");

    fin >> screenWidth >> screenHeight;

    fin.close();
    fin.open("stage3.txt");

    double centerTop, centerBottom, centerLeft, centerRight;
    double dx, dy;
    
    dx = 2.0 / screenWidth;
    dy = 2.0 / screenHeight;

    centerTop = 1.0 - dy/2.0; 
    centerBottom = -1.0 + dy/2.0;
    centerRight = 1.0 - dx/2.0;
    centerLeft= -1.0 + dx/2.0;

    vector<Vector> zBuffer(screenHeight, Vector(screenWidth, zMax));

    bitmap_image image = createBlackBMP();

    int n = no_of_triangles;
    while (n--) {
        Point p1, p2, p3;
        fin >> p1 >> p2 >> p3;

        Triangle triangle(p1, p2, p3);

        minX = max(getMin(p1.x, p2.x, p3.x), centerLeft);
        maxX = min(getMax(p1.x, p2.x, p3.x), centerRight);

        minY = max(getMin(p1.y, p2.y, p3.y), centerBottom);
        maxY = min(getMax(p1.y, p2.y, p3.y), centerTop);

        // range of y-values to iterate within
        int  yStart, yEnd;
        yStart = round((centerTop - maxY) / dy);
        yEnd = round((centerTop - minY) / dy);

        if (yStart < 0 || yEnd >= screenWidth) {
            cout << "Unexpected error: row coordinate out of range\n";
            return;
        }

        for(int i = yStart; i <= yEnd; i++) {
            double ys, xa, xb, za, zb;
            ys = centerTop - i*dy;

            Vector xz(4); // {xa, xb, za, zb}
            xz = calculateXandZ(ys, p1, p2, p3);    

            xa = xz[0]; xb = xz[1];
            za = xz[2]; zb = xz[3];

            // range of x-values to iterate within
            int xStart = round((xa - centerLeft) / dx);
            int xEnd = round((xb - centerLeft) / dx);

            if (xStart < 0 || xEnd >= screenWidth) {
                cout << "Unexpected error: column coordinate out of range\n";
                return;
            }

            // project the triangle-portions
            for(int j=xStart; j<=xEnd; j++) {
                double xp = centerLeft + j*dx;
                double zp = zb - (zb-za) * ((xb-xp) / (xb -xa));

                if (zp < zMin) {
                    continue;
                }

                if(zp < zBuffer[i][j]) {
                    zBuffer[i][j] = zp;
                    image.set_pixel(j, i, triangle.red, triangle.green, triangle.blue);
                }
            }
        }
    }

    fin.close();
    ofstream fout("z_buffer.txt");

    for (int i = 0; i < screenHeight; i++) {
        for (int j = 0; j < screenWidth; j++) {
            if (zBuffer[i][j] < zMax) {
                fout << setprecision(6) << fixed << zBuffer[i][j] << "\t";
            }
        }
        fout << "\n";
    }

    image.save_image("out.bmp");

    zBuffer.clear();
    zBuffer.shrink_to_fit();

    fout.close();
}

bitmap_image createBlackBMP() {
    bitmap_image image(screenWidth, screenHeight);

    for (int i=0; i<screenHeight; i++) {
        for (int j=0; j<screenWidth; j++) {
            image.set_pixel(i, j, 0, 0, 0);
        }
    }

    return image;
}

double getMin(double n1, double n2, double n3) {
    if (n1 < n2) {
        if (n1 < n3) {
            return n1;
        } else {
            return n3;
        }
    } else {
        if (n2 < n3) {
            return n2;
        } else {
            return n3;
        }
    }
}

double getMax(double n1, double n2, double n3) {
    if (n1 > n2) {
        if (n1 > n3) {
            return n1;
        } else {
            return n3;
        }
    } else {
        if (n2 > n3) {
            return n2;
        } else {
            return n3;
        }
    }
}

Vector calculateXandZ(double ys, Point p1, Point p2, Point p3) {
    double xa, xb, za, zb, ta, tb;
    xa = xb = za = zb = 0;

    double x[3] = {p1.x, p2.x, p3.x};
    double y[3] = {p1.y, p2.y, p3.y};
    double z[3] = {p1.z, p2.z, p3.z};

    bool flag = true;

    for(int i=0; i<3; i++) {
        int j = (i+1) % 3;

        if(y[i] == y[j]) {
            continue;
        }

        double multiplier = (y[i] - ys) / (y[i] - y[j]);

        if (ys >= min(y[i], y[j]) && ys <= max(y[i], y[j])) {
            if (flag) {
                ta = xa = x[i] - (x[i] - x[j]) * multiplier;
                za = z[i] - (z[i] - z[j]) * multiplier;

                flag = false;
            } else {
                tb = xb = x[i] - (x[i] - x[j]) * multiplier;
                zb = z[i] - (z[i] - z[j]) * multiplier;
            }
        }
    }

    if (xa < minX) {
        xa = minX;
    } else if (xa > maxX) {
        xa = maxX;
    }

    if (xb < minX) {
        xb = minX;
    } else if (xb > maxX) {
        xb = maxX;
    }

    za = zb - (zb - za) * (tb - xa) / (tb - ta);
    zb = zb - (zb - za) * (tb - xb) / (tb - ta);

    if (xa >= xb) {
        swap(xa, xb);
        swap(za, zb);
    }

    return Vector{xa, xb, za, zb};
}

