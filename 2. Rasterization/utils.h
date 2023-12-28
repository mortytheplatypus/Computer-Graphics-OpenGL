#include <bits/stdc++.h>
using namespace std;

#define PI acos(-1.0)
#define DegreeToRadian(theta) (theta * PI / 180.0)
#define DIM 4

static unsigned long int g_seed = 1;
inline int randomRGB() {
    g_seed = (214013 * g_seed + 2531011);
    return ((g_seed >> 16) & 0x7FFF) % 256;
}

typedef vector<double> Vector;
typedef vector<Vector> Matrix;

struct Point {
    double x, y, z, w;

    Point() {
        x = 0;
        y = 0;
        z = 0;
        w = 1.0;
    }

    Point(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = 1.0;
    }

    Point(double x, double y, double z, double w) {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = w;
    }

    void normalize() {
        double len = sqrt(x*x + y*y + z*z);

        x /= len;
        y /= len;
        z /= len;
    }

    void scaleDown() {
        x /= w;
        y /= w;
        z /= w;
        
        w = 1.0;
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

    friend istream& operator>>(istream &in, Point &p) {
        in >> p.x >> p.y >> p.z;

        return in;
    }

    friend ostream& operator<<(ostream &out, Point &p) {
        out << fixed << setprecision(7) << p.x << " " << p.y << " " << p.z;

        return out;
    }
};

struct Triangle {
    Point points[3];
    int red, green, blue;

    Triangle() {
        points[0] = Point(0, 0, 0);
        points[1] = Point(0, 0, 0);
        points[2] = Point(0, 0, 0);

        red = green = blue = 0;
    }

    Triangle(Point p1, Point p2, Point p3) {
        points[0] = p1;
        points[1] = p2;
        points[2] = p3;

        red = randomRGB();
        green = randomRGB();
        blue = randomRGB();
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

Matrix createNullMatrix() {
    Matrix matrix;
    matrix.resize(DIM, Vector(DIM, 0));

    return matrix;
}

Matrix createIdentityMatrix() {
    Matrix identityMatrix = createNullMatrix();

    for (int i=0; i<DIM; i++) {
        identityMatrix[i][i] = 1;
    }

    return identityMatrix;
}

Matrix multiplyMatrices(Matrix A, Matrix B) {
    Matrix P = createNullMatrix();

    for (int i=0; i<DIM; i++) {
        for (int j=0; j<DIM; j++) {
            for (int k=0; k<DIM; k++) {
                P[i][j] += (A[i][k] * B[k][j]);
            }
        }
    }

    return P;
}

Point multiplyMatrixVector(Matrix A, Point P) {
    vector<Vector> V(DIM, vector<double>(DIM));
    Vector product(DIM, 0);

    V[0][0] = P.x; V[1][0] = P.y; V[2][0] = P.z; V[3][0] = 1;

    for (int i=0; i<DIM; i++) {
        for (int k=0; k<DIM; k++) {
            product[i] += (A[i][k] * V[k][0]);
        }
    }

    Point point(product[0], product[1], product[2], product[3]);
    point.scaleDown();
    
    return point;
}

Matrix createTranslationMatrix(Point p) {
    Matrix T = createIdentityMatrix();

    T[0][3] = p.x;
    T[1][3] = p.y;
    T[2][3] = p.z;

    return T;
}

Matrix createScaleMatrix(Point p) {
    Matrix S = createIdentityMatrix();
    
    S[0][0] = p.x;
    S[1][1] = p.y;
    S[2][2] = p.z;

    return S;
}

Point applyRodrigues(Point x, Point a, double theta) {
    double th = DegreeToRadian(theta);

    return x*cos(th) + a*(1-cos(th))*dotProduct(a, x) + crossProduct(a, x)*sin(th);
}

Matrix createRotationMatrix(double theta, Point p) {
    Matrix R = createIdentityMatrix();

    Point a(p.x, p.y, p.z);
    a.normalize();

    Point x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);

    Point c1, c2, c3;

    c1 = applyRodrigues(x, a, theta);
    c2 = applyRodrigues(y, a, theta);
    c3 = applyRodrigues(z, a, theta);

    R[0][0] = c1.x;
    R[1][0] = c1.y;
    R[2][0] = c1.z;

    R[0][1] = c2.x;
    R[1][1] = c2.y;
    R[2][1] = c2.z;

    R[0][2] = c3.x;
    R[1][2] = c3.y;
    R[2][2] = c3.z;

    return R;
}

Matrix createViewMatrix(Point eye, Point look, Point up) {
    Matrix V = createIdentityMatrix();

    Point l, r, u;

    l = look - eye;
    l.normalize();

    r = crossProduct(l, up);
    r.normalize();

    u = crossProduct(r, l);
    u.normalize();

    Matrix T = createTranslationMatrix(eye*(-1));

    Matrix R = createIdentityMatrix();

    R[0][0] = r.x;
    R[0][1] = r.y;
    R[0][2] = r.z;

    R[1][0] = u.x;
    R[1][1] = u.y;
    R[1][2] = u.z;

    R[2][0] = -l.x;
    R[2][1] = -l.y;
    R[2][2] = -l.z;

    V = multiplyMatrices(R, T);

    return V; 
}

Matrix createProjectionMatrix(double fovY, double aspect, double near, double far) {
    Matrix P = createNullMatrix();
    double fovX = fovY * aspect, t, r;
    double th_t = DegreeToRadian(fovY/2.0), th_r = DegreeToRadian(fovX/2.0);

    t = near * tan(th_t);
    r = near * tan(th_r);

    P[0][0] = near/r;
    P[1][1] = near/t;
    P[2][2] = -(far+near)/(far-near);
    P[2][3] = -2.0*far*near/(far-near);
    P[3][2] = -1;

    return P;
}

void printMatrix(Matrix matrix) {
    for (int i=0; i<DIM; i++) {
        for (int j=0; j<DIM; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}
