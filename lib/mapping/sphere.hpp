#ifndef SPHERE_HPP
#define SPHERE_HPP

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <kdtree/geometry.h>
#include <kdtree/euclideanMetrics.h>

using kdtree::Point;
using kdtree::EuclideanMetrics;

class Sphere {
public:
    typedef Point<3, double> point_t;
    typedef EuclideanMetrics<3, double> metrics_t;

    Sphere(double X_, double Y_, double Z_, double R_, int typeID_);
    Sphere operator=(const Sphere& val);
//    static double dist2(std::vector<double> posA, std::vector<double> posB) {
//        double result = 0;
//        for (int i = 0; i < 3; i++) {
//            result += (posA[i] - posB[i]) * (posA[i] - posB[i]);
//        }
//        return result;
//    }

    std::vector<double> getBoundingBox() {
        return {X - R, Y - R, Z - R, X + R, Y + R, Z + R};
    }
    double getSize() {
        return R;
    }
    bool pointIsInside(point_t pos) {
        return metrics.distSqr(this->pos, pos) <= R*R;
    }
    bool isInsideBox(std::vector<double> pos, double dx) {
        return ((abs(this->pos[0] - pos[0]) < 0.5 * dx)
                && (abs(this->pos[1] - pos[1]) < 0.5 *dx)
                && (abs(this->pos[2] - pos[2]) < 0.5 * dx));
    }

    bool isInsideBox(point_t pos, double dx) {
        return ((abs(this->pos[0] - pos[0]) < 0.5 * dx)
                && (abs(this->pos[1] - pos[1]) < 0.5 *dx)
                && (abs(this->pos[2] - pos[2]) < 0.5 * dx));
    }

    void addMappingValue(double value) {
        sum += value;
        numOfSamples++;
        min = min > value ? value : min;
        max = max < value ? value : max;
    }
    double getMappingAverage() {
        if (numOfSamples == 0) {
            return 0;
        }
        return sum / numOfSamples;
    }
    double getMappingMax() {
        return max;
    }

    double getMappingMin() {
        return min;
    }
    void resetMapping() {
        sum = 0.0;
        numOfSamples = 0;
    }

    bool equals(Sphere& other) {
        return other.X == X
            && other.Y == Y
            && other.Z == Z
            && other.R == R;
    }

    const double X;
    const double Y;
    const double Z;
    const double R;

private:
    metrics_t metrics;

    const point_t pos;
    const int typeID;

    double sum;
    double numOfSamples;
    double min, max;

};

Sphere::Sphere(double X_, double Y_, double Z_, double R_, int typeID_) :
    X (X_), Y (Y_), Z(Z_), R(R_), pos(X_, Y_, Z_), typeID(typeID_) {
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::lowest();
    sum = 0;
    numOfSamples = 0;
}

Sphere Sphere::operator=(const Sphere& val) {
    Sphere res (val.X, val.Y, val.Z, val.R, val.typeID);
    return res;
}

#endif
