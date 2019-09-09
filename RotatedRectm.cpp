//
// Created by ubuntu on 3/22/19.
//

#include "RotatedRectm.h"

float widthRef[] = {8.0}; // Reference width of a cable

RotatedRectm::RotatedRectm(const RotatedRect &rRectIn, int X, int Y, int X_pixel, int Y_pixel,
        const Mat & Mask, Vec3f BGRCenter)
        : x(X), y(Y), x_pixel(X_pixel), y_pixel(Y_pixel), mask(Mask.clone()), BGRcenter(BGRCenter) {

    this->center = rRectIn.center;
    if (rRectIn.size.width >= rRectIn.size.height){ // copy as usual
        this->size = rRectIn.size;
        this->angle = rRectIn.angle;
    }
    else { // swap width (x) and height (y)
        this->size = Size2f(rRectIn.size.height, rRectIn.size.width);
        this->angle = rRectIn.angle + 90;
    }
    // make angle in range of [0, 180)
    while (this->angle < 0) this->angle += 180;
    while (this->angle >= 180) this->angle -= 180;

    this->calcConfidence();
    this->calcPoints8();
}

void RotatedRectm::calcConfidence() {
    double relativeDiff = abs(this->size.height - widthRef[0]) / widthRef[0];
    if (relativeDiff > 1)
        relativeDiff = 1;
    this->confidence = 1 - relativeDiff;
}

void RotatedRectm::calcPoints8() {
    // get the 4 vertices using the OpenCV function
    this->points(pts);
    // calculate the 4 mid-points
    for ( int i = 0; i < 4; i++){
        pts[i+4] = (pts[i] + pts[i + 1]) / 2;
    }
}

Point2f * RotatedRectm::getPoints8() {
    return this->pts;
}

void RotatedRectm::overlapAbstractArrow(Mat &dst, const Scalar& color,
        int thickness, int line_type, int shift, double tipLength) const {
    arrowedLine(dst, pts[4], pts[6], color, thickness, line_type, shift, tipLength);
}
