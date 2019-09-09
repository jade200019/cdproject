//
// Created by ubuntu on 3/22/19.
//

#ifndef TRYENV_ROTATEDRECTM_H
#define TRYENV_ROTATEDRECTM_H

#include "Config.h"

class RotatedRectm : public RotatedRect {
    // REQUIRES: this object is initialized using a RotatedRect
    // Called in the constructor
    // EFFECTS: calculate the confidence score, evaluating how likely the rectangle belongs to a cable
    // according to this->height (comparing to reference width of a cable)
    void calcConfidence();

    // Called in the constructor
    Point2f pts[8];
    void calcPoints8();

public:
    double confidence = 0.0; // between 0.0 and 1.0
    int x = -1, y = -1, r = -1;
    int x_pixel = 0, y_pixel = 0;
    Mat mask;   // k-means mask result, given by Mat Patch::getClosestColorMask() const
    Vec3f BGRcenter; // center color of the mask in k-means, given by Point3f Patch::getkClosestColorMaskColor() const

    RotatedRectm() = default;

    // REQUIRES: a RotatedRect object
    // EFFECTS: transform rotated rectangle to this object,
    // such that .size = [long side x short side], .angle within [0, 180)
    // i.e. .size.width is the long side, .size.height is the short side
    RotatedRectm(const RotatedRect & rRectIn, int X, int Y, int X_pixel, int Y_pixel,
            const Mat & Mask, Vec3f BGRCenter);

    // REQUIRES: instance constructed, so the confidence value is calculated
    // EFFECTS: set the index of rRect in the vector of rRect in a patch
    void setRIdx(const int rIdx){ r = rIdx; }

    /** returns 4 vertices and 4 mid-points of the rectangle
    @param pts The points array for storing rectangle vertices.
     The order is bottomLeft, topLeft, topRight, bottomRight.
     Then the mid-points of left, top, right, bottom.
     Allocated space should be 8 * sizeof(Point2f).
     Non-const output to avoid recalculation.
    */
    Point2f * getPoints8();

    // REQUIRES: dst size is same as the original image; BGR color
    // EFFECTS: overlap the arrow on dst
    void overlapAbstractArrow(Mat & dst, const Scalar& color,
            int thickness=1, int line_type=8, int shift=0, double tipLength=0.2) const;

};


#endif //TRYENV_ROTATEDRECTM_H
