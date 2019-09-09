//
// Created by jade on 2/27/19.
//

#ifndef TRYENV_PATCH_H
#define TRYENV_PATCH_H

#include "Config.h"
#include "RotatedRectm.h"

class Patch {
    Mat im;
    Mat grad_x, grad_y;
    Mat grad_rho, grad_theta;
    int histSize = 256; // number of bins, theta_hist size: 256 x 1
    Mat theta_hist;
    Mat top_k_index;
    Mat top_theta_mask;

    // k-means
    int k;
    Mat kLabels, kCenters;
    Mat kClosestColorMask;
    Point3f kClosestColorMaskColor; // BGR color
    Mat kMeansSeg; // optional, saved for display purpose
    vector<RotatedRectm> vecRRect; // Rotated min area rectangle

    // contour
    Mat kClosestColorMaskContour; // optional, saved for display purpose

    // REQUIRES: a rRect, and width of the desired ROI for canny
    // EFFECTS: return a ROI for canny, using pts[0] and pts[3]
    RotatedRect getRRect03(RotatedRectm & rRect, float cannyROIWidth);

    // REQUIRES: a rRect, and width of the desired ROI for canny
    // EFFECTS: return a ROI for canny, using pts[1] and pts[2]
    RotatedRect getRRect12(RotatedRectm & rRect, float cannyROIWidth);

    // REQUIRES: the source matrix, and a rotated rectangle ROI
    // EFFECTS: return the cropped area
    Mat extractRotatedrectArea(const Mat & src, RotatedRect rect);

    // REQUIRES: a slice (rotated rectangle ROI) of canny edges
    // EFFECTS: maximum canny edge "arc length" (i.e. area of one connected component in canny edges)
    float findMaxLenCannySlice(const Mat & cannySlice);

    // REQUIRES: a triangle, start pixel of the patch
    // EFFECTS: if all points are on the boundary of the patch, and if they can form a right angle, return the index
    // of the vertex of the right angle
    // if any of the two conditions cannot be satisfied, return -1
    int check_on_boundary_and_corner(int x_pixel, int y_pixel, const vector<Point2f> & triangle);

    float getCosOfVector(Point2f A1, Point2f A2, Point2f B1, Point2f B2);

    // REQUIRES: a triangle
    // EFFECTS: return the area of it
    float triangle_area(const vector<Point2f> & triangle);

    // REQUIRES: AB is a straight line, M is a point outside it
    // EFFECTS: return the projection of M on the line AB
    Point2f point_to_line_projection(Point2f & A, Point2f & B, Point2f & M);

public:
    int x, y; // Patch position in the original image, by index, start from 0
    int x_pixel, y_pixel; // Patch position by pixel in the original image, start from 0


    Patch(Mat patch) {
        im = patch.clone();
    }

    Patch(const Mat & patch, int X, int Y, bool shifted = false) : im(patch.clone()), x(X), y(Y), k(3) {
        if (!shifted){ // for not shifted patches
            x_pixel = X * PATCHX;
            y_pixel = Y * PATCHY;
        }
        else { // for patches shifted by half of the size in both direction
            x_pixel = X * PATCHX + PATCHX / 2;
            y_pixel = Y * PATCHY + PATCHY / 2;
        }
    }

    // Mat can deallocate its own memory
    ~Patch() = default;

    Mat getIm() const { return this->im; }

    // REQUIRES: member im
    // EFFECTS: Return the abs value of Sobel gradient
    // EFFECTS: Save results to member grad_x, grad_y
    Mat getSobelGradient();

    // REQUIRES: member grad_x, grad_y
    // EFFECTS: member grad_rho, grad_theta
    void getGradPolarCoordinate();

    // REQUIRES: this->grad_theta, argument mask represent the counted elements
    // by default, a mask is not used
    // EFFECTS: return histogram in this->theta_hist
    // default argument only defined once in .h file, not in .cpp file
    Mat getGradThetaHist(Mat mask = Mat(), int magThreshold = 0);

    // REQUIRES: this->theta_hist, argument int k
    // EFFECTS: get top k-ranked theta orientation, return bin index, size k x 1
    Mat getGradThetaHistTopK(int k);

    // REQUIRES: this->grad_theta, argument Mat bin_index
    // EFFECTS: this->top_theta_mask, plot according to the given theta bin
    Mat ThetaHistBin2Mask(Mat bin_index);

    // REQUIRES: argument threshold on gradient magnitude
    // EFFECTS: return a mask representing large enough gradient magnitude
    // Mat thresholdGradMag(int threshold);

    // REQUIRES: this->im is a BGR image, CV_8UC3
    // EFFECTS: k-means segmentation, save to this->kLabels, this->kCenters, this->kMeansSeg
    // return this->kMeansSeg
    Mat kMeans(int K);

    // REQUIRES: this->kLabels, this->kCenters, Vec3f RefBGR[]
    // EFFECTS: this->closestColorMask, choose the closest kCenter color to one of RefBGR[],
    // and it should be in the margin
    Mat closeBGRMask(const int margin, const vector<Vec3f> & vecRefBGR);

    // REQUIRES: this->kLabels, this->kCenters, & vecRefBGR
    // EFFECTS: choose in range mask, according to k-centers, not pixels
    Mat BGRInRangeMask(const int margin, const vector<Vec3f> & vecRefBGR);

    // REQUIRES: this->closestColorMask, 8UC1
    // EFFECTS: return it, later as a function variable to assemble back an image to display
    Mat getClosestColorMask() const { return this->kClosestColorMask; }

    // REQUIRES: this->kClosestColorMaskColor
    // EFFECTS: return it
    Point3f getkClosestColorMaskColor() const { return this->kClosestColorMaskColor; }

    // REQUIRES: this->kMeansSeg, 8U3C
    // EFFECTS: return it
    Mat getKMeansSeg() const {return this->kMeansSeg; }

    // REQUIRES: this->closestColorMask, better not using this->kMeansSeg
    // EFFECTS: find contour
    Mat closestBGRMaskContour();

    // REQUIRES: this->kClosestColorMaskContour, 8U3C
    // EFFECTS: return it
    Mat getkClosestColorMaskContour() const { return this->kClosestColorMaskContour; }

    // REQUIRES: this->closestColorMask; info about size of the image: imRows, imCols; size filter for each connected
    // components: minAreaPercentage, maxAreaPercentage, ranging from 0.0 to 1.0
    // EFFECTS: find connected components, filter with size, pad with 0, find minAreaRect, save the Rect
    // save to this->vecRRect
    vector<RotatedRectm> & findCompMinAreaRectImageLevel(int imRows, int imCols, float minAreaPer = 0.0, float maxAreaPer = 1.0, float minConfidence = 0.0);

    // REQUIRES: this->vecRRect
    // EFFECTS: return its reference
    vector<RotatedRectm> & getVecRRect() { return this->vecRRect; }

    // REQUIRES: this->closestColorMask; info about size of the image: imRows, imCols;
    // canny edges of the entire image: canny
    // min area percentage of the connected component in a mask over fitted rotated rectangle: minAreaPercentage
    // max relative deviation of the width of the rRect compared to the reference: maxWidthDev
    // EFFECTS: this->vecRRect
    vector<RotatedRectm> & findMinAreaRect(int imRows, int imCols, Mat & canny, vector<float> & vecRefWidth,
            float minRectWidth = 0.0, float minRectLength = 0.0,
            float minAreaPercentage = 0.0, float maxWidthDev = 10.0,
            float minCannyLeftConf = 0.0, float minCannyRightConf = 0.0, float minConfidence = 0.0);


};


#endif //TRYENV_PATCH_H
