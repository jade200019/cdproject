//
// Created by ubuntu on 3/28/19.
//

#ifndef TRYENV_COST_H
#define TRYENV_COST_H

#include "RotatedRectm.h"
#include "MergedComp.h"

class Cost {
    // Helper function for this->circle_from
    Point2f circle_center_from(float bx, float by, float cx, float cy);

    // REQUIRES: 3 points A, B, C
    // EFFECTS: return center and radius by reference
    void circle_from(Point2f & center, float & radius, Point2f A, Point2f B, Point2f C);

    // Helper function for Cost constructor
    float calc_angle_cost_more_patch(Point2f A1, Point2f B1, Point2f A2, Point2f B2, float angle1, float angle2);

    // REQUIRES: 4 points
    // EFFECTS: return the point of intersection of two straight lines
    Point2f lineLineIntersection(Point2f A, Point2f B, Point2f C, Point2f D);

    // REQUIRES: A, B are the two end point of Basel curve, with A1, B1 give information about direction
    // EFFECTS: return the length of the Basel curve
    float BezierLength(Point2f A, Point2f A1, Point2f B, Point2f B1);

    // same functionality as float BezierLength(Point2f A, Point2f A1, Point2f B, Point2f B1)
    // display procedure on dst
    float BezierLength(Mat & dst, Point2f A, Point2f A1, Point2f B, Point2f B1, Scalar color);

    // helper function
    float getPt( float n1 , float n2 , float perc );

    // REQUIRES: start point A, destination point B, intermediate point C, percentage of the curve perc
    // EFFECTS: return the corresponding Bezier point on the curve
    Point2f getBezierPoint(Point2f A, Point2f C, Point2f B, float perc);

    // REQUIRES: A, B are the two end point of Basel curve, with A1, B1 give information about direction
    // EFFECTS: sum of the cosine values: cosine of the angle between A1A and APt(0.01),
    // cosine of the angle between B1B and BPt(0.01)
    // can be optimized here
    float BezierAngleCosSum(Point2f A, Point2f A1, Point2f B, Point2f B1);

    // to replace BezierAngleCosSum(Point2f A, Point2f A1, Point2f B, Point2f B1)
    // angle between AC and CB
    float BezierAngleCos(Point2f A, Point2f A1, Point2f B, Point2f B1);

    // return cosine of angle between A1A1 and B1B2
    float getCosOfVector(Point2f A1, Point2f A2, Point2f B1, Point2f B2);

public:
    float dist;
    float angle;

    // REQUIRES: 2 rectangles (non-const to avoid recalculation of pts[8], nothing is changed)
    // displacement is to determine whether the two rRects are in the overlapped patches (HalfPatch),
    // neighboring patches (OnePatch), or displaced patches (MorePatch)
    // EFFECTS: .dist and .angle are calculated
    Cost(RotatedRectm & rRect1, RotatedRectm & rRect2, int displacement);

    // REQUIRES: coefficients for cost in distance and cost in angle, respectively
    // EFFECTS: return the result
    float calcTotalCost(float coeffDist, float coeffAngle, int displacement);


    float val;
    bool toMergeTail1, toMergeTail2;

    Cost(MergedComp & comp1, MergedComp & comp2);

    // EFFECTS: print the procedure of getting cost. print the bezier curve on dst.
    // same functionality as Cost(MergedComp & comp1, MergedComp & comp2);
    Cost(Mat & dst, MergedComp & comp1, MergedComp & comp2);

    float getVal(){
        return val;
    }

    // for debug Cost(), this function should only be called when debugging
    Cost(vector<RotatedRectm> & vecRect, MergedComp & comp1, MergedComp & comp2);
    void showPresentTwoComp(Mat & dst, vector<RotatedRectm> & vecRect, MergedComp & comp1, MergedComp & comp2);

    // Extend longer components, phase 2
    Cost(MergedComp & comp1, MergedComp & comp2, bool second);
    // to print
    Cost(Mat & dst, MergedComp & comp1, MergedComp & comp2, bool second);

    float lengthVal = 0, angleVal = 0;
    // common calculation procedure
    // write to this->lengthVal and this->angleVal
    void calculateCost(Mat & dst, MergedComp & comp1, MergedComp & comp2, bool printInfo);

    void calculateCost_v0(Mat &dst, MergedComp &comp1, MergedComp &comp2, bool printInfo);

    void calculateCost_v1(Mat &dst, MergedComp &comp1, MergedComp &comp2, bool printInfo);
};


#endif //TRYENV_COST_H
