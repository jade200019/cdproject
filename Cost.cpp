//
// Created by ubuntu on 3/28/19.
//

#include "Cost.h"

Cost::Cost(RotatedRectm &rRect1, RotatedRectm &rRect2, int displacement) {
    // Take the maximum long side of the two rectangles
    float Lmax = max(rRect1.size.width, rRect2.size.width);
    //Point2f pts1[8], pts2[8];
    Point2f * pts1, * pts2;
    pts1 = rRect1.getPoints8();
    pts2 = rRect2.getPoints8();
    float angle1, angle2;
    angle1 = rRect1.angle;
    angle2 = rRect2.angle;

    if (displacement == HalfPatch){ // rRect 2 is shifted by half of PATCHX and half of PATCHY
        this->dist = float(norm(pts1[7] - pts2[0])
                            + norm(pts1[5] - pts2[1])
                              + norm(pts1[3] - pts2[7])
                                + norm(pts1[2] - pts2[5])
                     ) / Lmax;
        this->angle = float(1 - cos((angle1 - angle2) * PI / 180));
    }
    else if (displacement == OnePatch){
        this->dist = float(norm(pts1[3] - pts2[0])
                           + norm(pts1[2] - pts2[1])
                            ) / Lmax;
        this->angle = float(1 - cos((angle1 - angle2) * PI / 180));
    }
    else { // MorePatch
        this->dist = float(norm(pts1[6] - pts2[4])
                            ) / Lmax;
        this->angle = calc_angle_cost_more_patch(pts1[4], pts1[6], pts2[4], pts2[6], angle1, angle2);
    }
//    cout << "this->dist: " << this->dist << " this->angle: " << this->angle << endl;
}

Point2f Cost::circle_center_from(float bx, float by, float cx, float cy) {
    // https://codeforces.com/blog/entry/17313
    float B=bx*bx+by*by, C=cx*cx+cy*cy, D=bx*cy-by*cx;
    //cout <<"bx:"<<bx<<" by:" << by << " cx:" << cx << " cy:" << cy << endl;
    //cout << "B:"<<B<<" C:" << C << " D:" << D << endl;
    return cv::Point2f((cy*B-by*C)/(2*D), (bx*C-cx*B)/(2*D));
}

void Cost::circle_from(Point2f &center, float &radius, Point2f A, Point2f B, Point2f C) {
    center = this->circle_center_from(B.x - A.x, B.y - A.y, C.x - A.x, C.y - A.y);
    radius = float(norm(center));
    center = center + A;
}

float Cost::calc_angle_cost_more_patch(Point2f A1, Point2f B1, Point2f A2, Point2f B2, float angle1, float angle2) {
    if (abs(angle1 - angle2) < EPSILON_ANGLE){
        // if A1, B1, A2 align on the same straight line, only consider the angle difference
        return float(1 - cos((angle1 - angle2) * PI / 180));
    }
    else {
        Point2f center;
        float radius;
        this->circle_from(center, radius, A1, B1, A2);

        Point2f OA2 = A2 - center;
        double OA2Angle = atan(OA2.y / OA2.x);
//        cout << "1 - abs(cos(OA2Angle + 90 - angle2)): " << 1 - abs(cos(OA2Angle + 90 - angle2))
//            << " atan(OA2.y / OA2.x): " << atan(OA2.y / OA2.x)
//            << " OA2.y:" << OA2.y << " OA2.x:" << OA2.x << endl
//            << " center:" << center.x << " " << center.y << endl;
        return float(1 - abs(cos(OA2Angle + 90 - angle2)));
    }
}

float Cost::calcTotalCost(float coeffDist, float coeffAngle, int displacement) {
    return coeffDist * this->dist + coeffAngle * this->angle;
}

Cost::Cost(MergedComp &comp1, MergedComp &comp2) {
    Mat dummy;
    calculateCost(dummy, comp1, comp2, false);
    float val = lengthVal + angleVal;
    this->val = val;

#if 0
    cout << "\t\tlengthARBL:" << lengthARBL << " lengthALBR:" << lengthALBR << endl
        << "\t\tcosARBL:" << cosARBL << " cosALBR:" << cosALBR << endl
        << "\t\tthis->val:" << this->val << endl
        << "--------------------------" << endl;
#endif

}



Point2f Cost::lineLineIntersection(Point2f A, Point2f B, Point2f C, Point2f D)
{
    // https://www.geeksforgeeks.org/program-for-point-of-intersection-of-two-lines/
    // Line AB represented as a1x + b1y = c1
    double a1 = B.y - A.y;
    double b1 = A.x - B.x;
    double c1 = a1*(A.x) + b1*(A.y);

    // Line CD represented as a2x + b2y = c2
    double a2 = D.y - C.y;
    double b2 = C.x - D.x;
    double c2 = a2*(C.x)+ b2*(C.y);

    double determinant = a1*b2 - a2*b1;

    if (determinant == 0)
    {
        // The lines are parallel. This is simplified
        // by returning a pair of FLT_MAX
        return Point2f(FLT_MAX, FLT_MAX);
    }
    else
    {
        double x = (b2*c1 - b1*c2)/determinant;
        double y = (a1*c2 - a2*c1)/determinant;
        return Point2f(x, y);
    }
}

float Cost::BezierLength(Point2f A, Point2f A1, Point2f B, Point2f B1) {
    Point2f C = lineLineIntersection(A, A1, B, B1); // intersection C

#if 0
    cout << "A:" << A << " A1:" << A1 << " B:" << B << " B1:" << B1 << endl
        << "\tC:" << C<<endl;
    Mat dst(Size(400, 400), CV_8UC3, Scalar(0));
    circle(dst, A, 4, Scalar(255, 0, 0), 3);       // blue
    circle(dst, A1, 2, Scalar(255, 153, 153));  // light blue
    line(dst, A, A1, Scalar(204, 51, 0));
    circle(dst, B, 4, Scalar(204, 51, 255), 3);    // pink
    circle(dst, B1, 2, Scalar(230, 153, 255));  // light pink
    line(dst, B, B1, Scalar(153, 0, 204));
    circle(dst, C, 6, Scalar(255, 255, 0));
#endif


    float length = 0;
    Point2f pt_prev = A;
    // https://stackoverflow.com/questions/785097/how-do-i-implement-a-b%C3%A9zier-curve-in-c/11435243#11435243
    for( float i = 0 ; i < 1 ; i += 0.01 )
    {
        // The Green Line
        //        xa = getPt( x1 , x2 , i );
        //        ya = getPt( y1 , y2 , i );
        //        xb = getPt( x2 , x3 , i );
        //        yb = getPt( y2 , y3 , i );
        Point2f pt = getBezierPoint(A, C, B, i);

        length += norm(pt - pt_prev);
        pt_prev = pt;

        //drawPixel( x , y , COLOR_RED );
#if 0
        circle(dst, pt, 1, Scalar(0, 255, 0)); // green
#endif
    }

#if 0
    cout << "\tlength:" << length << endl;

    imshow("intersection", dst);
    char c = waitKey(0);
#endif

    return length;
}

float Cost::getPt(float n1, float n2, float perc) {
    float diff = n2 - n1;
    return n1 + ( diff * perc );
}

Point2f Cost::getBezierPoint(Point2f A, Point2f C, Point2f B, float perc) {
    float xa = getPt(A.x, C.x, perc);
    float ya = getPt(A.y, C.y, perc);
    float xb = getPt(C.x, B.x, perc);
    float yb = getPt(C.y, B.y, perc);

    // The Black Dot
    float x = getPt( xa , xb , perc );
    float y = getPt( ya , yb , perc );
    return cv::Point2f(x, y);
}

float Cost::BezierAngleCosSum(Point2f A, Point2f A1, Point2f B, Point2f B1) {
    Point2f C = lineLineIntersection(A, A1, B, B1); // intersection C
    float cosA = getCosOfVector(A1, A, A, getBezierPoint(A, C, B, 0.01));
    float cosB = getCosOfVector(B1, B, B, getBezierPoint(A, C, B, 0.99));
    return cosA + cosB;
}

float Cost::BezierAngleCos(Point2f A, Point2f A1, Point2f B, Point2f B1) {
    Point2f C = lineLineIntersection(A, A1, B, B1); // intersection C
    return getCosOfVector(A1, A, B, B1);
}

float Cost::getCosOfVector(Point2f A1, Point2f A2, Point2f B1, Point2f B2) {
    Point2f A1A2 = A2 - A1;
    Point2f B1B2 = B2 - B1;
    return float((A1A2.x * B1B2.x + A1A2.y * B1B2.y) / norm(A1A2) / norm(B1B2));
}

Cost::Cost(Mat &dst, MergedComp &comp1, MergedComp &comp2) {
    calculateCost(dst, comp1, comp2, true);
    float val = lengthVal + angleVal;
    this->val = val;

    cout << "\t\tthis->val:" << this->val << endl
        << "--------------------------" << endl;
}

float Cost::BezierLength(Mat &dst, Point2f A, Point2f A1, Point2f B, Point2f B1, Scalar color) {
    Point2f C = lineLineIntersection(A, A1, B, B1); // intersection C

#if 0
    cout << "A:" << A << " A1:" << A1 << " B:" << B << " B1:" << B1 << endl
        << "\tC:" << C<<endl;
    circle(dst, A, 4, Scalar(255, 0, 0), 1);       // blue
    circle(dst, A1, 2, Scalar(255, 153, 153));  // light blue
    line(dst, A, A1, Scalar(204, 51, 0));
    circle(dst, B, 4, Scalar(204, 51, 255), 1);    // pink
    circle(dst, B1, 2, Scalar(230, 153, 255));  // light pink
    line(dst, B, B1, Scalar(153, 0, 204));
    circle(dst, C, 6, Scalar(255, 255, 0));
#endif


    float length = 0;
    Point2f pt_prev = A;
    // https://stackoverflow.com/questions/785097/how-do-i-implement-a-b%C3%A9zier-curve-in-c/11435243#11435243
    for( float i = 0 ; i < 1 ; i += 0.01 )
    {
        // The Green Line
        //        xa = getPt( x1 , x2 , i );
        //        ya = getPt( y1 , y2 , i );
        //        xb = getPt( x2 , x3 , i );
        //        yb = getPt( y2 , y3 , i );
        Point2f pt = getBezierPoint(A, C, B, i);

        length += norm(pt - pt_prev);
        pt_prev = pt;

        //drawPixel( x , y , COLOR_RED );
#if 1
        circle(dst, pt, 1, color);
#endif
    }

#if 0
    cout << "\tlength:" << length << endl;

    imshow("intersection", dst);
    char c = waitKey(0);
#endif

    return length;
}

void Cost::showPresentTwoComp(Mat & dst, vector<RotatedRectm> & vecRect, MergedComp &comp1, MergedComp &comp2) {


    Scalar color1 = Scalar( 0, 255, 0 );

    MergedComp & comp = comp1;
    for (int i : comp.listIdx){
        RotatedRectm & rRect = vecRect[i];
        Point2f vtx[4];
        rRect.points(vtx);

        // Draw the bounding box
        for (int i = 0; i < 4; i++)
            line(dst, vtx[i], vtx[(i + 1) % 4], color1, 1, LINE_AA);
    }


    Scalar color2 = Scalar( 255, 0, 0 );

    for (int i : comp2.listIdx){
        RotatedRectm & rRect = vecRect[i];
        Point2f vtx[4];
        rRect.points(vtx);

        // Draw the bounding box
        for (int i = 0; i < 4; i++)
            line(dst, vtx[i], vtx[(i + 1) % 4], color2, 1, LINE_AA);
    }

}

Cost::Cost(vector<RotatedRectm> &vecRect, MergedComp &comp1, MergedComp &comp2) {

    Mat dst(640, 640, CV_8UC3, Scalar(0));

    calculateCost(dst, comp1, comp2, true);

    float val = lengthVal + angleVal;
    this->val = val;

    showPresentTwoComp(dst, vecRect, comp1, comp2);

    //imshow("dst", dst);
    //char c = waitKey(0);
    imwrite("temp/" + std::to_string(comp1.headIdx) + "_" + std::to_string(comp2.headIdx) + "_NonShifted.jpg", dst);
}

Cost::Cost(MergedComp &comp1, MergedComp &comp2, bool second) {
    Mat dummy;
    calculateCost(dummy, comp1, comp2, false);
    float lengthScale = lengthVal * COST_VALUE_SCALE_LENGTH_ALPHA / (COST_VALUE_SCALE_LENGTH_ALPHA + comp1.listIdx.size() + comp2.listIdx.size());
    float angleScale = angleVal * COST_VALUE_SCALE_ANGLE_ALPHA / (COST_VALUE_SCALE_ANGLE_ALPHA + comp1.listIdx.size() + comp2.listIdx.size());
    float val = lengthScale + angleScale;
    this->val = val;
}

Cost::Cost(Mat &dst, MergedComp &comp1, MergedComp &comp2, bool second) {

    calculateCost(dst, comp1, comp2, true);

    float lengthScale = lengthVal * COST_VALUE_SCALE_LENGTH_ALPHA / (COST_VALUE_SCALE_LENGTH_ALPHA + comp1.listIdx.size() + comp2.listIdx.size());
    float angleScale = angleVal * COST_VALUE_SCALE_ANGLE_ALPHA / (COST_VALUE_SCALE_ANGLE_ALPHA + comp1.listIdx.size() + comp2.listIdx.size());
    float val = lengthScale + angleScale;
    this->val = val;

    cout << "\t\tthis->val:" << lengthVal +  angleVal << endl;
    cout << "\t\tmodified this->val:" << this->val << endl
         << "--------------------------" << endl;
}

void Cost::calculateCost(Mat &dst, MergedComp &comp1, MergedComp &comp2, bool printInfo) {
    calculateCost_v1(dst, comp1, comp2, printInfo);
}

void Cost::calculateCost_v0(Mat &dst, MergedComp &comp1, MergedComp &comp2, bool printInfo) {
    Point2f * ptsA = comp1.pts; // [8]
    Point2f * ptsB = comp2.pts; // [8]
    double minDist = MAX_LONG_SIDE * 2;
    int idxA = -1, idxB = -1;                                       // find the pair of outer points which has min distance
    for (int i = 0; i < 4; ++i){                                    // among (HL, HR, TL, TR) for comp1 and comp2
        for (int j = 0; j < 4; ++j){
            double dist = norm(ptsA[i] - ptsB[j]);
            if (dist < minDist){
                minDist = dist;
                idxA = i;
                idxB = j;
            }
        }
    }
    if (idxA == -1 || idxB == -1){
        std::cerr << "Cost::Cost idxA or idxB not correctly assigned.\n";
    }

    toMergeTail1 = (idxA > 1);                                      // if nearest point is 2/3, the tail is to be merged
    toMergeTail2 = (idxB > 1);

    Point2f AL = ptsA[(idxA / 2) * 2];                                    // Two pair of points used to calculate Basel
    Point2f AR = ptsA[(idxA / 2) * 2 + 1];
    Point2f BL = ptsB[(idxB / 2) * 2];
    Point2f BR = ptsB[(idxB / 2) * 2 + 1];

    Point2f AL1 = ptsA[(idxA / 2) * 2 + 4];                               // Points indicating the direction of head and tail rectangle
    Point2f AR1 = ptsA[(idxA / 2) * 2 + 5];
    Point2f BL1 = ptsB[(idxB / 2) * 2 + 4];
    Point2f BR1 = ptsB[(idxB / 2) * 2 + 5];


    float lengthARBL = float(norm(AR - BL));
    float lengthALBR = float(norm(AL - BR));
    float cosARBL = BezierAngleCos(AR, AR1, BL, BL1);
    float cosALBR = BezierAngleCos(AL, AL1, BR, BR1);
    float cosBezier = 1;

    if ( (lengthARBL > BEZIER_ESCAPE_SMALL_LENGTH                                // distance is large
          || lengthALBR > BEZIER_ESCAPE_SMALL_LENGTH) &&
         ! ((BezierAngleCos(AL1, AL, AL1, BR1) > COS_NEAR_PARALLEL && BezierAngleCos(BR1, BR, BR1, AL1) > COS_NEAR_PARALLEL) // AL1, AL, BR, BR1 on the same line
            ||
            (BezierAngleCos(AR1, AR, AR1, BL1) > COS_NEAR_PARALLEL && BezierAngleCos(BL1, BL, BL1, AR1) > COS_NEAR_PARALLEL) // AR1, AR, BL, BL1 on the same line
            ||
            (BezierAngleCos(AL1, AL, AL1, BL1) > COS_NEAR_PARALLEL && BezierAngleCos(BL1, BL, BL1, AL1) > COS_NEAR_PARALLEL) // AL1, AL, BL, BL1 on the same line
            ||
            (BezierAngleCos(AR1, AR, AR1, BR1) > COS_NEAR_PARALLEL && BezierAngleCos(BR1, BR, BR1, AR1) > COS_NEAR_PARALLEL) // AR1, AR, BR, BR1 on the same line
         )
        // no any of these four cases of same lines
            )
    {
        if (!dst.empty()) {
            lengthARBL = BezierLength(dst, AR, AR1, BL, BL1, Scalar(66, 244, 125)); // green
            lengthALBR = BezierLength(dst, AL, AL1, BR, BR1, Scalar(66, 244, 244)); // yellow
        }
        Point2f AM = (AL + AR) / 2;
        Point2f AM1 = (AL1 + AR1) / 2;
        Point2f BM = (BL + BR) / 2;
        Point2f BM1 = (BL1 + BR1) / 2;
        Point2f C = lineLineIntersection(AM1, AM, BM1, BM); // intersection C
        cosBezier = getCosOfVector(AM, C, C, BM);
        if (printInfo) {
            cout << "\t\t\t\t\t\t\t\t\t\t\t\t\tBezier\n";
        }
    }

    // connecting angle cost
    float cosAR1ARBL = (lengthARBL > ANGLE_ESCAPE_SMALL_LENGTH || lengthALBR > ANGLE_ESCAPE_SMALL_LENGTH) ? getCosOfVector(AR1, AR, AR, BL) : 1;
    float cosAL1ALBR = (lengthARBL > ANGLE_ESCAPE_SMALL_LENGTH || lengthALBR > ANGLE_ESCAPE_SMALL_LENGTH) ? getCosOfVector(AL1, AL, AL, BR) : 1;
    float cosBR1BRAL = (lengthARBL > ANGLE_ESCAPE_SMALL_LENGTH || lengthALBR > ANGLE_ESCAPE_SMALL_LENGTH) ? getCosOfVector(BR1, BR, BR, AL) : 1;
    float cosBL1BLAR = (lengthARBL > ANGLE_ESCAPE_SMALL_LENGTH || lengthALBR > ANGLE_ESCAPE_SMALL_LENGTH) ? getCosOfVector(BL1, BL, BL, AR) : 1;

    this->lengthVal = (lengthARBL + lengthALBR) ;

    this->angleVal = COSCOEFFICIENT * (2 - cosARBL - cosALBR) + COSCOEFFICIENT_1 * (4 - cosAR1ARBL - cosAL1ALBR - cosBR1BRAL - cosBL1BLAR)
                     + COS_COEFFICIENT_BEZIER * (1 - cosBezier);
    if (printInfo){
        cout << "\t\tlengthARBL:" << lengthARBL << " lengthALBR:" << lengthALBR << endl
             << "\t\tcosARBL:" << cosARBL << " cosALBR:" << cosALBR << " (2 - cosARBL - cosALBR):" << (2 - cosARBL - cosALBR) << endl
             << "\t\tIdx AL:" << (idxA / 2) * 2 << " AR:" << (idxA / 2) * 2 + 1 <<" BL:" << (idxB / 2) * 2 << " BR:"<< (idxB / 2) * 2 + 1 << endl
             << "--------------------------" << endl;
    }
}

void Cost::calculateCost_v1(Mat &dst, MergedComp &comp1, MergedComp &comp2, bool printInfo) {

    Point2f * ptsA = comp1.pts; // [8]
    Point2f * ptsB = comp2.pts; // [8]
    double minDist = MAX_LONG_SIDE * 2;
    int idxA = -1, idxB = -1;                                       // find the pair of outer points which has min distance
    for (int i = 0; i < 4; ++i){                                    // among (HL, HR, TL, TR) for comp1 and comp2
        for (int j = 0; j < 4; ++j){
            double dist = norm(ptsA[i] - ptsB[j]);
            if (dist < minDist){
                minDist = dist;
                idxA = i;
                idxB = j;
            }
        }
    }
    if (idxA == -1 || idxB == -1){
        std::cerr << "Cost::Cost idxA or idxB not correctly assigned.\n";
    }

    toMergeTail1 = (idxA > 1);                                      // if nearest point is 2/3, the tail is to be merged
    toMergeTail2 = (idxB > 1);

    Point2f AL = ptsA[(idxA / 2) * 2];                                    // Two pair of points used to calculate Basel
    Point2f AR = ptsA[(idxA / 2) * 2 + 1];
    Point2f BL = ptsB[(idxB / 2) * 2];
    Point2f BR = ptsB[(idxB / 2) * 2 + 1];

    Point2f AL1 = ptsA[(idxA / 2) * 2 + 4];                               // Points indicating the direction of head and tail rectangle
    Point2f AR1 = ptsA[(idxA / 2) * 2 + 5];
    Point2f BL1 = ptsB[(idxB / 2) * 2 + 4];
    Point2f BR1 = ptsB[(idxB / 2) * 2 + 5];

    float J_dist, J_angle, J_connect_angle;
    float lengthARBL = float(norm(AR - BL));
    float lengthALBR = float(norm(AL - BR));
    J_dist = lengthARBL + lengthALBR;


    Point2f C1 = lineLineIntersection(AL1, AL, BR1, BR); // intersection C1
    Point2f C2 = lineLineIntersection(AR1, AR, BL1, BL); // intersection C2
    float angle1 = (lengthALBR > ANGLE_ESCAPE_SMALL_LENGTH) ? getCosOfVector(AL, C1, C1, BR) : 1;
    float angle2 = (lengthARBL > ANGLE_ESCAPE_SMALL_LENGTH) ? getCosOfVector(AR, C2, C2, BL) : 1;
    J_angle = 2 - angle1 - angle2;

    // connecting angle cost
    float cosAR1ARBL = (lengthARBL > ANGLE_ESCAPE_SMALL_LENGTH || lengthALBR > ANGLE_ESCAPE_SMALL_LENGTH) ? getCosOfVector(AR1, AR, AR, BL) : 1;
    float cosAL1ALBR = (lengthARBL > ANGLE_ESCAPE_SMALL_LENGTH || lengthALBR > ANGLE_ESCAPE_SMALL_LENGTH) ? getCosOfVector(AL1, AL, AL, BR) : 1;
    float cosBR1BRAL = (lengthARBL > ANGLE_ESCAPE_SMALL_LENGTH || lengthALBR > ANGLE_ESCAPE_SMALL_LENGTH) ? getCosOfVector(BR1, BR, BR, AL) : 1;
    float cosBL1BLAR = (lengthARBL > ANGLE_ESCAPE_SMALL_LENGTH || lengthALBR > ANGLE_ESCAPE_SMALL_LENGTH) ? getCosOfVector(BL1, BL, BL, AR) : 1;
    J_connect_angle = 4 - cosAR1ARBL - cosAL1ALBR - cosBR1BRAL - cosBL1BLAR;

    this->lengthVal = J_dist;

    this->angleVal = 100 * J_connect_angle;
    //this->angleVal = 1000 * J_angle;
    //this->angleVal = 20 * J_connect_angle + 20 * J_angle;

    if (printInfo){
        cout << "\t\tlengthARBL:" << lengthARBL << " lengthALBR:" << lengthALBR << " J_dist:" << J_dist << endl
             //<< "\t\tangle1:" << angle1 << " angle2:" << angle2 << " J_angle:" << J_angle << endl
             << "\t\tcosAR1ARBL:" << cosAR1ARBL << " cosAL1ALBR:" << cosAL1ALBR << " cosBR1BRAL:" << cosBR1BRAL << " cosBL1BLAR:" << cosBL1BLAR << " J_connect_angle:" << J_connect_angle << endl
             << "\t\tIdx AL:" << (idxA / 2) * 2 << " AR:" << (idxA / 2) * 2 + 1 <<" BL:" << (idxB / 2) * 2 << " BR:"<< (idxB / 2) * 2 + 1 << endl
             << "--------------------------" << endl;
    }
}






