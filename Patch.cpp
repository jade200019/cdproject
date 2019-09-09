//
// Created by jade on 2/27/19.
//

#include "Patch.h"

Vec3f RefBGR[] = { {6, 5, 7}, {10, 11, 12} };

const int MaxArea = 1024; // 32 * 32

Mat Patch::getSobelGradient() {
    // Based on https://docs.opencv.org/3.4.5/d5/db5/tutorial_laplace_operator.html
    Mat image, src = im, src_gray, grad;
    int ksize = 1;
    int scale = 1;
    int delta = 0;
    int ddepth = CV_16S;
    // Remove noise by blurring with a Gaussian filter ( kernel size = 3 )
    // GaussianBlur(image, src, Size(3, 3), 0, 0, BORDER_DEFAULT); // TODO receive signal SIGABRT: raise.c no such file
    // Convert the image to grayscale
    cvtColor(src, src_gray, COLOR_BGR2GRAY);

    Mat abs_grad_x, abs_grad_y;
    Sobel(src_gray, grad_x, ddepth, 1, 0, ksize, scale, delta, BORDER_DEFAULT);
    Sobel(src_gray, grad_y, ddepth, 0, 1, ksize, scale, delta, BORDER_DEFAULT);
    // cout << "grad_x = "<< endl << " "  << grad_x << endl << endl;

    // converting back to CV_8U
    convertScaleAbs(grad_x, abs_grad_x);
    convertScaleAbs(grad_y, abs_grad_y);
    addWeighted(abs_grad_x, 0.5, abs_grad_y, 0.5, 0, grad);

    // const String window_name = "Sobel Demo - Simple Edge Detector";
    //    imshow(window_name, grad);
    //    char key = (char)waitKey(0);

    return grad;
}


void Patch::getGradPolarCoordinate() {
    Mat grad_x_float, grad_y_float;
    // cartToPolar requires parameter type in CV_32F or CV_64F
    grad_x.convertTo(grad_x_float, CV_32F);
    grad_y.convertTo(grad_y_float, CV_32F);
    cartToPolar(grad_x_float, grad_y_float, grad_rho, grad_theta, false); // false: using angles in radian
    //cout << "grad_theta = "<< endl << " "  << grad_theta << endl << endl;
}

Mat Patch::getGradThetaHist(Mat mask, int magThreshold) {
    // https://docs.opencv.org/3.4.5/d8/dbc/tutorial_histogram_calculation.html
    // https://docs.opencv.org/3.4.5/d6/dc7/group__imgproc__hist.html#ga4b2b5fd75503ff9e6844cc4dcdaed35d

    float range[] = { 0, 2 * PI }; //the upper boundary is exclusive
    const float* histRange = { range };

    bool uniform = true, accumulate = false;

    if (mask.empty())
        mask = (grad_rho > magThreshold);
    cout << "grad_rho=\n" << grad_rho << endl << endl;
    cout << "mask=\n" << mask << endl << endl;
    calcHist( &grad_theta, 1, 0, mask, theta_hist, 1, &histSize, &histRange, uniform, accumulate );

    // draw the hist
    int hist_w = 512, hist_h = 400;
    int bin_w = cvRound( (double) hist_w/histSize );
    Mat histImage( hist_h, hist_w, CV_8UC3, Scalar( 0,0,0) );

    for( int i = 1; i < histSize; i++ )
    {
        line( histImage, Point( bin_w*(i-1), hist_h - cvRound(theta_hist.at<float>(i-1)) ),
              Point( bin_w*(i), hist_h - cvRound(theta_hist.at<float>(i)) ),
              Scalar( 255, 0, 0), 2, 8, 0  );

    }
    //    imshow("calcHist Demo", histImage );
    //    waitKey();
    //cout << "theta_hist = "<< endl << " "  << theta_hist << endl << endl;
    //cout << "size" << theta_hist.size << endl;
    return histImage;
}

Mat Patch::getGradThetaHistTopK(int k){
    // https://stackoverflow.com/questions/17831753/sorting-cvmat-in-opencv
    Mat dst;
    cv::sortIdx(theta_hist, dst, SORT_EVERY_COLUMN + SORT_DESCENDING );
    //cout << "theta_hist = "<< endl << " "  << theta_hist << endl << endl;
    //cout << "dst = "<< endl << " "  << dst << endl << endl;
    top_k_index = dst( Rect(0, 0, 1, k) ).clone();
    cout << "top_k_index = "<< endl << " "  << top_k_index << endl << endl;
    return top_k_index;
}

Mat Patch::ThetaHistBin2Mask(Mat bin_index){
    // must be CV_8U to fit the type of OR operation
    Mat mask(grad_theta.size(), CV_8U, Scalar(0));
    for (int i = 0; i < bin_index.rows; i++){
        int ind = bin_index.at<int>(i, 0);
        double left = ind * 2 * PI / histSize;
        float right = (ind + 1) * 2 * PI / histSize;
        Mat m = (this->grad_theta >= left) & (this->grad_theta < right);
        //cout << "m=" << endl << m << endl << endl;
        mask |= m;
    }
    top_theta_mask = mask;
    return mask;
}

Mat Patch::kMeans(int K){
    // http://answers.opencv.org/question/182006/opencv-c-k-means-color-clustering/
    // https://docs.opencv.org/3.1.0/d5/d38/group__core__cluster.html#ga9a34dc06c6ec9460e90860f15bcd2f88
    k = K;
    Mat ocv = this->im;

    // convert to float & reshape to a [3 x W*H] Mat
    //  (so every pixel is on a row of it's own)
    Mat data;
    ocv.convertTo(data,CV_32F);
    data = data.reshape(1,data.total());

    // do kmeans
    kmeans(data, k, kLabels, TermCriteria(CV_TERMCRIT_ITER, 10, 1.0), 3,
           KMEANS_PP_CENTERS, kCenters);

    // reshape both to a single row of Vec3f pixels:
    //cout << "kCenters=\n" << kCenters << endl<<endl;
    kCenters = kCenters.reshape(3, kCenters.rows);
    //cout << "kCenters=\n" << kCenters << endl<<endl;
    data = data.reshape(3,data.rows);

    // replace pixel values with their center value:
    Vec3f *p = data.ptr<Vec3f>();
    for (size_t i=0; i<data.rows; i++) {
        int center_id = kLabels.at<int>(i);
        p[i] = kCenters.at<Vec3f>(center_id);
    }

    // back to 2d, and uchar:
    ocv = data.reshape(3, ocv.rows);
    ocv.convertTo(ocv, CV_8U);

    this->kMeansSeg = ocv; // optional, for display purpose

    return ocv;
}

Mat Patch::closeBGRMask(const int margin, const vector<Vec3f> & vecRefBGR){
    int minDistLabel = -1;
    double minDist = BGR_MAX_DISTANCE_CONST;
    kClosestColorMaskColor.x = 0; // B
    kClosestColorMaskColor.y = 0; // G
    kClosestColorMaskColor.z = 0; // R

    for (int i = 0; i < this->k; i++){
        //double dist = cv::norm(this->kCenters.at<Vec3f>(i) - vecRefBGR[0]);
        double dist = DBL_MAX;
        for (const Vec3f & ref : vecRefBGR){
            if (norm(this->kCenters.at<Vec3f>(i) - ref) < dist){
                dist = norm(this->kCenters.at<Vec3f>(i) - ref);     // choose the minimum distance to one of the ref BGR
            }
        }


        if (dist < minDist && dist < margin){
            minDist = dist;
            minDistLabel = i;
            kClosestColorMaskColor = this->kCenters.at<Vec3f>(i);
        }
    }
    Mat kLabels2D = this->kLabels.reshape(1, PATCHY);
    this->kClosestColorMask = (kLabels2D == minDistLabel);
    this->kClosestColorMask.convertTo(this->kClosestColorMask, CV_8U);
    return this->kClosestColorMask;
}

Mat Patch::BGRInRangeMask(const int margin, const vector<Vec3f> &vecRefBGR) {
    vector<int> inRangeLabel;
    for (int i = 0; i < this->k; i++){
        double dist = DBL_MAX;
        for (const Vec3f & ref : vecRefBGR){
            if (norm(this->kCenters.at<Vec3f>(i) - ref) < dist){
                dist = norm(this->kCenters.at<Vec3f>(i) - ref);     // choose the minimum distance to one of the ref BGR
            }
        }
        if (dist < margin){
            inRangeLabel.push_back(i);
        }
    }
    Mat kLabels2D = this->kLabels.reshape(1, PATCHY);
    if (!inRangeLabel.empty()){
        this->kClosestColorMask = (kLabels2D == inRangeLabel[0]);
        for (int j : inRangeLabel){
            this->kClosestColorMask |= (kLabels2D == inRangeLabel[j]);
        }
    }

    this->kClosestColorMask.convertTo(this->kClosestColorMask, CV_8U);


    if (kClosestColorMask.dims == 0){       // prevent exception in morphologyEx()
        //cout << "kClosestColorMask.dims == 0" << endl;
        return this->kClosestColorMask;
    }

    // https://docs.opencv.org/3.4/d4/d86/group__imgproc__filter.html#ga67493776e3ad1a3df63883829375201f
    Mat element = getStructuringElement( MORPH_ELLIPSE, Size( 2* MASK_MORPH_SIZE + 1, 2 * MASK_MORPH_SIZE + 1 ),
            Point( MASK_MORPH_SIZE, MASK_MORPH_SIZE ) );
    Mat dst = this->kClosestColorMask.clone();
    morphologyEx( this->kClosestColorMask, dst, MORPH_CLOSE, element );
//#define SHOW_MORPH_RESULT
#ifdef SHOW_MORPH_RESULT
    imshow("mask", kClosestColorMask);
    imshow("dst", dst);
    static int i = 0;
    imwrite("bin/morph" + std::to_string(i++) + "_" + "mask" + ".jpg", dst);
    imwrite("bin/morph" + std::to_string(i++) + "_" + "dst" + ".jpg", dst);
    //waitKey(0);
#endif
    this->kClosestColorMask = dst.clone();

    return this->kClosestColorMask;
}


Mat Patch::closestBGRMaskContour() {
    // https://docs.opencv.org/3.4.5/df/d0d/tutorial_find_contours.html
    RNG rng(12345); // random number generator
    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;
    findContours( this->kClosestColorMask, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE );

    // save the display
    Mat drawing = Mat::zeros( this->kClosestColorMask.size(), CV_8UC3 );
    for( size_t i = 0; i< contours.size(); i++ )
    {
        Scalar color = Scalar( rng.uniform(0, 256), rng.uniform(0,256), rng.uniform(0,256) );
        drawContours( drawing, contours, (int)i, color, 1, LINE_8, hierarchy, 0 );
    }
    kClosestColorMaskContour = drawing.clone();
    return drawing.clone();
}

vector<RotatedRectm> & Patch::findCompMinAreaRectImageLevel(int imRows, int imCols, float minAreaPer, float maxAreaPer, float minConfidence) {
    if (this->kClosestColorMask.empty() )
        return this->vecRRect; // an empty vector
    Mat labels, stats, centroids;
    vector<Mat> vecComp; // store individual connected components, as 8U1C Mat;
                            // unnecessary if usage limited to this function

    // find connected components
    connectedComponentsWithStats(this->kClosestColorMask, labels, stats, centroids, 8, CV_16U);

    // "band pass" according to area threshold
    for (int i = 1; i < stats.rows; i++){ // ignore stats[0], it is for background
        int area = stats.at<int>(i, CC_STAT_AREA);
        if (area >= minAreaPer * MaxArea && area <= maxAreaPer * MaxArea){
            Mat binaryComp;
            inRange(labels, i, i, binaryComp); // find Mat binaryComp = (labels == i) in Python
            vecComp.push_back(binaryComp);
            //cout << "binaryComp=\n" << binaryComp << endl << endl;
        }
    }

    int rIdx = 0;
    for (Mat & binaryComp : vecComp){
        // pad with zero and extend to the size of the image
        Mat padded = Mat::zeros(imRows, imCols, CV_8U);

        Mat dst_roi = padded(Rect(this->x_pixel, this->y_pixel, PATCHX, PATCHY));
        binaryComp.copyTo(dst_roi);

        //imwrite("bin/Patch_" + std::to_string(this->x) + "_" + std::to_string(this->y) + "_padded.jpg", padded);

        // find rectangles
        vector<vector<Point> > contours;
        vector<Vec4i> hierarchy;
        findContours( padded, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE );


        RotatedRect rRect = minAreaRect(contours[0]);
        RotatedRectm transformedRRect(rRect, this->x, this->y, this->x_pixel, this->y_pixel,
                                      binaryComp, this->kClosestColorMaskColor);

        if (transformedRRect.confidence < minConfidence)
            continue; // does not add this RRect to the vector in the patch, TODO Modify the mechanism of adding and deleting RRects

        transformedRRect.setRIdx(rIdx++);               // set the index in vecRRect to the rRect instance
        this->vecRRect.push_back(transformedRRect);

    }

    return this->vecRRect;
}

vector<RotatedRectm> &
Patch::findMinAreaRect(int imRows, int imCols, Mat &canny, vector<float> & vecRefWidth,
                       float minRectWidth, float minRectLength,
                       float minAreaPercentage, float maxWidthDev,
                       float minCannyLeftConf, float minCannyRightConf, float minConfidence) {
    if (this->kClosestColorMask.empty() )
        return this->vecRRect; // an empty vector
    Mat labels, stats, centroids;

    // find connected components
    connectedComponentsWithStats(this->kClosestColorMask, labels, stats, centroids, 8, CV_16U);

    int rIdx = 0;
    for (int i = 1; i < stats.rows; i++) { // ignore stats[0], it is for background
        int area = stats.at<int>(i, CC_STAT_AREA);
        Mat binaryComp;
        inRange(labels, i, i, binaryComp); // find Mat binaryComp = (labels == i) in Python

        // pad with zero and extend to the size of the image
        Mat padded = Mat::zeros(imRows, imCols, CV_8U);
        Mat dst_roi = padded(Rect(this->x_pixel, this->y_pixel, PATCHX, PATCHY));
        binaryComp.copyTo(dst_roi);

        // find rectangles
        vector<vector<Point> > contours;
        vector<Vec4i> hierarchy;
        findContours( padded, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE );

        if (contours.size() < 1) { continue; } // may not exist
        RotatedRect rRectCV = minAreaRect(contours[0]);

        float areaPercentage = area / (rRectCV.size.width * rRectCV.size.height);      // area percentage initial
        if (area / minRectWidth / minRectLength < minAreaPercentage / 3) { continue; }  // filter out too small ones

        vector<Point2f> triangle1;
        minEnclosingTriangle(contours[0], triangle1);
        if (triangle1.size() != 3){
            std::cout << "Patch::findMinAreaRect triangle1.size() != 3, triangle1.size(): " << triangle1.size() << endl;
            continue;
            //std::cerr << "Patch::findMinAreaRect triangle1.size() != 3, triangle1.size(): " << triangle1.size() << endl;
        }
        int indM = check_on_boundary_and_corner(this->x_pixel, this->y_pixel, triangle1); // return -1 if not on boundary and corner

        bool useTriangle = false;
        if (indM != -1){
            Point2f & A = triangle1[(indM + 1) % 3];
            Point2f & B = triangle1[(indM + 2) % 3];
            Point2f & M = triangle1[indM];
            float triangle1_area = float(norm(A - M) * norm(B - M) / 2);    // right angle
            float triangle_perc = area / triangle1_area;

            if (triangle_perc > MIN_TRIANGLE_AREA_PERCENTAGE && triangle_perc > areaPercentage) {
                Point2f P = point_to_line_projection(A, B, M);
                Point2f A1 = A + (M - P); // move A along the vector PM
                RotatedRect rRectCV1(A1, A, B);     // construct a rotated rectangle using A1, A, B
                rRectCV.size = rRectCV1.size;       // update the rectangle
                rRectCV.center = rRectCV1.center;
                rRectCV.angle = rRectCV1.angle;
                areaPercentage = triangle_perc;
                useTriangle = true;
            }
        }

        RotatedRectm rRect(rRectCV, this->x, this->y, this->x_pixel, this->y_pixel,
                           binaryComp, this->kClosestColorMaskColor);

        float length = rRect.size.width, width = rRect.size.height;
        if (length <= 0 || width <= 0                                                   // skip empty or too small rectangles
            || length < minRectLength || width < minRectWidth){                         // filter out
            continue;
        }

        //float areaPercentage = area / (length * width);                               // area percentage
        if (areaPercentage < minAreaPercentage) { continue; }                           // filter out

        if (vecRefWidth.empty()){
            std::cerr << "Patch::findMinAreaRect vecRefWidth empty.\n";
        }
        //float widthDev = abs(width - vecRefWidth[0]) / vecRefWidth[0];                  // min width deviation
        float widthDev = abs(width - vecRefWidth[0]);                  // min width deviation
        for (float & ref : vecRefWidth){                                                // compared to reference width
            //float temp = abs(width - ref) / ref;
            float temp = abs(width - ref);
            if (temp < widthDev){
                widthDev = temp;
            }
        }
        if (widthDev > maxWidthDev) { continue; }                                       // filter out
#if ENABLE_CANNY
        RotatedRect rROI03 = getRRect03(rRect, CANNY_ROTATED_RECTANGLE_ROI_WIDTH);     // construct rotated rectangle ROI
        RotatedRect rROI12 = getRRect12(rRect, CANNY_ROTATED_RECTANGLE_ROI_WIDTH);     // to clip canny edges

        Mat cannySlice03 = extractRotatedrectArea(canny, rROI03);                      // clipped canny edges around pts[0] and pts[3]
        if (cannySlice03.empty()){
            std::cerr << "Patch::findMinAreaRect cannySlice03 empty.\n";
            std::cerr << "\trROI03: " << "rROI03.center: " << rROI03.center << " rROI03.size: " << rROI03.size
            << " rROI03.angle: " << rROI03.angle << endl;
        }
        Mat cannySlice12 = extractRotatedrectArea(canny, rROI12);                      // clipped canny edges around pts[1] and pts[2]
        if (cannySlice12.empty()){
            std::cerr << "Patch::findMinAreaRect cannySlice12 empty.\n";
            std::cerr << "\trROI12: " << "rROI12.center: " << rROI12.center << " rROI12.size: " << rROI12.size
                      << " rROI12.angle: " << rROI12.angle << endl;
        }

        float maxLenCannySlice03 = findMaxLenCannySlice(cannySlice03);
        float maxLenCannySlice12 = findMaxLenCannySlice(cannySlice12);


        float CannyLeftConf = (maxLenCannySlice03 / length > 1) ? 1 : maxLenCannySlice03 / length;
        float CannyRightConf = (maxLenCannySlice12 / length > 1) ? 1 : maxLenCannySlice12 / length;
        if (CannyLeftConf < minCannyLeftConf) { continue; }                             // filter out
        if (CannyRightConf < minCannyRightConf) { continue; }
#endif

#if 0
        cout << "length: " << length << " width: "<< width  << endl
            << " areaPercentage: " << areaPercentage << " widthDev: " << widthDev  << endl
            << " maxLenCannySlice03: " << maxLenCannySlice03 << " maxLenCannySlice12: " << maxLenCannySlice12 << endl
            << " CannyLeftConf: " << CannyLeftConf << " CannyRightConf: " << CannyRightConf << endl;
        imshow("canny", canny);
        imshow("patch im", this->im);
        imshow("patch mask", this->kClosestColorMask);
        imshow("canny03", cannySlice03);
        imshow("canny12", cannySlice12);
        char c = waitKey(0);
#endif

        rRect.setRIdx(rIdx++);                                                          // set the index in vecRRect to the rRect instance
        this->vecRRect.push_back(rRect);
//#define SHOW_TRIANGLE_EFFECT
#ifdef SHOW_TRIANGLE_EFFECT
        if (useTriangle){
            cout << "triangle rRectCV: " << rRectCV.size << " " << rRectCV.center << " " << rRectCV.angle << endl;
            Mat dst(padded.size(), CV_8UC3, Scalar(0));
            cvtColor(padded, dst, COLOR_GRAY2BGR);
            rRect.overlapAbstractArrow(dst, Scalar(0, 255, 0));
            static int i = 0;
            imwrite("bin/" + std::to_string(i++) + "_" + "show effec" + ".jpg", dst);
            imshow("show effec", dst);
            //char c = waitKey(0);
        }
#endif
    }

    return this->vecRRect;
}

RotatedRect Patch::getRRect03(RotatedRectm & rRect, float cannyROIWidth) {
    Point2f * pts;                          //Point2f pts[8];
    pts = rRect.getPoints8();

    // RotatedRect(const Point2f& center, const Size2f& size, float angle);
    Point2f center = (pts[0] + pts[3]) / 2;
    Size2f size;
    size.width = rRect.size.width;          // long side of the rRect
    size.height = cannyROIWidth;
    RotatedRect ret(center, size, rRect.angle);

    return ret;
}

RotatedRect Patch::getRRect12(RotatedRectm & rRect, float cannyROIWidth) {
    Point2f * pts;                          //Point2f pts[8];
    pts = rRect.getPoints8();

    // RotatedRect(const Point2f& center, const Size2f& size, float angle);
    Point2f center = (pts[1] + pts[2]) / 2;
    Size2f size;
    size.width = rRect.size.width;          // long side of the rRect
    size.height = cannyROIWidth;
    RotatedRect ret(center, size, rRect.angle);

    return ret;
}

Mat Patch::extractRotatedrectArea(const Mat & src, RotatedRect rect) {
    // http://answers.opencv.org/question/497/extract-a-rotatedrect-area/

    // rect is the RotatedRect (I got it from a contour...)
    //RotatedRect rect;
    // matrices we'll use
    Mat M, rotated, cropped;
    // get angle and size from the bounding box
    float angle = rect.angle;
    Size rect_size = rect.size;
    // thanks to http://felix.abecassis.me/2011/10/opencv-rotation-deskewing/
    if (rect.angle < -45.) {
        angle += 90.0;
        swap(rect_size.width, rect_size.height);
    }
    // get the rotation matrix
    M = getRotationMatrix2D(rect.center, angle, 1.0);
    // perform the affine transformation
    warpAffine(src, rotated, M, src.size(), INTER_CUBIC);
    // crop the resulting image
    getRectSubPix(rotated, rect_size, rect.center, cropped);

    return cropped;
}

float Patch::findMaxLenCannySlice(const Mat &cannySlice) {
//    Mat labels, stats, centroids;
//
//    // find connected components
//    connectedComponentsWithStats(cannySlice, labels, stats, centroids, 8, CV_16U);
//
//    //    int maxArea = 0;
//    for (int i = 1; i < stats.rows; i++) { // ignore stats[0], it is for background
//        //        int area = stats.at<int>(i, CC_STAT_AREA);  // value is too large to represent arc length of the clipped canny edge
//        //        if (area > maxArea){
//        //            maxArea = area;
//        //        }
//        Mat binaryComp;
//        inRange(labels, i, i, binaryComp); // find Mat binaryComp = (labels == i) in Python
//
//
//
//    }
//    //return maxArea;

    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;
    findContours( cannySlice, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE );

    double maxContourLen = 0;
    for( size_t i = 0; i < contours.size(); i++ )
    {
        double contourLen = arcLength(contours[i], true); // closed curve
        if (contourLen > maxContourLen){
            maxContourLen = contourLen;
        }
    }
    return float(maxContourLen / 2);
}

int Patch::check_on_boundary_and_corner(int x_pixel, int y_pixel, const vector<Point2f> &triangle) {
    for (int i = 0; i < 3; ++i){
        const Point2f & pt = triangle[i];
        if (((abs(pt.x - x_pixel) < EPSILON_DISTANCE || abs(pt.x - (x_pixel + PATCHX - 1)) < EPSILON_DISTANCE)
                && (pt.y > y_pixel - EPSILON_DISTANCE && pt.y < y_pixel + PATCHY - 1 + EPSILON_DISTANCE) )
            || ((abs(pt.y - y_pixel) < EPSILON_DISTANCE || abs(pt.y - (y_pixel + PATCHY - 1)) < EPSILON_DISTANCE)
                && (pt.x > x_pixel - EPSILON_DISTANCE && pt.x < x_pixel + PATCHX - 1 + EPSILON_DISTANCE) )
            ){
            continue;   // ok
        }
        else { // coordinate not in the range of boundary of the patch
            return -1;
        }
    }
    for (int i = 0; i < 3; ++i){
        const Point2f & M = triangle[i];
        const Point2f & A = triangle[(i + 1) % 3];
        const Point2f & B = triangle[(i + 2) % 3];
        float cosVal = getCosOfVector(M, A, M, B);
        if (abs(cosVal) < EPSILON_COS){    // return the index of the vertex M
            return i;
        }
    }
    return -1;
}

float Patch::getCosOfVector(const Point2f A1, const Point2f A2, const Point2f B1, const Point2f B2) {
    Point2f A1A2 = A2 - A1;
    Point2f B1B2 = B2 - B1;
    return float((A1A2.x * B1B2.x + A1A2.y * B1B2.y) / norm(A1A2) / norm(B1B2));
}

float Patch::triangle_area(const vector<Point2f> &triangle) {
    return 1;
}

Point2f Patch::point_to_line_projection(Point2f &A, Point2f &B, Point2f &M) {
    // https://stackoverflow.com/questions/22087193/optimized-functions-to-compute-projection-of-a-point-on-a-line
    // modified
    Point2f AB = B - A;
    Point2f AM = M - A;
    float d1 = AB.dot(AM) / norm(AB);
    Point2f AP = AB / norm(AB) * d1;
    Point2f P = AP + A;
    return P;
}


