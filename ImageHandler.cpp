//
// Created by ubuntu on 3/13/19.
//

#include "ImageHandler.h"

#include "Config.h"

Mat ImageHandler::Resize(){
    int maxSide = max(this->im.rows, this->im.cols);
    if (maxSide < MAX_LONG_SIDE)
        return this->im;
    Mat dst;
    resize(this->im, dst, cv::Size(), MAX_LONG_SIDE/maxSide, MAX_LONG_SIDE/maxSide);
    this->im = dst;
    return dst;
}

vector<vector<Patch>> & ImageHandler::SplitToPatch(){
    this->vecVecPatch.clear();
    int i = 0, j = 0;
    cv::Mat patchMat, image = this->im;
    for (; i < (image.cols / PATCHX); i++) {
        j = 0;
        vector<Patch> vecPatch;
        for (; j < image.rows / PATCHY; j++) {
            cv::Rect rangeRect(i * PATCHX, j * PATCHY, PATCHX, PATCHY);
            patchMat = image(rangeRect); // shallow copy, reference to the original image
            Patch patch(patchMat, i, j);
            //this->vecPatch.push_back(patch);
            vecPatch.push_back(patch);
//#ifdef DEBUG_SAVE_PATCH_TO_FILE
#if 1
            cv::Mat patchLocation;
            patchLocation = image.clone(); // deep copy, another instance
            cv::rectangle(patchLocation, rangeRect, cv::Scalar(0, 255, 0));

            if (patchMat.empty()) {
                std::cerr << "Cannot get the patch." << std::endl;
            }
            // Save the frame into a file
            imwrite("bin/Patch_" + std::to_string(i) + "_" + std::to_string(j) + ".jpg", patchMat); // A JPG FILE IS BEING SAVED
            imwrite("bin/Patch_" + std::to_string(i) + "_" + std::to_string(j) +  "_Location.jpg", patchLocation);
#endif
        }
        this->vecVecPatch.push_back(vecPatch);
    }

    this->patchColNum = image.cols / PATCHX;
    this->patchRowNum = image.rows / PATCHY;

    cout << "Patch Num Col: "<<patchColNum << " Row: " <<patchRowNum << endl;

    return this->vecVecPatch;
}

vector<vector<Patch>> &ImageHandler::SplitToShiftedPatch() {
    this->vecVecShiftedPatch.clear();
    int x = PATCHX / 2, y = PATCHY / 2;
    cv::Mat patchMat, image = this->im;
    for (; x <= image.cols - PATCHX; x += PATCHX){
        y = PATCHY / 2;
        vector<Patch> vecShiftedPatch;
        for (; y <= image.rows - PATCHY; y += PATCHY){
            Rect rangeRect(x, y, PATCHX, PATCHY);
            patchMat = image(rangeRect);
            Patch patch(patchMat, x / PATCHX, y / PATCHY, true); // true: shifted patch constructor
            vecShiftedPatch.push_back(patch);
        }
        this->vecVecShiftedPatch.push_back(vecShiftedPatch);
    }
    this->shiftedPatchColNum = x / PATCHX;
    this->shiftedPatchRowNum = y / PATCHY;

    cout << "Shifted Patch Num Col: "<<shiftedPatchColNum << " Row: " <<shiftedPatchRowNum << endl;

    return this->vecVecShiftedPatch;
}


void ImageHandler::DisplayMaskInImage(const vector<vector<Patch>> &vecVecPatch, String FileName, Mat (Patch::*getMask)() const, bool shifted){
    Mat dst(this->im.size(), CV_8U, Scalar(0));
    for (auto &vecPatch : vecVecPatch) {
        for (auto &patch : vecPatch) {
            int left = patch.x_pixel;
            int top = patch.y_pixel;
            Mat dst_roi = dst(Rect(left, top, PATCHX, PATCHY));

            //getMask().copyTo(dst_roi); // not recognized as pointer-to-function

            Mat src = patch.getClosestColorMask();
            src.copyTo(dst_roi);

        }
    }
    imwrite("bin/" + FileName + "Mask.jpg", dst);

    // Draw boundaries of patches on the mask
    Mat maskWithBox(this->im.size(), this->im.type(), Scalar(0));
    cvtColor(dst, maskWithBox, COLOR_GRAY2BGR);

    // Draw boundaries of patches on the mask
    if (!shifted) {
        this->overlayAllPatchBoundary(maskWithBox, Scalar(50, 50, 50), shifted);
        imwrite("bin/" + FileName + "Mask_withBox_NonShifted.jpg", maskWithBox);
    }
    else {
        this->overlayAllPatchBoundary(maskWithBox, Scalar(50, 50, 50), shifted);
        imwrite("bin/" + FileName + "Mask_withBox_Shifted.jpg", maskWithBox);
    }
}


void ImageHandler::DisplayImageWithBox(String FileName){
    // Draw boundaries of patches on the mask
    Mat maskWithBox = this->im.clone();

    this->overlayAllPatchBoundary(maskWithBox, Scalar(0, 255, 0), false);
    imwrite("bin/" + FileName + "ImageWithBox_NonShifted.jpg", maskWithBox);

    maskWithBox = this->im.clone();
    this->overlayAllPatchBoundary(maskWithBox, Scalar(255, 0, 0), true);
    imwrite("bin/" + FileName + "ImageWithBox_Shifted.jpg", maskWithBox);

}

void ImageHandler::DisplayBGRPatchInImage(const vector<vector<Patch>> &vecVecPatch, String IndexPrefix, String FileName, Mat (Patch::*memberFuncPtr)() const, bool shifted){
    Mat dst(this->im.size(), this->im.type(), Scalar(0));
    for (auto &vecPatch : vecVecPatch) {
        for (auto &patch : vecPatch) {
            int left = patch.x_pixel;
            int top = patch.y_pixel;
            Mat dst_roi = dst(Rect(left, top, PATCHX, PATCHY));

            //getMask().copyTo(dst_roi); // not recognized as pointer-to-function

            //Mat src = patch.getKMeansSeg();
            Mat src = (patch.*memberFuncPtr)();
            src.copyTo(dst_roi);

        }
    }

    // Draw boundaries of patches on the mask
    Mat maskWithBox = dst;
    if (!shifted) {
        imwrite("bin/" + IndexPrefix + FileName + "_NonShifted.jpg", dst);
        this->overlayAllPatchBoundary(maskWithBox, Scalar(50, 50, 50), shifted);
        imwrite("bin/" + IndexPrefix + FileName + "_withBox_NonShifted.jpg", maskWithBox);
    }
    else {
        imwrite("bin/" + IndexPrefix + FileName + "_Shifted.jpg", dst);
        this->overlayAllPatchBoundary(maskWithBox, Scalar(50, 50, 50), shifted);
        imwrite("bin/" + IndexPrefix + FileName + "_withBox_Shifted.jpg", maskWithBox);
    }

}

void ImageHandler::overlayAllPatchBoundary(Mat &dst, const Scalar &color, bool shifted) {
    cv::Mat & image = this->im;
    Rect rangeRect(0, 0, im.cols, im.rows);
    cv::rectangle(dst, rangeRect, color);
    int x = 0, y = 0;
    if (shifted) {
        x = PATCHX / 2;
    }
    for (; x <= image.cols - PATCHX; x += PATCHX){
        y = 0;
        if (shifted) {
            y = PATCHY / 2;
        }
        for (; y <= image.rows - PATCHY; y += PATCHY){
            Rect rangeRect(x, y, PATCHX, PATCHY);
            cv::rectangle(dst, rangeRect, color);

            if (y == 0 || y == PATCHY / 2){                                // put text of x index in the top row
                putText(dst, std::to_string(x / PATCHX + 1), Point(x + PATCHX / 2 - 2, y + 8), FONT_HERSHEY_SIMPLEX, 0.3, Scalar(0, 255, 0));
            }
            if (x == 0 || x == PATCHY / 2){                                // put text of x index in the top row
                putText(dst, std::to_string(y / PATCHX + 1), Point(x + 6, y + PATCHY / 2 - 2), FONT_HERSHEY_SIMPLEX, 0.3, Scalar(0, 255, 0));
            }
        }
    }
}


void ImageHandler::DisplayRRectInImage(String IndexPrefix, String FileName) {
    Mat dst(this->im.size(), this->im.type(), Scalar(0));
    for (auto & vecPatch : this->vecVecPatch) {
        for (auto & patch : vecPatch) {
            // for each rectangle in the patch
            for (auto &rRect : patch.getVecRRect()) {
                Point2f vtx[4];
                rRect.points(vtx);

                // Draw the bounding box
                for (int i = 0; i < 4; i++)
                    line(dst, vtx[i], vtx[(i + 1) % 4], Scalar(0, 255, 0), 1, LINE_AA);

                //            cout << "vtx=" << endl;
                //            for( int i = 0; i < 4; i++ )
                //                cout << vtx[i] << endl;
                //            cout << "size = " << rRect.size.height << " " << rRect.size.width << endl;
                //            cout << "angle = " << rRect.angle << endl;
                //            imshow("show", dst);
                //            char c = waitKey(); // show and pause during execution
            }
        }
    }

    // for the rectangles in shifted patches
    for (auto & vecShiftedPatch : this->vecVecShiftedPatch) {
        for (auto &patch : vecShiftedPatch) {
            // for each rectangle in the patch
            for (auto &rRect : patch.getVecRRect()) {
                Point2f vtx[4];
                rRect.points(vtx);

                // Draw the bounding box
                for (int i = 0; i < 4; i++)
                    line(dst, vtx[i], vtx[(i + 1) % 4], Scalar(255, 0, 0), 1, LINE_AA);

            }
        }
    }

    this->overlayStartRRectInImage(dst, Scalar(0, 255, 0), false);
    this->overlayStartRRectInImage(dst, Scalar(255, 0, 0), true);


    imwrite("bin/" + IndexPrefix + FileName + ".jpg", dst);
    imshow("RRect", dst);
}

void ImageHandler::DisplayRRectInImage(String IndexPrefix, String FileName, const Scalar& color, bool shifted) {
    Mat dst(this->im.size(), this->im.type(), Scalar(0));
    this->overlayAllPatchBoundary(dst, Scalar(50, 50, 50), shifted);

    vector<vector<Patch>> vecvecPatch;
    if (!shifted){
        vecvecPatch = this->vecVecPatch;
    }
    else {
        vecvecPatch = this->vecVecShiftedPatch;
    }
    for (auto & vecPatch : vecvecPatch) {
        for (auto & patch : vecPatch) {
            // for each rectangle in the patch
            for (auto &rRect : patch.getVecRRect()) {
                Point2f vtx[4];
                rRect.points(vtx);

                // Draw the bounding box
                for (int i = 0; i < 4; i++)
                    line(dst, vtx[i], vtx[(i + 1) % 4], Scalar(0, 255, 0), 1, LINE_AA);

                circle( dst, vtx[2], 2, Scalar( 0, 84, 211 ), FILLED, LINE_8 ); // red, temp for display

                rRect.overlapAbstractArrow(dst, Scalar(244, 208, 63)); // mint, temp for display
            }
        }
    }


    if (!shifted) {
        //this->overlayStartRRectInImage(dst, Scalar(0, 255, 0), false);
        imwrite("bin/" + IndexPrefix + FileName + "_NonShifted.jpg", dst);
    }
    else {
        //this->overlayStartRRectInImage(dst, Scalar(255, 0, 0), true);
        imwrite("bin/" + IndexPrefix + FileName + "_Shifted.jpg", dst);
    }
}

void ImageHandler::calcStartRRect(double confThreshold, bool shifted) {
    vector<vector<Patch>> vecvecPatch;
    if (!shifted){
        vecvecPatch = this->vecVecPatch;
    }
    else {
        vecvecPatch = this->vecVecShiftedPatch;
    }
    if (!shifted){
        this->vecStartRRectIndex.clear();
    }
    else {
        this->vecShiftedStartRRectIndex.clear();
    }
    vector<vector<int>> vecStart;

    for (int m = 0; m < vecvecPatch.size(); m++) {
        vector<Patch> &vecPatch = vecvecPatch[m];
        for (int i = 0; i < vecPatch.size(); i++) {
            Patch &patch = vecPatch[i];
            for (int j = 0; j < patch.getVecRRect().size(); j++) {
                RotatedRectm &rRect = patch.getVecRRect()[j];
                if (rRect.confidence >= confThreshold) {
                    // add to index vector
                    vector<int> vecInt;
                    vecInt.push_back(m); // index in x direction
                    vecInt.push_back(i); // index in y direction
                    vecInt.push_back(j); // index of RRect within the patch
                    vecStart.push_back(vecInt);
                }
            }
        }
    }
    if (!shifted){
        this->vecStartRRectIndex = vecStart;
    }
    else {
        this->vecShiftedStartRRectIndex = vecStart;
    }
}

void ImageHandler::overlayStartRRectInImage(Mat & dst, const Scalar& color, bool shifted) {
    vector<vector<int>> vecStart;
    if (!shifted){
        vecStart = this->vecStartRRectIndex;
    }
    else {
        vecStart = this->vecShiftedStartRRectIndex;
    }

    if (vecStart.empty() ) return;


    vector<vector<Patch>> vecvecPatch;
    if (!shifted){
        vecvecPatch = this->vecVecPatch;
    }
    else {
        vecvecPatch = this->vecVecShiftedPatch;
    }

    for (int i = 0; i < vecStart.size(); i++){
        int patchIndx = vecStart[i][0];
        int patchIndy = vecStart[i][1];
        int rRectInd = vecStart[i][2];
        RotatedRectm & rRect = vecvecPatch[patchIndx][patchIndy].getVecRRect()[rRectInd];
            Point2f vtx[4];
            rRect.points(vtx);

            // Draw the bounding box
            for( int i = 0; i < 4; i++ )
                line(dst, vtx[i], vtx[(i+1)%4], color, 2, LINE_AA);

//                imshow("show", dst);
//                char c = waitKey(); // show and pause during execution

    }
}

void ImageHandler::calcCableFromStartRRect() {
    vector<vector<vector<vector<bool>>>> isUsed (2,
            vector<vector<vector<bool>>>(im.cols,
                    vector<vector<bool>>(im.rows,
                            vector<bool>(MAX_RRECT_NUM_IN_PATCH, false))));
    // Deal with non-shifted start points first
    bool shifted = false;
    vector<vector<int>> vecStart = this->vecStartRRectIndex;
    vector<vector<Patch>> vecvecPatch = this->vecVecPatch;
    for (const vector<int> & index : vecStart){
        int x = index[0];
        int y = index[1];
        int k = index[2];
        RotatedRectm & rRect1 = vecvecPatch[x][y].getVecRRect()[k];
        float angle = rRect1.angle;
        vector<vector<int>> vecVecNeighborIndex;
        //! [NONSHIFTED/SHIFTED, x, y, displacement in enum]
        calc_neighbor_patch_index(vecVecNeighborIndex, NONSHIFTED, x, y, angle);

        //TODO delete out-of-range index
        cout << "x = " << x << " y = " << y << " k = " << k <<" angle = " << angle << endl;
//        cout << "vecVecNeighborIndex=" << endl;
//        for (vector<int> vec : vecVecNeighborIndex){
//            for (int i  : vec){
//                cout << i << " ";
//            }
//            cout << endl;
//        }


        auto it = vecVecNeighborIndex.begin();
        while (it != vecVecNeighborIndex.end() ){
            if (this->is_out_of_range(*it)){
                it = vecVecNeighborIndex.erase(it);
                continue;
            }
            ++it;
        }

        cout << "after deleting" << endl;
        cout << "vecVecNeighborIndex=" << endl;
        for (vector<int> & vec : vecVecNeighborIndex){
            for (int i  : vec){
                cout << i << " ";
            }
            cout << endl;
        }

//        //TODO calculate Cost and find min
//        // [NONSHIFTED/SHIFTED, x, y, displacement in enum, k], cost
//        // also delete index for patch with no rRect
//        // cause segmentation fault here
//        auto minIt = vecVecNeighborIndex.end();
//        int minCost = INT_MAX;
//        it = vecVecNeighborIndex.begin();
//        while (it != vecVecNeighborIndex.end() ){ // for each vector<int>
//            auto vec = *it;
//            int x2 = vec[1];
//            int y2 = vec[2];
//            int displacement = vec[3];
//            vector<vector<Patch>> vecvecPatch2; // choose the correct vector to find rRect2
//            if (vec[0] == NONSHIFTED){
//                vecvecPatch2 = this->vecVecPatch;
//            }
//            else {
//                vecvecPatch2 = this->vecVecShiftedPatch;
//            }
//            auto & vecRect = vecvecPatch[x2][y2].getVecRRect();
//            if (vecRect.empty() ){ // no rRects in this patch
//                it = vecVecNeighborIndex.erase(it);
//                continue;
//            }
//            for (RotatedRectm & rRect2 : vecRect){
//                Cost cost(rRect1, rRect2, displacement);
//                float costVal = cost.calcTotalCost(1, 1, 0);
//                // it->push_back(costVal); //TODO type is float
//            }
//
//
//
//
//            ++it;
//        }




        //TODO when add the RRect, set the flag
    }

}

void ImageHandler::calc_neighbor_patch_index(vector<vector<int>> & ret, int shifted, int x, int y, float angle) {
    // return a collection of neighboring patch index [NONSHIFTED/SHIFTED, x, y, displacement in enum]
    //vector<vector<int>> ret;
    if (shifted == NONSHIFTED){
        if (angle < 90){ // to the right bottom
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+1, y, OnePatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x, y+1, OnePatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+1, y+1, OnePatch));

            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x, y+2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+1, y+2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+2, y+2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+2, y+1, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+2, y, TwoPatch));

            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-1, y, HalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x, y-1, HalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x, y, HalfPatch));

            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-1, y+1, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x,   y+1, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x+1, y+1, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x+1, y, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x+1, y-1, OneAndHalfPatch));

        }
        else { // to the left bottom
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x-1, y, OnePatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x, y-1, OnePatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x-1, y-1, OnePatch));

            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x, y-2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x-1, y-2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x-2, y-2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x-2, y-1, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x-2, y, TwoPatch));

            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-1, y, HalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x, y-1, HalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-1, y-1, HalfPatch));

            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x, y+1, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-1, y+1, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-2, y+1, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-2, y, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-2, y-1, OneAndHalfPatch));

        }
    }
    else { // SHIFTED
        if (angle < 90){
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x+1, y, OnePatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x, y+1, OnePatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x+1, y+1, OnePatch));

            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x, y+2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x+1, y+2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x+2, y+2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x+2, y+1, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x+2, y, TwoPatch));

            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+1, y, HalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x, y+1, HalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+1, y+1, HalfPatch));

            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x, y+2, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+1, y+2, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+2, y+2, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+2, y+1, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+2, y, OneAndHalfPatch));
        }
        else {
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-1, y, OnePatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x, y-1, OnePatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-1, y-1, OnePatch));

            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x, y-2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-1, y-2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-2, y-2, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-2, y-1, TwoPatch));
            ret.push_back(helper_calc_neighbor_patch_index(SHIFTED, x-2, y, TwoPatch));

            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+1, y, HalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x, y+1, HalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x, y, HalfPatch));

            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x+1, y+2, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x, y+2, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x-1, y+2, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x-1, y+1, OneAndHalfPatch));
            ret.push_back(helper_calc_neighbor_patch_index(NONSHIFTED, x-1, y, OneAndHalfPatch));
        }
    }

    auto it = ret.begin();
    while (it != ret.end() ){
        if (this->is_out_of_range(*it)){
            it = ret.erase(it);
            continue;
        }
        ++it;
    }
}

vector<int> ImageHandler::helper_calc_neighbor_patch_index(int shifted, int x, int y, int displacement) {
    vector<int> ret = {shifted, x, y, displacement};
    return ret;
}

bool ImageHandler::is_out_of_range(const vector<int> & vec) {
    if (vec[0] == NONSHIFTED) {
        if (vec[1] < 0 || vec[1] >= this->patchColNum) return true; // x
        if (vec[2] < 0 || vec[2] >= this->patchRowNum) return true; // y
    }
    if (vec[0] == SHIFTED) {
        if (vec[1] < 0 || vec[1] >= this->shiftedPatchColNum) return true; // x
        if (vec[2] < 0 || vec[2] >= this->shiftedPatchRowNum) return true; // y
    }
    return false;
}

void ImageHandler::saveToAllPatch() {
    this->allPatch.push_back(this->vecVecPatch);
    this->allPatch.push_back(this->vecVecShiftedPatch);
}

void ImageHandler::rRect2Node() {
    for (int s1 = 0; s1 <= 1; ++s1){
        vector<vector<Patch>> & vecvecPatch = allPatch[s1];
        for (int x1 = 0; x1 < vecvecPatch.size(); ++x1){
            vector<Patch> & vP = vecvecPatch[x1];
            for (int y1 = 0; y1 < vP.size(); ++y1){
                Patch & patch1 = vP[y1];                                            // for each patch1 [s1, x1, y1]
                for (int r1 = 0; r1 < patch1.getVecRRect().size(); ++r1){
                    RotatedRectm & rRect1 = patch1.getVecRRect()[r1];               // for each rRect1 in patch1
                    vector<vector<int>> vecVecNeighborIndex;
                    calc_neighbor_patch_index(vecVecNeighborIndex, s1, x1, y1, rRect1.angle);
                    for (vector<int> vec : vecVecNeighborIndex){
                        int s2 = vec[0], x2 = vec[1], y2 = vec[2], displacement = vec[3];
//                        if (s2 > 1 || x2 >= allPatch[s2].size() || y2 >= allPatch[s2][x2].size()){
//                            cout << "neighbor index out of range" << endl;
//                            continue;
//                        }
                        Patch & patch2 = allPatch[s2][x2][y2];                      // for each neighboring patch2 (risky)
//                        cout << "x1: " << x1 << " y1:" << y1<<endl;
                        for (int r2 = 0; r2 < patch2.getVecRRect().size(); ++r2){
                            RotatedRectm & rRect2 = patch2.getVecRRect()[r2];       // for each rRect2 in patch2
                            Cost cost(rRect1, rRect2, displacement);                // calculate cost with 1.0 1.0
                            float costVal = cost.calcTotalCost(1.0, 1.0, displacement);
                            int v1 = this->graphHandler.sxyrIndex2vIndex(s1, x1, y1, r1);
                            int v2 = this->graphHandler.sxyrIndex2vIndex(s2, x2, y2, r2);
                            //graphHandler.addEdge(v1, v2, costVal);
                            graphHandler.addEdge(v1, v2, costVal,
                                    rRect1.getPoints8()[4].x, rRect1.getPoints8()[4].y,
                                                 rRect2.getPoints8()[4].x, rRect2.getPoints8()[4].y);
#ifdef DEBUG_PRINT_COSTVAL
                            cout << "costVal: "<<costVal<<endl;
#endif
//                            cout << "s1, x1, y1, r1: " << s1 << " " << x1 << " " << y1 << " " << r1 << endl;
//                            cout << "s2, x2, y2, r2: " << s2 << " " << x2 << " " << y2 << " " << r2 << endl;
//                            auto vec1 = this->graphHandler.vIndex2sxyrIndex(v1);
//                            auto vec2 = this->graphHandler.vIndex2sxyrIndex(v2);
//                            cout << "vec1: ";
//                            for (auto i : vec1) cout << i << " ";
//                            cout << "vec2: " << endl;
//                            for (auto i : vec2) cout << i << " ";
//                            cout << endl;
                        }
                    }
                }
            }
        }
    }
    graphHandler.minimumSpanningTree();
    graphHandler.saveDotImage(graphHandler.g, "bin/graph_" + std::to_string(id) + ".dot");
    graphHandler.visualizeGraph(graphHandler.g, "bin/Graph_" + std::to_string(id) + ".jpg");
    graphHandler.cutBranch("bin/gMST_" + std::to_string(id) + ".jpg", "bin/gCut_" + std::to_string(id) + ".jpg");
}

vector<RotatedRectm> &ImageHandler::patch2rRect(vector<vector<Patch>> &vecVecPatch, const int K, const int margin) {
    for (auto &vecPatch : vecVecPatch) {
        for (auto &patch : vecPatch) { // range based for-loop, use reference for modifying objects in vector
            Mat kMeansSeg = patch.kMeans(K);
            //imwrite("bin/Patch_" + std::to_string(patch.x) + "_" + std::to_string(patch.y) + "_kMeansSeg.jpg",
            //kMeansSeg);
            //Mat kMask = patch.closeBGRMask(margin, vecRefBGR);
            Mat kMask = patch.BGRInRangeMask(margin, vecRefBGR);
            //imwrite("bin/Patch_" + std::to_string(patch.x) + "_" + std::to_string(patch.y) + "_kMask.jpg", kMask);
            Mat contourDrawing = patch.closestBGRMaskContour();
            //imwrite("bin/" +  std::to_string(i)+"_" + "Patch_" + std::to_string(patch.x) + "_" + std::to_string(patch.y) + "_contourDrawing.jpg", contourDrawing);

            //vector<RotatedRectm> vecR = patch.findCompMinAreaRectImageLevel(this->im.rows, this->im.cols, 0.1, 0.5, 0.1); // TODO change mechanism
            vector<RotatedRectm> vecR = patch.findMinAreaRect(this->im.rows, this->im.cols, this->canny, this->vecRefWidth,
                                                              MIN_RECTANGLE_WIDTH, MIN_RECTANGLE_LENGTH,
                                                              MIN_AREA_PERCENTAGE, MAX_WIDTH_DEVIATION,
                                                              MIN_CANNY_LEFT_CONFIDENCE, MIN_CANNY_RIGHT_CONFIDENCE, 0.0);

            for (auto rRect : vecR){
                this->vecRRect.push_back(rRect);
            }
        }
    }
    return this->vecRRect;
}

vector<RotatedRectm> &
ImageHandler::shiftedPatch2rRect(vector<vector<Patch>> &vecVecPatch, const int K, const int margin) {
    for (auto & vecShiftedPatch : vecVecPatch) {
        for (auto &patch : vecShiftedPatch) { // range based for-loop, use reference for modifying objects in vector
            Mat kMeansSeg = patch.kMeans(K);
            //imwrite("bin/Patch_" + std::to_string(patch.x) + "_" + std::to_string(patch.y) + "_kMeansSeg.jpg",
            //kMeansSeg);
            //Mat kMask = patch.closeBGRMask(margin, vecRefBGR);
            Mat kMask = patch.BGRInRangeMask(margin, vecRefBGR);
            //imwrite("bin/Patch_" + std::to_string(patch.x) + "_" + std::to_string(patch.y) + "_kMask.jpg", kMask);
            Mat contourDrawing = patch.closestBGRMaskContour();
            //imwrite("bin/" +  std::to_string(i)+"_" + "Patch_" + std::to_string(patch.x) + "_" + std::to_string(patch.y) + "_contourDrawing.jpg", contourDrawing);

            //vector<RotatedRectm> vecR = patch.findCompMinAreaRectImageLevel(this->im.rows, this->im.cols, 0.1, 0.5, 0.1);
            vector<RotatedRectm> vecR = patch.findMinAreaRect(this->im.rows, this->im.cols, this->canny, this->vecRefWidth,
                                                              MIN_RECTANGLE_WIDTH, MIN_RECTANGLE_LENGTH,
                                                              MIN_AREA_PERCENTAGE, MAX_WIDTH_DEVIATION,
                                                              MIN_CANNY_LEFT_CONFIDENCE, MIN_CANNY_RIGHT_CONFIDENCE, 0.0);

            for (auto rRect : vecR){
                this->vecShiftedRRect.push_back(rRect);
            }
        }
    }
    return this->vecShiftedRRect;
}

Mat &ImageHandler::imCanny(double thres1, double thres2) {
    Canny(this->im, this->canny, thres1, thres2, 3, true);  // true for L2 norm gradient
    imshow("Canny", canny);
    return this->canny;
}

vector<MergedComp> &ImageHandler::rRect2InitialMergedComp(vector<RotatedRectm> &vecRRect) {
    for (int i = 0; i < vecRRect.size(); i++){
        MergedComp comp(i, vecRRect[i].getPoints8());
        this->vecMergedComp.push_back(comp);
    }
    return this->vecMergedComp;
}

vector<MergedComp> &ImageHandler::shiftedRRect2InitialMergedComp(vector<RotatedRectm> &vecRRect) {
    for (int i = 0; i < vecRRect.size(); i++){
        MergedComp comp(i, vecRRect[i].getPoints8());
        this->vecShiftedMergedComp.push_back(comp);
    }
    return this->vecShiftedMergedComp;
}

void ImageHandler::mergeComp(vector<MergedComp> &vecMergedComp, int iterationN, float maxCostToMerge, bool shifted) {
    int iteration = 0;
    float minCost = -1;
    bool minToMergeTail1 = false, minToMergeTail2 = false;

    Mat dst(this->im.size(), this->im.type(), Scalar(0)); // to display
    overlayAllRRect(dst, shifted, 1);

    std::string winNameMergeString = (shifted) ? winNameMergeShifted : winNameMerge;
    imshow(winNameMergeString, dst);
    //char c = waitKey(0); // uncomment to show rectangles before merge


    while (iteration < iterationN && minCost <= maxCostToMerge){

        if (vecMergedComp.size() < 2) {
            cout << "vecMergedComp.size() < 2, break.\n";
            break;
        }

        minCost = -1;
        minToMergeTail1 = false, minToMergeTail2 = false;

        cout << "iteration:" << iteration << endl;

        int idx1 = 0, idx2 = 1;
        for (int i = 0; i < vecMergedComp.size(); i++){
            for (int j = 0; j < i; j++){                                    // for each pair of merged components

                if (vecMergedComp.size() < 2) {
                    cout << "vecMergedComp.size() < 2, break.\n";
                    break;
                }

#if 1
                Cost cost(vecMergedComp[j], vecMergedComp[i]);
#endif
#if 0
                // debug
                vector<RotatedRectm> & vecRect = (shifted) ? this->vecShiftedRRect : this->vecRRect;
                Cost cost( vecRect , vecMergedComp[j], vecMergedComp[i]);
#endif
                float costVal = cost.getVal();
                if (costVal < minCost
                        || minCost == -1){                                  // first visit
                    minCost = costVal;                                      // find the minimal cost
                    idx1 = j;   // the smaller index
                    idx2 = i;
                    minToMergeTail1 = cost.toMergeTail1;
                    minToMergeTail2 = cost.toMergeTail2;

#if 0
                    cout << "iteration:" << iteration << " i:" << i << " j:" << j << " vecMergedComp.size():" << vecMergedComp.size()<< endl
                        << "\tminCost:" << minCost << " idx1:" << idx1 <<" idx2:" << idx2 << endl;
#endif
                }
            }
        }


        if (vecMergedComp.size() < 2) {
            cout << "vecMergedComp.size() < 2, break.\n";
            break;
        }

        cout << "iteration:" << iteration << " vecMergedComp.size():" << vecMergedComp.size()<< endl
             << "\tminCost:" << minCost << " idx1:" << idx1 <<" idx2:" << idx2 << endl;



        if (minCost > maxCostToMerge) break;                                // exit the loop

        if (idx1 >= idx2){
            std::cerr << "ImageHandler::mergeComp idx1, idx2 wrong value. idx1: "<< idx1 << " idx2: " << idx2 << endl;
        }

        Cost cost(dst, vecMergedComp[idx1], vecMergedComp[idx2]);   // display the Bezier curve

        vecMergedComp[idx1].merge(vecMergedComp[idx2],
                minToMergeTail1, minToMergeTail2);                          // merge to vecMergedComp[idx1]
        vecMergedComp.erase(vecMergedComp.begin() + idx2,
                vecMergedComp.begin() + idx2 + 1);                          // delete vecMergedComp[idx2] from the vector


        overlayComp(dst, idx1, shifted, 2);   // display the merged one
        //imshow(winNameMergeString, dst);
        //char c = waitKey(0);  // uncomment to see adding rectangles step-by-step

        iteration ++;
    }

    imshow(winNameMergeString, dst);
    imwrite("bin/" + std::to_string(this->id) + "_" + winNameMergeString + ".jpg", dst);
    //char c1 = waitKey(0); // uncomment to see instance result specified by rectangles

}

void ImageHandler::mergeComp(vector<MergedComp> &vecMergedComp, int iterationN, float maxCostToMerge, bool shifted,
                             bool second) {
    int iteration = 0;
    float minCost = -1;
    bool minToMergeTail1 = false, minToMergeTail2 = false;

    Mat dst(this->im.size(), this->im.type(), Scalar(0)); // to display
    overlayAllRRect(dst, shifted, true, 1); // no text version
    overlayCompIndex(dst, vecMergedComp, shifted);

    std::string winNameMergeString = (shifted) ? winNameMergeShifted : winNameMerge;
    //imshow(winNameMergeString, dst);
    //char c = waitKey(0); // uncomment to show rectangles before merge


    while (iteration < iterationN && minCost <= maxCostToMerge){

        minCost = -1;
        minToMergeTail1 = false, minToMergeTail2 = false;

        cout << "iteration:" << iteration << endl;

        int idx1 = 0, idx2 = 1;
        for (int i = 0; i < vecMergedComp.size(); i++){
            for (int j = 0; j < i; j++){                                    // for each pair of merged components

                if (vecMergedComp.size() < 2) {
                    cout << "vecMergedComp.size() < 2, break.\n";
                    break;
                }
#if 1
                Cost cost(vecMergedComp[j], vecMergedComp[i], second);
#endif
#if 0
                // debug
                cout << "\t\ti:" << i << " j:" << j << endl;
                vector<RotatedRectm> & vecRect = (shifted) ? this->vecShiftedRRect : this->vecRRect;
                Cost cost( vecRect , vecMergedComp[j], vecMergedComp[i]);
#endif
                float costVal = cost.getVal();
                if (costVal < minCost
                    || minCost == -1){                                  // first visit
                    minCost = costVal;                                      // find the minimal cost
                    idx1 = j;   // the smaller index
                    idx2 = i;
                    minToMergeTail1 = cost.toMergeTail1;
                    minToMergeTail2 = cost.toMergeTail2;

#if 0
                    cout << "iteration:" << iteration << " i:" << i << " j:" << j << " vecMergedComp.size():" << vecMergedComp.size()<< endl
                        << "\tminCost:" << minCost << " idx1:" << idx1 <<" idx2:" << idx2 << endl;
#endif
                }
            }
        }


        if (vecMergedComp.size() < 2) {
            cout << "vecMergedComp.size() < 2, break.\n";
            break;
        }

        cout << "iteration:" << iteration << " vecMergedComp.size():" << vecMergedComp.size()<< endl
             << "\tminCost:" << minCost << " idx1:" << idx1 <<" idx2:" << idx2 << endl;

        if (minCost > maxCostToMerge) break;                                // exit the loop

        if (idx1 >= idx2){
            std::cerr << "ImageHandler::mergeComp idx1, idx2 wrong value. idx1: "<< idx1 << " idx2: " << idx2 << endl;
        }

        Cost cost(dst, vecMergedComp[idx1], vecMergedComp[idx2], second);   // display the Bezier curve

        vecMergedComp[idx1].merge(vecMergedComp[idx2],
                                  minToMergeTail1, minToMergeTail2);                          // merge to vecMergedComp[idx1]
        vecMergedComp.erase(vecMergedComp.begin() + idx2,
                            vecMergedComp.begin() + idx2 + 1);                          // delete vecMergedComp[idx2] from the vector


        overlayComp(dst, idx1, shifted, 2);   // display the merged one
        //imshow(winNameMergeString + " extend", dst);
        //char c1 = waitKey(0);   // uncomment to see adding rectangles step-by-step

        iteration ++;
    }

    imshow(winNameMergeString + " extend", dst);
    imwrite("bin/" + std::to_string(this->id) + "_" + winNameMergeString + "_extend" + ".jpg", dst);
    //char c1 = waitKey(0); // uncomment to see instance result specified by rectangles
}

void ImageHandler::overlayAllRRect(Mat &dst, bool shifted, const int thickness, const Scalar &color) {
#if 1
    vector<RotatedRectm> vecRect;
    if (!shifted){
        vecRect = this->vecRRect;
    }
    else {
        vecRect = this->vecShiftedRRect;
    }
    for (int i = 0; i < vecRect.size(); ++i) {
        RotatedRectm & rRect = vecRect[i];
        Point2f vtx[4];
        rRect.points(vtx);

        // Draw the bounding box
        for (int i = 0; i < 4; i++)
            line(dst, vtx[i], vtx[(i + 1) % 4], color, thickness, LINE_AA);

        //circle( dst, vtx[2], 2, Scalar( 0, 84, 211 ), FILLED, LINE_8 ); // red, temp for display

        rRect.overlapAbstractArrow(dst, Scalar(102, 51, 0)); // dark blue, temp for display

        putText(dst, std::to_string(i), vtx[0], FONT_HERSHEY_COMPLEX, 0.5, Scalar(255, 255, 255));
    }

#endif
}

void ImageHandler::overlayAllRRect(Mat &dst, bool shifted, bool noText, const int thickness, const Scalar &color) {
    vector<vector<Patch>> vecvecPatch;
    if (!shifted){
        vecvecPatch = this->vecVecPatch;
    }
    else {
        vecvecPatch = this->vecVecShiftedPatch;
    }
    for (auto & vecPatch : vecvecPatch) {
        for (auto & patch : vecPatch) {
            // for each rectangle in the patch
            for (auto &rRect : patch.getVecRRect()) {
                Point2f vtx[4];
                rRect.points(vtx);

                // Draw the bounding box
                for (int i = 0; i < 4; i++)
                    line(dst, vtx[i], vtx[(i + 1) % 4], color, thickness, LINE_AA);

                //circle( dst, vtx[2], 2, Scalar( 0, 84, 211 ), FILLED, LINE_8 ); // red, temp for display

                rRect.overlapAbstractArrow(dst, Scalar(102, 51, 0)); // dark blue, temp for display
            }
        }
    }
}

void ImageHandler::overlayCompIndex(Mat &dst, vector<MergedComp> &vecMergedComp, bool shifted) {
    vector<RotatedRectm> vecRect;
    if (!shifted){
        vecRect = this->vecRRect;
    }
    else {
        vecRect = this->vecShiftedRRect;
    }
    for (int i = 0; i < vecMergedComp.size(); ++i){
        MergedComp & comp = vecMergedComp[i];
        if (!comp.listIdx.empty()){
            int j = comp.listIdx.front();
            int k = comp.listIdx.back();
            {
                RotatedRectm &rRect = vecRect[j];
                Point2f vtx[4];
                rRect.points(vtx);
                putText(dst, std::to_string(i)+'f', vtx[0], FONT_HERSHEY_COMPLEX, 0.5, Scalar(255, 255, 255));
            }
            {
                RotatedRectm &rRect = vecRect[k];
                Point2f vtx[4];
                rRect.points(vtx);
                putText(dst, std::to_string(i)+'b', vtx[0], FONT_HERSHEY_COMPLEX, 0.5, Scalar(255, 255, 255));
            }

        }
    }
}

RNG rng(12345);
void ImageHandler::overlayComp(Mat &dst, int compIdx, bool shifted, const int thickness) {

    Scalar color = Scalar( rng.uniform(10, 256), rng.uniform(10, 256), rng.uniform(10, 256) );

    vector<MergedComp> vecComp;
    vector<RotatedRectm> vecRect;
    if (!shifted){
        vecComp = this->vecMergedComp;
        vecRect = this->vecRRect;
    }
    else {
        vecComp = this->vecShiftedMergedComp;
        vecRect = this->vecShiftedRRect;
    }
    if (vecComp.size() <= compIdx){
        return; // TODO
        std::cerr << "ImageHandler::overlayComp compIdx " << compIdx << " exceeds vecComp.size() " << vecComp.size() << endl;
    }
    MergedComp & comp = vecComp[compIdx];
    for (int i : comp.listIdx){
        RotatedRectm & rRect = vecRect[i];
        Point2f vtx[4];
        rRect.points(vtx);

        // Draw the bounding box
        for (int i = 0; i < 4; i++)
            line(dst, vtx[i], vtx[(i + 1) % 4], color, thickness, LINE_AA);
    }
}

void
ImageHandler::mergedComp2InstanceMask(vector<MergedComp> &vecMergedComp, vector<RotatedRectm> &vecRRect, bool shifted,
        int minNumRect) {
    int instanceIndex = 0;
    //Mat dst(this->im.size(), this->im.type(), Scalar(0)); // to display
    Mat dst = this->im.clone();

    std::string winNameMergeString = (shifted) ? winNameInstanceShifted : winNameInstance;

    for (MergedComp & comp : vecMergedComp){
        if (comp.listIdx.empty()){
            std::cerr << "ImageHandler::mergedComp2InstanceMask comp.listIdx empty.\n";
        }
        if (comp.listIdx.size() < minNumRect){                  // skip too short instances
            continue;
        }
        Mat instanceMask(this->im.size(), CV_8UC1, Scalar(0));
        for (int i : comp.listIdx){
            RotatedRectm & rRect = vecRRect[i];
            Mat instanceMask_roi = instanceMask(Rect(rRect.x_pixel, rRect.y_pixel, PATCHX, PATCHY));
            bitwise_or(instanceMask_roi, rRect.mask, instanceMask_roi);

        }
        comp.mask = instanceMask.clone();
        overlayInstanceMaskWithIndex(dst, instanceMask, ++instanceIndex);

//        imshow(winNameMergeString, dst);
//        char c1 = waitKey(0);         // uncomment to see step-by-step result
    }
    imshow(winNameMergeString, dst);
    imwrite("bin/" + std::to_string(this->id) + "_" + winNameMergeString + ".jpg", dst);
    //char c1 = waitKey(0);         // uncomment to pause here
}

void ImageHandler::overlayInstanceMaskWithIndex(Mat &dst, Mat &instanceMask, int instanceIndex) {
    ////https://docs.opencv.org/3.4.3/df/d0d/tutorial_find_contours.html
    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;
    findContours( instanceMask, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE );

    Scalar color = Scalar( rng.uniform(10, 256), rng.uniform(10, 256), rng.uniform(10, 256) ); // lighter color

    for( size_t i = 0; i< contours.size(); i++ )
    {
        drawContours( dst, contours, (int)i, color, FILLED, LINE_8, hierarchy, 0 );
        putText(dst, std::to_string(instanceIndex), contours[i][0], FONT_HERSHEY_COMPLEX, 0.5, Scalar(255, 255, 255));
    }
}

vector<CombinedComp> &ImageHandler::mergedComp2InitialCombinedComp(vector<MergedComp> &vecMergedComp,
                                                                   vector<MergedComp> &vecMergedCompShifted) {
    for (int i = 0; i < vecMergedComp.size(); ++i){
        MergedComp & mComp = vecMergedComp[i];
        if (!mComp.mask.empty()) {
            CombinedComp cComp(i, false, mComp.mask);
            this->vecCombinedComp.push_back(cComp);
        }
    }
    for (int i = 0; i < vecMergedCompShifted.size(); ++i){
        MergedComp & mComp = vecMergedCompShifted[i];
        if (!mComp.mask.empty()) {
            CombinedComp cComp(i, true, mComp.mask);
            this->vecCombinedComp.push_back(cComp);
        }
    }
    return this->vecCombinedComp;
}

void ImageHandler::combineComp(vector<CombinedComp> &vecCombinedComp, int iterationN, float minAreaToCombine) {
    int iteration = 0;
    float maxANDArea = FLT_MAX;
    int ind1, ind2;

    while (iteration < iterationN && maxANDArea >= minAreaToCombine){
        maxANDArea = 0;
        ind1 = -1;
        ind2 = -1;
        for (int i = 0; i < vecCombinedComp.size(); ++i){
            for (int j = 0; j < i; ++j){
                CombinedComp & comp1 = vecCombinedComp[i];
                CombinedComp & comp2 = vecCombinedComp[j];
                Mat Mask_AND_result;
                bitwise_and(comp1.mask, comp2.mask, Mask_AND_result);


                Mat labels, stats, centroids;

                // find connected components
                connectedComponentsWithStats(Mask_AND_result, labels, stats, centroids, 8, CV_16U);

//                int maxArea = 0; // max area of this AND result
                int maxSumArea = 0;
                for (int m = 1; m < stats.rows; m++) { // ignore stats[0], it is for background
                    int area = stats.at<int>(m, CC_STAT_AREA);
//                    if (area > maxArea){
//                        maxArea = area;
//                    }
                    maxSumArea += area;
                }

//                if (maxArea > maxANDArea || ind1 == -1){
//                    ind1 = j;
//                    ind2 = i;
//                    maxANDArea = maxArea;
//                }
                if (maxSumArea > maxANDArea || ind1 == -1){
                    ind1 = j;
                    ind2 = i;
                    maxANDArea = maxSumArea;
                }

            }
        }

        if (maxANDArea < minAreaToCombine){         // too small to combine
            break;
        }

        vecCombinedComp[ind1].combine(vecCombinedComp[ind2]);
        vecCombinedComp.erase(vecCombinedComp.begin() + ind2,
                              vecCombinedComp.begin() + ind2 + 1);

        iteration++;
    }
}

void ImageHandler::combinedComp2InstanceMask(vector<CombinedComp> &vecCombinedComp) {
    int instanceIndex = 0;
    //Mat dst(this->im.size(), this->im.type(), Scalar(0)); // to display
    Mat dst = this->im.clone();

    for (CombinedComp & cComp : vecCombinedComp){
        overlayInstanceMaskWithIndex(dst, cComp.mask, ++instanceIndex);
        vecPredictionMask.push_back(cComp.mask);
    }
    imwrite("bin/" + std::to_string(this->id) + "_" + winNameCombineString + ".jpg", dst);

#ifdef PAUSE_AND_SHOW_COMBINED
    imshow(winNameCombineString, dst);
    char c1 = waitKey(0);
#endif
}





