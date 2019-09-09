//
// Created by ubuntu on 3/13/19.
//

#ifndef TRYENV_IMAGEHANDLER_H
#define TRYENV_IMAGEHANDLER_H

#include "Config.h"
#include "Patch.h"
#include "Cost.h"
#include "GraphHandler.h"
#include "MergedComp.h"
#include "CombinedComp.h"

class ImageHandler {
    Mat im;
    int id;
    vector<vector<Patch>> vecVecPatch;
    int patchColNum = 0, patchRowNum = 0;
    vector<vector<Patch>> vecVecShiftedPatch;
    int shiftedPatchColNum = 0, shiftedPatchRowNum = 0;

    vector<RotatedRectm> vecRRect;  // index by integer
    vector<RotatedRectm> vecShiftedRRect;

    vector<MergedComp> vecMergedComp;
    vector<MergedComp> vecShiftedMergedComp;

    Mat canny;
    vector<float> vecRefWidth;
    vector<Vec3f> vecRefBGR;

    vector<CombinedComp> vecCombinedComp;

    vector< vector<vector<Patch>> > allPatch; // index by [s, x, y]
    GraphHandler graphHandler;

    vector<vector<int>> vecStartRRectIndex; // [x, y, rRectInd]
    vector<vector<int>> vecShiftedStartRRectIndex;

    vector<int> helper_calc_neighbor_patch_index(int shifted, int x, int y, int displacement);

    // Helper function for calcCableFromStartRRect
    // REQUIRES: info about the patch where rRect1 locates: NONSHIFTED/SHIFTED, x, y, and angle of rRect1
    // EFFECTS: return a collection of neighboring patch index [NONSHIFTED/SHIFTED, x, y, displacement in enum]
    void calc_neighbor_patch_index(vector<vector<int>> & ret, int shifted, int x, int y, float angle);

    bool is_out_of_range(const vector<int> & vec);

public:
    vector<Mat> vecPredictionMask;

    ImageHandler(Mat &image) : im(image) {}
    ImageHandler(const Mat &image, const int i) : im(image.clone()), id(i) {
        vecRefWidth.push_back(DEFAULT_CABLE_WIDTH_1); // set a default value
        vecRefWidth.push_back(DEFAULT_CABLE_WIDTH_2);
        vecRefBGR.push_back(DEFAULT_CABLE_COLOR_1);
        vecRefBGR.push_back(DEFAULT_CABLE_COLOR_2);
    }
    ImageHandler(const Mat &image, const int i, const vector<float> & vecWidth, const vector<Vec3f> & vecBGR) : im(image.clone()), id(i) {
        vecRefWidth = vecWidth; // to be called by user, or inherit from the previous frame
        if (vecRefWidth.empty()){
            vecRefWidth.push_back(DEFAULT_CABLE_WIDTH_1); // set a default value
            vecRefWidth.push_back(DEFAULT_CABLE_WIDTH_2);
        }
        vecRefBGR = vecBGR;
        if (vecRefBGR.empty()){
            vecRefBGR.push_back(DEFAULT_CABLE_COLOR_1);
            vecRefBGR.push_back(DEFAULT_CABLE_COLOR_2);
        }
    }
    ImageHandler(const Mat &image, const int i, vector<float> para) : im(image.clone()), id(i) {
        if (vecRefWidth.empty()){
            vecRefWidth.push_back(para[1]);
            vecRefWidth.push_back(para[2]);
        }
        if (vecRefBGR.empty()){
            vecRefBGR.push_back(Vec3f({para[3], para[4], para[5]}));
            vecRefBGR.push_back(Vec3f({para[6], para[7], para[8]}));
        }
    }
    ~ImageHandler() = default;

    // REQUIRES: this->im
    // EFFECTS: resize this->im to MAX_LONG_SIDE
    Mat Resize();

    // REQUIRES: this->im resized
    // EFFECTS: this->vecPatch, split the image into patches, save and return the vector
    // this->patchColNum, this->patchRowNum: number of patches in the column / row
    vector<vector<Patch>> & SplitToPatch();

    // REQUIRES: this->im resized
    // EFFECTS: this->vecShiftedPatch, split the image into patches, save and return the vector
    // this->shiftedPatchColNum, this->shiftedPatchRowNum: number of patches in the column / row
    // splitting of patches is shifted PATCHX/2 to the right, PATCHY/2 to the bottom
    vector<vector<Patch>> & SplitToShiftedPatch();

    // REQUIRES: processed masks in the vector
    // EFFECTS: assemble the masks to an image to display, save to a file
    void DisplayMaskInImage(const vector<vector<Patch>> &vecVecPatch, String FileName, Mat (Patch::*getMask)() const = &Patch::getClosestColorMask, bool shifted = false);

    // REQUIRES: this->im
    // EFFECTS: display all boundaries of patches, save to a file
    void DisplayImageWithBox(String FileName);

    // REQUIRES: processed k-means segmentation or other BGR patches in the vector
    // EFFECTS: assemble the patch segmentation or other BGR patches to an image to display, save to a file
    void DisplayBGRPatchInImage(const vector<vector<Patch>> &vecVecPatch, String IndexPrefix, String FileName, Mat (Patch::*memberFuncPtr)() const, bool shifted = false);

    // REQUIRES: Vector of rotated rectangles are found in patches
    // EFFECTS: Display rotated rectangles in the image, save to a file
    // Display rectangles from both non-shifted patches and shifted patches
    void DisplayRRectInImage(String IndexPrefix, String FileName);

    // REQUIRES: Vector of rotated rectangles are found in patches
    // EFFECTS: Display rotated rectangles in the image, save to a file
    // Display rectangles from either non-shifted patches or shifted patches
    void DisplayRRectInImage(String IndexPrefix, String FileName, const Scalar& color, bool shifted);

    // with index
    void overlayAllRRect(Mat & dst, bool shifted, const int thickness = 1, const Scalar& color = Scalar(102, 153, 51));
    // without index
    void overlayAllRRect(Mat & dst, bool shifted, bool noText, const int thickness = 1, const Scalar& color = Scalar(102, 153, 51));
    // overlay the component index
    void overlayCompIndex(Mat & dst, vector<MergedComp> &vecMergedComp, bool shifted);

    // REQUIRES: dst, 8U3C; color, e.g. green is cv::Scalar(0, 255, 0); non-shifted or shifted
    // EFFECTS: overlay all patch boundaries on dst
    // Called in previous Display functions
    void overlayAllPatchBoundary(Mat & dst, const Scalar& color, bool shifted = false);

    // REQUIRES: confidence threshold
    // EFFECTS: a vector of identifiers of starting point rectangles,
    // starting point rectangles are defined by large enough confidence score
    void calcStartRRect(double confThreshold, bool shifted = false);

    // REQUIRES: indices of starting point rectangles are found
    // EFFECTS: Overlay starting point rotated rectangles in the image, save to a file
    void overlayStartRRectInImage(Mat & dst, const Scalar& color, bool shifted = false);

    // REQUIRES: this->calcStartRRect is executed
    // EFFECTS: search around start points and add the best into cable, iteratively
    void calcCableFromStartRRect();

    // REQUIRES: all operations are done on vecVecPatch and vecVecShiftedPatch
    // EFFECTS: save them to allPatch
    void saveToAllPatch();

    // REQUIRES: patches with rRects are saved in allPatch
    // EFFECTS: add the cost between two neighboring rRects to the graph in graphHandler
    void rRect2Node();

    // REQUIRES: none
    // EFFECTS: this->canny
    Mat & imCanny(double thres1, double thres2);

    // REQUIRES: SplitToPatch() is executed, imCanny() is executed
    // EFFECTS: do processing to patches (vecVecPatch or vecVecShiftedPatch),
    // store the rotated rectangles in corresponding patches,
    // and store to vector of RRects this->vecRRect
    vector<RotatedRectm> & patch2rRect(vector<vector<Patch>> & vecVecPatch, const int K, const int margin);

    // REQUIRES: SplitToShiftedPatch() is executed, imCanny() is executed
    // EFFECTS: do processing to patches (vecVecPatch or vecVecShiftedPatch),
    // store the rotated rectangles in corresponding patches,
    // and store to vector of RRects this->vecShiftedRRect
    vector<RotatedRectm> & shiftedPatch2rRect(vector<vector<Patch>> & vecVecPatch, const int K, const int margin);

    // REQUIRES: patch2rRect() executed
    // EFFECTS: create initial this->vecMergedComp, and this->vecShiftedMergedComp
    // and return its reference
    vector<MergedComp> & rRect2InitialMergedComp(vector<RotatedRectm> & vecRRect);

    // REQUIRES: patch2rRect() executed
    // EFFECTS: create initial this->vecShiftedMergedComp,
    // and return its reference
    vector<MergedComp> & shiftedRRect2InitialMergedComp(vector<RotatedRectm> & vecRRect);

    // REQUIRES: rRect2InitialMergedComp() executed for & vecMergedComp, not the internal one
    // EFFECTS: merge component N times
    void mergeComp(vector<MergedComp> & vecMergedComp, int iterationN, float maxCostToMerge, bool shifted);

    // REQUIRES: rRect2InitialMergedComp() executed for & vecMergedComp, not the internal one
    // EFFECTS: merge component N times, to extend
    void mergeComp(vector<MergedComp> & vecMergedComp, int iterationN, float maxCostToMerge, bool shifted, bool second);

    // REQUIRES: shiftedRRect2InitialMergedComp() executed
    // EFFECTS: merge component N times
    void mergeShiftedComp(vector<MergedComp> & vecMergedComp, int iterationN, float maxCostToMerge);


    void overlayComp(Mat & dst, int compIdx, bool shifted, const int thickness);

    // REQUIRES: mergeComp() executed for & vecMergedComp
    // patch2rRect() or shiftedPatch2rRect() executed for & vecRRect
    // if shifted: shifted
    // number of rectangle filter for mask: minNumRect
    // vecMergedComp must correspond to vecRRect
    // EFFECTS: turn each merged component to instance mask
    // save the instance mask to comp.mask
    void mergedComp2InstanceMask(vector<MergedComp> & vecMergedComp, vector<RotatedRectm> & vecRRect, bool shifted,
            int minNumRect = 0);

    // called in mergedComp2InstanceMask()
    // index not implemented yet
    void overlayInstanceMaskWithIndex(Mat & dst, Mat & instanceMask, int instanceIndex);

    // REQUIRES: mergedComp2InstanceMask() executed for both & vecMergedComp and & vecMergedCompShifted
    // EFFECTS:
    vector<CombinedComp> & mergedComp2InitialCombinedComp(vector<MergedComp> & vecMergedComp, vector<MergedComp> & vecMergedCompShifted);

    // REQUIRES: mergedComp2InitialCombinedComp() executed
    // EFFECTS: combine based on mask overlapping
    void combineComp(vector<CombinedComp> &vecCombinedComp, int iterationN, float minAreaToCombine);

    // REQUIRES: combineComp() executed
    // EFFECTS: display instance masks
    void combinedComp2InstanceMask(vector<CombinedComp> &vecCombinedComp);

};


#endif //TRYENV_IMAGEHANDLER_H
