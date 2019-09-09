//
// Created by ubuntu on 5/14/19.
//

#ifndef TRYENV_EVALUATION_H
#define TRYENV_EVALUATION_H

#include "Config.h"

class Evaluation {
public:
    vector<int> vecTP; // at every IoU level, total number of TP in all images
    vector<int> vecFP; // number of remaining prediction masks that does not have a match
    vector<int> vecFN; // number of remaining ground truth masks that does not have a match

    vector<float> vecIoULevel;
    vector<float> vecPrecision;
    vector<float> vecRecall;
    vector<float> vecPrecPerIoU, vecRecallPerIoU;
    std::ofstream output_file;

    // initialize vecIoULevel
    Evaluation();
    ~Evaluation(){
        output_file.close();
    }

    // Evaluate the result from all images
    // for each IoU level, run the program and compare with ground truth for each of the images

    // Called before run the algorithm for one IoU level
    void resetCount(){ vecPrecision.clear(); vecRecall.clear(); vecTP.clear(); vecFP.clear(); vecFN.clear(); }

    // Called after run the algorithm for each image
    // REQUIRES: vector of ground truth masks, vector of prediction masks, a level of IoU
    // EFFECTS: calculate precision and recall for one image
    void accumulateCount(const Mat & im, const int id, const vector<Mat> & vecGroundTruth, const vector<Mat> & vecPrediction, float IoULevel);

    // Called after run the algorithm for all images at one IoU level
    void calculateMeanPrecisionRecall(float IoULevel);

    void showMatchedMasksInImage(const Mat & im, const int id, const vector<Mat> & vecPrediction, const vector<int> GIndForP);

    void printTable();

};


#endif //TRYENV_EVALUATION_H
