//
// Created by ubuntu on 5/14/19.
//

#include "Evaluation.h"

Evaluation::Evaluation() {
//    for (float i = 1; i > 0; i -= 0.05){
//        this->vecIoULevel.push_back(i);
//    }
    vecIoULevel.push_back(0.3);
    vecIoULevel.push_back(0.5);
    vecIoULevel.push_back(0.7);
    resetCount();
    output_file.open("result_evaluation.txt");
}

void Evaluation::accumulateCount(const Mat & im, const int id, const vector<Mat> &vecGroundTruth, const vector<Mat> &vecPrediction, float IoULevel) {
    /// initialize a vector of int with size of vecPrediction and content of -1
    vector<int> PIndForG(vecGroundTruth.size(), -1);
    vector<int> GIndForP(vecPrediction.size(), -1); // to display matched prediction

    vector< std::pair<float, std::pair<int,int>>> data;

    for (int i = 0; i < vecGroundTruth.size(); ++i){
        for (int j = 0; j < vecPrediction.size(); ++j){
            Mat BinaryG;
            cvtColor(vecGroundTruth[i], BinaryG, COLOR_BGR2GRAY);
            threshold(BinaryG, BinaryG, 1, 255, THRESH_BINARY);

            Mat I;
            bitwise_and(BinaryG, vecPrediction[j], I);
            Mat U;
            bitwise_or(BinaryG, vecPrediction[j], U);
            int iCount = countNonZero(I);
            int uCount = countNonZero(U);
            float IoU = float(iCount) / float(uCount);
            if (IoU > IoULevel){
                cout << "IoU: " << IoU << " i:" << i << " " << iCount << " j:" << j << " " << uCount << endl;
                data.push_back(std::make_pair(IoU, std::make_pair(i, j)));
            }

        }
    }

    /// sort according to the first entry
    std::sort(data.begin(), data.end());
    for (auto x = data.rbegin(); x != data.rend(); ++x){
        int i = x->second.first;
        int j = x->second.second;
        if (PIndForG[i] == -1 && GIndForP[j] == -1){
            PIndForG[i] = j;
            GIndForP[j] = i;
        }
    }

    int notMatchedP = int(std::count(GIndForP.begin(), GIndForP.end(), -1)); // FP
    int notMatchedG = int(std::count(PIndForG.begin(), PIndForG.end(), -1)); // FN
    int numTP = int(GIndForP.size()) - notMatchedP; // TP

    float precision = (GIndForP.size() == 0) ? 0 : float(numTP) / float(GIndForP.size());
    float recall = (PIndForG.size() == 0) ? 0 : float(numTP) / float(PIndForG.size());
    vecPrecision.push_back(precision);
    vecRecall.push_back(recall);
    vecTP.push_back(numTP);
    vecFP.push_back(notMatchedP);
    vecFN.push_back(notMatchedG);

    cout << "Image: " << id << "\tTP:" << numTP <<"\tFP:" << notMatchedP << "\tFN:" << notMatchedG
            << "\tTP+FN:" << numTP + notMatchedG << endl
            << "\t\tIoULevel:" << IoULevel << " precision:" << precision << " recall:" << recall << endl;
    output_file << "Image: " << id << "\tIoU:" << IoULevel << "\tTP:" << numTP <<"\tFP:" << notMatchedP << "\tFN:" << notMatchedG << "\tTP+FN:" << numTP + notMatchedG
                << "\tprecision:" << precision << "\trecall:" << recall << endl;

    this->showMatchedMasksInImage(im, id, vecPrediction, GIndForP);

}

static Vec3f PlotBGRMatched[] = {{0, 0, 0xff},      // red
                                  {0, 0xff, 0xff},  // yellow
                                  {0, 0xff, 0},     // green
                                  {0xff, 0xff, 0},  // light blue
                                  {0xff, 0, 0},     // blue
                                  {0xff, 0, 0xff},  // pink #ff00ff
                                  {0, 0x80, 0xff},  // orange #ff8000
                                  {0, 0xff, 0x80},  // light green #80ff00
                                  {0x80, 0xff, 0},  // blueish green #00ff80
                                  {0xff, 0x80, 0},  // faint blue #0080ff
                                  {0xff, 0, 0x80},  // purple #8000ff
                                  {0x80, 0, 0xff},  // dark pink #ff0080
                                  };                // size: 12

static Vec3f PlotBGRNotMatched[] = {{0xc0, 0x80, 0x80},       // dark blue
                                    {0xc0, 0x80, 0xc0},
                                    {0x80, 0x80, 0xc0},
//                                    {0xcc, 0xcc, 0xcc}, // light gray
//                                    {0x80, 0x80, 0x80},
//                                    {0x40, 0x40, 0x40}, // dark gray
                                    };


void Evaluation::showMatchedMasksInImage(const Mat &im, const int id, const vector<Mat> &vecPrediction, const vector<int> GIndForP) {
    Mat dst = Mat(im.size(), im.type(), Scalar(0));
    Mat dst1 = Mat(im.size(), im.type(), Scalar(0));
    int NotMatchedInd = 0;
    for (int j = 0; j < vecPrediction.size(); ++j) {
        const Mat &mask = vecPrediction[j];

        vector<vector<Point> > contours;
        vector<Vec4i> hierarchy;
        findContours(mask, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE);

        Scalar color;
        if (GIndForP[j] != -1) {
            color = PlotBGRMatched[GIndForP[j] % 12];
        }
        else {
            color = PlotBGRNotMatched[(NotMatchedInd++) % 3];
        }
        for( size_t i = 0; i< contours.size(); i++ )
        {
            drawContours( dst, contours, (int)i, color, FILLED, LINE_8, hierarchy, 0 );
            putText(dst, std::to_string(j + 1), contours[i][0], FONT_HERSHEY_COMPLEX, 0.5, Scalar(255, 255, 255));

            drawContours( dst1, contours, (int)i, color, FILLED, LINE_8, hierarchy, 0 );
        }
    }

    imwrite("bin/" + std::to_string(id) + "_matched_prediction.png", dst);  // black background

    addWeighted(im, 1, dst1, 0.4, 0, dst1);


    for (int j = 0; j < vecPrediction.size(); ++j) {
        const Mat &mask = vecPrediction[j];

        vector<vector<Point> > contours;
        vector<Vec4i> hierarchy;
        findContours(mask, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE);

        Scalar color;
        if (GIndForP[j] != -1) {
            color = PlotBGRMatched[GIndForP[j] % 12];
        }
        else {
            color = PlotBGRNotMatched[(NotMatchedInd++) % 3];
        }
        for( size_t i = 0; i< contours.size(); i++ )
        {
            if (GIndForP[j] != -1)
                drawContours( dst1, contours, (int)i, color, 2, LINE_8, hierarchy, 0 );
            else
                drawContours( dst1, contours, (int)i, color, 1, LINE_8, hierarchy, 0 );
            putText(dst1, std::to_string(j + 1), contours[i][0], FONT_HERSHEY_SIMPLEX, 0.5, Scalar(255, 255, 255));
        }
    }

    imwrite("bin/" + std::to_string(id) + "_matched_prediction_transparency.png", dst1);  // overlay on image
    imwrite("result/" + std::to_string(id) + "_matched.png", dst1);  // overlay on image
#ifdef PAUSE_AND_SHOW_MATCHED_PREDICTION
    imshow("matched_prediction", dst1);
    waitKey(0);
#endif
}

void Evaluation::calculateMeanPrecisionRecall(float IoULevel) {
    double meanPrecision = std::accumulate(vecPrecision.begin(), vecPrecision.end(), 0.0) / vecPrecision.size();
    double meanRecall = std::accumulate(vecRecall.begin(), vecRecall.end(), 0.0) / vecRecall.size();
    cout << "IoULevel:" << IoULevel<< " meanPrecision: " << meanPrecision << " meanRecall:" << meanRecall << endl;
    output_file << "Summary of IoULevel: " << IoULevel << " \t\t" << meanPrecision << "\t\t" << meanRecall
                << " \t\t" << accumulate(vecTP.begin(), vecTP.end(),0)
                << " \t\t" << accumulate(vecFP.begin(), vecFP.end(),0)
                << " \t\t" << accumulate(vecFN.begin(), vecFN.end(),0) << endl;
    output_file << "total ground truth: " << accumulate(vecTP.begin(), vecTP.end(),0)  + accumulate(vecFN.begin(), vecFN.end(),0) << endl;
    output_file << endl << endl;

    vecPrecPerIoU.push_back(meanPrecision);
    vecRecallPerIoU.push_back(meanRecall);
}

void Evaluation::printTable() {
    output_file << "Table line:" << endl;
    for (int i = 0; i < vecPrecPerIoU.size(); ++i) {
        output_file << vecPrecPerIoU[i] << "\t" << vecRecallPerIoU[i] << "\t";
    }
    output_file << endl << endl;
    cout << "\nTable line:" << endl;
    for (int i = 0; i < vecPrecPerIoU.size(); ++i) {
        cout << vecPrecPerIoU[i] << "\t" << vecRecallPerIoU[i] << "\t";
    }
    cout << endl << endl;
}
