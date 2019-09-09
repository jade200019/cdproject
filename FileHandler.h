//
// Created by jade on 3/11/19.
//

#ifndef TRYENV_FILEHANDLER_H
#define TRYENV_FILEHANDLER_H

#include "Config.h"

class FileHandler {
public:
    vector<cv::String> fn;

    FileHandler() = default;
    ~FileHandler() = default;

    void fileList2VecMat(const char fileList[], vector<Mat> &images);

    // to use '+' add of strings in the caller
    void fileList2VecMat(const String fileList, vector<Mat> &images);

    void readParameters(const String fileList, vector<float> &parameters);

    void readParameters(const String file, vector<vector<float>> &allParameters);
};


#endif //TRYENV_FILEHANDLER_H
