//
// Created by jade on 3/11/19.
//

#include "FileHandler.h"

void FileHandler::fileList2VecMat(const char fileList[], vector<Mat> &images){
    //vector<cv::String> fn;
    glob(fileList, fn, false);

    size_t count = fn.size(); //number of png files in images folder
    for (size_t i=0; i<count; i++){
        images.push_back(imread(fn[i], cv::IMREAD_LOAD_GDAL | cv::IMREAD_COLOR ));
        cout << "Exist file: " << fn[i] << endl;
    }
}

void FileHandler::fileList2VecMat(const String fileList, vector<Mat> &images) {

    glob(fileList, fn, false);

    size_t count = fn.size(); //number of png files in images folder
    for (size_t i=0; i<count; i++){
        images.push_back(imread(fn[i], cv::IMREAD_LOAD_GDAL | cv::IMREAD_COLOR ));
        cout << "Read file: " << fn[i] << endl;
    }
}

void FileHandler::readParameters(const String fileList, vector<float> &parameters) {
    vector<cv::String> fn1;

    glob(fileList, fn, false);

    size_t count = fn.size(); //number of png files in images folder
    for (size_t i=0; i<count; i++){
        std::ifstream input_file;
        input_file.open(fn[i]);
        if (!input_file){
            std::cerr << "Parameter file " << fn[i] << " is not open.\n";
            continue;
        }
        int num;
        input_file >> num;
        float para;
        for (int j = 0; j < num; ++j){
            input_file >> para;
            parameters.push_back(para);
        }
        cout << "Read parameter: " << fn[i] << endl;
    }
}

void FileHandler::readParameters(const String fileStr, vector<vector<float>> &allParameters) {
    std::ifstream input_file;
    input_file.open(fileStr);

    if (!input_file){
        std::cerr << "Parameter file " << fileStr << " is not open.\n";
        return;
    }

    int numIm, numPara;
    float dummyIndex;
    input_file >> numIm >> numPara;
    for (int i = 0; i < numIm; ++i){
        vector<float> vec;
        input_file >> dummyIndex;
        vec.push_back(dummyIndex);
        for (int j = 0; j < numPara; ++j){
            float para;
            input_file >> para;
            vec.push_back(para);
        }
        allParameters.push_back(vec);
    }

}
