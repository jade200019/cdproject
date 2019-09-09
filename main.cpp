//
// Created by jade on 2/26/19.
//

#include "Config.h"
#include "Patch.h"
#include "FileHandler.h"
#include "ImageHandler.h"
#include "GraphHandler.h"
#include "GUIApplication.h"
#include "MergedComp.h"
#include "Evaluation.h"

#include <ctime>

GUIApplication guiapp;                                                      // global definition
static void on_mouse( int event, int x, int y, int flags, void* param )
{
    guiapp.mouseClick( event, x, y, flags, param );
}

int main( int argc, char* argv[] ){
//    GraphHandler graphHandler;
//    //graphHandler.test_image();
//    graphHandler.test_random_image();
//    graphHandler.cutBranch();


    namedWindow( winName, WINDOW_AUTOSIZE );
    setMouseCallback( winName, on_mouse, 0 );
    namedWindow( "RRect", WINDOW_AUTOSIZE );
    setMouseCallback( "RRect", on_mouse, 0 );
    namedWindow( "Canny", WINDOW_AUTOSIZE );
    setMouseCallback( "Canny", on_mouse, 0 );

    //namedWindow( winNameMerge, WINDOW_AUTOSIZE );

    FileHandler fileHandler;
    //vector<Mat> images;
    //fileHandler.fileList2VecMat("input/*.jpg", images); // not used
    vector<vector<float>> allParameters;
    fileHandler.readParameters("input/parameters.txt", allParameters);

    Evaluation evaluation;
    std::ofstream output_runtime;
    output_runtime.open("result_runtime.txt");

#ifdef DO_EVALUATION_FOR_MULTIPLE_IOU
    for (float IoULevel : evaluation.vecIoULevel) {
#else
        float IoULevel = 0.3;
#endif
        evaluation.resetCount();
        int i = 0;
        for (vector<float> & parameters : allParameters) {


            Mat image = imread("input/" + std::to_string(int(parameters[0])) + ".jpg", cv::IMREAD_LOAD_GDAL | cv::IMREAD_COLOR );
            cout << "Read Image" << "input/" + std::to_string(int(parameters[0])) + ".jpg" << endl;
            vector<Mat> MatGroundTruth;
            FileHandler fileHandlerG;
            fileHandlerG.fileList2VecMat("input/" + std::to_string(int(parameters[0])) + "/*.png", MatGroundTruth);
            if (MatGroundTruth.empty()) {
                std::cerr << "main: MatGroundTruth.empty()\n";
            }
            //vector<float> parameters = allParameters[i];

            cout << "image size: " << image.cols << " " << image.rows << endl;

            //ImageHandler imageHandler(image);
            ImageHandler imageHandlerTemp(image, i); // to resize before create Window
            image = imageHandlerTemp.Resize();
            //imwrite("resized/" + std::to_string(i) + ".png", image);

            guiapp.reset();
            guiapp.setImage(image);
            Mat dst = image.clone();
            imshow(winName, dst);
#ifndef DISABLE_USER_INTERACTION
            char c = (char)waitKey(0);
#else
            char c = '\x1b';
#endif

            while (c != '\x1b') {        // ESC to exit

                if (c == 'c') {          // c to clear
                    guiapp.reset();
                }

                guiapp.showImage();

                if (c == 'r'){
                    ImageHandler imageHandler(image, i, parameters);
                    image = imageHandler.Resize();

                    Mat &canny = imageHandler.imCanny(CANNY_THRESHOLD_1, CANNY_THRESHOLD_2);

                    cout << "image size: " << image.cols << " " << image.rows << endl;

                    imageHandler.DisplayImageWithBox(std::to_string(i) + "_");

                    vector<vector<Patch>> &vecVecPatch = imageHandler.SplitToPatch();
                    vector<vector<Patch>> &vecVecShiftedPatch = imageHandler.SplitToShiftedPatch();

                    vector<RotatedRectm> &vecRRect = imageHandler.patch2rRect(vecVecPatch, K_MEANS_K, REF_BGR_MARGIN);
                    imageHandler.DisplayMaskInImage(vecVecPatch, std::to_string(i) + "_");

                    vector<RotatedRectm> &vecShiftedRRect = imageHandler.shiftedPatch2rRect(vecVecShiftedPatch, K_MEANS_K,
                                                                                            REF_BGR_MARGIN);

                    vector<MergedComp> &vecMergedComp = imageHandler.rRect2InitialMergedComp(vecRRect);
                    vector<MergedComp> &vecShiftedMergedComp = imageHandler.shiftedRRect2InitialMergedComp(vecShiftedRRect);

                    imageHandler.mergeComp(vecMergedComp, MERGE_NUM_ITERATION, MERGE_MAX_COST, false);
                    imageHandler.mergeComp(vecMergedComp, MERGE_NUM_ITERATION, MERGE_MAX_COST_PHASE2, false, true);
                    imageHandler.mergedComp2InstanceMask(vecMergedComp, vecRRect, false, INSTANCE_MIN_NUM_RRECT);

                    imageHandler.mergeComp(vecShiftedMergedComp, MERGE_NUM_ITERATION, MERGE_MAX_COST, true); // shifted
                    imageHandler.mergeComp(vecShiftedMergedComp, MERGE_NUM_ITERATION, MERGE_MAX_COST_PHASE2, true,
                                           true); // shifted
                    imageHandler.mergedComp2InstanceMask(vecShiftedMergedComp, vecShiftedRRect, true, INSTANCE_MIN_NUM_RRECT);

                    vector<CombinedComp> &vecCombinedComp = imageHandler.mergedComp2InitialCombinedComp(vecMergedComp,
                                                                                                        vecShiftedMergedComp);
                    imageHandler.combineComp(vecCombinedComp, 100, COMBINE_MIN_AREA);
                    imageHandler.combinedComp2InstanceMask(vecCombinedComp);

                    imageHandler.DisplayRRectInImage(std::to_string(i) + "_", "RRect");

                    if (MatGroundTruth.empty()) {
                        std::cerr << "main: MatGroundTruth.empty()\n";
                    }
                    for (Mat &mask : MatGroundTruth) {
                        ImageHandler imageHandlerTempMask(mask); // use the constructor that does not clone, but use reference
                        mask = imageHandlerTempMask.Resize().clone();
                    }

                    evaluation.accumulateCount(image, i, MatGroundTruth, imageHandler.vecPredictionMask, IoULevel);


                    imageHandler.saveToAllPatch();
                    imageHandler.DisplayMaskInImage(vecVecShiftedPatch, std::to_string(i) + "_", &Patch::getClosestColorMask,
                                                    true);

                    //imageHandler.DisplayBGRPatchInImage(vecPatch, std::to_string(i)+"_", &Patch::getClosestColorMask); // only one channel, should not use this display function
                    imageHandler.DisplayBGRPatchInImage(vecVecPatch, std::to_string(i) + "_", "KMeans", &Patch::getKMeansSeg);
                    imageHandler.DisplayBGRPatchInImage(vecVecPatch, std::to_string(i) + "_", "Contour",
                                                        &Patch::getkClosestColorMaskContour);


                    imageHandler.DisplayBGRPatchInImage(vecVecShiftedPatch, std::to_string(i) + "_", "KMeans",
                                                        &Patch::getKMeansSeg, true);
                    imageHandler.DisplayRRectInImage(std::to_string(i) + "_", "RRect", Scalar(0, 255, 0), false);
                    imageHandler.DisplayRRectInImage(std::to_string(i) + "_", "RRect", Scalar(255, 0, 0), true);



                }
#ifndef DISABLE_USER_INTERACTION
                c = (char)waitKey(0);
#else
                c = '\x1b';
#endif
            }


            std::clock_t start;
            start = std::clock();

            //ImageHandler imageHandler(image, i, guiapp.vecRefWidth, guiapp.vecRefBGR); // to process
            ImageHandler imageHandler(image, int(parameters[0]), parameters);
            image = imageHandler.Resize();

            Mat &canny = imageHandler.imCanny(CANNY_THRESHOLD_1, CANNY_THRESHOLD_2);

            cout << "image size: " << image.cols << " " << image.rows << endl;

#if SAVE_TO_FILE
            imageHandler.DisplayImageWithBox(std::to_string(int(parameters[0])) + "_");
#endif
            vector<vector<Patch>> &vecVecPatch = imageHandler.SplitToPatch();
            vector<vector<Patch>> &vecVecShiftedPatch = imageHandler.SplitToShiftedPatch();

            vector<RotatedRectm> &vecRRect = imageHandler.patch2rRect(vecVecPatch, K_MEANS_K, REF_BGR_MARGIN);

#if SAVE_TO_FILE
            imageHandler.DisplayMaskInImage(vecVecPatch, std::to_string(int(parameters[0])) + "_");
#endif

            vector<RotatedRectm> &vecShiftedRRect = imageHandler.shiftedPatch2rRect(vecVecShiftedPatch, K_MEANS_K,
                                                                                    REF_BGR_MARGIN);

            vector<MergedComp> &vecMergedComp = imageHandler.rRect2InitialMergedComp(vecRRect);
            vector<MergedComp> &vecShiftedMergedComp = imageHandler.shiftedRRect2InitialMergedComp(vecShiftedRRect);

            imageHandler.mergeComp(vecMergedComp, MERGE_NUM_ITERATION, MERGE_MAX_COST, false);
            imageHandler.mergeComp(vecMergedComp, MERGE_NUM_ITERATION, MERGE_MAX_COST_PHASE2, false, true);
            imageHandler.mergedComp2InstanceMask(vecMergedComp, vecRRect, false, INSTANCE_MIN_NUM_RRECT);

            imageHandler.mergeComp(vecShiftedMergedComp, MERGE_NUM_ITERATION, MERGE_MAX_COST, true); // shifted
            imageHandler.mergeComp(vecShiftedMergedComp, MERGE_NUM_ITERATION, MERGE_MAX_COST_PHASE2, true,
                                   true); // shifted
            imageHandler.mergedComp2InstanceMask(vecShiftedMergedComp, vecShiftedRRect, true, INSTANCE_MIN_NUM_RRECT);

            vector<CombinedComp> &vecCombinedComp = imageHandler.mergedComp2InitialCombinedComp(vecMergedComp,
                                                                                                vecShiftedMergedComp);
            imageHandler.combineComp(vecCombinedComp, 100, COMBINE_MIN_AREA);
            imageHandler.combinedComp2InstanceMask(vecCombinedComp);
            output_runtime << (float)(std::clock() - start) / CLOCKS_PER_SEC << endl;

            #if SAVE_TO_FILE
            imageHandler.DisplayRRectInImage(std::to_string(int(parameters[0])) + "_", "RRect");
            #endif

            if (MatGroundTruth.empty()) {
                std::cerr << "main: MatGroundTruth.empty()\n";
            }
            for (Mat &mask : MatGroundTruth) {
                ImageHandler imageHandlerTempMask(mask); // use the constructor that does not clone, but use reference
                mask = imageHandlerTempMask.Resize().clone();
                //imshow("mask", mask);
                //char c = waitKey(0);
            }

            evaluation.accumulateCount(image, int(parameters[0]), MatGroundTruth, imageHandler.vecPredictionMask, IoULevel);

            #if SAVE_TO_FILE
            imageHandler.saveToAllPatch();
            //imageHandler.rRect2Node();


            imageHandler.DisplayMaskInImage(vecVecShiftedPatch, std::to_string(int(parameters[0])) + "_", &Patch::getClosestColorMask,
                                            true);

            //imageHandler.DisplayBGRPatchInImage(vecPatch, std::to_string(i)+"_", &Patch::getClosestColorMask); // only one channel, should not use this display function
            imageHandler.DisplayBGRPatchInImage(vecVecPatch, std::to_string(int(parameters[0])) + "_", "KMeans", &Patch::getKMeansSeg);
            imageHandler.DisplayBGRPatchInImage(vecVecPatch, std::to_string(int(parameters[0])) + "_", "Contour",
                                                &Patch::getkClosestColorMaskContour);


            imageHandler.DisplayBGRPatchInImage(vecVecShiftedPatch, std::to_string(int(parameters[0])) + "_", "KMeans",
                                                &Patch::getKMeansSeg, true);
            //imageHandler.DisplayBGRPatchInImage(vecVecShiftedPatch, std::to_string(i) + "_", "Contour",
            //                                    &Patch::getkClosestColorMaskContour, true);


            //imageHandler.calcStartRRect(0.8, false);
            //imageHandler.calcStartRRect(0.8, true);
            imageHandler.DisplayRRectInImage(std::to_string(int(parameters[0])) + "_", "RRect", Scalar(0, 255, 0), false);
            imageHandler.DisplayRRectInImage(std::to_string(int(parameters[0])) + "_", "RRect", Scalar(255, 0, 0), true);
            //imageHandler.DisplayRRectInImage(std::to_string(i) + "_",
            //                                 "RRectOverlay"); // display both shifted and non-shifted on the same image

            //imageHandler.calcCableFromStartRRect();
            #endif

            i++;
        }
        evaluation.calculateMeanPrecisionRecall(IoULevel);
#ifdef DO_EVALUATION_FOR_MULTIPLE_IOU
    }
#endif
    evaluation.printTable();
    output_runtime.close();

#if 0
    int i = 0, j = 0;
    cv::Mat patchMat;
    cv::Mat patchLocation;
    cout << image.cols / PATCHX << " " << image.rows / PATCHY << endl;

    for (; i < (image.cols / PATCHX); i++) {
        j = 0;
        for (; j < image.rows / PATCHY; j++) {
            cv::Rect rangeRect(i * PATCHX, j * PATCHY, PATCHX, PATCHY);
            patchMat = image(rangeRect); // shallow copy, reference to the original image
            patchLocation = image.clone(); // deep copy, another instance
            cv::rectangle(patchLocation, rangeRect, cv::Scalar(0, 255, 0));

            if (patchMat.empty()) {
                std::cerr << "Cannot get the patch." << std::endl;
            }
            // Save the frame into a file
            imwrite("bin/Patch_" + std::to_string(i) + "_" + std::to_string(j) + ".jpg", patchMat); // A JPG FILE IS BEING SAVED
            imwrite("bin/Patch_" + std::to_string(i) + "_" + std::to_string(j) +  "_Location.jpg", patchLocation);

            // Create Patch instance
            Patch patch(patchMat);

#ifdef TEST_GRADIENT_METHOD
            // Calculate Sobel gradient
            Mat grad = patch.getSobelGradient();
            imwrite("bin/patch_" + std::to_string(i) + "_" + std::to_string(j) +  "_Grad.jpg", grad);

            patch.getGradPolarCoordinate();
            Mat gradThetaHist = patch.getGradThetaHist(Mat(), 5);
            imwrite("bin/patch_" + std::to_string(i) + "_" + std::to_string(j) +  "_GradThetaHist.jpg", gradThetaHist);

            Mat topK = patch.getGradThetaHistTopK(3);
            Mat topKMask = patch.ThetaHistBin2Mask(topK);
            imwrite("bin/patch_" + std::to_string(i) + "_" + std::to_string(j) +  "_topKMask.jpg", topKMask);
#endif

            // K-means
            Mat kMeansSeg = patch.kMeans(3);
            imwrite("bin/Patch_" + std::to_string(i) + "_" + std::to_string(j) +  "_kMeansSeg.jpg", kMeansSeg);
        }
        // cout << i << " " << image.cols / PATCHX << endl;
    }
#endif

    return(0);
}
