//
// Created by ubuntu on 4/15/19.
//

#ifndef TRYENV_GUIAPPLICATION_H
#define TRYENV_GUIAPPLICATION_H

#include "Config.h"

class GUIApplication {
public:
    int clickCount = 0;
    vector<Point2i> vecClick; // at most [2]
    vector<float> vecRefWidth;
    vector<Vec3f> vecRefBGR;
    Mat im;
    Mat dst;

    // REQUIRES: image
    // EFFECTS: copy the reference here
    // must be called when ever a new image is to be used
    void setImage(const Mat & Im){
        im = Im;
        dst = Im.clone();
    }

    void mouseClick( int event, int x, int y, int flags, void* param );

    void resetClick(){
        clickCount = 0;
        vecClick.clear();
    }

    void reset(){
        clickCount = 0;
        vecClick.clear();
        vecRefWidth.clear();
        vecRefBGR.clear();
        dst = im.clone();
        cout << "reset\n";
    }

    void showImage() {
        for (Point2i & pt : vecClick){
            drawMarker(dst, pt, Scalar(0, 0, 255, 0.5), MARKER_CROSS, 5, 2);
            //circle(dst, pt, 2, Scalar(0, 0, 255), FILLED);
        }
        imshow(winName, dst);
    }
};


#endif //TRYENV_GUIAPPLICATION_H
