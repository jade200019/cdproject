//
// Created by ubuntu on 4/15/19.
//

#include "GUIApplication.h"

void GUIApplication::mouseClick(int event, int x, int y, int flags, void *param) {
    switch (event){
        case EVENT_MOUSEMOVE: {

        }
            break;
        case EVENT_LBUTTONDOWN: {
            cout << "x:" << x <<" y:" << y << endl;

            if (clickCount % 2 == 0) {                      // odd click
                vecClick.emplace_back(Point2i(x, y));


                cout << "vecClick = ";
                for (Point2i & pt : vecClick){
                    cout << pt << " ";
                }
                cout << endl;

                clickCount++;
            }
            else if (clickCount % 2 == 1){                  // even click
                vecClick.emplace_back(Point2i(x, y));
                if (norm(vecClick[vecClick.size() - 1]
                             - vecClick[vecClick.size() - 2]) > 2)    // double click, only add color, does not add width
                    vecRefWidth.emplace_back(norm(vecClick[vecClick.size() - 1] - vecClick[vecClick.size() - 2]));
                Point2i midPt = (vecClick[vecClick.size() - 1] + vecClick[vecClick.size() - 2]) / 2;
                if (im.cols < midPt.x || im.rows < midPt.y){
                    cout << "image is too small to select mid point:" << midPt << endl;
                }
                else {
                    ////https://answers.opencv.org/question/1870/find-pixel-color-out-of-cvmat-on-specific-position/
                    vecRefBGR.emplace_back(im.at<Vec3b>(midPt));
                }

                cout << "vecClick = ";
                for (Point2i &pt : vecClick) {
                    cout << pt << " ";
                }
                cout << endl;
                cout << "vecRefWidth = ";
                for (float i : vecRefWidth) {
                    cout << i << " ";
                }
                cout << endl;
                cout << "vecRefBGR = ";
                for (Vec3f &color : vecRefBGR) {
                    cout << color << " ";
                }
                cout << endl;


                clickCount++;
            }

        }
            break;
    }
    this->showImage();
}
