//
// Created by ubuntu on 4/23/19.
//

#ifndef TRYENV_MERGEDCOMP_H
#define TRYENV_MERGEDCOMP_H

#include "Config.h"
#include "RotatedRectm.h"

class MergedComp {
public:
    list<int> listIdx;
    int headIdx;
    int tailIdx;
    Point2f pts[8]; // HL, HR, TL, TR; HL1, HR1, TL1, TR1
    Mat mask;       // check if empty
    list<RotatedRectm> listRRect;   // TODO include the mask in the process of merging; TODO use RRect to combine

    MergedComp(int i, Point2f * pts8) : headIdx(i), tailIdx(i) {
        listIdx.push_back(i);
        if (pts8 == NULL){
            std::cerr << "MergedComp::MergedComp() failed. pts8 is null pointer.\n";
        }
        else {
            for (int j = 0; j < 4; ++j){
                pts[j] = pts8[j];
            }
            pts[4] = pts8[3]; // HL1 = pts8[3] when there's only one rectangle in the comp
            pts[5] = pts8[2]; // HR1
            pts[6] = pts8[1]; // TL1
            pts[7] = pts8[0]; // TR1
        }
    }

    void merge(MergedComp & comp2, bool toMergeTail1, bool toMergeTail2);
};


#endif //TRYENV_MERGEDCOMP_H
