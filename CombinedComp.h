//
// Created by ubuntu on 5/1/19.
//

#ifndef TRYENV_COMBINEDCOMP_H
#define TRYENV_COMBINEDCOMP_H

#include "Config.h"

class CombinedComp {
public:
    list<int> listIdx;          // index of merged component in vecMergedComp
    list<int> listIdxShifted;   // index of shifted merged component in vecShiftedMergedComp
    Mat mask;

    CombinedComp(int ind, bool shifted, Mat & Mask) : mask(Mask.clone()){
        if (shifted){
            listIdxShifted.push_back(ind);
        }
        else {
            listIdx.push_back(ind);
        }
    }

    void combine(CombinedComp & comp2){
        while ( !comp2.listIdx.empty() ){
            int temp = comp2.listIdx.front();                           // comp2 pop from front
            comp2.listIdx.pop_front();
            this->listIdx.push_back(temp);                              // this push from back
        }
        while ( !comp2.listIdxShifted.empty() ){
            int temp = comp2.listIdxShifted.front();                           // comp2 pop from front
            comp2.listIdxShifted.pop_front();
            this->listIdxShifted.push_back(temp);                              // this push from back
        }
        bitwise_or(this->mask, comp2.mask, this->mask);
    }
};


#endif //TRYENV_COMBINEDCOMP_H
