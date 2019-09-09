//
// Created by ubuntu on 3/29/19.
//

#include "Cable.h"

int Cable::addTrack(vector<int> index, bool Up) {
    if (!Up){
        this->trackDown.push_back(index);
    }
    else {
        this->trackUp.push_back(index);
    }
    this->trackAll.push_back(index);
    this->length++;
    return this->length;
}
