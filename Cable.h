//
// Created by ubuntu on 3/29/19.
//

#ifndef TRYENV_CABLE_H
#define TRYENV_CABLE_H

#include "Config.h"

class Cable {
    // Vector to store index of RRects. [NONSHIFTED/SHIFTED, x of the patch, y of the patch, k-th rectangle in the patch]
    vector<int> start;
    vector<vector<int>> trackDown; // in the direction of angle
    vector<vector<int>> trackUp; // against the direction of angle
    vector<vector<int>> trackAll;
    int length;
public:
    // REQUIRES: the index of the start point
    // EFFECTS: initialize the object
    Cable(vector<int> Start) {
        if (Start.size() == 4){
            this->start = std::move(Start);
            this->trackAll.push_back(this->start);
            this->length = 1;
        }
        else {
            this->length = 0;
        }
    }

    // REQUIRES: the index vector to add; false: down, true: up
    // EFFECTS: return the length after added
    int addTrack(vector<int> index, bool Up);

    // REQUIRES: this->trackAll
    // EFFECTS: return it
    vector<vector<int>> getTrackAll() const {
        return this->trackAll;
    }
};


#endif //TRYENV_CABLE_H
