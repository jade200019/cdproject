//
// Created by ubuntu on 4/23/19.
//

#include "MergedComp.h"

void MergedComp::merge(MergedComp &comp2, bool toMergeTail1, bool toMergeTail2) {
    if (toMergeTail1 == false && toMergeTail2 == false){                // [T----COMP----H] [H----COMP2----T]
        while ( !comp2.listIdx.empty() ){
            int temp = comp2.listIdx.front();                           // comp2 pop from front
            comp2.listIdx.pop_front();
            this->listIdx.push_front(temp);                             // this push from front
        }

        this->headIdx = comp2.tailIdx;                                  // assign tail info of comp2 to head info of this
        this->pts[HL] = comp2.pts[TL];
        this->pts[HR] = comp2.pts[TR];
        this->pts[HL1] = comp2.pts[TL1];
        this->pts[HR1] = comp2.pts[TR1];

    }
    else if (toMergeTail1 == false && toMergeTail2 == true){            // [T----COMP----H] [T----COMP2----H]
        while ( !comp2.listIdx.empty() ){
            int temp = comp2.listIdx.back();                            // comp2 pop from back (tail)
            comp2.listIdx.pop_back();
            this->listIdx.push_front(temp);                             // this push from front
        }

        this->headIdx = comp2.headIdx;                                  // assign head info of comp2 to head info of this
        this->pts[HL] = comp2.pts[HL];
        this->pts[HR] = comp2.pts[HR];
        this->pts[HL1] = comp2.pts[HL1];
        this->pts[HR1] = comp2.pts[HR1];

    }
    else if (toMergeTail1 == true && toMergeTail2 == false){            // [H----COMP----T] [H----COMP2----T]
        while ( !comp2.listIdx.empty() ){
            int temp = comp2.listIdx.front();                           // comp2 pop from front
            comp2.listIdx.pop_front();
            this->listIdx.push_back(temp);                              // this push from back (tail)
        }

        this->tailIdx = comp2.tailIdx;                                  // assign tail info of comp2 to tail info of this
        this->pts[TL] = comp2.pts[TL];
        this->pts[TR] = comp2.pts[TR];
        this->pts[TL1] = comp2.pts[TL1];
        this->pts[TR1] = comp2.pts[TR1];

    }
    else {                                                              // [H----COMP----T] [T----COMP2----H]
        while ( !comp2.listIdx.empty() ){
            int temp = comp2.listIdx.back();                            // comp2 pop from back (tail)
            comp2.listIdx.pop_back();
            this->listIdx.push_back(temp);                              // this push from back (tail)
        }

        this->tailIdx = comp2.headIdx;                                  // assign head info of comp2 to tail info of this
        this->pts[TL] = comp2.pts[HL];
        this->pts[TR] = comp2.pts[HR];
        this->pts[TL1] = comp2.pts[HL1];
        this->pts[TR1] = comp2.pts[HR1];

    }
}
