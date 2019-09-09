//
// Created by jade on 2/27/19.
//

#ifndef TRYENV_CONFIG_H
#define TRYENV_CONFIG_H

#include <iostream>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include <vector>
#include <list>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/connected_components.hpp>
#include <fstream>
#include <sstream> // for debug

#include <math.h>
#include <numeric> // for accumulate
#include <fstream> // print to file

using std::cout;
using std::endl;
using namespace cv;

using std::vector;
using std::list;

// Math definition
#define PI 3.14159265

// Work resolution
#define MAX_LONG_SIDE 640.0

// Patch size
#define PATCHX 32 // from left to right
#define PATCHY 32 // from top to bottomS

// Displacement to calculate Cost
// enum Displacement {HalfPatch, OnePatch, OneAndHalfPatch, TwoPatch};
#define HalfPatch 1
#define OnePatch 2
#define OneAndHalfPatch 3
#define TwoPatch 4

// Angle precision
#define EPSILON_ANGLE 0.5

// Shifted macro definition, use in the index vector
#define NONSHIFTED 0
#define SHIFTED 1

// Max number of RotatedRects in one patch
#define MAX_RRECT_NUM_IN_PATCH 5

// Graph
// number of nodes = 2 x 20 x 20 x 3 = 2400
#define MAX_NUM_NODE 2400
//#define GRAPH_DISPLAY_PATCH_SIZE 128
#define GRAPH_DISPLAY_PATCH_SIZE 32
// #define DEBUG_PRINT_COSTVAL                              // in ImageHandler::rRect2Node()
// #define DEBUG_PRINT_VERTICES_AND_OUT_DEGREE              // in GraphHandler::cutBranch
// #define DEBUG_PRINT_OUTEDGE_OF_THE_VERTEX
// #define DEBUG_PRINT_PUSH_EDGE_TO_REMOVE
// #define DEBUG_PRINT_REMOVED_EDGES
// #define DEBUG_PRINT_CONNECTED_COMPONENTS

// Debug flag
//#define DEBUG_SAVE_PATCH_TO_FILE

// const maximum distance in BGR
#define BGR_MAX_DISTANCE_CONST 445

// Triangle vertex precision
#define EPSILON_DISTANCE    1.0
#define EPSILON_COS         0.1

// GUI
const std::string winName = "image";
const std::string winNameGraph = "graph";
const std::string winNameMerge = "merge_non-shifted";
const std::string winNameMergeShifted = "merge_shifted";
const std::string winNameInstance = "instance_non-shifted";
const std::string winNameInstanceShifted = "instance_shifted";
const std::string winNameCombineString = "combine";

// MergeComp constants
#define HL 0
#define HR 1
#define TL 2
#define TR 3
#define HL1 4
#define HR1 5
#define TL1 6
#define TR1 7

// --------------- Options -------------------
//#define PAUSE_AND_SHOW_COMBINED
#define DISABLE_USER_INTERACTION
//#define DO_EVALUATION_FOR_MULTIPLE_IOU
//#define PAUSE_AND_SHOW_MATCHED_PREDICTION
#define SAVE_TO_FILE 1
#define PRINT_MESSAGES_TO_TERMINAL 0

// parameters
// default width and color, change in input, not used later
#define DEFAULT_CABLE_WIDTH_1         13.0
#define DEFAULT_CABLE_WIDTH_2         15.0
#define DEFAULT_CABLE_COLOR_1       Vec3f({5, 5, 5})
#define DEFAULT_CABLE_COLOR_2       Vec3f({20, 20, 20})

// canny parameters
#define ENABLE_CANNY                1
#define CANNY_ROTATED_RECTANGLE_ROI_WIDTH 4
#define CANNY_THRESHOLD_1           150
#define CANNY_THRESHOLD_2           200

// patch to initial masks
#define K_MEANS_K                   3
#define REF_BGR_MARGIN              50

#define MASK_MORPH_SIZE             1

// initial masks to rectangle parameters
#define MIN_RECTANGLE_WIDTH         DEFAULT_CABLE_WIDTH_1 / 2
#define MIN_RECTANGLE_LENGTH        10.0
#define MIN_AREA_PERCENTAGE         0.5
#define MAX_WIDTH_DEVIATION         3.0
#define MIN_CANNY_LEFT_CONFIDENCE   0.2
#define MIN_CANNY_RIGHT_CONFIDENCE  0.2

#define MIN_TRIANGLE_AREA_PERCENTAGE    0.90

// MergedComp cost parameters
#define COSCOEFFICIENT              50
#define COSCOEFFICIENT_1            20
#define COS_COEFFICIENT_BEZIER      100
#define BEZIER_ESCAPE_SMALL_LENGTH  5
#define ANGLE_ESCAPE_SMALL_LENGTH   15
#define COS_NEAR_PARALLEL           0.99
#define COST_VALUE_SCALE_LENGTH_ALPHA   20.0
#define COST_VALUE_SCALE_ANGLE_ALPHA    30.0

// merging process parameters
#define MERGE_NUM_ITERATION         100
#define MERGE_MAX_COST              150.0
#define MERGE_MAX_COST_PHASE2       150.0

// instance mask
#define INSTANCE_MIN_NUM_RRECT      3

// combine
#define COMBINE_MIN_AREA            MIN_RECTANGLE_WIDTH * MIN_RECTANGLE_LENGTH * 3







#endif //TRYENV_CONFIG_H
