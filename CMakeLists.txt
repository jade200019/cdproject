cmake_minimum_required(VERSION 3.13)
project(tryEnv)

set(CMAKE_CXX_STANDARD 14)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#set(OpenCV_LIBS "opencv_core opencv_imgproc opencv_calib3d opencv_video opencv_features2d opencv_ml opencv_highgui opencv_objdetect opencv_contrib opencv_legacy opencv_gpu")
set(OpenCV_LIBS "/usr/local/include/opencv2")
find_package( OpenCV REQUIRED )

find_package(Boost 1.67.0 REQUIRED)

if(NOT Boost_FOUND)
    message(FATAL_ERROR "Could not find boost!")
endif()

include_directories(${Boost_INCLUDE_DIR})

add_executable(tryEnv main.cpp Patch.cpp Patch.h Config.h FileHandler.cpp FileHandler.h ImageHandler.cpp ImageHandler.h RotatedRectm.cpp RotatedRectm.h Cost.h Cost.cpp Cable.cpp Cable.h GraphHandler.cpp GraphHandler.h GUIApplication.cpp GUIApplication.h MergedComp.cpp MergedComp.h CombinedComp.cpp CombinedComp.h Evaluation.cpp Evaluation.h)
#add_executable(tryEnv samples/cpp/kmeans.cpp)

include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries(tryEnv ${OpenCV_LIBS})


file(REMOVE_RECURSE "./bin")
file(MAKE_DIRECTORY "./bin")