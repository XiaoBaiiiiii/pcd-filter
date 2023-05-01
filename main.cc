#include <string>
#include <iostream>
#include "PointCloudFilter.h"

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " input.ply output.ply\n";
        return 0;
    }
    const std::string input_filename = argv[1];
    const std::string output_filename = argv[2];
    std::cout<<"start processing"<<std::endl;
    if (pointcloudfilter::pointCloudFilter(input_filename, output_filename)) {
        std::cout<< "Point cloud filtered successfully."<<std::endl;
    }

    return 0;
}