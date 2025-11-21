#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include "OpenFOAMMeshReader.h"
#include "ExodusWriter.h"

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " <OpenFOAM_case_dir> <output.exo>\n";
    std::cout << "\nConverts OpenFOAM mesh to Exodus II format with boundary patches as sidesets\n";
    std::cout << "\nArguments:\n";
    std::cout << "  OpenFOAM_case_dir : Path to OpenFOAM case directory containing polyMesh\n";
    std::cout << "  output.exo        : Output Exodus II file path\n";
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        printUsage(argv[0]);
        return 1;
    }

    std::string casePath = argv[1];
    std::string outputPath = argv[2];

    try {
        std::cout << "Reading OpenFOAM mesh from: " << casePath << std::endl;
        
        OpenFOAMMeshReader reader(casePath);
        reader.readMesh();
        
        std::cout << "Mesh statistics:" << std::endl;
        std::cout << "  Points: " << reader.getNumPoints() << std::endl;
        std::cout << "  Cells: " << reader.getNumCells() << std::endl;
        std::cout << "  Boundary patches: " << reader.getNumBoundaryPatches() << std::endl;
        
        std::cout << "\nWriting Exodus II file: " << outputPath << std::endl;
        
        ExodusWriter writer(outputPath);
        writer.writeMesh(reader);
        
        std::cout << "Conversion completed successfully!" << std::endl;
        std::cout << "Boundary patches written as sidesets." << std::endl;
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
