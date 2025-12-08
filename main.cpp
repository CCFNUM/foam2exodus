#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include "OpenFOAMMeshReader.h"
#include "ExodusWriter.h"
#include "MergedMeshReader.h"

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " <OpenFOAM_case_dir> <output.exo>\n";
    std::cout << "   or: " << programName << " --multiple <case_dir1> <case_dir2> ... <output.exo>\n";
    std::cout << "\nConverts OpenFOAM mesh(es) to Exodus II format with boundary patches as sidesets\n";
    std::cout << "\nArguments:\n";
    std::cout << "  OpenFOAM_case_dir : Path to OpenFOAM case directory containing polyMesh\n";
    std::cout << "  output.exo        : Output Exodus II file path\n";
    std::cout << "\nMulti-mesh mode (--multiple):\n";
    std::cout << "  Merges multiple OpenFOAM meshes into a single Exodus file\n";
    std::cout << "  Each mesh's element blocks and sidesets remain separate\n";
    std::cout << "  Example: " << programName << " --multiple case1 case2 case3 merged.exo\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        printUsage(argv[0]);
        return 1;
    }

    try {
        // Check if --multiple flag is present
        std::string firstArg = argv[1];
        bool multipleMode = (firstArg == "--multiple");

        if (multipleMode) {
            // Multi-mesh mode: --multiple case1 case2 ... output.exo
            if (argc < 4) {
                std::cerr << "Error: --multiple requires at least one case directory and an output file\n";
                printUsage(argv[0]);
                return 1;
            }

            // All arguments except first (--multiple) and last (output.exo) are case directories
            std::vector<std::string> casePaths;
            for (int i = 2; i < argc - 1; ++i) {
                casePaths.push_back(argv[i]);
            }
            std::string outputPath = argv[argc - 1];

            std::cout << "=== Multi-mesh merge mode ===" << std::endl;
            std::cout << "Merging " << casePaths.size() << " OpenFOAM meshes into: " << outputPath << std::endl;
            std::cout << std::endl;

            // Read and merge all meshes
            MergedMeshReader mergedReader(casePaths);
            mergedReader.readMeshes();

            std::cout << "\nWriting merged Exodus II file: " << outputPath << std::endl;

            // Write merged mesh with custom naming to keep element blocks separate
            ExodusWriter writer(outputPath);
            writer.writeMesh(mergedReader,
                           mergedReader.getElementBlockNames(),
                           mergedReader.getSidesetNames());

            std::cout << "\nMerge completed successfully!" << std::endl;
            std::cout << "Element blocks from different meshes remain separate." << std::endl;

        } else {
            // Single mesh mode (backward compatible)
            if (argc != 3) {
                std::cerr << "Error: Single mesh mode requires exactly 2 arguments\n";
                printUsage(argv[0]);
                return 1;
            }

            std::string casePath = argv[1];
            std::string outputPath = argv[2];

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
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
