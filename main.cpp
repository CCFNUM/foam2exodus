#include "ExodusWriter.h"
#include "MergedMeshReader.h"
#include "OpenFOAMFieldReader.h"
#include "OpenFOAMMeshReader.h"
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

void printUsage(const char* programName)
{
    std::cout << "Usage: " << programName
              << " [--fields all|U,p,...]"
                 " [--field-names velocity,pressure,...]"
                 " [--time latest|TIME]"
                 " <OpenFOAM_case_dir> <output.exo>\n";
    std::cout << "   or: " << programName
              << " --multiple <case_dir1> <case_dir2> ... <output.exo>\n";
    std::cout << "\nConverts OpenFOAM mesh(es) to Exodus II format with "
                 "boundary patches as sidesets.\n";
    std::cout << "Selected volume fields are conservatively projected to CVFEM "
                 "nodes.\n";
    std::cout << "\nArguments:\n";
    std::cout << "  OpenFOAM_case_dir : Path to OpenFOAM case directory "
                 "containing polyMesh\n";
    std::cout << "  output.exo        : Output Exodus II file path\n";
    std::cout << "  --fields all      : Convert every supported volume field\n";
    std::cout << "  --fields U,p,k    : Convert only the comma/space-separated "
                 "fields\n";
    std::cout << "  --field-names ... : Exodus base names in the same order as "
                 "--fields\n";
    std::cout
        << "  --time TIME       : Field time directory (default: latest)\n";
    std::cout << "\nMulti-mesh mode (--multiple):\n";
    std::cout
        << "  Merges multiple OpenFOAM meshes into a single Exodus file\n";
    std::cout << "  Each mesh's element blocks and sidesets remain separate\n";
    std::cout << "  Example: " << programName
              << " --multiple case1 case2 case3 merged.exo\n";
}

std::vector<std::string> parseFieldNames(std::string text)
{
    for (char& c : text)
        if (c == ',' || c == '(' || c == ')' || c == ';')
            c = ' ';
    std::istringstream input(text);
    std::vector<std::string> names;
    std::string name;
    while (input >> name)
        names.push_back(name);
    return names;
}

int main(int argc, char* argv[])
{
    if (argc == 1)
    {
        printUsage(argv[0]);
        return 1;
    }

    try
    {
        bool multipleMode = false;
        bool fieldsRequested = false;
        bool allFields = false;
        std::vector<std::string> fieldNames;
        std::vector<std::string> outputFieldNames;
        bool outputNamesRequested = false;
        std::string timeSelection = "latest";
        std::vector<std::string> positional;

        for (int i = 1; i < argc; ++i)
        {
            const std::string argument = argv[i];
            if (argument == "-h" || argument == "--help")
            {
                printUsage(argv[0]);
                return 0;
            }
            if (argument == "--multiple")
            {
                multipleMode = true;
                continue;
            }
            if (argument == "--time")
            {
                if (++i >= argc)
                    throw std::runtime_error("--time requires a value");
                timeSelection = argv[i];
                continue;
            }
            if (argument.rfind("--time=", 0) == 0)
            {
                timeSelection = argument.substr(7);
                continue;
            }
            std::string outputNamesArgument;
            if (argument == "--field-names")
            {
                if (++i >= argc)
                    throw std::runtime_error("--field-names requires a value");
                outputNamesArgument = argv[i];
            }
            else if (argument.rfind("--field-names=", 0) == 0)
            {
                outputNamesArgument = argument.substr(14);
            }
            if (!outputNamesArgument.empty())
            {
                outputNamesRequested = true;
                std::vector<std::string> parsed =
                    parseFieldNames(outputNamesArgument);
                outputFieldNames.insert(
                    outputFieldNames.end(), parsed.begin(), parsed.end());
                continue;
            }
            std::string fieldsArgument;
            if (argument == "--fields")
            {
                if (++i >= argc)
                    throw std::runtime_error("--fields requires a value");
                fieldsArgument = argv[i];
            }
            else if (argument.rfind("--fields=", 0) == 0)
            {
                fieldsArgument = argument.substr(9);
            }
            if (!fieldsArgument.empty())
            {
                fieldsRequested = true;
                std::string lower = fieldsArgument;
                std::transform(lower.begin(),
                               lower.end(),
                               lower.begin(),
                               [](unsigned char c)
                { return static_cast<char>(std::tolower(c)); });
                allFields = (lower == "all");
                if (!allFields)
                {
                    std::vector<std::string> parsed =
                        parseFieldNames(fieldsArgument);
                    fieldNames.insert(
                        fieldNames.end(), parsed.begin(), parsed.end());
                }
                continue;
            }
            if (!argument.empty() && argument[0] == '-')
                throw std::runtime_error("Unknown option: " + argument);
            positional.push_back(argument);
        }

        if (fieldsRequested && !allFields && fieldNames.empty())
            throw std::runtime_error("--fields did not contain a field name");
        if (outputNamesRequested && !fieldsRequested)
            throw std::runtime_error(
                "--field-names requires an explicit --fields list");
        if (outputNamesRequested && allFields)
            throw std::runtime_error(
                "--field-names cannot be combined with --fields all; list the "
                "source fields explicitly");
        if (outputNamesRequested &&
            outputFieldNames.size() != fieldNames.size())
            throw std::runtime_error(
                "--field-names must contain exactly one name for each entry "
                "in --fields");

        if (multipleMode)
        {
            if (fieldsRequested)
                throw std::runtime_error(
                    "--fields is currently supported only for a single case");
            // Multi-mesh mode: --multiple case1 case2 ... output.exo
            if (positional.size() < 2)
            {
                std::cerr << "Error: --multiple requires at least one case "
                             "directory and an output file\n";
                printUsage(argv[0]);
                return 1;
            }

            std::vector<std::string> casePaths;
            casePaths.assign(positional.begin(), positional.end() - 1);
            std::string outputPath = positional.back();

            std::cout << "=== Multi-mesh merge mode ===" << std::endl;
            std::cout << "Merging " << casePaths.size()
                      << " OpenFOAM meshes into: " << outputPath << std::endl;
            std::cout << std::endl;

            // Read and merge all meshes
            MergedMeshReader mergedReader(casePaths);
            mergedReader.readMeshes();

            std::cout << "\nWriting merged Exodus II file: " << outputPath
                      << std::endl;

            // Write merged mesh with custom naming to keep element blocks
            // separate
            ExodusWriter writer(outputPath);
            writer.writeMesh(mergedReader,
                             mergedReader.getElementBlockNames(),
                             mergedReader.getSidesetNames());

            std::cout << "\nMerge completed successfully!" << std::endl;
            std::cout << "Element blocks from different meshes remain separate."
                      << std::endl;
        }
        else
        {
            // Single mesh mode (backward compatible)
            if (positional.size() != 2)
            {
                std::cerr
                    << "Error: Single mesh mode requires exactly 2 arguments\n";
                printUsage(argv[0]);
                return 1;
            }

            std::string casePath = positional[0];
            std::string outputPath = positional[1];

            std::cout << "Reading OpenFOAM mesh from: " << casePath
                      << std::endl;

            OpenFOAMMeshReader reader(casePath);
            reader.readMesh();

            std::cout << "Mesh statistics:" << std::endl;
            std::cout << "  Points: " << reader.getNumPoints() << std::endl;
            std::cout << "  Cells: " << reader.getNumCells() << std::endl;
            std::cout << "  Boundary patches: "
                      << reader.getNumBoundaryPatches() << std::endl;

            std::cout << "\nWriting Exodus II file: " << outputPath
                      << std::endl;

            ExodusWriter writer(outputPath);
            writer.writeMesh(reader);

            if (fieldsRequested)
            {
                OpenFOAMFieldReader fieldReader(
                    casePath, timeSelection, reader);
                std::vector<VolFieldInfo> fields =
                    fieldReader.selectFields(fieldNames, allFields);
                for (size_t i = 0; i < outputFieldNames.size(); ++i)
                    fields[i].outputName = outputFieldNames[i];
                std::cout << "\nSelected OpenFOAM time "
                          << fieldReader.getTimeName() << " and "
                          << fields.size() << " field(s)" << std::endl;
                writer.writeNodalFields(
                    reader, fieldReader, fields, fieldReader.getTimeValue());
            }

            std::cout << "Conversion completed successfully!" << std::endl;
            std::cout << "Boundary patches written as sidesets." << std::endl;
        }

        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
