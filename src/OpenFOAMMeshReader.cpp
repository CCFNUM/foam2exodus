// File       : OpenFOAMMeshReader.cpp
// Created    : Thu Mar 19 2026
// Author     : Mhamad Mahdi Alloush
// Description: 
// Copyright (c) 2026 CCFNUM, Lucerne University of Applied Sciences and
// Arts.
// SPDX-License-Identifier: BSD-3-Clause

#include "OpenFOAMMeshReader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <cctype>

OpenFOAMMeshReader::OpenFOAMMeshReader(const std::string& casePath) 
    : meshPath(casePath + "/constant/polyMesh") {}

void OpenFOAMMeshReader::readMesh() {
    readPoints();
    readFaces();
    readOwner();
    readNeighbour();
    readBoundary();
    readCellZones();
    constructCells();
}

std::string OpenFOAMMeshReader::stripComments(const std::string& line) {
    size_t pos = line.find("//");
    if (pos != std::string::npos) {
        return line.substr(0, pos);
    }
    return line;
}

bool OpenFOAMMeshReader::skipHeaderAndDetectFormat(std::ifstream& file) {
    std::string line;
    bool foundFoamFile = false;
    bool isBinary = false;
    int braceCount = 0;
    
    while (std::getline(file, line)) {
        line = stripComments(line);
        
        if (line.find("FoamFile") != std::string::npos) {
            foundFoamFile = true;
        }
        
        if (foundFoamFile) {
            if (line.find("format") != std::string::npos) {
                if (line.find("binary") != std::string::npos) {
                    isBinary = true;
                }
            }
            
            for (char c : line) {
                if (c == '{') braceCount++;
                if (c == '}') braceCount--;
            }
            
            if (braceCount == 0 && line.find('}') != std::string::npos) {
                return isBinary;
            }
        }
    }
    return isBinary;
}

template<typename T>
T OpenFOAMMeshReader::readBinary(std::ifstream& file) {
    T value;
    file.read(reinterpret_cast<char*>(&value), sizeof(T));
    return value;
}

std::string getNextNonEmptyLineIncludingDelimiters(std::ifstream& file) {
    std::string line;
    while (std::getline(file, line)) {
        size_t firstNonSpace = line.find_first_not_of(" \t\r\n");
        if (firstNonSpace != std::string::npos) {
            line.erase(0, firstNonSpace);
        } else {
            line.clear();
        }
        
        size_t lastNonSpace = line.find_last_not_of(" \t\r\n");
        if (lastNonSpace != std::string::npos) {
            line.erase(lastNonSpace + 1);
        } else {
            line.clear();
        }
        
        size_t commentPos = line.find("//");
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
            lastNonSpace = line.find_last_not_of(" \t\r\n");
            if (lastNonSpace != std::string::npos) {
                line.erase(lastNonSpace + 1);
            } else {
                line.clear();
            }
        }
        
        if (!line.empty()) {
            return line;
        }
    }
    return "";
}

std::string getNextNonEmptyLine(std::ifstream& file) {
    std::string line;
    while (std::getline(file, line)) {
        size_t firstNonSpace = line.find_first_not_of(" \t\r\n");
        if (firstNonSpace != std::string::npos) {
            line.erase(0, firstNonSpace);
        } else {
            line.clear();
        }
        
        size_t lastNonSpace = line.find_last_not_of(" \t\r\n");
        if (lastNonSpace != std::string::npos) {
            line.erase(lastNonSpace + 1);
        } else {
            line.clear();
        }
        
        size_t commentPos = line.find("//");
        if (commentPos != std::string::npos) {
            line = line.substr(0, commentPos);
            lastNonSpace = line.find_last_not_of(" \t\r\n");
            if (lastNonSpace != std::string::npos) {
                line.erase(lastNonSpace + 1);
            } else {
                line.clear();
            }
        }
        
        if (!line.empty() && line != "(" && line != ")" && line != "{" && line != "}" && line != ";") {
            return line;
        }
    }
    return "";
}

int OpenFOAMMeshReader::parseLeadingInt(const std::string& line, size_t* endPos) {
    size_t idx = 0;
    while (idx < line.size() && std::isspace(static_cast<unsigned char>(line[idx]))) {
        idx++;
    }
    
    size_t start = idx;
    if (idx < line.size() && (line[idx] == '+' || line[idx] == '-')) {
        idx++;
    }
    
    size_t digitsStart = idx;
    while (idx < line.size() && std::isdigit(static_cast<unsigned char>(line[idx]))) {
        idx++;
    }
    
    if (digitsStart == idx) {
        throw std::runtime_error("Error parsing integer from line: '" + line + "'");
    }
    
    if (endPos) {
        *endPos = idx;
    }
    
    return std::stoi(line.substr(start, idx - start));
}

std::vector<int> OpenFOAMMeshReader::parseIntListFromString(const std::string& line) {
    std::string cleaned = line;
    for (char& c : cleaned) {
        if (c == '(' || c == ')' || c == ';') {
            c = ' ';
        }
    }
    
    std::istringstream iss(cleaned);
    int value;
    std::vector<int> values;
    while (iss >> value) {
        values.push_back(value);
    }
    return values;
}

void OpenFOAMMeshReader::readPointsBinary(std::ifstream& file, int nPoints) {
    // Skip whitespace and read opening parenthesis
    char c;
    while (file.get(c)) {
        if (c == '(') {
            break;
        } else if (!std::isspace(c)) {
            throw std::runtime_error("Expected '(' before binary points data, got: '" + std::string(1, c) + "'");
        }
    }

    // For vectorField binary format, there's no list size prefix
    // The data starts immediately with the doubles
    points.reserve(nPoints);
    for (int i = 0; i < nPoints; ++i) {
        Point p;
        p.x = readBinary<double>(file);
        p.y = readBinary<double>(file);
        p.z = readBinary<double>(file);
        if (!file.good()) {
            throw std::runtime_error("Error reading binary point data at index " + std::to_string(i));
        }
        points.push_back(p);
    }

    // Skip whitespace and read closing parenthesis
    while (file.get(c)) {
        if (c == ')') {
            break;
        } else if (!std::isspace(c)) {
            throw std::runtime_error("Expected ')' after binary points data, got: '" + std::string(1, c) + "'");
        }
    }
}

void OpenFOAMMeshReader::readFacesBinary(std::ifstream& file, int nFaces) {
    // For faceCompactList binary format:
    // nFaces (ASCII)
    // ( offset_0 offset_1 ... offset_nFaces-1 ) (binary ints, no list size prefix)
    // Note: OpenFOAM stores only nFaces offsets, not nFaces+1
    // The last offset value equals totalData
    // totalData (ASCII)
    // ( point_0 point_1 ... point_totalData-1 ) (binary ints, no list size prefix)

    // Skip whitespace and read opening parenthesis for indices
    char c;
    while (file.get(c)) {
        if (c == '(') {
            break;
        } else if (!std::isspace(c)) {
            throw std::runtime_error("Expected '(' before binary faces indices, got: '" + std::string(1, c) + "'");
        }
    }

    // Read offset indices array - only nFaces offsets (not nFaces+1)
    std::vector<int> indices(nFaces);
    for (int i = 0; i < nFaces; ++i) {
        indices[i] = readBinary<int>(file);
        if (!file.good()) {
            throw std::runtime_error("Error reading binary face index at position " + std::to_string(i));
        }
    }

    // Skip whitespace and read closing parenthesis for indices
    while (file.get(c)) {
        if (c == ')') {
            break;
        } else if (!std::isspace(c)) {
            throw std::runtime_error("Expected ')' after binary faces indices, got: '" + std::string(1, c) + "'");
        }
    }

    // Read totalData count (ASCII)
    std::string line = getNextNonEmptyLine(file);
    int totalData;
    try {
        totalData = parseLeadingInt(line);
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing total face data count from line: '" + line + "' - " + e.what());
    }

    // Verify totalData matches last offset
    if (nFaces > 0 && totalData != indices[nFaces - 1]) {
        throw std::runtime_error("Face data count mismatch: expected " + std::to_string(indices[nFaces - 1]) + " but got " + std::to_string(totalData));
    }

    // Add implicit last offset for easier processing
    indices.push_back(totalData);

    // Skip whitespace and read opening parenthesis for data
    while (file.get(c)) {
        if (c == '(') {
            break;
        } else if (!std::isspace(c)) {
            throw std::runtime_error("Expected '(' before binary faces data, got: '" + std::string(1, c) + "'");
        }
    }

    // Read point indices data array (no list size prefix)
    std::vector<int> data(totalData);
    for (int i = 0; i < totalData; ++i) {
        data[i] = readBinary<int>(file);
        if (!file.good()) {
            throw std::runtime_error("Error reading binary face data at position " + std::to_string(i));
        }
    }

    // Skip whitespace and read closing parenthesis for data
    while (file.get(c)) {
        if (c == ')') {
            break;
        } else if (!std::isspace(c)) {
            throw std::runtime_error("Expected ')' after binary faces data, got: '" + std::string(1, c) + "'");
        }
    }

    // Construct faces from indices and data
    faces.reserve(nFaces);
    for (int i = 0; i < nFaces; ++i) {
        Face face;
        int start = indices[i];
        int end = indices[i + 1];
        for (int j = start; j < end; ++j) {
            face.pointIndices.push_back(data[j]);
        }
        faces.push_back(face);
    }
}

void OpenFOAMMeshReader::readOwnerBinary(std::ifstream& file, int nOwner) {
    // For labelList binary format: no list size prefix
    // Skip whitespace and read opening parenthesis
    char c;
    while (file.get(c)) {
        if (c == '(') {
            break;
        } else if (!std::isspace(c)) {
            throw std::runtime_error("Expected '(' before binary owner data, got: '" + std::string(1, c) + "'");
        }
    }

    // Read owner data array (no list size prefix)
    owner.reserve(nOwner);
    for (int i = 0; i < nOwner; ++i) {
        owner.push_back(readBinary<int>(file));
        if (!file.good()) {
            throw std::runtime_error("Error reading binary owner data at index " + std::to_string(i));
        }
    }

    // Skip whitespace and read closing parenthesis
    while (file.get(c)) {
        if (c == ')') {
            break;
        } else if (!std::isspace(c)) {
            throw std::runtime_error("Expected ')' after binary owner data, got: '" + std::string(1, c) + "'");
        }
    }
}

void OpenFOAMMeshReader::readNeighbourBinary(std::ifstream& file, int nNeighbour) {
    // For labelList binary format: no list size prefix
    // Skip whitespace and read opening parenthesis
    char c;
    while (file.get(c)) {
        if (c == '(') {
            break;
        } else if (!std::isspace(c)) {
            throw std::runtime_error("Expected '(' before binary neighbour data, got: '" + std::string(1, c) + "'");
        }
    }

    // Read neighbour data array (no list size prefix)
    neighbour.reserve(nNeighbour);
    for (int i = 0; i < nNeighbour; ++i) {
        neighbour.push_back(readBinary<int>(file));
        if (!file.good()) {
            throw std::runtime_error("Error reading binary neighbour data at index " + std::to_string(i));
        }
    }

    // Skip whitespace and read closing parenthesis
    while (file.get(c)) {
        if (c == ')') {
            break;
        } else if (!std::isspace(c)) {
            throw std::runtime_error("Expected ')' after binary neighbour data, got: '" + std::string(1, c) + "'");
        }
    }
}

void OpenFOAMMeshReader::readPoints() {
    std::ifstream file(meshPath + "/points", std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open points file: " + meshPath + "/points");
    }
    
    bool isBinary = skipHeaderAndDetectFormat(file);
    
    if (isBinary) {
        std::string line = getNextNonEmptyLine(file);
        int nPoints = parseLeadingInt(line);
        readPointsBinary(file, nPoints);
    } else {
        std::string line = getNextNonEmptyLine(file);
        int nPoints;
        try {
            nPoints = parseLeadingInt(line);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error parsing number of points from line: '" + line + "' - " + e.what());
        }
        
        points.reserve(nPoints);
        for (int i = 0; i < nPoints; ++i) {
            line = getNextNonEmptyLine(file);
            
            line.erase(std::remove(line.begin(), line.end(), '('), line.end());
            line.erase(std::remove(line.begin(), line.end(), ')'), line.end());
            
            std::istringstream iss(line);
            Point p;
            iss >> p.x >> p.y >> p.z;
            points.push_back(p);
        }
    }
    
    std::cout << "Read " << points.size() << " points" << std::endl;
}

void OpenFOAMMeshReader::readFaces() {
    std::ifstream file(meshPath + "/faces", std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open faces file: " + meshPath + "/faces");
    }
    
    bool isBinary = skipHeaderAndDetectFormat(file);
    
    std::string line = getNextNonEmptyLine(file);
    int nFaces;
    try {
        nFaces = parseLeadingInt(line);
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing number of faces from line: '" + line + "' - " + e.what());
    }
    
    if (isBinary) {
        readFacesBinary(file, nFaces);
    } else {
        faces.reserve(nFaces);
        for (int i = 0; i < nFaces; ++i) {
            line = getNextNonEmptyLine(file);
            
            size_t pos = line.find('(');
            if (pos != std::string::npos) {
                std::string numStr = line.substr(0, pos);
                int nPointsInFace = parseLeadingInt(numStr);
                
                std::string pointsStr = line.substr(pos + 1);
                pointsStr.erase(std::remove(pointsStr.begin(), pointsStr.end(), ')'), pointsStr.end());
                
                Face face;
                std::istringstream iss(pointsStr);
                int pointIdx;
                while (iss >> pointIdx) {
                    face.pointIndices.push_back(pointIdx);
                }
                faces.push_back(face);
            }
        }
    }
    
    std::cout << "Read " << faces.size() << " faces" << std::endl;
}

void OpenFOAMMeshReader::readOwner() {
    std::ifstream file(meshPath + "/owner", std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open owner file: " + meshPath + "/owner");
    }
    
    bool isBinary = skipHeaderAndDetectFormat(file);
    
    std::string line = getNextNonEmptyLine(file);
    int nFaces;
    size_t consumedChars = 0;
    try {
        nFaces = parseLeadingInt(line, &consumedChars);
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing number of faces from owner file, line: '" + line + "' - " + e.what());
    }
    
    if (isBinary) {
        readOwnerBinary(file, nFaces);
    } else {
        owner.reserve(nFaces);
        
        std::vector<int> pending = parseIntListFromString(line.substr(consumedChars));
        size_t pendingIdx = 0;
        while (pendingIdx < pending.size() && owner.size() < nFaces) {
            owner.push_back(pending[pendingIdx++]);
        }
        
        while (owner.size() < nFaces) {
            line = getNextNonEmptyLine(file);
            auto values = parseIntListFromString(line);
            for (int value : values) {
                owner.push_back(value);
                if (owner.size() == nFaces) break;
            }
        }
    }
    
    std::cout << "Read owner information for " << owner.size() << " faces" << std::endl;
}

void OpenFOAMMeshReader::readNeighbour() {
    std::ifstream file(meshPath + "/neighbour", std::ios::binary);
    if (!file.is_open()) {
        std::cout << "Note: neighbour file not found (may be empty for some meshes)" << std::endl;
        return;
    }
    
    bool isBinary = skipHeaderAndDetectFormat(file);
    
    std::string line = getNextNonEmptyLine(file);
    int nInternalFaces;
    size_t consumedChars = 0;
    try {
        nInternalFaces = parseLeadingInt(line, &consumedChars);
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing number of internal faces from neighbour file, line: '" + line + "' - " + e.what());
    }
    
    if (isBinary) {
        readNeighbourBinary(file, nInternalFaces);
    } else {
        neighbour.reserve(nInternalFaces);
        
        std::vector<int> pending = parseIntListFromString(line.substr(consumedChars));
        size_t pendingIdx = 0;
        while (pendingIdx < pending.size() && neighbour.size() < nInternalFaces) {
            neighbour.push_back(pending[pendingIdx++]);
        }
        
        while (neighbour.size() < nInternalFaces) {
            line = getNextNonEmptyLine(file);
            auto values = parseIntListFromString(line);
            for (int value : values) {
                neighbour.push_back(value);
                if (neighbour.size() == nInternalFaces) break;
            }
        }
    }
    
    std::cout << "Read neighbour information for " << neighbour.size() << " internal faces" << std::endl;
}

void OpenFOAMMeshReader::readBoundary() {
    std::ifstream file(meshPath + "/boundary");
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open boundary file: " + meshPath + "/boundary");
    }
    
    skipHeaderAndDetectFormat(file);
    
    std::string line = getNextNonEmptyLine(file);
    int nPatches;
    try {
        nPatches = parseLeadingInt(line);
    } catch (const std::exception& e) {
        throw std::runtime_error("Error parsing number of patches from boundary file, line: '" + line + "' - " + e.what());
    }
    
    for (int i = 0; i < nPatches; ++i) {
        BoundaryPatch patch;
        
        line = getNextNonEmptyLine(file);
        patch.name = line;
        patch.name.erase(std::remove(patch.name.begin(), patch.name.end(), ' '), patch.name.end());
        patch.name.erase(std::remove(patch.name.begin(), patch.name.end(), '\t'), patch.name.end());
        
        std::getline(file, line);
        
        while (std::getline(file, line)) {
            line = stripComments(line);
            if (line.find("type") != std::string::npos) {
                size_t pos = line.find("type");
                std::string typeStr = line.substr(pos + 4);
                typeStr.erase(std::remove(typeStr.begin(), typeStr.end(), ' '), typeStr.end());
                typeStr.erase(std::remove(typeStr.begin(), typeStr.end(), ';'), typeStr.end());
                patch.type = typeStr;
            } else if (line.find("nFaces") != std::string::npos) {
                size_t pos = line.find("nFaces");
                std::string nFacesStr = line.substr(pos + 6);
                nFacesStr.erase(std::remove(nFacesStr.begin(), nFacesStr.end(), ' '), nFacesStr.end());
                nFacesStr.erase(std::remove(nFacesStr.begin(), nFacesStr.end(), ';'), nFacesStr.end());
                patch.nFaces = parseLeadingInt(nFacesStr);
            } else if (line.find("startFace") != std::string::npos) {
                size_t pos = line.find("startFace");
                std::string startFaceStr = line.substr(pos + 9);
                startFaceStr.erase(std::remove(startFaceStr.begin(), startFaceStr.end(), ' '), startFaceStr.end());
                startFaceStr.erase(std::remove(startFaceStr.begin(), startFaceStr.end(), ';'), startFaceStr.end());
                patch.startFace = parseLeadingInt(startFaceStr);
            } else if (line.find('}') != std::string::npos) {
                break;
            }
        }
        
        boundaryPatches.push_back(patch);
        std::cout << "  Boundary patch: " << patch.name << " (type: " << patch.type 
                  << ", faces: " << patch.nFaces << ", start: " << patch.startFace << ")" << std::endl;
    }
}

void OpenFOAMMeshReader::readCellZones() {
    std::ifstream file(meshPath + "/cellZones");
    if (!file.is_open()) {
        std::cout << "No cellZones file found (optional)" << std::endl;
        return;
    }
    
    skipHeaderAndDetectFormat(file);
    
    std::string line = getNextNonEmptyLine(file);
    int nZones;
    try {
        nZones = parseLeadingInt(line);
    } catch (const std::exception& e) {
        std::cout << "Could not parse number of zones, assuming 0 zones" << std::endl;
        return;
    }
    
    if (nZones == 0) {
        std::cout << "No cell zones defined" << std::endl;
        return;
    }
    
    std::getline(file, line);
    
    for (int i = 0; i < nZones; ++i) {
        CellZone zone;
        
        line = getNextNonEmptyLine(file);
        zone.name = line;
        zone.name.erase(std::remove(zone.name.begin(), zone.name.end(), ' '), zone.name.end());
        zone.name.erase(std::remove(zone.name.begin(), zone.name.end(), '\t'), zone.name.end());
        
        while (std::getline(file, line)) {
            line = stripComments(line);
            
            size_t firstNonSpace = line.find_first_not_of(" \t\r\n");
            if (firstNonSpace != std::string::npos) {
                line.erase(0, firstNonSpace);
            } else {
                line.clear();
            }
            
            size_t lastNonSpace = line.find_last_not_of(" \t\r\n");
            if (lastNonSpace != std::string::npos) {
                line.erase(lastNonSpace + 1);
            } else {
                line.clear();
            }
            
            if (line.empty() || line == "{" || line == "(") {
                continue;
            }
            
            if (line.find("cellLabels") != std::string::npos) {
                size_t parenPos = line.find('(');
                
                if (parenPos != std::string::npos) {
                    std::string beforeParen = line.substr(0, parenPos);
                    std::istringstream iss(beforeParen);
                    std::string word;
                    int nCells = 0;
                    while (iss >> word) {
                        if (word != "cellLabels") {
                            try {
                                nCells = parseLeadingInt(word);
                            } catch (...) {}
                        }
                    }
                    
                    std::string cellData = line.substr(parenPos + 1);
                    cellData.erase(std::remove(cellData.begin(), cellData.end(), ')'), cellData.end());
                    cellData.erase(std::remove(cellData.begin(), cellData.end(), ';'), cellData.end());
                    
                    std::istringstream cellStream(cellData);
                    int cellIdx;
                    while (cellStream >> cellIdx) {
                        zone.cellIndices.push_back(cellIdx);
                    }
                } else {
                    line = getNextNonEmptyLine(file);
                    int nCells;
                    try {
                        nCells = parseLeadingInt(line);
                    } catch (const std::exception& e) {
                        throw std::runtime_error("Error parsing number of cells in zone, line: '" + line + "' - " + e.what());
                    }
                    
                    std::getline(file, line);
                    
                    zone.cellIndices.reserve(nCells);
                    for (int j = 0; j < nCells; ++j) {
                        line = getNextNonEmptyLine(file);
                        zone.cellIndices.push_back(parseLeadingInt(line));
                    }
                }
                
                break;
            } else if (line.find('}') != std::string::npos) {
                break;
            }
        }
        
        cellZones.push_back(zone);
        std::cout << "  Cell zone: " << zone.name << " (" << zone.cellIndices.size() << " cells)" << std::endl;
    }
}

void OpenFOAMMeshReader::constructCells() {
    if (owner.empty()) {
        throw std::runtime_error("Owner data not available");
    }
    
    int nCells = *std::max_element(owner.begin(), owner.end()) + 1;
    cells.resize(nCells);
    
    for (size_t faceIdx = 0; faceIdx < owner.size(); ++faceIdx) {
        int cellIdx = owner[faceIdx];
        cells[cellIdx].faceIndices.push_back(faceIdx);
    }
    
    for (size_t faceIdx = 0; faceIdx < neighbour.size(); ++faceIdx) {
        int cellIdx = neighbour[faceIdx];
        if (cellIdx >= 0 && cellIdx < nCells) {
            cells[cellIdx].faceIndices.push_back(faceIdx);
        }
    }
    
    for (auto& cell : cells) {
        int nFaces = cell.faceIndices.size();
        int nTri = 0, nQuad = 0;
        for (int fi : cell.faceIndices) {
            if (fi >= 0 && fi < (int)faces.size()) {
                int nPts = faces[fi].pointIndices.size();
                if (nPts == 3) ++nTri;
                else if (nPts == 4) ++nQuad;
            }
        }
        if      (nFaces == 4 && nTri == 4 && nQuad == 0) cell.type = "tet";
        else if (nFaces == 6 && nQuad == 6 && nTri == 0) cell.type = "hex";
        else if (nFaces == 5 && nQuad == 1 && nTri == 4) cell.type = "pyr";
        else if (nFaces == 5 && nQuad == 3 && nTri == 2) cell.type = "wedge";
        else                                              cell.type = "unknown";
    }
    
    std::cout << "Constructed " << cells.size() << " cells" << std::endl;
}
