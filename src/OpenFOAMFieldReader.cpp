// File       : OpenFOAMFieldReader.cpp
// Description: Standalone ASCII/binary OpenFOAM volume-field reader.
// SPDX-License-Identifier: BSD-3-Clause

#include "OpenFOAMFieldReader.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace
{

std::vector<char> readFile(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input)
        throw std::runtime_error("Cannot open OpenFOAM field: " + path);
    input.seekg(0, std::ios::end);
    const std::streamoff length = input.tellg();
    input.seekg(0, std::ios::beg);
    std::vector<char> data(static_cast<size_t>(length));
    if (length > 0)
        input.read(data.data(), length);
    if (!input && length > 0)
        throw std::runtime_error("Failed while reading OpenFOAM field: " +
                                 path);
    return data;
}

std::string headerValue(const std::string& text, const std::string& key)
{
    const std::regex expression("\\b" + key + R"(\s+("?[^";\s]+"?)\s*;)");
    std::smatch match;
    if (!std::regex_search(text, match, expression))
        return "";
    std::string value = match[1].str();
    if (value.size() >= 2 && value.front() == '"' && value.back() == '"')
        value = value.substr(1, value.size() - 2);
    return value;
}

std::string headerText(const std::vector<char>& bytes)
{
    const std::string all(bytes.begin(), bytes.end());
    const size_t foam = all.find("FoamFile");
    if (foam == std::string::npos)
        throw std::runtime_error("Missing FoamFile header");
    const size_t open = all.find('{', foam);
    if (open == std::string::npos)
        throw std::runtime_error("Malformed FoamFile header");
    int depth = 0;
    for (size_t i = open; i < all.size(); ++i)
    {
        if (all[i] == '{')
            ++depth;
        else if (all[i] == '}' && --depth == 0)
            return all.substr(foam, i - foam + 1);
    }
    throw std::runtime_error("Unterminated FoamFile header");
}

bool hostIsLittleEndian()
{
    const unsigned short one = 1;
    return *reinterpret_cast<const unsigned char*>(&one) == 1;
}

class FoamLexer
{
public:
    FoamLexer(const std::vector<char>& bytes,
              bool binary,
              bool fileLittleEndian,
              int scalarBytes)
        : bytes(bytes), binary(binary), fileLittleEndian(fileLittleEndian),
          scalarBytes(scalarBytes)
    {
    }

    bool eof()
    {
        skipSpaceAndComments();
        return pos >= bytes.size();
    }

    std::string next()
    {
        skipSpaceAndComments();
        if (pos >= bytes.size())
            return "";

        const char c = bytes[pos];
        if (std::strchr("{}();[]", c))
        {
            ++pos;
            return std::string(1, c);
        }
        if (c == '"')
        {
            const size_t begin = ++pos;
            while (pos < bytes.size() && bytes[pos] != '"')
                ++pos;
            const std::string token(bytes.data() + begin, pos - begin);
            if (pos < bytes.size())
                ++pos;
            return token;
        }

        const size_t begin = pos;
        while (pos < bytes.size())
        {
            const char value = bytes[pos];
            if (std::isspace(static_cast<unsigned char>(value)) ||
                std::strchr("{}();[]\"", value))
                break;
            ++pos;
        }
        return std::string(bytes.data() + begin, pos - begin);
    }

    std::vector<double> readEntry(const std::string& modifier,
                                  int expectedComponents,
                                  int expectedCount,
                                  bool expandUniform)
    {
        if (modifier == "uniform")
        {
            std::vector<double> tuple;
            while ((int)tuple.size() < expectedComponents)
            {
                const std::string token = next();
                if (token.empty())
                    throw std::runtime_error(
                        "Unexpected end of uniform field value");
                if (token == "(" || token == ")" || token == ";")
                    continue;
                tuple.push_back(parseDouble(token));
            }
            if (!expandUniform)
                return tuple;
            std::vector<double> result(static_cast<size_t>(expectedCount) *
                                       expectedComponents);
            for (int i = 0; i < expectedCount; ++i)
                std::copy(tuple.begin(),
                          tuple.end(),
                          result.begin() +
                              static_cast<size_t>(i) * expectedComponents);
            return result;
        }

        if (modifier != "nonuniform")
            throw std::runtime_error("Unsupported OpenFOAM field value '" +
                                     modifier + "'");

        const std::string listType = next();
        const int storedComponents = componentsFromListType(listType);
        const std::string countToken = next();
        const int count = parseInt(countToken);
        if (next() != "(")
            throw std::runtime_error(
                "Expected '(' before nonuniform field data");

        std::vector<double> stored;
        stored.reserve(static_cast<size_t>(count) * storedComponents);
        if (binary)
        {
            const size_t valueCount =
                static_cast<size_t>(count) * storedComponents;
            const size_t byteCount = valueCount * scalarBytes;
            if (pos + byteCount > bytes.size())
                throw std::runtime_error(
                    "Truncated binary OpenFOAM field data");
            for (size_t i = 0; i < valueCount; ++i)
            {
                stored.push_back(readBinaryScalar(pos));
                pos += scalarBytes;
            }
        }
        else
        {
            int listDepth = 1;
            while ((int)stored.size() < count * storedComponents)
            {
                const std::string token = next();
                if (token.empty())
                    throw std::runtime_error(
                        "Unexpected end of nonuniform field data");
                if (token == "(")
                {
                    ++listDepth;
                    continue;
                }
                if (token == ")")
                {
                    --listDepth;
                    continue;
                }
                if (token == ";")
                    continue;
                stored.push_back(parseDouble(token));
            }
            while (listDepth > 0)
            {
                const std::string token = next();
                if (token.empty())
                    throw std::runtime_error(
                        "Unexpected end of nonuniform field data");
                if (token == "(")
                    ++listDepth;
                else if (token == ")")
                    --listDepth;
            }
            return finishList(stored,
                              storedComponents,
                              expectedComponents,
                              count,
                              expectedCount);
        }

        // Consume the list closing delimiter. Binary data ends exactly before
        // it.
        std::string close = next();
        while (!close.empty() && close != ")")
            close = next();
        if (close != ")")
            throw std::runtime_error(
                "Expected ')' after nonuniform field data");

        return finishList(
            stored, storedComponents, expectedComponents, count, expectedCount);
    }

    // Skip an unknown nonuniform entry without interpreting its values.
    void skipNonuniform()
    {
        const std::string listType = next();
        const int components = componentsFromListType(listType);
        const int count = parseInt(next());
        if (next() != "(")
            throw std::runtime_error(
                "Expected '(' before nonuniform dictionary data");
        if (binary)
        {
            const size_t byteCount =
                static_cast<size_t>(count) * components * scalarBytes;
            if (pos + byteCount > bytes.size())
                throw std::runtime_error(
                    "Truncated binary OpenFOAM dictionary data");
            pos += byteCount;
            std::string close = next();
            while (!close.empty() && close != ")")
                close = next();
            if (close != ")")
                throw std::runtime_error(
                    "Expected ')' after nonuniform dictionary data");
            return;
        }

        int values = 0;
        int listDepth = 1;
        while (values < count * components)
        {
            const std::string token = next();
            if (token.empty())
                throw std::runtime_error(
                    "Unexpected end of nonuniform dictionary data");
            if (token == "(")
            {
                ++listDepth;
                continue;
            }
            if (token == ")")
            {
                --listDepth;
                continue;
            }
            if (token == ";")
                continue;
            (void)parseDouble(token);
            ++values;
        }
        while (listDepth > 0)
        {
            const std::string token = next();
            if (token.empty())
                throw std::runtime_error(
                    "Unexpected end of nonuniform dictionary data");
            if (token == "(")
                ++listDepth;
            else if (token == ")")
                --listDepth;
        }
    }

private:
    static std::vector<double> finishList(std::vector<double> stored,
                                          int storedComponents,
                                          int expectedComponents,
                                          int count,
                                          int expectedCount)
    {
        if (count != expectedCount)
            throw std::runtime_error(
                "OpenFOAM field list has " + std::to_string(count) +
                " entries; expected " + std::to_string(expectedCount));
        if (storedComponents != expectedComponents)
            throw std::runtime_error(
                "OpenFOAM field component count does not match its class");
        return stored;
    }

    const std::vector<char>& bytes;
    size_t pos = 0;
    bool binary;
    bool fileLittleEndian;
    int scalarBytes;

    void skipSpaceAndComments()
    {
        for (;;)
        {
            while (pos < bytes.size() &&
                   std::isspace(static_cast<unsigned char>(bytes[pos])))
                ++pos;
            if (pos + 1 >= bytes.size() || bytes[pos] != '/')
                return;
            if (bytes[pos + 1] == '/')
            {
                pos += 2;
                while (pos < bytes.size() && bytes[pos] != '\n')
                    ++pos;
            }
            else if (bytes[pos + 1] == '*')
            {
                pos += 2;
                while (pos + 1 < bytes.size() &&
                       !(bytes[pos] == '*' && bytes[pos + 1] == '/'))
                    ++pos;
                if (pos + 1 < bytes.size())
                    pos += 2;
            }
            else
                return;
        }
    }

    static int parseInt(const std::string& token)
    {
        size_t used = 0;
        const long value = std::stol(token, &used);
        if (used != token.size() || value < 0 ||
            value > std::numeric_limits<int>::max())
            throw std::runtime_error("Invalid OpenFOAM list size '" + token +
                                     "'");
        return static_cast<int>(value);
    }

    static double parseDouble(const std::string& token)
    {
        size_t used = 0;
        const double value = std::stod(token, &used);
        if (used != token.size())
            throw std::runtime_error("Invalid OpenFOAM scalar '" + token + "'");
        return value;
    }

    static int componentsFromListType(const std::string& type)
    {
        if (type.find("sphericalTensor") != std::string::npos)
            return 1;
        if (type.find("symmTensor") != std::string::npos)
            return 6;
        if (type.find("tensor") != std::string::npos)
            return 9;
        if (type.find("vector") != std::string::npos)
            return 3;
        if (type.find("scalar") != std::string::npos)
            return 1;
        throw std::runtime_error("Unsupported OpenFOAM list type '" + type +
                                 "'");
    }

    double readBinaryScalar(size_t offset) const
    {
        if (scalarBytes != 4 && scalarBytes != 8)
            throw std::runtime_error(
                "Only 32-bit and 64-bit OpenFOAM scalars are supported");
        unsigned char raw[8] = {};
        std::memcpy(raw, bytes.data() + offset, scalarBytes);
        if (fileLittleEndian != hostIsLittleEndian())
            std::reverse(raw, raw + scalarBytes);
        if (scalarBytes == 8)
        {
            double value;
            std::memcpy(&value, raw, 8);
            return value;
        }
        float value;
        std::memcpy(&value, raw, 4);
        return value;
    }
};

bool isNumericDirectory(const std::string& name, double& value)
{
    try
    {
        size_t used = 0;
        value = std::stod(name, &used);
        return used == name.size() && std::isfinite(value);
    }
    catch (...)
    {
        return false;
    }
}

} // namespace

OpenFOAMFieldReader::OpenFOAMFieldReader(const std::string& casePath,
                                         const std::string& timeSelection,
                                         const OpenFOAMMeshReader& mesh)
    : casePath(casePath), timeName(chooseTime(casePath, timeSelection)),
      mesh(mesh)
{
}

double OpenFOAMFieldReader::getTimeValue() const
{
    double value = 0.0;
    if (!isNumericDirectory(timeName, value))
        throw std::runtime_error("Selected OpenFOAM time '" + timeName +
                                 "' is not numeric");
    return value;
}

int OpenFOAMFieldReader::componentCount(const std::string& foamClass)
{
    if (foamClass == "volScalarField" || foamClass == "volSphericalTensorField")
        return 1;
    if (foamClass == "volVectorField")
        return 3;
    if (foamClass == "volSymmTensorField")
        return 6;
    if (foamClass == "volTensorField")
        return 9;
    return 0;
}

std::vector<std::string> OpenFOAMFieldReader::componentSuffixes(int count)
{
    if (count == 1)
        return {""};
    if (count == 3)
        return {"x", "y", "z"};
    if (count == 6)
        return {"xx", "xy", "xz", "yy", "yz", "zz"};
    if (count == 9)
        return {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"};
    return {};
}

VolFieldInfo OpenFOAMFieldReader::inspectField(const std::string& path)
{
    const std::vector<char> bytes = readFile(path);
    const std::string header = headerText(bytes);
    VolFieldInfo info;
    info.path = path;
    info.foamClass = headerValue(header, "class");
    info.numComponents = componentCount(info.foamClass);
    info.name = headerValue(header, "object");
    if (info.name.empty())
        info.name = std::filesystem::path(path).filename().string();
    info.outputName = info.name;
    info.componentNames = componentSuffixes(info.numComponents);
    return info;
}

std::string OpenFOAMFieldReader::chooseTime(const std::string& casePath,
                                            const std::string& selection)
{
    namespace fs = std::filesystem;
    if (selection != "latest" && selection != "latestTime")
    {
        const fs::path selected = fs::path(casePath) / selection;
        if (!fs::is_directory(selected))
            throw std::runtime_error(
                "OpenFOAM time directory does not exist: " + selected.string());
        return selection;
    }

    double latest = -std::numeric_limits<double>::infinity();
    std::string latestName;
    for (const fs::directory_entry& entry :
         fs::directory_iterator(fs::path(casePath)))
    {
        if (!entry.is_directory())
            continue;
        double value = 0.0;
        const std::string name = entry.path().filename().string();
        if (isNumericDirectory(name, value) && value > latest)
        {
            latest = value;
            latestName = name;
        }
    }
    if (latestName.empty())
        throw std::runtime_error(
            "No numeric OpenFOAM time directories found in " + casePath);
    return latestName;
}

std::vector<VolFieldInfo>
OpenFOAMFieldReader::selectFields(const std::vector<std::string>& names,
                                  bool all) const
{
    namespace fs = std::filesystem;
    const fs::path directory = fs::path(casePath) / timeName;
    std::vector<VolFieldInfo> result;

    if (all)
    {
        for (const fs::directory_entry& entry :
             fs::directory_iterator(directory))
        {
            if (!entry.is_regular_file())
                continue;
            VolFieldInfo info = inspectField(entry.path().string());
            if (info.numComponents > 0)
                result.push_back(std::move(info));
        }
        std::sort(result.begin(),
                  result.end(),
                  [](const VolFieldInfo& a, const VolFieldInfo& b)
        { return a.name < b.name; });
    }
    else
    {
        for (const std::string& name : names)
        {
            const fs::path path = directory / name;
            if (!fs::is_regular_file(path))
                throw std::runtime_error(
                    "Requested OpenFOAM field not found: " + path.string());
            VolFieldInfo info = inspectField(path.string());
            if (info.numComponents == 0)
                throw std::runtime_error(
                    "Requested object '" + name +
                    "' is not a supported volume field (class=" +
                    info.foamClass + ")");
            result.push_back(std::move(info));
        }
    }

    if (result.empty())
        throw std::runtime_error(
            "No supported OpenFOAM volume fields selected");
    return result;
}

VolField OpenFOAMFieldReader::readField(const VolFieldInfo& info) const
{
    const std::vector<char> bytes = readFile(info.path);
    const std::string header = headerText(bytes);
    const bool binary = headerValue(header, "format") == "binary";
    std::string arch;
    std::smatch archMatch;
    if (std::regex_search(
            header, archMatch, std::regex(R"foam(\barch\s+"([^"]+)"\s*;)foam")))
        arch = archMatch[1].str();
    const bool littleEndian = arch.find("MSB") == std::string::npos;
    int scalarBytes = 8;
    std::smatch scalarMatch;
    if (std::regex_search(
            arch, scalarMatch, std::regex(R"(scalar\s*=\s*(32|64))")))
        scalarBytes = std::stoi(scalarMatch[1].str()) / 8;

    FoamLexer lexer(bytes, binary, littleEndian, scalarBytes);
    VolField field;
    field.info = info;

    std::string token;
    while (!(token = lexer.next()).empty() && token != "internalField")
    {
    }
    if (token.empty())
        throw std::runtime_error("Missing internalField in " + info.path);
    field.internalValues = lexer.readEntry(
        lexer.next(), info.numComponents, mesh.getNumCells(), true);

    while (!(token = lexer.next()).empty() && token != "boundaryField")
    {
    }
    if (token.empty())
        throw std::runtime_error("Missing boundaryField in " + info.path);
    if (lexer.next() != "{")
        throw std::runtime_error("Malformed boundaryField in " + info.path);

    std::map<std::string, int> patchSizes;
    for (const BoundaryPatch& patch : mesh.getBoundaryPatches())
        patchSizes[patch.name] = patch.nFaces;

    while (!lexer.eof())
    {
        const std::string patchName = lexer.next();
        if (patchName == "}")
            break;
        if (patchName.empty())
            break;
        if (lexer.next() != "{")
            throw std::runtime_error("Malformed boundary patch '" + patchName +
                                     "' in " + info.path);

        BoundaryFieldValues patchValues;
        int depth = 1;
        while (depth > 0 && !lexer.eof())
        {
            const std::string item = lexer.next();
            if (item == "{")
            {
                ++depth;
                continue;
            }
            if (item == "}")
            {
                --depth;
                continue;
            }
            if (item == "nonuniform")
            {
                lexer.skipNonuniform();
                continue;
            }
            if (depth != 1)
                continue;
            if (item == "type")
            {
                patchValues.type = lexer.next();
            }
            else if (item == "value")
            {
                const auto sizeIt = patchSizes.find(patchName);
                if (sizeIt == patchSizes.end())
                    throw std::runtime_error(
                        "Field contains unknown boundary patch '" + patchName +
                        "'");
                const std::string modifier = lexer.next();
                if (modifier == "uniform" || modifier == "nonuniform")
                {
                    patchValues.values = lexer.readEntry(
                        modifier, info.numComponents, sizeIt->second, true);
                    patchValues.hasValues = true;
                }
            }
        }
        field.boundary[patchName] = std::move(patchValues);
    }

    return field;
}
