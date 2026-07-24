// File       : OpenFOAMFieldReader.h
// Description: Standalone reader for OpenFOAM volume fields.
// SPDX-License-Identifier: BSD-3-Clause

#ifndef OPENFOAM_FIELD_READER_H
#define OPENFOAM_FIELD_READER_H

#include "OpenFOAMMeshReader.h"
#include <map>
#include <string>
#include <vector>

struct VolFieldInfo
{
    std::string path;
    std::string name;
    std::string outputName;
    std::string foamClass;
    int numComponents = 0;
    std::vector<std::string> componentNames;
};

struct BoundaryFieldValues
{
    std::string type;
    bool hasValues = false;
    // Face-major storage: values[face * numComponents + component].
    std::vector<double> values;
};

struct VolField
{
    VolFieldInfo info;
    // Cell-major storage: internalValues[cell * numComponents + component].
    std::vector<double> internalValues;
    std::map<std::string, BoundaryFieldValues> boundary;
};

class OpenFOAMFieldReader
{
public:
    OpenFOAMFieldReader(const std::string& casePath,
                        const std::string& timeSelection,
                        const OpenFOAMMeshReader& mesh);

    const std::string& getTimeName() const
    {
        return timeName;
    }

    double getTimeValue() const;

    // "all" selects every supported volume field in the selected time
    // directory. Otherwise, names are OpenFOAM object/file names.
    std::vector<VolFieldInfo>
    selectFields(const std::vector<std::string>& names, bool all) const;

    VolField readField(const VolFieldInfo& info) const;

private:
    std::string casePath;
    std::string timeName;
    const OpenFOAMMeshReader& mesh;

    static VolFieldInfo inspectField(const std::string& path);
    static int componentCount(const std::string& foamClass);
    static std::vector<std::string> componentSuffixes(int count);
    static std::string chooseTime(const std::string& casePath,
                                  const std::string& selection);
};

#endif
