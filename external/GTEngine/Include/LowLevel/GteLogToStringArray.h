// Geometric Tools LLC, Redmond WA 98052
// Copyright (c) 1998-2015
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 2.0.0 (2015/09/23)

#pragma once

#include <LowLevel/GteLogger.h>
#include <vector>

namespace gte
{

class GTE_IMPEXP LogToStringArray : public Logger::Listener
{
public:
    LogToStringArray(std::string const& name, int flags);

    std::string const& GetName() const;
    std::vector<std::string> const& GetMessages() const;
    std::vector<std::string>& GetMessages();

private:
    virtual void Report(std::string const& message);

    std::string mName;
    std::vector<std::string> mMessages;
};

}
