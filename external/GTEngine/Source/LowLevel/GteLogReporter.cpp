// David Eberly, Geometric Tools, Redmond WA 98052
// Copyright (c) 1998-2017
// Distributed under the Boost Software License, Version 1.0.
// http://www.boost.org/LICENSE_1_0.txt
// http://www.geometrictools.com/License/Boost/LICENSE_1_0.txt
// File Version: 3.0.0 (2016/06/19)

#include <GTEnginePCH.h>
#include <LowLevel/GteLogReporter.h>
using namespace gte;

LogReporter::~LogReporter()
{
    if (mLogToStdout)
    {
        Logger::Unsubscribe(mLogToStdout.get());
    }

    if (mLogToFile)
    {
        Logger::Unsubscribe(mLogToFile.get());
    }

#if defined(__MSWINDOWS__)
    if (mLogToOutputWindow)
    {
        Logger::Unsubscribe(mLogToOutputWindow.get());
    }

    if (mLogToMessageBox)
    {
        Logger::Unsubscribe(mLogToMessageBox.get());
    }
#endif
}

LogReporter::LogReporter(std::string const& logFile, int logFileFlags,
    int logStdoutFlags, int logMessageBoxFlags, int logOutputWindowFlags)
    :
    mLogToFile(nullptr),
    mLogToStdout(nullptr)
#if defined(__MSWINDOWS__)
    ,
    mLogToMessageBox(nullptr),
    mLogToOutputWindow(nullptr)
#endif
{
    if (logFileFlags != Logger::Listener::LISTEN_FOR_NOTHING)
    {
        mLogToFile = std::make_unique<LogToFile>(logFile, logFileFlags);
        Logger::Subscribe(mLogToFile.get());
    }

    if (logStdoutFlags != Logger::Listener::LISTEN_FOR_NOTHING)
    {
        mLogToStdout = std::make_unique<LogToStdout>(logStdoutFlags);
        Logger::Subscribe(mLogToStdout.get());
    }

#if defined(__MSWINDOWS__)
    if (logMessageBoxFlags != Logger::Listener::LISTEN_FOR_NOTHING)
    {
        mLogToMessageBox = std::make_unique<LogToMessageBox>(logMessageBoxFlags);
        Logger::Subscribe(mLogToMessageBox.get());
    }

    if (logOutputWindowFlags != Logger::Listener::LISTEN_FOR_NOTHING)
    {
        mLogToOutputWindow = std::make_unique<LogToOutputWindow>(logOutputWindowFlags);
        Logger::Subscribe(mLogToOutputWindow.get());
    }
#endif
}
