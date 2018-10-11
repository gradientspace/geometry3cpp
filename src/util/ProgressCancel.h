#pragma once

#include <g3types.h>

namespace g3
{



/// <summary>
/// interface that provides a cancel function
/// </summary>
class ICancelSource
{
public:
	virtual bool Cancelled() = 0;
};


/// <summary>
/// Just wraps a func<bool> as an ICancelSource
/// </summary>
class CancelFunction : public ICancelSource
{
public:
	const std::function<bool()> & CancelF;
    CancelFunction(const std::function<bool()> & cancelF) : CancelF(cancelF) {
    }
    bool Cancelled() { return CancelF(); }
};


/// <summary>
/// This class is intended to be passed to long-running computes to 
///  1) provide progress info back to caller (not implemented yet)
///  2) allow caller to cancel the computation
/// </summary>
class ProgressCancel
{
public:
	std::shared_ptr<ICancelSource> Source;

    bool WasCancelled = false;  // will be set to true if CancelF() ever returns true

    ProgressCancel(std::shared_ptr<ICancelSource> source)
    {
        Source = source;
    }
    ProgressCancel(const std::function<bool()> & cancelF)
    {
        Source = std::make_shared<CancelFunction>(cancelF);
    }

    /// <summary>
    /// Check if client would like to cancel
    /// </summary>
    bool Cancelled()
    {
        if (WasCancelled)
            return true;
        WasCancelled = Source->Cancelled();
        return WasCancelled;
    }
};



}