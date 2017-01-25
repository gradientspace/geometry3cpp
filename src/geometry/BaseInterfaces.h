#pragma once

#include <g3types.h>


namespace g3
{

//
// This class provides an incrementable counter that we
//   use for change tracking in various places. Many core
//   geometry classes implement this interface.
//
class g3External ITimeStamped
{
public:
	virtual ~ITimeStamped() {}

	virtual void updateTimeStamp() {
		m_nTimestamp++;
	}
	virtual unsigned int getTimeStamp() const {
		return m_nTimestamp;
	}

private:
	unsigned int m_nTimestamp = 1;
};




}

