#pragma once

#include <g3types.h>

namespace g3 
{

class Colors
{
public:
	//
	// Standard extreme-rgb colors
	//
	static const Color4b White;
	static const Color4b Black;
	static const Color4b Blue;
	static const Color4b Green;
	static const Color4b Red;
	static const Color4b Yellow;
	static const Color4b Cyan;
	static const Color4b Magenta;

	//
	// Video-codec-safe versions of standard colors.
	// These will not "shimmer" like full-saturation colors will under certain codecs.
	//
	static const Color4b VideoWhite;
	static const Color4b VideoBlack;
	static const Color4b VideoBlue;
	static const Color4b VideoGreen;
	static const Color4b VideoRed;
	static const Color4b VideoYellow;
	static const Color4b VideoCyan;
	static const Color4b VideoMagenta;

	//
	// Miscellaneous color set
	//
	static const Color4b Purple;
	static const Color4b Orange;
	static const Color4b Gold;
	static const Color4b DarkYellow;
	static const Color4b BlueMetal;

};


}

