#pragma once

#include <g3types.h>

namespace g3
{

class IPixelHitRadius
{
public:
	virtual ~IPixelHitRadius() {}
	virtual float GetWorldHitRadius( const g3::Vector3f & vHit ) const = 0;
};

class HitTestRay
{
public:
	g3::Vector3f vOrigin;
	g3::Vector3f vDirection;
	IPixelHitRadius * pHitThresh;

	HitTestRay() { 
		vOrigin = vDirection = g3::Vector3f::ZERO; pHitThresh = nullptr; 
	}
	HitTestRay( const g3::Vector3f & o, const g3::Vector3f & d ) {
		vOrigin = o; vDirection = d; pHitThresh = nullptr; 
	}
	HitTestRay( const g3::Vector3f & o, const g3::Vector3f & d, IPixelHitRadius * pPixelRadius ) {
		vOrigin = o; vDirection = d; pHitThresh = pPixelRadius; 
	}
};


}

