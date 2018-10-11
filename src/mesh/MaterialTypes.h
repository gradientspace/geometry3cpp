#pragma once

#include <g3types.h>

namespace g3
{

class GenericMaterial
{
public:
	static constexpr float Invalidf = std::numeric_limits<float>::max();
	static Vector3f InvalidColor() { return Vector3f(-1, -1, -1); }

	std::string name;
	int id;

	Vector3f diffuse_color;
	float alpha;

	virtual Vector3f DiffuseColor() = 0;
	virtual float Alpha() = 0;

	enum class KnownMaterialTypes
	{
		OBJ_MTL_Format
	};
	KnownMaterialTypes Type;


	GenericMaterial() {
		diffuse_color = InvalidColor();
		alpha = Invalidf;
	}
};
typedef std::shared_ptr<GenericMaterial> GenericMaterialPtr;





// details: http://www.fileformat.info/format/material/
// Note: if value is initialized to Invalid vector, -1, or NaN, it was not defined in material file
class OBJMaterial : public GenericMaterial
{
public:
	Vector3f Ka;     // rgb ambient reflectivity
	Vector3f Kd;     // rgb diffuse reflectivity 
	Vector3f Ks;     // rgb specular reflectivity
	Vector3f Ke;     // rgb emissive
	Vector3f Tf;        // rgb transmission filter
	int illum;          // illumination model 0-10
	float d;            // dissolve (alpha)
	float Ns;           // specular exponent (shininess)
	float sharpness;    // reflection sharpness
	float Ni;            // index of refraction / optical density

	std::string map_Ka;
	std::string map_Kd;
	std::string map_Ks;
	std::string map_Ke;
	std::string map_d;
	std::string map_Ns;

	std::string bump;
	std::string disp;
	std::string decal;
	std::string refl;

	// [TODO] texture materials


	OBJMaterial()
	{
		Type = KnownMaterialTypes::OBJ_MTL_Format;
		id = -1;
		name = "///INVALID_NAME";
		Ka = Kd = Ks = Ke = Tf = InvalidColor();
		illum = -1;
		d = Ns = sharpness = Ni = Invalidf;
	}

	virtual Vector3f DiffuseColor() {
		return (Kd == InvalidColor()) ? Vector3f(1, 1, 1) : Kd;;
	}
	virtual float Alpha() {
		return (d == Invalidf) ? 1.0f : d;
	}
};





}