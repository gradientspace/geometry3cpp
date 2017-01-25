#include <geometry3PCH.h>

#include <g3types.h>

std::ostream& operator<<(std::ostream& os, const g3::Vector2f & v)
{
	os << v[0] << "," << v[1];
	return os;
}
std::ostream& operator<<(std::ostream& os, const g3::Vector2d & v)
{
	os << v[0] << "," << v[1];
	return os;
}

std::ostream& operator<<(std::ostream& os, const g3::Vector3f & v)
{
	os << v[0] << "," << v[1] << "," << v[2];
	return os;
}
std::ostream& operator<<(std::ostream& os, const g3::Vector3d & v)
{
	os << v[0] << "," << v[1] << "," << v[2];
	return os;
}


std::ostream& operator<<(std::ostream& os, const g3::Vector4f & v)
{
	os << v[0] << "," << v[1] << "," << v[2] << "," << v[3];
	return os;
}
std::ostream& operator<<(std::ostream& os, const g3::Vector4d & v)
{
	os << v[0] << "," << v[1] << "," << v[2] << "," << v[3];
	return os;
}


std::ostream& operator<<(std::ostream& os, const g3::Vector2i & v)
{
	os << v[0] << "," << v[1];
	return os;
}
std::ostream& operator<<(std::ostream& os, const g3::Vector3i & v)
{
	os << v[0] << "," << v[1] << "," << v[2];
	return os;
}
std::ostream& operator<<(std::ostream& os, const g3::Vector4i & v)
{
	os << v[0] << "," << v[1] << "," << v[2] << "," << v[3];
	return os;
}


std::ostream& operator<<(std::ostream& os, const g3::AxisAlignedBox2f & b)
{
	os << b.Min[0] << "-" << b.Max[0] << "," << b.Min[1] << "-" << b.Max[1];
	return os;
}
std::ostream& operator<<(std::ostream& os, const g3::AxisAlignedBox2d & b)
{
	os << b.Min[0] << "-" << b.Max[0] << "," << b.Min[1] << "-" << b.Max[1];
	return os;
}


std::ostream& operator<<( std::ostream& os, const g3::AxisAlignedBox3f & b )
{
	os << b.Min[0] << "-" << b.Max[0] << "," << b.Min[1] << "-" << b.Max[1] << "," << b.Min[2] << "-" << b.Max[2]  ;
	return os;
}
std::ostream& operator<<( std::ostream& os, const g3::AxisAlignedBox3d & b )
{
	os << b.Min[0] << "-" << b.Max[0] << "," << b.Min[1] << "-" << b.Max[1] << "," << b.Min[2] << "-" << b.Max[2]  ;
	return os;
}
