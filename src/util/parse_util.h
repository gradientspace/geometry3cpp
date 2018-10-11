#pragma once

#include <vector>
#include <string>

namespace g3
{

class ParseUtil
{
public:
	ParseUtil() = delete;

	enum class Options {
		RemoveEmpty = 1
	};

	static int ToInt(const std::string & s) {
		return std::stoi(s, 0);
	}
	static float ToFloat(const std::string & s) {
		return std::stof(s, 0);
	}
	static double ToDouble(const std::string & s){
		return std::stod(s, 0);
	}

	static bool Contains(const std::string & s, char c)
	{
		return s.find(c, 0) != std::string::npos;
	}
	static bool Contains(const std::string & s, const std::string & sub)
	{
		return s.find(sub, 0) != std::string::npos;
	}


	static void Split(const std::string& s, char delimiter, std::vector<std::string> & tokens, Options options)
	{
		bool skip_empty = (((int)options & (int)Options::RemoveEmpty) != 0);

		tokens.clear();
		auto cur = s.begin(), end = s.end();
		auto cur_start = cur;
		cur++;
		while (cur != end) {
			if (*cur == delimiter) {
				bool is_empty = (cur_start == cur);
				if (is_empty == false || skip_empty == false) {
					tokens.push_back(std::string(cur_start, cur));
				}
				cur++;  // advance to next char
				cur_start = cur;
			} else {
				cur++;
			}
		}
		if (cur_start != cur)
			tokens.push_back(std::string(cur_start, end) );
	}

	static std::vector<std::string> Split(const std::string& s, char delimiter, Options options)
	{
		std::vector<std::string> tokens;
		std::string token;
		std::istringstream tokenStream(s);
		while (std::getline(tokenStream, token, delimiter))
		{
			if (((int)options & (int)Options::RemoveEmpty) != 0 && token.size() == 0)
				continue;
			tokens.push_back(token);
		}
		return tokens;
	}
};


}