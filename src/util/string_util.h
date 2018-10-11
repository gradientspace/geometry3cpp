#pragma once

#include <string>
#include <regex>
#include <algorithm>
#include <utility>
#include <tuple>

namespace g3 {

class StringUtil
{
public:
	StringUtil() = delete;

	std::string ToString(int i) {
		return std::to_string(i);
	}
	std::string ToString(float f) {
		return std::to_string(f);
	}
	std::string ToString(double d) {
		return std::to_string(d);
	}




	/*Generalized toStdString*/
	template <typename T>
	static std::string toStdString(const T &rhs)
	{
		std::stringstream stringStream{};
		stringStream << rhs;
		return stringStream.str();
	}

	/*Template specializations for toStdString*/
	static std::string toStdString(const char *rhs) { return std::string{ rhs }; }
	static std::string toStdString(const std::string &rhs) { return rhs; }
	static std::string toStdString(char *rhs) { return std::string{ rhs }; }
	static std::string toStdString(char rhs) { return std::string{ 1, rhs }; }

	/*Base case to break recursion*/
	static std::string Format(const char *formatting)
	{
		return std::string{ formatting };
	}


	/*C# style String.Format()*/
	template <typename First, typename ... Args>
	static std::string Format(const char *formatting, const First& first, const Args& ... args)
	{
		/* Match exactly one opening brace, one or more numeric digit,
		* then exactly one closing brace, identifying a token
		* Ex: {0} will match, {-1} will not */
		static const std::regex targetRegex{ "\\{[0-9]+\\}" };
		std::smatch match;

		/* Copy the formatting string to a std::string, to
		* make for easier processing, which will eventually
		* be used (the .c_str() method) to pass the remainder
		* of the formatting recursively */
		std::string returnString{ formatting };

		/* Copy the formatting string to another std::string, which
		* will get modified in the regex matching loop, to remove the
		* current match from the string and find the next match */
		std::string copyString{ formatting };

		/* std::tuple to hold the current smallest valued brace token,
		* wrapped in a std::vector because there can be multiple brace
		* tokens with the same value. For example, in the following format string:
		* "There were {0} books found matching the title {1}, {0}/{2}",
		* this pass will save the locations of the first and second {0} */
		using TokenInformation = std::tuple<int, size_t, size_t>;
		std::vector<TokenInformation> smallestValueInformation{ std::make_tuple(-1, 0, 0) };

		/*Iterate through string, finding position and lengths of all matches {x}*/
		while (std::regex_search(copyString, match, targetRegex)) {
			/*Get the absolute position of the match in the original return string*/
			size_t foundPosition{ match.position() + (returnString.length() - copyString.length()) };
			int regexMatchNumericValue{ 0 };
			try {
				/*Convert the integer value between the opening and closing braces to an int to compare */
				regexMatchNumericValue = std::stoi(returnString.substr(foundPosition + 1, (foundPosition + match.str().length())));
			}
			catch (std::exception &e) {
				/*The value between the braces was not an integer value, so throw an
				* exception. Is this check actually necessary? I am not familiar enough with
				* regex to know whether it is even possible to have a non-numeric value picked up */
				throw std::runtime_error(Format("ERROR: In Format() - Formatted string contains non-numeric token(formatting = {0}, ex = {1})", formatting, e.what()));
			}
			/*Do not allow negative numbers, although this should never get picked up the regex anyway*/
			if (regexMatchNumericValue < 0) {
				throw std::runtime_error(Format("ERROR: In Format() - Formatted string is invalid (formatting = {0})", formatting));
			}
			/* If the numeric value in the curly brace token is smaller than
			* the current smallest (or if the smallest value has not yet been set,
			* ie it is the first match), set the corresponding smallestX variables
			* and wrap them up into a TokenInformation and add it to the std::vector */
			int smallestValue{ std::get<0>(smallestValueInformation.at(0)) };
			if ((smallestValue == -1) || (regexMatchNumericValue < smallestValue)) {
				smallestValueInformation.clear();
				smallestValueInformation.push_back(std::make_tuple(regexMatchNumericValue,
					foundPosition,
					match.str().length()));
			}
			else if (regexMatchNumericValue == smallestValue) {
				smallestValueInformation.push_back(std::make_tuple(regexMatchNumericValue,
					foundPosition,
					match.str().length()));
			}
			/*Set the copy string to just past the match token, so we don't accidentally
			* match it again on the next iteration */
			copyString = match.suffix();
		}
		int smallestValue{ std::get<0>(smallestValueInformation.at(0)) };
		/*If the smallest value is still the initialized value of -1, no other value was
		*found. This is also checking whether not we found ANY tokens (the smallestValueInformation
		*will receive whatever the first match is by default) */
		if (smallestValue == -1) {
			throw std::runtime_error(Format("ERROR: In Format() - Formatted string is invalid (formatting = {0})", formatting));
		}
		/* Set the returnString to be up to the brace token, then the string
		* representation of current argument in line (first), then the remainder
		* of the format string, effectively removing the token and replacing it
		* with the requested item in the final string, then pass it off recursively */

		std::string firstString{ toStdString(first) };
		int index{ 0 };
		for (const auto &it : smallestValueInformation) {
			size_t smallestValueLength{ std::get<2>(it) };


			/* Since the original string will be modified, the adjusted position must be
			calculated for any repeated brace tokens, kept track of by index.
			The length of string representation of first mutiplied by which the iterationn count
			is added, and the length of the brace token multiplied by the iteration count is
			subtracted, resulting in the correct starting position of the current brace token */
			size_t lengthOfTokenBracesRemoved{ index * smallestValueLength };
			size_t lengthOfStringAdded{ index * firstString.length() };
			size_t smallestValueAdjustedPosition{ std::get<1>(it) + lengthOfStringAdded - lengthOfTokenBracesRemoved };
			returnString = returnString.substr(0, smallestValueAdjustedPosition)
				+ firstString
				+ returnString.substr(smallestValueAdjustedPosition + smallestValueLength);
			index++;
		}
		/*Call template recursively with newly replaced string*/
		return Format(returnString.c_str(), args...);
	}


};


}