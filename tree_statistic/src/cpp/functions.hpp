#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

template<typename T> std::string toString(T v)
{
	std::stringstream ss;
	ss << v;
	return ss.str();
}

template<typename T> std::string toString(const std::vector<T>& v)
{
	std::stringstream ss;
	typename std::vector<T>::const_iterator it = v.begin();
	ss << (*it);
	it++;
	while(it != v.end())
	{
		ss << "," << *it;
		it++;
	}
	return ss.str();
}

#endif
