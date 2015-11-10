#ifndef CORE_HPP
#define CORE_HPP


/**
* @brief Different link functions
*/
enum {
  IDENTITY,
  LOG,
  LOGIT
};

/**
* @brief Convert a value to a string
*/
template<typename T> std::string toString(T v) {
  std::stringstream ss;

  ss << v;

  return ss.str();
}

template<typename T> std::string toString(const std::vector<T>& v)  {
  unsigned int i;
  std::string res="";
  for(i = 0; i < v.size()-1; i++) {
	res += toString(v[i]) +",";
  }
  res += toString(v[i]);
  return res;
}

/**
* @ brief Convert an array of values in an array of std::string
*/
template<typename T> std::vector<std::string> InString(const std::vector<T>& v)  {
  unsigned int i;
  std::vector<std::string> res;
  for(i = 0; i < v.size(); i++) {
	res.push_back(toString(v[i]));
  }
  return res;
}

/**
* @ brief Convert a matrix of values in a matrix of std::string
*/
template<typename T> std::vector< std::vector<std::string> > InString(const std::vector< std::vector<T> >& M)  {
  unsigned int i,j;
  std::vector< std::vector<std::string> > res(M.size());
  for(i = 0; i < M.size(); i++) {
	res[i].clear();
	for(j = 0; j < M[i].size(); j++) {
	  res[i].push_back(toString(M[i][j]));
	}
  }
  return res;
}

/**
* @brief Transpose a std::vector< std::vector<T> >
*/
template<typename T> std::vector< std::vector<T> > t(const std::vector< std::vector<T> >& M) {
  unsigned int i,j;
  std::vector< std::vector<T> > P(M[0].size(), std::vector<T>(M.size()));

  for(i = 0; i < M.size(); i++) {
	for(j = 0; j < M[i].size(); j++) {
	  P[j][i] = M[i][j];
	}
  }


  return P;
}

/**
* @brief Gerate arandom number between 0 and 1
*/
double generator() {
  return (double)rand()/((double)RAND_MAX + 1);
}

/**
* @brief Display vector in a stream
*/
template<typename T> std::ostream& display(const std::vector<T>& V, const std::vector<std::string>& labels, std::ostream& os) {
  unsigned int i;

  if(labels.empty()) {
	for(i = 0; i < V.size(); i++) {
	  os << toString(V[i]) << " ";
	}
	os << std::endl;
  } else {
	std::vector< std::vector<T> > M;
	M.push_back(V);
	display(M, labels, std::vector<std::string>(), os);
  }

  return os;
}

/**
* @brief Display a matrix in a stream
*/
template<typename T> std::ostream& display(const std::vector< std::vector<T> >& M, std::ostream& os) {
  unsigned int i,j,k;
  std::vector<unsigned int> max(M[0].size(),0);
  std::vector< std::vector<std::string> > S(M.size(), std::vector<std::string>(M[0].size()));

  for(i = 0; i < M.size(); i++) {
	for(j = 0; j < M[i].size(); j++) {
	  S[i][j] = toString(M[i][j]);
	  max[j] = MAX(max[j], S[i][j].size());
	}
  }
  for(i = 0; i < S.size(); i++) {
	for(j = 0; j < S[i].size(); j++) {
	  for(k = 0; k < (max[j]-S[i][j].size()); k++) {
		os << " ";
	  }
	  os << S[i][j] << " ";
	}
	os << std::endl;
  }

  return os;
}

/**
* @brief Display a matrix in a stream
*/
template<typename T> std::ostream& display(const std::vector< std::vector<T> >& M, const std::vector<std::string>& clabels, const std::vector<std::string>& rlabels, std::ostream& os) {
  unsigned int i,j,k;
  std::vector<unsigned int> max(M[0].size()+1,0);
  std::vector< std::vector<std::string> > S(M.size()+1, std::vector<std::string>(M[0].size()+1));

  S[0][0] = " ";
  for(i = 0; i < clabels.size(); i++) {
	S[0][i+1] = clabels[i];
	max[i+1] = MAX(max[i+1], S[0][i+1].size());
  }
  if(rlabels.size() > 0) {
	for(i = 0; i < rlabels.size(); i++) {
	  S[i+1][0] = rlabels[i];
	  max[0] = MAX(max[0], S[i+1][0].size());
	}
  }
  for(i = 0; i < M.size(); i++) {
	for(j = 0; j < M[i].size(); j++) {
	  S[i+1][j+1] = toString(M[i][j]);
	  max[j+1] = MAX(max[j+1], S[i+1][j+1].size());
	}
  }
  for(i = 0; i < S.size(); i++) {
	for(j = 0; j < S[i].size(); j++) {
	  if(max[j] > 0) {
		for(k = 0; k < (max[j]-S[i][j].size()); k++){
		  os << " ";
		}
		os << S[i][j] << " ";
	  }
	}
	os << std::endl;
  }

  return os;
}

/**
* @brief Split a string using a separator
*/
std::vector<std::string> split_string(const std::string& line, const std::string& separator) {
  boost::tokenizer< boost::char_separator<char> >::const_iterator i;
  std::vector<std::string> splitted_line;
  boost::char_separator<char> sep(separator.c_str());
  boost::tokenizer< boost::char_separator<char> > token(line, sep);

  for (i = token.begin(); i != token.end(); i++){
	splitted_line.push_back(*i);
  }

  return splitted_line;
}

template<typename T> std::vector< std::set<T> > compute_categories(const std::vector< std::vector<T> >& values) {
  unsigned int i,j;
  std::vector< std::set<T> > icategories(values.size());
  typename std::set<T>::iterator it;

  for(i = 0; i < values.size(); i++) {
	for(j = 0; j < values[i].size(); j++) {
	  it = icategories[i].find(values[i][j]);
	  if(it == icategories[i].end()) {
		icategories[i].insert(values[i][j]);
	  }
	}
  }

  return icategories;
}

template<typename T> std::vector< std::vector<T> > compute_combinations(const std::vector< std::set<T> >& icategories) {
  unsigned int i;
  std::vector< std::vector<T> > icombinations;
  std::vector< std::set<T> > iicategories;
  typename std::set<T>::iterator it1;
  typename std::vector< std::vector<T> >::iterator it2;

  if(icategories.size() > 1) {
	iicategories = icategories;
	iicategories.erase(iicategories.begin());
	std::vector< std::vector<T> > iicombinations = compute_combinations(iicategories);
	it1 = icategories[0].begin();
	while(it1 != icategories[0].end()) {
	  for(it2 = iicombinations.begin(); it2 != iicombinations.end(); it2++) {
		icombinations.push_back(*it2);
		icombinations.back().insert(icombinations.back().begin(), *it1);
	  }
	  it1++;
	}
  } else {
	for(it1 = icategories[0].begin() ; it1 != icategories[0].end(); it1++) {
	  icombinations.push_back(std::vector<T>(1,*it1));
	}
  }

  return icombinations;
}

std::vector< std::vector<unsigned int> > compute_powerset(const std::vector< std::vector<unsigned int> >& set, unsigned int invariables) {
  unsigned int i,j;
  std::vector< std::vector<unsigned int> > ipowerset;

  for(i = 0; i < set.size(); i++){
	if(set[i].size() == 1) {
	  ipowerset.push_back(set[i]);
	}
	for(j = set[i].back()+1; j < invariables; j++){
	  ipowerset.push_back(set[i]);
	  ipowerset.back().push_back(j);
	}
  }
  if(set.size() != pow(2,invariables)-1) {
	ipowerset = compute_powerset(ipowerset, invariables);
  }
  return ipowerset;
}

std::vector<double> Eigen2std(const Eigen::VectorXd& v) {
	unsigned int i;
	std::vector<double> res;
	for(i = 0; i < v.size(); i++) {
		res.push_back(v(i));
	}
	return res;
}

std::vector< std::vector<double> > Eigen2std(const Eigen::MatrixXd& M) {
	unsigned int i,j;
	std::vector< std::vector<double> > res(M.rows());
	for(i = 0; i < M.rows(); i++) {
		res[i].resize(M.cols());
		for(j = 0; j < M.cols(); j++) {
			res[i][j] = M(i,j);
		}
	}
	return res;
}

template<typename T> Eigen::VectorXd std2Eigen(const std::vector<T>& v) {
	unsigned int i;
	Eigen::VectorXd res(v.size());
	for(i = 0; i < v.size(); i++) {
		res(i) = v[i];
	}
	return res;
}

template<typename T> Eigen::MatrixXd std2Eigen(const std::vector< std::vector<T> >& M) {
	unsigned int i,j;
	Eigen::MatrixXd res(M.size(), M[0].size());
	for(i = 0; i < M.size(); i++) {
		for(j = 0; j < M[i].size(); j++) {
			res(i,j) = M[i][j];
		}
	}
	return res;
}

template<typename T> std::vector< std::vector<T> > v2M(const std::vector<T>& v) {
	unsigned int i,j,loc=0;
	unsigned int size = (unsigned int)((sqrt(v.size()*8+1)-1)/2.);
	std::vector< std::vector<T> > M(size, std::vector<T>(size,0));
	for(i = 0; i < size; i++) {
		M[i][i] = v[loc];
		loc++;
		for(j = i+1; j < size; j++) {
			M[i][j] = v[loc];
			M[j][i] = v[loc];
			loc++;
		}
	}
	return M;
}

#endif
