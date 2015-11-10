#ifndef CONTINGENCY_TABLE_HPP
#define CONTINGENCY_TABLE_HPP

template<typename T> ContingencyTable<T>::ContingencyTable() {
}

template<typename T> ContingencyTable<T>::ContingencyTable(const std::vector< std::vector<T> >& values) {
  try {
	unsigned int i,j;
	std::vector< std::vector<T> > tvalues(values[0].size(), std::vector<T>(values.size()));
	for(i = 0; i < values.size(); i++) {
		for(j = 0; j < values[i].size(); j++) {
			tvalues[j][i] = values[i][j];
		}
	}
	nvariables = tvalues.size();
	if(nvariables == 0) {
	  throw std::string("Error: no variables in input !");
	}
	nobs = tvalues[0].size();
	for(i = 1; i < tvalues.size(); i++) {
	  if(tvalues[i].size() != nobs)
		throw std::string("Error: variable "+ toString(i) + " have "+ toString(tvalues[i].size()) + " observations but expecting "+ toString(nobs) +"!");
	}
	std::vector< std::set<T> > categories = compute_categories(tvalues);
	parameters_names = compute_parametersnames(categories);
	nparameters = parameters_names.size()+1;
	std::vector< std::vector<T> > dvalues = tvalues;
	i = 0;
	while(i < categories.size()) {
	  if(categories[i].size() == 1) {
		std::cerr << "Warning: variable "+ toString(i) + " will not be considered because there is only one category (considered as a determinist variable)!" << std::endl;
		categories.erase(categories.begin()+i);
		dvalues.erase(dvalues.begin()+i);
	  } else {
		i++;
	  }
	}
	glm = GLM_init(categories,dvalues);
	varcovar = Get_covariances(tvalues);
  } catch(const std::string& error) {
		std::cerr << error << std::endl;
  }
}

template<typename T> ContingencyTable<T>::ContingencyTable(const ContingencyTable<T>& ct) {
  nobs = ct.nobs;
  ncases = ct.ncases;
  nvariables = ct.nvariables;
  nparameters = ct.nparameters;
  parameters_names = ct.parameters_names;
  glm = ct.glm;
	varcovar = ct.varcovar;
}

template<typename T> ContingencyTable<T>::~ContingencyTable() {
  parameters_names.clear();
  //delete glm;
}

template<typename T> std::vector< std::set<T> > ContingencyTable<T>::compute_categories(const std::vector< std::vector<T> >& values) const{
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

template<typename T> std::vector< std::vector<T> > ContingencyTable<T>::compute_combinations(const std::vector< std::set<T> >& icategories) const {
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

template<typename T> std::vector< std::vector<unsigned int> > ContingencyTable<T>::compute_powerset(const std::vector< std::vector<unsigned int> >& set, unsigned int invariables) const {
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

template<typename T> std::vector<double> ContingencyTable<T>::Get_means(const std::vector< std::vector<T> >& values) const {
  unsigned int i,j;
  std::vector<double> means(nvariables, 0);

  for(i = 0; i < nobs; i++){
	for(j = 0; j < nvariables; j++){
	  means[j] += values[j][i]/((double)nobs);
	}
  }

  return means;
}

template<typename T> std::vector< std::vector<double> > ContingencyTable<T>::Get_covariances(const std::vector< std::vector<T> >& values) const {
  unsigned int i,j,k;
  std::vector<double> means = Get_means(values);
  std::vector< std::vector<double> > varcovari(nvariables, std::vector<double>(nvariables,0));

  for(i = 0; i < nobs; i++){
	for(j = 0; j < nvariables; j++){
	  varcovari[j][j] += pow(values[j][i]-means[j],2)/((double)nobs-1);
	  for(k = j+1; k < nvariables; k++){
		varcovari[j][k] += values[j][i]*values[k][i]/((double)nobs);
	  }
	}
  }

  for(i = 0; i < nvariables-1; i++){
	for(j = i+1; j < nvariables; j++){
	  varcovari[i][j] -= means[i]*means[j];
	  varcovari[j][i] = varcovari[i][j];
	}
  }

  return varcovari;
}

template<typename T> std::vector<std::string> ContingencyTable<T>::compute_parametersnames(const std::vector< std::set<T> >& icategories) const {
  unsigned int i,j,k;
  bool determinist;
  std::string name;
  std::vector<std::string> names;
  std::vector< std::vector<unsigned int> > set;
  for(i = 0; i < nvariables; i++) {
	set.push_back(std::vector<unsigned int>(1,i));
  }
  std::vector< std::vector<unsigned int> > powerset = compute_powerset(set, nvariables);
  for(i = 0; i < powerset.size(); i++) {
	determinist = false;
	name.clear();
	j = 0;
	k = 1;
	while(j < powerset[i].size()-1 && !determinist) {
	  if(icategories[powerset[i][j]].size() == 1)
		determinist = true;
	  else
		name += toString(powerset[i][j]) +"--";
		k *= icategories[powerset[i][j]].size()-1;
		j++;
	}
	if(icategories[powerset[i][j]].size() == 1)
	  determinist = true;
	else {
	  name += toString(powerset[i][j]);
	  k *= icategories[powerset[i][j]].size()-1;
	}
	if(!determinist)
	  for(j = 0; j < k; j++) {
		names.push_back(name);
	  }
  }
  return names;
}

template<typename T> void ContingencyTable<T>::Estimate() {
  try {
	std::string comment="";
	if(!SILENT)
	  std::cout << "# Estimation using GLM procedures..." << std::endl;
	glm->Estimate(comment);
	if(comment.size() != 0) {
	  throw comment;
	}
  } catch(const std::string& error) {
	std::cerr << error << std::endl;
  }
}

template<typename T> ContingencyTable<T>* ContingencyTable<T>::Drop_parameter(const std::string& parameter) const {
  ContingencyTable<T>* nCT = new ContingencyTable<T>(*this);
  std::vector<std::string> sparameter = split_string(parameter, "-");
  bool allfinded,end;
  unsigned int i,j,pos;
  std::vector<bool> todel;
  for(i = 0; i < nCT->parameters_names.size(); i++) {
	allfinded = true;
	end = false;
	j = 0;
	while(allfinded && !end) {
	  pos = nCT->parameters_names[i].find(sparameter[j]);
	  if(pos >= nCT->parameters_names[i].size())
		allfinded = false;
	  j++;
	  if(j == sparameter.size())
		end = true;
	}
	if(allfinded)
	  todel.push_back(true);
	else
	  todel.push_back(false);
  }
  for(i = nCT->parameters_names.size(); i > 1; i--){
	if(todel[i-1])
	  nCT->parameters_names.erase(nCT->parameters_names.begin()+i-1);
  }
  todel.insert(todel.begin(), false);
  nCT->nparameters = nCT->parameters_names.size()+1;
  nCT->glm = glm->Drop_parameter(todel);
  return nCT;
}

template<typename T> GLM* ContingencyTable<T>::GLM_init(const std::vector< std::set<T> >& icategories, std::vector< std::vector<T> > dvalues) {
  unsigned int i;
  std::vector< std::vector<T> > combinations = compute_combinations(icategories);
  ncases = combinations.size();
	if(ncases < 1000) {
	  i = 0;
		typename std::set<T>::iterator it;
		std::vector< std::set<T> > tcategories = icategories;
		for(i = 0; i < tcategories.size(); i++) {
			it = tcategories[i].end();
			it--;
			tcategories[i].erase(it);
		}
		Eigen::VectorXd y(ncases);
		Eigen::MatrixXd X(ncases, nparameters);
		std::vector< std::vector<unsigned int> > set;
		for(i = 0; i < icategories.size(); i++) {
			set.push_back(std::vector<unsigned int>(1,i));
		}
		std::vector< std::vector<unsigned int> > powerset = compute_powerset(set, icategories.size());
		std::vector< std::set<T> > varscat;
		std::vector< std::vector<T> > tcombinations, valuescombinations, variablescombinations;
		std::vector<T> ivalue;
		unsigned int j;
		typename std::vector< std::vector<T> >::iterator itc;
		for(i = 0; i < powerset.size(); i++) {
			varscat.clear();
			for(j = 0; j < powerset[i].size(); j++) {
				varscat.push_back(tcategories[powerset[i][j]]);
			}
			tcombinations = compute_combinations(varscat);
			for(itc = tcombinations.begin(); itc != tcombinations.end(); itc++) {
				valuescombinations.push_back(*itc);
				variablescombinations.push_back(powerset[i]);
			}
		}
		unsigned int c=0, rep;
		typename std::vector< std::vector<T> >::iterator itv;
		dvalues = t(dvalues);

		for(itc = combinations.begin(); itc != combinations.end(); itc++) {
			rep = 0;
			itv = dvalues.begin();
			while(itv != dvalues.end()) {
				if(std::equal((*itv).begin(), (*itv).end(), (*itc).begin())) {
					dvalues.erase(itv);
					rep++;
				} else {
					itv++;
				}
			}
			y(c) = rep;
			X(c,0) = 1;
			for(i = 0; i < valuescombinations.size(); i++) {
				ivalue.clear();
				for(j = 0; j < variablescombinations[i].size(); j++) {
					ivalue.push_back((*itc)[variablescombinations[i][j]]);
				}
				if(std::equal(valuescombinations[i].begin(), valuescombinations[i].end(), ivalue.begin())) {
					X(c,i+1) = 1;
				} else {
					X(c,i+1) = 0;
				}
			}
			c++;
		}
		GLM* iglm = new GLM(y,X,stat_tool::POISSON,LOG);
		return iglm;
	} else {
		return NULL;
	}
}

template<typename T> std::ostream& ContingencyTable<T>::Summary(std::ostream& os) const {
  os << "# Number of observations:" << std::endl << nobs << std::endl;
  os << "# Number of cases:" << std::endl << ncases << std::endl;

  os << "# Number of parameters:" << std::endl << nparameters << std::endl;
	os << "# Variance-covariance matrix:" << std::endl;
	display(varcovar);
  os << "# Associated edges:" << std::endl;
  std::set<std::string> edges = extract_edges();
	std::set<std::string>::iterator it;
	for(it = edges.begin(); it != edges.end(); it++) {
		os << *it << std::endl;
  }
	os << "# Associated decomposition:" << std::endl;
	edges = extract_cliques();
	edges = get_decomposition(edges);
	for(it = edges.begin(); it != edges.end(); it++) {
		os << *it << std::endl;
  }
}

/*template<typename T> Gnuplot<double, double, double>* ContingencyTable<T>::Plot_residuals() const {
  Gnuplot<double, double, double> *G = new Gnuplot<double, double, double>();
  std::vector<std::string> with(1,"points");
  G->set_with(with);
  std::vector<std::string> labels;
  labels.push_back("Predicted");
  labels.push_back("Residuals");
  G->set_label(labels);
  G->set_data(glm->Get_residuals());
  G->plot();
  return G;
}*/

template<typename T> std::set<std::string> ContingencyTable<T>::extract_edges() const {
  unsigned int i;
  std::set<std::string> edges;
  std::set<std::string>::iterator it;
  for(i = 0; i < parameters_names.size(); i++) {
		it = edges.find(parameters_names[i]);
		if(it == edges.end() && split_string(parameters_names[i],"-").size() == 2) {
			edges.insert(parameters_names[i]);
		}
  }
  return edges;
}

template<typename T> std::set<std::string> ContingencyTable<T>::extract_cliques() const {
  unsigned int i;
  std::set<std::string> edges;
  std::set<std::string>::iterator it;
  for(i = 0; i < parameters_names.size(); i++) {
		it = edges.find(parameters_names[i]);
		if(it == edges.end()) {
			if(samecovsign(split_string(parameters_names[i],"-"))) {
				edges.insert(parameters_names[i]);
			}
		}
  }
	edges = get_decomposition(edges);
  return edges;
}

template<typename T> std::set<std::string> ContingencyTable<T>::get_decomposition(const std::set<std::string>& iedges) const {
	std::set<std::string>::iterator it, it2, itr;
	std::set<std::string> edges = iedges;
	unsigned int i;
	bool loop=true,allfinded,end
	;
	it = edges.begin();
	while(it != edges.end()) {
		itr = edges.begin();
		loop = true;
		while(loop) {
			if(split_string(*itr,"-").size() > split_string(*it,"-").size()) {
				allfinded = true;
				end = false;
				std::vector<std::string> vars = split_string(*it,"-");
				i = 0;
				while(allfinded && !end) {
					allfinded *= ((*itr).find(vars[i])<(*itr).size());
					i++;
					if(i == vars.size())
						end = true;
				}
				if(allfinded) {
					edges.erase(it);
					it++;
					loop = false;
				} else {
					itr++;
				}
			} else {
				itr++;
			}
			if(itr == edges.end()) {
				loop = false;
				it++;
			}
		}
	}

	return edges;
}

template<typename T> unsigned int ContingencyTable<T>::toInt(const std::string& s) const {
	std::stringstream strStream(s);
	unsigned int res;
	strStream >> res;
	return res;
}
template<typename T> bool ContingencyTable<T>::samecovsign(const std::vector<std::string>& vars) const {
	if(vars.size() <=2)
		return true;
	else {
		unsigned int i,j,s=0;
		for(i = 0; i < vars.size()-1; i++) {
			for(j = i+1; j < vars.size(); j++) {
				s += (int)(varcovar[toInt(vars[i])][toInt(vars[j])]>0);
			}
		}
		if(s == 0 || s == (vars.size()*(vars.size()-1))/2.)
			return true;
		else
			return false;
	}
}

template<typename T> ContingencyTable<T>* ContingencyTable<T>::Parsimoniest() const {
  ContingencyTable<T>* bestCT = new ContingencyTable<T>(*this);
  ContingencyTable<T> *CTtemp = NULL, *bestCTtemp = NULL;
  std::set<std::string> edges = bestCT->extract_edges();
  std::set<std::string>::iterator it;
  bool change=true;
  std::string comment="",bedge;
	if(ncases < 1000) {
		if(!SILENT)
			std::cout << "# Step: estimation of the less parsimoniest model..." << std::endl;
		bestCT->glm->Estimate(comment);
		if(comment.size() > 0) {
		if(!SILENT)
			std::cout << comment << std::endl;
		change = false;
		}
		while(change) {
			change = false;
			for(it = edges.begin(); it != edges.end(); it++) {
				if(!SILENT)
					std::cout << "# Step: removing edge (" << *it << ")..." << std::endl;
				comment = "";
				CTtemp = new ContingencyTable<T>(*bestCT->Drop_parameter(*it));
				CTtemp->glm->Estimate(comment);
				if(comment.size() == 0) {
					if(!SILENT) {
						if(CTtemp->accept_against(*bestCT))
							std::cout << "Could be removed !" << std::endl;
						else
							std::cout << "Could not be removed !" << std::endl;
					}
					if(bestCTtemp == NULL && CTtemp->accept_against(*bestCT)) {
						change = true;
						bedge = *it;
						bestCTtemp = new ContingencyTable<T>(*CTtemp);
					} else {
						if(bestCTtemp != NULL && CTtemp->accept_against(*bestCT) && (CTtemp->Get_loglikelihood() > bestCTtemp->Get_loglikelihood())) {
							bestCTtemp = new ContingencyTable<T>(*CTtemp);
							bedge = *it;
						}
					}
				} else {
					if(!SILENT)
						std::cout << comment << std::endl;
				}
				delete CTtemp;
			}
			if(change) {
				if(!SILENT)
					std::cout << "# Update: dropped the edge ("+ bedge +")!" << std::endl;
				delete bestCT;
				bestCT = new ContingencyTable<T>(*bestCTtemp);
				delete bestCTtemp;
				bestCTtemp = NULL;
				edges = bestCT->extract_edges();
			}
		}
		if(!SILENT) {
			std::cout << "# End: the previous update was the last one selected !" << std::endl;
			bestCT->Summary();
		}
	} else {
		if(!SILENT) {
			std::cout << "# Error: no estimation performed (too many parameters) !" << std::endl;
			bestCT->Summary();
		}
	}

  return bestCT;
}

template<typename T> double ContingencyTable<T>::Get_loglikelihood() const {
  return glm->Get_loglikelihood();
}

template<typename T> double ContingencyTable<T>::accept_against(const ContingencyTable<T>& ct) const {
  bool accept;
  //std::cout << ct.Get_loglikelihood() << " " << Get_loglikelihood() << " " << ct.nparameters-nparameters << " " << quantile(complement(boost::math::chi_squared(ct.nparameters-nparameters), 0.05)) << std::endl;
  double deviance = 2*(ct.Get_loglikelihood() - Get_loglikelihood());

  return (deviance < (quantile(complement(boost::math::chi_squared(ct.nparameters-nparameters), 0.05))));
}
#endif
