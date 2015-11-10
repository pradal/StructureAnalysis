#ifndef DISCRETE_UNIVARIATE_CONDITIONAL_DISTRIBUTION_HPP
#define DISCRETE_UNIVARIATE_CONDITIONAL_DISTRIBUTION_HPP

template<typename T> DiscreteUnivariateConditionalReestimation<T>::DiscreteUnivariateConditionalReestimation()
{
	nb_covariables = 0;
	nb_elements = 0;
}


template<typename T> DiscreteUnivariateConditionalReestimation<T>::DiscreteUnivariateConditionalReestimation(const DiscreteUnivariateConditionalReestimation<T>& duce)
{
	conditional_histogram = duce.conditional_histogram;
	nb_covariables = duce.nb_covariables;
	nb_elements = duce.nb_elements;
	nb_different_elements = duce.nb_different_elements;
}

template<typename T> DiscreteUnivariateConditionalReestimation<T>::DiscreteUnivariateConditionalReestimation(unsigned int inb_covariables)
{
	nb_covariables = inb_covariables;
	nb_elements = 0;
	nb_different_elements = 0;
}
/*template<> DiscreteUnivariateConditionalReestimation<unsigned int>::DiscreteUnivariateConditionalReestimation(const Vectors& vectors, int response_pos)
{
	int response;
	nb_covariables = vectors.get_nb_variable()-1;
	std::vector<int> explanatory(nb_covariables,0);
	nb_elements = vectors.get_nb_vector();
	std::map< std::vector<int>, std::pair<unsigned int, std::map<int, unsigned int> > >::iterator it0;
	std::map<int, unsigned int> tmap;
	std::map<int, unsigned int>::iterator it1;
	for(unsigned int i = 0; i < nb_elements; ++i)
	{
		for(unsigned int j = 0; j < nb_covariables+1; ++j)
		{
			if(j == response_pos)
				response = vectors.get_int_vector(i,j);
			else
				explanatory[j-(int)(j > response_pos)] = vectors.get_int_vector(i, j);
		}
		it0 = conditional_histogram.find(explanatory);
		if(it0 == conditional_histogram.end())
		{
			tmap.clear();
			tmap.insert(std::pair<int, unsigned int>(response, 1));
			conditional_histogram.insert(std::pair<std::vector<int>, std::pair<unsigned int, std::map<int, unsigned int> > >(explanatory, std::pair<unsigned int, std::map<int, unsigned int> >(1, tmap)));
		} else {
			it1 = (((*it0).second).second).find(response);
			if(it1 == (((*it0).second).second).end())
			{
				(((*it0).second).second).insert(std::pair<int, unsigned int>(response, 1));
				(((*it0).second).first) += 1;
			} else {
				(*it1).second += 1;
				(((*it0).second).first) += 1;
			}
		}
	}
	nb_different_elements = 0;
	for(std::map<std::vector<int>, std::pair<unsigned int, std::map<int, unsigned int> > >::iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		nb_different_elements += (((*it).second).second).size();
	}
}*/

template<typename T> DiscreteUnivariateConditionalReestimation<T>::DiscreteUnivariateConditionalReestimation(const std::map<std::vector<int>, T>& histogram, const std::set<unsigned int>& covariables, unsigned int variable)
{
	nb_covariables = covariables.size();
	std::vector<int> covariates(covariables.size());
	typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::iterator ith;
	typename std::map<int, T>::iterator itmh;
	std::map<int, T> tmap;
	for(typename std::map<std::vector<int>, T>::const_iterator it = histogram.begin(); it != histogram.end(); ++it)
	{
		for(std::set<unsigned int>::const_iterator itc = covariables.begin(); itc != covariables.end(); ++itc)
			covariates[distance(covariables.begin(), itc)] = ((*it).first)[*itc];
		ith = conditional_histogram.find(covariates);
		if(ith == conditional_histogram.end())
		{
			tmap.clear();
			tmap.insert(std::pair<int, T>(((*it).first)[variable], (*it).second));
			conditional_histogram.insert(std::pair<std::vector<int>, std::pair<T, std::map<int, T> > >(covariates, std::pair<T, std::map<int, T> >((*it).second, tmap)));
		} else {
			itmh = (((*ith).second).second).find(((*it).first)[variable]);
			if(itmh == (((*ith).second).second).end())
				(((*ith).second).second).insert(std::pair<int, T>(((*it).first)[variable], (*it).second));
			else
				(*itmh).second += (*it).second;
			((*ith).second).first += (*it).second;
		}
	}
	nb_elements = 0;
	nb_different_elements = 0;
	for(typename std::map<std::vector<int>, std::pair<T, std::map<int, T> > >::iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		nb_elements += ((*it).second).first;
		nb_different_elements += (((*it).second).second).size();
	}
}

template<typename T> DiscreteUnivariateConditionalReestimation<T>::~DiscreteUnivariateConditionalReestimation()
{
}

template<typename T> void DiscreteUnivariateConditionalReestimation<T>::update(const std::vector<int>& covariates, const std::map<int, T>& marginal_histogram)
{
	assert(covariates.size() == nb_covariables);
	typename std::map<std::vector<int>, std::pair<T, std::map<int, T> > >::iterator itcm;
	typename std::map<int, T>::iterator itm;
	itcm = conditional_histogram.find(covariates);
	if(itcm == conditional_histogram.end())
	{
		conditional_histogram.insert(std::pair<std::vector<int>, std::pair<T, std::map<int, T> > >(covariates, std::pair<T, std::map<int, T> >(0, marginal_histogram)));
		itcm = conditional_histogram.find(covariates);
		for(typename std::map<int, T>::iterator it = (((*itcm).second).second).begin(); it != (((*itcm).second).second).end(); ++it)
			((*itcm).second).first += (*it).second;
	} else {
		for(typename std::map<int, T>::const_iterator it = marginal_histogram.begin(); it != marginal_histogram.end(); ++it)
		{
			itm = (((*itcm).second).second).find((*it).first);
			if(itm == (((*itcm).second).second).end())
				(((*itcm).second).second).insert(*it);
			else
				(*itm).second += (*it).second;
			((*itcm).second).first += (*it).second;
		}
	}
}

template<typename T> void DiscreteUnivariateConditionalReestimation<T>::nb_value_computation()
{
	nb_elements = 0;
	nb_different_elements = 0;
	for(typename std::map<std::vector<int>, std::pair<T, std::map<int, T> > >::iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		nb_elements += ((*it).second).first;
		nb_different_elements += (((*it).second).second).size();
	}
}

template<typename T> std::ostream& DiscreteUnivariateConditionalReestimation<T>::display(std::ostream& os) const
{
	std::vector< std::vector<std::string> > Output(1, std::vector<std::string>(3,""));
	Output[0][0] = "Covariates";
	Output[0][1] = "Response";
	Output[0][2] = "Occurancy";
	unsigned int p = 1;
	bool not_setted;
	for(typename std::map< std::vector<int>, std::pair<unsigned int, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		not_setted = true;
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
		{
			Output.push_back(std::vector<std::string>(3,""));
			Output[p][1] = toString((*itb).first);
			Output[p][2] = toString((*itb).second);
			if(distance((((*it).second).second).begin(), itb) >= (int)(distance((((*it).second).second).begin(), (((*it).second).second).end())/2.) && not_setted)
			{
				Output[p][0] = toString((*it).first);
				not_setted = false;
			}
			p++;
		}
	}
	for(std::vector< std::vector<std::string> >::iterator it = Output.begin(); it != Output.end(); ++it)
	{
		for(std::vector<std::string>::iterator itb = (*it).begin(); itb != (*it).end(); ++itb)
			os << " " << *itb << " ";
		os << std::endl;
	}
	return os;
}

template<typename T> double DiscreteUnivariateConditionalReestimation<T>::entropy_computation() const
{
	double e = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
			e -= (*itb).second * (log((*itb).second) - log(((*it).second).first));
	}
	return e;
}

template<typename T> DiscreteUnivariateConditionalDistribution* DiscreteUnivariateConditionalReestimation<T>::distribution_estimation() const 
{
	std::map< std::vector<int>, std::map<int, double> > cmass;
	std::map<int, double> mass;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		mass.clear();
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
			mass.insert(std::pair<int, double>((*itb).first, (double)((*itb).second)));
		cmass.insert(std::pair< std::vector<int>, std::map<int, double> >((*it).first, mass));
	}
	return new DiscreteUnivariateConditionalDistribution(nb_covariables, cmass);
}

template<typename T> DiscreteUnivariateConditionalParametric* DiscreteUnivariateConditionalReestimation<T>::type_parametric_estimation(unsigned int maxits) const
{
	DiscreteUnivariateConditionalParametric *best = NULL, *current = NULL;
	best = new DiscreteUnivariateConditionalParametric(*(parametric_estimation(stat_tool::BINOMIAL, DiscreteUnivariateConditionalParametric::LOGIT, maxits)));
	current = new DiscreteUnivariateConditionalParametric(*(parametric_estimation(stat_tool::BINOMIAL, DiscreteUnivariateConditionalParametric::LOGLOG, maxits)));
	if(likelihood_computation(*best) < likelihood_computation(*current))
	{
		delete best;
		best = new DiscreteUnivariateConditionalParametric(*current);
	}
	delete current;
	current = new DiscreteUnivariateConditionalParametric(*(parametric_estimation(stat_tool::POISSON, DiscreteUnivariateConditionalParametric::LOG, maxits)));
	if(likelihood_computation(*best) < likelihood_computation(*current))
	{
		delete best;
		best = new DiscreteUnivariateConditionalParametric(*current);
	}
	delete current;
	current = new DiscreteUnivariateConditionalParametric(*(parametric_estimation(stat_tool::POISSON, DiscreteUnivariateConditionalParametric::IDENTITY, maxits)));
	if(likelihood_computation(*best) < likelihood_computation(*current))
	{
		delete best;
		best = new DiscreteUnivariateConditionalParametric(*current);
	}
	delete current;
	return best;
}

template<typename T> DiscreteUnivariateConditionalParametric* DiscreteUnivariateConditionalReestimation<T>::parametric_estimation(int family, int link, unsigned int maxits) const
{
	Eigen::VectorXd y, z, beta;
	Eigen::MatrixXd X, W;
	y = y_init();
	if(family == stat_tool::BINOMIAL)
	{
		bool success = make_binomial(y);
		if(!success)
			return new DiscreteUnivariateConditionalParametric();
	}
	X = X_init();
	W = W_init(y);
	switch(link)
	{
		case DiscreteUnivariateConditionalParametric::IDENTITY :
			z = identity_link_z_init(y);
			break;
		case DiscreteUnivariateConditionalParametric::LOG :
			z = log_link_z_init(y);
			break;
		case DiscreteUnivariateConditionalParametric::LOGIT :
			z = logit_link_z_init(y);
			break;
		case DiscreteUnivariateConditionalParametric::LOGLOG :
			z = loglog_link_z_init(y);
			break;
	}
	beta = beta_compute(X, W, z);
	DiscreteUnivariateConditionalParametric *previous = NULL, *current = NULL;
	previous = new DiscreteUnivariateConditionalParametric(nb_covariables, family, link, to_stdvector(beta));
	double plh = likelihood_computation(*previous), clh;
	int iterations = 0;
	bool end = false;
	while(!end)
	{
		switch(link)
		{
			case DiscreteUnivariateConditionalParametric::IDENTITY :
				identity_link_z_update(y, z, X, beta);
				poisson_identity_link_W_update(W, X, beta);
				break;
			case DiscreteUnivariateConditionalParametric::LOG :
				log_link_z_update(y, z, X, beta);
				poisson_log_link_W_update(W, X, beta);
				break;
			case DiscreteUnivariateConditionalParametric::LOGIT :
				logit_link_z_update(y, z, X, beta);
				binomial_logit_link_W_update(W, X, beta);
				break;
			case DiscreteUnivariateConditionalParametric::LOGLOG :
				loglog_link_z_update(y, z, X, beta);
				binomial_loglog_link_W_update(W, X, beta);
				break;
		}
		beta = beta_compute(X, W, z);
		current = new DiscreteUnivariateConditionalParametric(nb_covariables, family, link, to_stdvector(beta));
		clh = likelihood_computation(*current);
		if((clh-plh)/abs(plh) < 10e-4)
			end = true;
		plh = clh;
		delete previous;
		previous = new DiscreteUnivariateConditionalParametric(*current);
		delete current;
		iterations++;
		if(iterations > maxits)
			end = true;
	}
	return previous;
}

template<typename T> double DiscreteUnivariateConditionalReestimation<T>::likelihood_computation(DiscreteUnivariateConditionalDistribution& ducd, bool log_computation) const
{
	double lh;
	if(log_computation)
		lh = 0;
	else
		lh = 1;
	std::vector<int> event;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		event = (*it).first;
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
		{
			if(log_computation)
				lh += (*itb).second * ducd.get_mass(event, (*itb).first, true);
			else
				lh *= pow(ducd.get_mass(event, (*itb).first, false), (*itb).second);
		}
	}
	return lh;
}

template<typename T> double DiscreteUnivariateConditionalReestimation<T>::kullback_distance_computation(DiscreteUnivariateConditionalDistribution& ducd) const
{
	double d = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
			d += (*itb).second/((double)nb_elements) * (log((*itb).second) - log(((*it).second).first) - ducd.get_mass((*it).first, (*itb).first, true));
	}
	return d;
}

template<typename T> Eigen::VectorXd DiscreteUnivariateConditionalReestimation<T>::y_init() const
{
	Eigen::VectorXd y(nb_different_elements);
	unsigned int i_element = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
		{
			y(i_element) = (*itb).first;
			i_element++;
		}
	}
	return y;
}

template<typename T> bool DiscreteUnivariateConditionalReestimation<T>::make_binomial(Eigen::VectorXd& y) const
{
	std::set<int> categories;
	for(unsigned int i = 0; i < y.size(); ++i)
		categories.insert(y(i));
	if(categories.size() != 2)
	{
		return false;
	} else {
		for(unsigned int i = 0; i < y.size(); ++i)
			y(i) = distance(categories.begin(), categories.find(y(i)));
		return true;
	}
}

template<typename T> Eigen::VectorXd DiscreteUnivariateConditionalReestimation<T>::identity_link_z_init(const Eigen::VectorXd& y) const
{
	return y;
}

template<typename T> Eigen::VectorXd DiscreteUnivariateConditionalReestimation<T>::log_link_z_init(const Eigen::VectorXd& y) const
{
	Eigen::VectorXd z(y.size());
	for(unsigned int i = 0; i < z.size(); ++i)
		z(i) = log(MAX(y(i),0.5));
	return z;
}

template<typename T> Eigen::VectorXd DiscreteUnivariateConditionalReestimation<T>::logit_link_z_init(const Eigen::VectorXd& y) const
{
	Eigen::VectorXd z(y.size());
	for(unsigned int i = 0; i < z.size(); ++i)
		z(i) = log(0.5)/(1-log(0.5));
	return z;
}

template<typename T> Eigen::VectorXd DiscreteUnivariateConditionalReestimation<T>::loglog_link_z_init(const Eigen::VectorXd& y) const
{
	Eigen::VectorXd z(y.size());
	for(unsigned int i = 0; i < z.size(); ++i)
		z(i) = log(-log(0.5));
	return z;
}

template<typename T> void DiscreteUnivariateConditionalReestimation<T>::identity_link_z_update(const Eigen::VectorXd& y, Eigen::VectorXd& z, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const
{
	for(unsigned int i = 0; i < z.size(); ++i)
		z(i) = y(i);
}

template<typename T> void DiscreteUnivariateConditionalReestimation<T>::log_link_z_update(const Eigen::VectorXd& y, Eigen::VectorXd& z, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const
{
	z = X * beta;
	for(unsigned int i = 0; i < z.size(); ++i)
	 z(i) += (y(i)-exp(z(i)))/exp(z(i));
}

template<typename T> void DiscreteUnivariateConditionalReestimation<T>::logit_link_z_update(const Eigen::VectorXd& y, Eigen::VectorXd& z, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const
{
	z = X * beta;
	for(unsigned int i = 0; i < z.size(); ++i)
	 z(i) += (y(i) - exp(z(i))/((1+exp(z(i)))))/(exp(z(i))/(1+exp(z(i)))*(1-exp(z(i))/(1+exp(z(i)))));
}

template<typename T> void DiscreteUnivariateConditionalReestimation<T>::loglog_link_z_update(const Eigen::VectorXd& y, Eigen::VectorXd& z, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const
{
	z = X * beta;
	for(unsigned int i = 0; i < z.size(); ++i)
	 z(i) += (y(i) - (1- exp(-exp(z(i))))/((1-(1-exp(-exp(z(i)))))*log(1-(1-exp(-exp(z(i)))))));
}

template<typename T> Eigen::MatrixXd DiscreteUnivariateConditionalReestimation<T>::X_init() const
{
	Eigen::MatrixXd X(nb_different_elements, nb_covariables+1);
	unsigned int p = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(unsigned int i = 0; i < (((*it).second).second).size(); ++i)
		{
			X(p, 0) = 1;
			for(unsigned int j = 0; j < nb_covariables; ++j)
				X(p, j+1) = ((*it).first)[j];
			p++;
		}
	}
	return X;
}

template<typename T> Eigen::MatrixXd DiscreteUnivariateConditionalReestimation<T>::W_init(const Eigen::VectorXd& y) const
{
	Eigen::MatrixXd W(y.size(), y.size());//; = Eigen::MatrixXd::Zero(y.size(), y.size());
	for(unsigned int i = 0; i < y.size(); ++i)
	{
		for(unsigned int j = 0; j < y.size(); ++j)
			W(i,j) = 0;
		W(i, i) = MAX(y(i), 0.5);
	}
	unsigned int p = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
		{
			W(p, p) = W(p, p)*(*itb).second/((*it).second).first;
			p++;
		}
	}
	return W;
}

template<typename T> void DiscreteUnivariateConditionalReestimation<T>::poisson_identity_link_W_update(Eigen::MatrixXd& W, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const
{
	Eigen::VectorXd Xb = X * beta;
	for(unsigned int i = 0; i < W.cols(); ++i)
		W(i, i) = 1/Xb(i);
		unsigned int p = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
		{
			W(p, p) = W(p, p)*(*itb).second/((*it).second).first;
			p++;
		}
	}
}

template<typename T> void DiscreteUnivariateConditionalReestimation<T>::poisson_log_link_W_update(Eigen::MatrixXd& W, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const
{
	Eigen::VectorXd Xb = X * beta;
	for(unsigned int i = 0; i < W.cols(); ++i)
		W(i, i) = exp(Xb(i));
	unsigned int p = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
		{
			W(p, p) = W(p, p)*(*itb).second/((*it).second).first;
			p++;
		}
	}
}

template<typename T> void DiscreteUnivariateConditionalReestimation<T>::binomial_logit_link_W_update(Eigen::MatrixXd& W, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const
{
	Eigen::VectorXd Xb = X * beta;
	for(unsigned int i = 0; i < W.cols(); ++i)
		W(i, i) = exp(Xb(i))/(1+exp(Xb(i)))*(1-exp(Xb(i))/(1+exp(Xb(i))));
	unsigned int p = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
		{
			W(p, p) = W(p, p)*(*itb).second/((*it).second).first;
			p++;
		}
	}
}

template<typename T> void DiscreteUnivariateConditionalReestimation<T>::binomial_loglog_link_W_update(Eigen::MatrixXd& W, const Eigen::MatrixXd& X, const Eigen::VectorXd& beta) const
{
	Eigen::VectorXd Xb = X * beta;
	for(unsigned int i = 0; i < W.cols(); ++i)
		W(i, i) = (1-(1-exp(-exp(Xb(i)))))*pow((1-(1-exp(-exp(Xb(i))))),2)/((1-exp(-exp(Xb(i)))));//1/(Xb(i)*(1-Xb(i)));
	unsigned int p = 0;
	for(typename std::map< std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
		{
			W(p, p) = W(p, p)*(*itb).second/((*it).second).first;
			p++;
		}
	}
}

template<typename T> Eigen::VectorXd DiscreteUnivariateConditionalReestimation<T>::beta_compute(const Eigen::MatrixXd& X, const Eigen::MatrixXd& W, const Eigen::VectorXd& z) const
{
	return (X.transpose() * W * X).colPivHouseholderQr().solve(X.transpose() * W * z);
}

template<typename T> Parameters DiscreteUnivariateConditionalReestimation<T>::to_stdvector(const Eigen::VectorXd& beta) const
{
	Parameters parameters(beta.size());
	for(unsigned int i = 0; i < beta.size(); ++i)
		parameters[i] = beta(i);
	return parameters;
}

template<typename T> std::multimap<double, std::pair<double, double> > DiscreteUnivariateConditionalReestimation<T>::get_pearson_residuals(DiscreteUnivariateConditionalDistribution& ducd) const
{
	std::multimap<double, std::pair<double, double> > pearson_residuals;
	double m, s;
	for(typename std::map<std::vector<int>, std::pair<T, std::map<int, T> > >::const_iterator it = conditional_histogram.begin(); it != conditional_histogram.end(); ++it)
	{
		m = ducd.get_mean((*it).first);
		s = pow(ducd.get_variance((*it).first), 0.5);
		for(typename std::map<int, T>::const_iterator itb = (((*it).second).second).begin(); itb != (((*it).second).second).end(); ++itb)
			pearson_residuals.insert(std::pair<double, std::pair<double, double> >(m, std::pair<double, double>(((*itb).first-m)/s, (*itb).second/((double)(((*it).second).first)))));
	}
	return pearson_residuals;
}
#endif
