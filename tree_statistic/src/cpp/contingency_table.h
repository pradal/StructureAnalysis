#ifndef CONTINGENCY_TABLE_H
#define CONTINGENCY_TABLE_H

template<typename T> class ContingencyTable {

  public:
	ContingencyTable();
	ContingencyTable(const std::vector< std::vector<T> >& values);
	ContingencyTable(const ContingencyTable<T>& ct);
	~ContingencyTable();

	std::ostream& Summary(std::ostream& os=std::cout) const;
	//Gnuplot<double, double, double>* Plot_residuals() const;
	//Display(const std::vector< std::vector<T> >& values);
	void Estimate();
	double Get_loglikelihood() const;
	ContingencyTable<T>* Drop_parameter(const std::string& parameter) const;
	ContingencyTable<T>* Parsimoniest() const;

  private:
	unsigned int nobs;
	unsigned int ncases;
	unsigned int nvariables;
	unsigned int nparameters;
	std::vector<std::string> parameters_names;
	std::vector< std::vector<double> > varcovar;
	GLM* glm;

	std::vector<std::string> compute_parametersnames(const std::vector< std::set<T> >& icategories) const;
	GLM* GLM_init(const std::vector< std::set<T> >& icategories, std::vector< std::vector<T> > ivalues);
	std::set<std::string> extract_edges() const;
	double accept_against(const ContingencyTable<T>& ct) const;
	std::vector< std::set<T> > compute_categories(const std::vector< std::vector<T> >& values) const;
	std::vector< std::vector<T> > compute_combinations(const std::vector< std::set<T> >& icategories) const;
	std::vector< std::vector<unsigned int> > compute_powerset(const std::vector< std::vector<unsigned int> >& set, unsigned int invariables) const;
	 std::vector<double> Get_means(const std::vector< std::vector<T> >& values) const;
	 std::vector< std::vector<double> > Get_covariances(const std::vector< std::vector<T> >& values) const;
	 unsigned int toInt(const std::string& s) const;
	 bool samecovsign(const std::vector<std::string>& vars) const;
	 std::set<std::string> extract_cliques() const;
	 std::set<std::string> get_decomposition(const std::set<std::string>& iedges) const;
};

#include "./contingency_table.hpp"
#endif
