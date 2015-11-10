#ifndef GLM_H
#define GLM_H

class GLM {
  public:
	GLM(const GLM& glm);
	GLM(const Eigen::VectorXd& yi, const Eigen::MatrixXd& Xi, unsigned int ifamily, unsigned int ilink);
	~GLM();

	bool Estimate(std::string& comment);
	std::vector<double> Get_parameters() const;
	unsigned int Get_nobs() const;
	unsigned int Get_family() const;
	unsigned int Get_link() const;
	double Get_loglikelihood() const;
	std::multimap<double, double> Get_residuals() const;
	GLM* Drop_parameter(const std::vector<bool>& todel) const;
	std::vector<double> Wald_test() const;

  private:
	void copy(const GLM& glm);

	unsigned int family;
	unsigned int link;

	Eigen::VectorXd y;
	Eigen::VectorXd beta;
	Eigen::VectorXd z;
	Eigen::MatrixXd X;
	Eigen::MatrixXd W;
	Eigen::MatrixXd I;

	bool compute_beta();
	Eigen::VectorXd z_init() const;
	void z_update();
	Eigen::MatrixXd W_init() const;
	void W_update();

	double compute_loglikelihood() const;
};

#include "./glm.hpp"
#endif
