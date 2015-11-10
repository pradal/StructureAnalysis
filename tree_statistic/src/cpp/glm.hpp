#ifndef GLM_HPP
#define GLM_HPP

GLM::GLM(const Eigen::VectorXd& yi, const Eigen::MatrixXd& Xi, unsigned int ifamily, unsigned int ilink) {
  y = yi;
  X = Xi;
  family = ifamily;
  link = ilink;

  z = z_init();
  W = W_init();
}

GLM::GLM(const GLM& glm) {
  copy(glm);
}

void GLM::copy(const GLM& glm) {
  y = glm.y;
  X = glm.X;
  family = glm.family;
  link = glm.link;
  z = glm.z;
  W = glm.W;
  I = glm.I;
  beta = glm.beta;
}

GLM::~GLM() {
}

GLM* GLM::Drop_parameter(const std::vector<bool>& todel) const {
  unsigned int i,j,d,s;
  s = std::accumulate(todel.begin(), todel.end(),0);
  std::vector< std::vector<double> > X0(X.cols()-s, std::vector<double>(X.rows()));

  for(i = 0; i < X.rows(); i++) {
	d = 0;
	for(j = 0; j < X.cols(); j++) {
	  if(todel[j])
		d++;
	  else
	   X0[j-d][i] = X(i,j);
	}
  }

  Eigen::MatrixXd Xnew(X0[0].size(), X0.size());
  for(i = 0; i < Xnew.rows(); i++) {
	for(j = 0; j < Xnew.cols(); j++) {
	  Xnew(i,j) = X0[j][i];
	}
  }
  GLM* glmt = new GLM((*this).y,Xnew,(*this).family,(*this).link);
  return glmt;
}

Eigen::VectorXd GLM::z_init() const {
  unsigned int i;
  Eigen::VectorXd z0 = y;

  switch(link) {
	case IDENTITY : {
	  break;
	}
	case LOG : {
	  for(i = 0; i < z0.size(); i++) {
		z0(i) = std::log(MAX(z0(i),0.5));
	  }
	  break;
	}
	case LOGIT : {
	  for(i = 0; i < z0.size(); i++) {
		z0(i) = y(i);//std::log(MAX(0.5,z0(i))/(1-MIN(0.5,y(i))));
	  }
	  break;
	}
  }

  return z0;
}

void GLM::z_update() {
  unsigned int i;
  z = X * beta;

  switch(link) {
	case IDENTITY : {
	  break;
	}
	case LOG : {
	  for(i = 0; i < z.size(); i++) {
		z(i) = z(i) + (y(i)-std::exp(z(i)))/(std::exp(z(i)));
	  }
	  break;
	}
	case LOGIT : {
	  for(i = 0; i < z.size(); i++) {
		z(i) = z(i) + (y(i)-std::exp(z(i))/(1+std::exp(z(i))))/(std::exp(z(i))/(1+std::pow(1+std::exp(z(i)),2)));
	  }
	  break;
	}
  }
}

Eigen::MatrixXd GLM::W_init() const {
  unsigned int i,j;
  Eigen::MatrixXd W0(y.size(), y.size());

  for(i = 0; i < W0.rows(); i++) {
	W0(i,i) = MAX(0.5,y(i));
	for(j = i+1; j < W0.cols(); j++) {
	  W0(i,j) = 0;
	  W0(j,i) = 0;
	}
  }

  return W0;
}

void GLM::W_update() {
  unsigned int i;
  Eigen::VectorXd Xb = X * beta;

  switch(link) {
	case IDENTITY : {
	  break;
	}
	case LOG : {
	  for(i = 0; i < W.rows(); i++) {
		W(i,i) = std::exp(Xb(i));
	  }
	  break;
	}
	case LOGIT : {
	  for(i = 0; i < W.rows(); i++) {
		W(i,i) = std::exp(Xb(i))/(std::pow((1+std::exp(Xb(i))),2));
	  }
	  break;
	}
  }
}

bool GLM::compute_beta() {
	beta = (X.transpose() * W * X).colPivHouseholderQr().solve(X.transpose() * W * z);
	double relative_error = ((X.transpose() * W * X) * beta - (X.transpose() * W * z)).norm()/((X.transpose() * W * z).norm());
	if(relative_error > EPSILON) {
		beta = (X.transpose() * W * X).fullPivHouseholderQr().solve(X.transpose() * W * z);
		relative_error = ((X.transpose() * W * X) * beta - (X.transpose() * W * z)).norm()/((X.transpose() * W * z).norm());
		if(relative_error > EPSILON)
		  return false;
		else
		  return true;
	} else {
	  return true;
	}
}

double GLM::compute_loglikelihood() const {
  unsigned int i;
  Eigen::VectorXd Xb = X * beta;
  //Eigen::VectorXd b(y.size());
  Eigen::VectorXd a(y.size());
  Eigen::VectorXd b(y.size());
  Eigen::VectorXd c(y.size());

  switch(link) {
	case IDENTITY : {
	  break;
	}
	case LOG : {
	  for(i = 0; i < b.size(); i++) {
		b(i) = std::exp(Xb(i));
	  }
	  break;
	}
	case LOGIT : {
	  for(i = 0; i < b.size(); i++) {
		b(i) = std::log(1+std::exp(Xb(i)));
	  }
	  break;
	}
  }

  switch(family) {
        case stat_tool::POISSON : {
	  for(i = 0; i < y.size(); i++) {
		a(i) = 1;
		c(i) = -boost::math::lgamma(MAX(y(i),1)+1);
	  }
	  break;
	}
	case stat_tool::BINOMIAL : {
	  for(i = 0; i < y.size(); i++) {
		a(i) = 1;
		c(i) = 0;
	  }
	}
  }

  double loglikelihood = 0;
  for(i = 0; i < y.size(); i++) {
	loglikelihood += (y(i)*Xb(i)-b(i))/a(i)+c(i);
  }

  return loglikelihood;
}

bool GLM::Estimate(std::string& comment) {
  bool success = compute_beta(), conv = false;
  unsigned int it = 1;

  if(success) {
	double prev = compute_loglikelihood(), cur;
	if(boost::math::isnan(prev))
	  success = false;
	while(it < MAXITS && success && !conv) {
	  z_update();
	  W_update();
	  success = compute_beta();
	  if(success) {
		cur = compute_loglikelihood();
		if(boost::math::isnan(cur))
		  success = false;
		if((cur-prev)/std::abs(prev) < EPSILON) {
		  conv = true;
		} else {
		  prev = cur;
		  it++;
		}
	  } else {
		if(comment.size() != 0)
		  comment += "\n";
		comment += "Error: update of parameters failed !";
	  }
	}
	if(it > MAXITS) {
	  if(comment.size() != 0)
		comment += "\n";
	  comment += "Error: no convergence of estimation algorithm !";
	}
  }
  if(conv && !SILENT) {
	std::cout << "Estimated in " + toString(it) + " iteration(s)." << std::endl;
  }
  if(!conv && !SILENT) {
	std::cout << comment << std::endl;
	comment.clear();
  }

  return (success && conv);
}

std::vector<double> GLM::Get_parameters() const {
  unsigned int i;
  std::vector<double> vparameters;

  for(i = 0; i < beta.size(); i++) {
	vparameters.push_back(beta(i));
  }

  return vparameters;
}

unsigned int GLM::Get_family() const {
  return family;
}

unsigned int GLM::Get_link() const {
  return link;
}

std::multimap<double, double> GLM::Get_residuals() const {
  unsigned int i;
  double yi;
  std::multimap<double, double> residuals;
  Eigen::VectorXd Xb = X * beta;
  Eigen::MatrixXd H = X * (X.transpose()* W * X).inverse() * X.transpose();
  double s = 0;

  switch(link) {
	case LOG : {
	  for(i = 0; i < y.size(); i++) {
		residuals.insert(std::pair<double, double>(std::exp(Xb(i)), (y(i)-std::exp(Xb(i)))/(sqrt(std::exp(Xb(i))))));
		s +=  pow(y(i)-std::exp(Xb(i)),2);
	  }
	  break;
	}
	case LOGIT : {
	  for(i = 0; i < y.size(); i++) {
		residuals.insert(std::pair<double, double>((std::exp(Xb(i)))/(1+std::exp(Xb(i))),(y(i)-(std::exp(Xb(i)))/(1+std::exp(Xb(i))))/(sqrt(std::exp(Xb(i))*(1-std::exp(Xb(i)))))));
		s += pow(y(i)-(std::exp(Xb(i)))/(1+std::exp(Xb(i))),2);
	  }
	  break;
	}
  }
  if(y.size() != beta.size()) {
	s /= (y.size()-beta.size());
	/*std::multimap<double, double>::iterator it;
	for(it = residuals.begin(); it != residuals.end(); it++) {
	  (*it).second /= s;
	}*/
  }

  return residuals;
}

unsigned int GLM::Get_nobs() const {
  return y.size();
}

double GLM::Get_loglikelihood() const {
  return compute_loglikelihood();
}

std::vector<double> GLM::Wald_test() const {
  unsigned int i;
  std::vector<double> wald_stat;
  Eigen::MatrixXd I = X.transpose() * W * X;
  Eigen::MatrixXd Iinv = I.inverse();
  for(i = 0; i < beta.size(); i++) {
	wald_stat.push_back(pow(beta(i),2)/Iinv(i,i));
  }
  return wald_stat;
}
#endif
