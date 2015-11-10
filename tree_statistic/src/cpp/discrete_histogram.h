#ifndef DISCRETE_HISTOGRAM_H
#define DISCRETE_HISTOGRAM_H

#ifndef I_MAX
#define I_MAX std::numeric_limits<int>::max()
#endif

#ifndef I_MIN
#define I_MIN std::numeric_limits<int>::min()
#endif

#include <numeric>
#include <limits>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>
#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <algorithm>
#include <set>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>

typedef double Moment;
typedef std::vector< Moment > Moments;
typedef std::vector< Moments > PairedMoments;

class DiscreteUnivariateHistogram
{
	public:
		DiscreteUnivariateHistogram();
		DiscreteUnivariateHistogram(const DiscreteUnivariateHistogram& hist);
		~DiscreteUnivariateHistogram();
		DiscreteUnivariateHistogram(const std::vector<int>& vector);
		double GetMean() const;
		double GetVariance() const;
		int GetMinimum() const;
		int GetMaximum() const;
		void AddEvent(int event);
		void AddEvents(const std::vector<int>& events);
		std::vector<int> GetEvents(bool exhaustive=false) const;
		unsigned int GetNbEvent(int event) const;
		unsigned int GetNbEvents() const;
		template<typename Y> std::vector< std::pair<Y, int> > GetPlotable() const;
	
	private:
		unsigned int total;
		std::map<int, unsigned int> histogram;
		Moment mean;
		Moment variance;
		void mean_computation();
		void variance_computation();
};

class DiscreteMultivariateHistogram
{
	public:
		DiscreteMultivariateHistogram();
		DiscreteMultivariateHistogram(const DiscreteMultivariateHistogram& hist);
		~DiscreteMultivariateHistogram();
		DiscreteMultivariateHistogram(const std::vector< std::vector<int> >& vectors);
		DiscreteMultivariateHistogram(unsigned int invariables);
		void AddEvent(std::vector<int> vector);
		void AddEvents(const std::vector< std::vector<int> >& events);
		DiscreteMultivariateHistogram* GetHistogram(const std::set<unsigned int>& variables) const;
		DiscreteUnivariateHistogram* GetMarginalHistogram(unsigned int with) const;
		DiscreteMultivariateHistogram* GetMarginalHistogram(const std::pair<unsigned int, unsigned int>& with) const;
		DiscreteMultivariateHistogram* GetMarginalHistogram(const std::set<unsigned int>& with) const;
		DiscreteMultivariateHistogram* GetMarginalizedHistogram(unsigned int with) const;
		DiscreteMultivariateHistogram* GetMarginalizedHistogram(const std::pair<unsigned int, unsigned int>& with) const;
		DiscreteMultivariateHistogram* GetMarginalizedHistogram(const std::set<unsigned int>& with) const;
		Moments GetMeans() const;
		PairedMoments GetCVM() const;
		std::vector<int> GetMinima() const;
		std::vector<int> GetMaxima() const;
		std::vector< std::vector<int> > GetEvents(bool exhaustive=false) const;
		unsigned int GetNbBreaks() const;
		unsigned int GetNbEvent(const std::vector<int>& event) const;
		unsigned int GetNbEvents() const;
		unsigned int GetNbVariables() const;
		std::vector< std::pair< unsigned int, unsigned int> > GetPairs() const;
		template<typename Y> std::vector< std::pair<Y, int> > GetPlotable(unsigned int variable) const;
		template<typename Z> std::vector< boost::tuple<Z, int, int> > GetPlotable(const std::pair<unsigned int, unsigned int>& variables) const;

	private:
		unsigned int total;
		std::set<unsigned int> toIgnore;
		unsigned int nvariables;
		std::map< std::vector<int>, unsigned int > histogram;
		Moments means;
		PairedMoments cvm;
		std::vector<int> minima;
		std::vector<int> maxima;
		void covariance_computation();
};

#endif
