//
// GibbsSamplerFromGMOM.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef GIBBSGMOM
#define GIBBSGMOM

#include<stdlib.h>
#include<math.h>
#include<cmath>
#include<iostream>
#include<vector>
#include<numeric>
#include<memory>
#include<random>
#include<iomanip>
#include<fstream>
#include<limits>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include"CsvFileParser.hpp"
#include"utils.hpp"


class GibbsSamplerFromGMOM{
protected:
    const std::vector<std::vector<double> > &_O, &_U;
    const unsigned int _N, _M, _G;
    const unsigned int _iterationNumber, _burnIn, _samplingInterval;
    const double _alpha, _beta, _A;
    std::vector<std::vector<unsigned int> > _V;
    std::vector<std::vector<double> > _P, _sumOfPSampled, _gamma;
    unsigned int _samplingCount;
    std::vector<double> _logLikelihood;
public:
    GibbsSamplerFromGMOM(const CsvFileParser<double> &orthologFile, const CsvFileParser<double> &microbeFile, double A, double alpha, double beta, unsigned int iterationNumber, unsigned int burnIn, unsigned int samplingInterval);
    virtual ~GibbsSamplerFromGMOM();
    virtual void initializeParameters();
    virtual void updateGamma(unsigned int j, unsigned int k, int deltaVjk);
    virtual void sampleV();
    virtual void sampleP();
    virtual void calculateLogLikelihood();
    virtual void writeParameters(std::string PFilename, std::string VFilename)const;
    virtual void writeLogLikelihood(std::string logLikelihoodFilename)const;
    virtual void storeSamples();
    virtual void runIteraions();
};

#endif
