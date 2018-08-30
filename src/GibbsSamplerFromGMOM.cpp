//
// GibbsSamplerFromGMOM.cpp
//
// Copyright (c) 2018 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/GibbsSamplerFromGMOM.hpp"

using namespace std;

double calculateLogBetaFunction(vector<double> alpha){//{{{
    double resultValue = 0;
    for(int i=0; i<alpha.size(); i++){
        if(alpha[i] == 0.0){
            resultValue += std::lgamma(1.0e-10);
        }else{
            resultValue += std::lgamma(alpha[i]);
        }
    }
    resultValue -= std::lgamma(sum(alpha));
    return resultValue;
}//}}}

double calculateDirichletLogPDF(vector<double> x, vector<double> alpha){//{{{
    if(sum(x) > 1.1 || sum(x) < 0.9) cout<<"Sum of sample is not equals to 1."<<endl;
    if(x.size() != alpha.size()) cout<<"Sample dimension do not match parameter dimension."<<endl;
    double term1 = calculateLogBetaFunction(alpha);
    double term2 = 0;
    for(int i=0; i<x.size(); i++){
        if(x[i] == 0.0){
            term2 += (alpha[i] - 1.0) * log(1.0e-10);
        }else{
            term2 += (alpha[i] - 1.0) * log(x[i]);
        }
    }
    return -term1 + term2;
}//}}}

unsigned int calculateFactorial(unsigned int x){//{{{
    unsigned int res = 1;
    for(int i=0; i<x; i++){
        res *= i + 1;
    }
    return res;
}//}}}

double calculatePoissonLogPMF(unsigned int x, double lambda){//{{{
    double res = pow(lambda, x) * exp(-lambda) + calculateFactorial(x);
    return res;
}//}}}

GibbsSamplerFromGMOM::GibbsSamplerFromGMOM(const CsvFileParser<double> &orthologFile, const CsvFileParser<double> &microbeFile, double A, double alpha, double beta, unsigned int iterationNumber, unsigned int burnIn, unsigned int samplingInterval)//{{{
    :_O(orthologFile.getExtractedMatrix()),
     _U(microbeFile.getExtractedMatrix()),
     _N(_O.size()),
     _M(_U[0].size()),
     _G(_O[0].size()),
     _V(_M),
     _P(_M),
     _A(A),
     _alpha(alpha),
     _beta(beta),
     _iterationNumber(iterationNumber),
     _burnIn(burnIn),
     _samplingInterval(samplingInterval),
     _samplingCount(0)
{
    this->initializeParameters();
}//}}}

GibbsSamplerFromGMOM::~GibbsSamplerFromGMOM(){//{{{

}//}}}

void GibbsSamplerFromGMOM::initializeParameters(){//{{{
    _V = vector<vector<unsigned int> >(_M, vector<unsigned int>(_G, 0));
    _P = vector<vector<double> >(_M, vector<double>(_G, 0.0));
    _sumOfPSampled = vector<vector<double> >(_M, vector<double>(_G, 0.0));
    // random_device rnd;
    // mt19937 mt(rnd());
    // uniform_real_distribution<> rand(0.0, 1.0);
    // for(int j=0; j<_M; j++){
    //     for(int k=0; k<_G; k++){
    //         boost::math::beta_distribution<> distribution(_alpha, _beta);
    //         _P[j][k] = boost::math::quantile(distribution, rand(mt));
    //     }
    // }
    random_device rnd;
    mt19937 mt(rnd());
    uniform_int_distribution<> rand(0, 1);
    for(int j=0; j<_M; j++){
        for(int k=0; k<_G; k++){
            _V[j][k] = rand(mt);
        }
    }
    _gamma = dot<double, double, unsigned int>(_U, _V);
    for(int i=0; i<_N; i++){
        for(int k=0; k<_G; k++){
            _gamma[i][k] *= _A;
            _gamma[i][k]++;
        }
    }
}//}}}

void GibbsSamplerFromGMOM::updateGamma(unsigned int j, unsigned int k, int deltaVjk){//{{{
    for(int i=0; i<_N; i++){
        _gamma[i][k] += deltaVjk * _U[i][j] * _A;
    }
}//}}}

void GibbsSamplerFromGMOM::sampleV(){//{{{
    for(int j=0; j<_M; j++){
        for(int k=0; k<_G; k++){
            unsigned int oldVjk = _V[j][k];
            double cumulativeProbability = 0.0;
            vector<double> logSamplingDistribution(2, 0);
            unsigned int v = 0;
            while(cumulativeProbability > 0.999){
                double thisPMF = calculatePoissonLogPMF(v, _P[k]);
                cumulativeProbability += thisPMF;
                logSamplingDistribution.push_back(thisPMF);
                v++;
            }
            for(int v=0; v<logSamplingDistribution.size(); v++){
                this->updateGamma(j, k, v - oldVjk);
                for(int i=0; i<_N; i++){
                    logSamplingDistribution[v] += calculateDirichletLogPDF(_O[i], _gamma[i]);
                }
            }
            _V[j][k] = sampleDiscreteValues(logSamplingDistribution, true);
            if(_V[j][k] != v){
                int deltaVjk = _V[j][k] - v;
                this->updateGamma(j, k, deltaVjk);
            }
        }
    }
}//}}}

void GibbsSamplerFromGMOM::sampleP(){//{{{
    random_device rnd;
    mt19937 mt(rnd());
    uniform_real_distribution<> rand(0.0, 1.0);
    for(int j=0; j<_M; j++){
        for(int k=0; k<_G; k++){
            boost::math::beta_distribution<> distribution(_V[j][k] + _alpha, 1 - _V[j][k] + _beta);
            _P[j][k] = boost::math::quantile(distribution, rand(mt));
        }
    }
}//}}}

void GibbsSamplerFromGMOM::calculateLogLikelihood(){//{{{
    double logLikelihood = 0;
    for(int i=0; i<_N; i++){
        logLikelihood += calculateDirichletLogPDF(_O[i], _gamma[i]);
    }
    cout<<setprecision(numeric_limits<double>::max_digits10);
    cout<<logLikelihood<<endl;
    _logLikelihood.push_back(logLikelihood);
}//}}}

void GibbsSamplerFromGMOM::writeParameters(string PFilename, string VFilename)const{//{{{
    outputVector(_sumOfPSampled, PFilename);
    outputVector(_V, VFilename);
}//}}}

void GibbsSamplerFromGMOM::writeLogLikelihood(string logLikelihoodFilename)const{//{{{
    outputVector(_logLikelihood, logLikelihoodFilename);
}//}}}

void GibbsSamplerFromGMOM::storeSamples(){//{{{
    _samplingCount++;
    for(int j=0; j<_M; j++){
        for(int k=0; k<_G; k++){
            _sumOfPSampled[j][k] += _P[j][k];
        }
    }
}//}}}

void GibbsSamplerFromGMOM::runIteraions(){//{{{
    for(int i=0; i<_iterationNumber; i++){
        this->sampleP();
        this->sampleV();
        this->calculateLogLikelihood();
        if((i > _burnIn - 1) && ((i + 1) % _samplingInterval == 0)){
            this->storeSamples();
        }
    }
    for(int j=0; j<_M; j++){
        for(int k=0; k<_G; k++){
            _sumOfPSampled[j][k] /= _samplingCount;
        }
    }
}//}}}
