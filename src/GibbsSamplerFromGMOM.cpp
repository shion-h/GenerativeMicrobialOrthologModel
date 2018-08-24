//
// GibbsSamplerFromGMOM.cpp
//
// Copyright (c) 2018 Shion Hosoda
//
// This software is released u_nder the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#include"include/GibbsSamplerFromGMOM.hpp"

using namespace std;

double calculateLogBetaFunction(vector<double> alpha){//{{{
    double resultValue = 0;
    for(int i=0; i<alpha.size(); i++){
        resultValue += boost::math::lgamma(alpha[i]);
    }
    resultValue -= boost::math::lgamma(sum(alpha));
    return resultValue;
}//}}}

double calculateDirichletLogPDF(vector<double> x, vector<double> alpha){//{{{
    if(sum(x) > 1.1 || sum(x) < 0.9) cout<<"Sum of sample is not equals to 1."<<endl;
    if(x.size() != alpha.size()) cout<<"Sample dimension do not match parameter dimension."<<endl;
    double term1 = calculateLogBetaFunction(alpha);
    double term2 = 0;
    for(int i=0; i<x.size(); i++){
        term2 += (alpha[i] - 1.0) * log(x[i]);
    }
    return -term1 + term2;
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
    _gamma = prod<double, double, unsigned int>(_U, _V);
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
            bool isFrom0 = true;
            if(oldVjk == 1) isFrom0 = false;
            vector<double> logSamplingDistribuion(2, 0);
            int step = 1;
            unsigned int initialIndex = 0;
            if(!isFrom0){
                step = -1;
                initialIndex = 1;
            }
            unsigned int v;
            for(int b=0; b<2; b++){
                v = (b - initialIndex) * step;
                if(b == 1){
                    this->updateGamma(j, k, step);
                }
                for(int i=0; i<_N; i++){
                    logSamplingDistribuion[v] += calculateDirichletLogPDF(_O[i], _gamma[i]);
                }
                if(v == 0){
                    logSamplingDistribuion[v] += log(1 - _P[j][k]);
                }else{
                    logSamplingDistribuion[v] += log(_P[j][k]);
                }
            }
            _V[j][k] = sampleDiscreteValues(logSamplingDistribuion, true);
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
