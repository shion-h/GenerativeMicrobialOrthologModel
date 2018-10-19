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

double calculateLogBetaFunction(vector<double> alpha, int i/*=-1*/){//{{{
    double resultValue = 0;
    if(i == -1){
        for(int i=0; i<alpha.size(); i++){
            if(alpha[i] == 0.0){
                resultValue += -230.26/* == lgamma(1.0e-100)*/;
            }else{
                resultValue += std::lgamma(alpha[i]);
            }
        }
        resultValue -= std::lgamma(sum(alpha));
        return resultValue;
    }else{
        if(alpha[i] == 0.0){
            resultValue += -230.26/* == lgamma(1.0e-100)*/;
        }else{
            resultValue += std::lgamma(alpha[i]);
        }
        resultValue -= std::lgamma(sum(alpha));
        return resultValue;
    }
}//}}}

double calculateDirichletLogPDF(vector<double> x, vector<double> alpha, int i/*=-1*/){//{{{
    if(sum(x) > 1.1 || sum(x) < 0.9) cout<<"Sum of sample is not equals to 1."<<endl;
    if(x.size() != alpha.size()) cout<<"Sample dimension do not match parameter dimension."<<endl;
    if(i == -1){
        double term1 = calculateLogBetaFunction(alpha);
        double term2 = 0;
        for(int i=0; i<x.size(); i++){
            if(x[i] == 0.0){
                term2 += (alpha[i] - 1.0) * (-230.26)/* == log(1.0e-100)*/;
            }else{
                term2 += (alpha[i] - 1.0) * log(x[i]);
            }
        }
        return -term1 + term2;
    }else{
        double term1 = calculateLogBetaFunction(alpha, i);
        double term2 = 0;
        if(x[i] == 0.0){
            term2 += (alpha[i] - 1.0) * (-230.26)/* == log(1.0e-100)*/;
        }else{
            term2 += (alpha[i] - 1.0) * log(x[i]);
        }
        return -term1 + term2;
    }
}//}}}

double calculateLogFactorial(unsigned int x){//{{{
    double res = 0;
    for(int i=0; i<x; i++){
        res += log(i + 1);
    }
    return res;
}//}}}

double calculatePoissonPMF(unsigned int x, double lambda){//{{{
    double res = x * log(lambda) - lambda - calculateLogFactorial(x);
    return exp(res);
}//}}}

GibbsSamplerFromGMOM::GibbsSamplerFromGMOM(const CsvFileParser<double> &orthologFile, const CsvFileParser<double> &microbeFile, double A, double k, double theta, unsigned int iterationNumber, unsigned int burnIn, unsigned int samplingInterval, bmpi::communicator &world)//{{{
    :_logO(orthologFile.getExtractedMatrix()),
     _U(microbeFile.getExtractedMatrix()),
     _N(_logO.size()),
     _M(_U[0].size()),
     _G(_logO[0].size()),
     _V(_M),
     _P(_M),
     _A(A),
     _k(k),
     _theta(theta),
     _iterationNumber(iterationNumber),
     _burnIn(burnIn),
     _samplingInterval(samplingInterval),
     _samplingCount(0),
     _world(world)
{
    this->initializeParameters();
}//}}}

GibbsSamplerFromGMOM::~GibbsSamplerFromGMOM(){//{{{

}//}}}

void GibbsSamplerFromGMOM::initializeParameters(){//{{{
    _V = vector<vector<unsigned int> >(_M, vector<unsigned int>(_G, 0));
    _sumOfVSampled = vector<vector<double> >(_M, vector<double>(_G, 0.0));
    _P = vector<vector<double> >(_M, vector<double>(_G, 0));
    _sumOfPSampled = vector<vector<double> >(_M, vector<double>(_G, 0));
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
    poisson_distribution<> rand(1.0);
    if(_world.rank()==0){
        for(int j=0; j<_M; j++){
            for(int k=0; k<_G; k++){
                _V[j][k] = rand(mt);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    bmpi::broadcast(_world, _V, 0);
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

double GibbsSamplerFromGMOM::calculateDirichletLogPDF(int k/*=-1*/){//{{{
    double term1 = 0;
    double term2 = 0;
    if(k == -1){
        for(int i=0; i<_N; i++){
            term1 += calculateLogBetaFunction(_gamma[i]);
            for(int k=0; k<_G; k++){
                term2 += (_gamma[i][k] - 1.0) * _logO[i][k];
            }
        }
    }else{
        for(int i=0; i<_N; i++){
            term1 += calculateLogBetaFunction(_gamma[i], k);
            term2 += (_gamma[i][k] - 1.0) * _logO[i][k];
        }
    }
    return -term1 + term2;
}//}}}

void GibbsSamplerFromGMOM::sampleV(){//{{{
    for(int j=0; j<_M; j++){
        for(int k=0; k<_G; k++){
            double cumulativeProbability = 0.0;
            vector<double> logSamplingDistribution;
            unsigned int v = 0;
            while(cumulativeProbability < 0.999){
                double thisPMF = calculatePoissonPMF(v, _P[j][k]);
                cumulativeProbability += thisPMF;
                logSamplingDistribution.push_back(log(thisPMF));
                v++;
            }
            this->updateGamma(j, k, 0 - _V[j][k]);
            // if(_world.rank() == 0)cout<<logSamplingDistribution.size()<<endl;
            for(v=0; v<logSamplingDistribution.size(); v++){
                if((v%_world.size()) == _world.rank()){
                    logSamplingDistribution[v] += this->calculateDirichletLogPDF(k);
                }
                this->updateGamma(j, k, 1);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            for(v=0; v<logSamplingDistribution.size(); v++){
                bmpi::broadcast(_world, logSamplingDistribution[v], v%_world.size());
            }
            if(_world.rank() == 0){
                _V[j][k] = sampleDiscreteValues(logSamplingDistribution, true);
            }
            bmpi::broadcast(_world, _V[j][k], 0);
            this->updateGamma(j, k, _V[j][k] - v);
        }
    }
}//}}}

void GibbsSamplerFromGMOM::sampleP(){//{{{
    random_device rnd;
    mt19937 mt(rnd());
    for(int j=0; j<_M; j++){
        for(int k=0; k<_G; k++){
            if(((j * _G + k)%_world.size()) == _world.rank()){
                std::gamma_distribution<> distribution(_V[j][k] + _k, 1 / (1 + 1/_theta));
                _P[j][k] = distribution(mt);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for(int j=0; j<_M; j++){
        for(int k=0; k<_G; k++){
            bmpi::broadcast(_world, _P[j][k], (j * _G + k)%_world.size());
        }
    }
}//}}}

void GibbsSamplerFromGMOM::calculateLogLikelihood(){//{{{
    double logLikelihood = this->calculateDirichletLogPDF();
    cout<<setprecision(numeric_limits<double>::max_digits10);
    cout<<logLikelihood<<endl;
    _logLikelihood.push_back(logLikelihood);
}//}}}

void GibbsSamplerFromGMOM::writeParameters(string PFilename, string VFilename)const{//{{{
    outputVector(_sumOfPSampled, PFilename);
    outputVector(_sumOfVSampled, VFilename);
}//}}}

void GibbsSamplerFromGMOM::writeLogLikelihood(string logLikelihoodFilename)const{//{{{
    outputVector(_logLikelihood, logLikelihoodFilename);
}//}}}

void GibbsSamplerFromGMOM::storeSamples(){//{{{
    _samplingCount++;
    for(int j=0; j<_M; j++){
        for(int k=0; k<_G; k++){
            _sumOfVSampled[j][k] += _V[j][k];
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
            _sumOfVSampled[j][k] /= _samplingCount;
            _sumOfPSampled[j][k] /= _samplingCount;
        }
    }
}//}}}
