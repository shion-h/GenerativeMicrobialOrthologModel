//
// GMOM.cpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

//include{{{
#include <stdlib.h>
#include <boost/program_options.hpp>
#include "include/CsvFileParser.hpp"
#include "include/GibbsSamplerFromGMOM.hpp"
//}}}

using namespace std;

int main(int argc, char *argv[]){

    //Options{{{
    boost::program_options::options_description opt("Options");
    opt.add_options()
    ("help,h", "show help")
    ("otpt,o", boost::program_options::value<string>()->default_value("./"), "directory name for output")
    ("itnm,n", boost::program_options::value<unsigned int>()->default_value(1000), "the number of iteration")
    ("intr,i", boost::program_options::value<unsigned int>()->default_value(5), "sampling interval")
    ("bnin,b", boost::program_options::value<unsigned int>()->default_value(500), "burn-in term")
    ("k,k", boost::program_options::value<double>()->default_value(1.0), "k value")
    ("theta,t", boost::program_options::value<double>()->default_value(1.0), "theta value")
    ("A,A", boost::program_options::value<double>()->default_value(1000.0), "A value")
    ;

    boost::program_options::positional_options_description pd;
    pd.add("orthologfile", 1);
    pd.add("microbefile", 2);

    boost::program_options::options_description hidden("hidden");
    hidden.add_options()
        ("orthologfile", boost::program_options::value<string>(), "hidden")
        ("microbefile", boost::program_options::value<string>(), "hidden")
        ;
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(opt).add(hidden);

    boost::program_options::variables_map vm;
    try{
        boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pd).run(), vm);
    }catch(const boost::program_options::error_with_option_name& e){
        cout<<e.what()<<endl;
    }
    boost::program_options::notify(vm);
    string orthologFilename;
    string microbeFilename;
    double k, theta, A;
    unsigned int iterationNumber;
    unsigned int samplingInterval;
    unsigned int burnIn;
    string outputDirectory;
    string PFilename;
    string VFilename;
    string logLikelihoodFilename;
    if (vm.count("help") || !vm.count("orthologfile") || !vm.count("microbefile")){
        cout<<"Usage:\n GMOM [ortholog file] [microbe file] [-options] "<<endl;
        cout<<endl;
        cout<<opt<<endl;
        exit(1);
    }else{
        orthologFilename = vm["orthologfile"].as<std::string>();
        microbeFilename = vm["microbefile"].as<std::string>();
        if(vm.count("k"))k = vm["k"].as<double>();
        if(vm.count("theta"))theta = vm["theta"].as<double>();
        if(vm.count("A"))A = vm["A"].as<double>();
        if(vm.count("itnm"))iterationNumber = vm["itnm"].as<unsigned int>();
        if(vm.count("intr"))samplingInterval = vm["intr"].as<unsigned int>();
        if(vm.count("bnin"))burnIn = vm["bnin"].as<unsigned int>();
        if(vm.count("otpt"))outputDirectory = vm["otpt"].as<std::string>();
        if(outputDirectory[outputDirectory.size()-1] != '/')outputDirectory.push_back('/');
        PFilename = outputDirectory + "P.csv";
        VFilename = outputDirectory + "V.csv";
        logLikelihoodFilename = outputDirectory + "LogLikelihood.csv";
    }

    CsvFileParser<double> orthologFile(orthologFilename);
    CsvFileParser<double> microbeFile(microbeFilename);
    orthologFile.convertLog();
    //}}}

//estimation{{{
    GibbsSamplerFromGMOM *estimator;
    estimator = new GibbsSamplerFromGMOM(orthologFile, microbeFile, A, k, theta, iterationNumber, burnIn, samplingInterval);
    estimator->runIteraions();
    estimator->writeParameters(PFilename, VFilename);
    estimator->writeLogLikelihood(logLikelihoodFilename);
    delete estimator;
//}}}
    return 0;
}
