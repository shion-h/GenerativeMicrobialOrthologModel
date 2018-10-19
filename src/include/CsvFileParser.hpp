//
// CsvFileParser.hpp
//
// Copyright (c) 2017 Shion Hosoda
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//

#ifndef CsvFILEPASER
#define CsvFILEPASER
#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include<cmath>

template<class T>
class CsvFileParser{
protected:
    std::ifstream _inputText;
    std::vector<std::string> _columnNameVector;
    std::vector<std::string> _rowNameVector;
    unsigned int _rowDimension;
    unsigned int _columnDimension;
    std::vector<std::vector<T> > _extractedMatrix;
public:
    CsvFileParser(std::string filename);
    ~CsvFileParser();
    void readCsvFile();
    void convertLog();
    const std::vector<std::vector<T> > &getExtractedMatrix()const{
        return _extractedMatrix;
    }
    const std::vector<std::string> &getColumnNameVector()const{
        return _columnNameVector;
    }
    const std::vector<std::string> &getRowNameVector()const{
        return _rowNameVector;
    }
    const unsigned int &getColumnDimension()const{
        return _columnDimension;
    }
    const unsigned int &getRowDimension()const{
        return _rowDimension;
    }
};

template<class T>
CsvFileParser<T>::CsvFileParser(std::string filename):_inputText(filename),_rowDimension(0),_columnDimension(0){//{{{
    std::cout<<filename<<" was opened."<<std::endl;
    this->readCsvFile();
}//}}}

template<class T>
CsvFileParser<T>::~CsvFileParser(){//{{{
    _inputText.close();
}//}}}

template<class T>
void CsvFileParser<T>::readCsvFile(){//{{{
    if(!_inputText){
        std::cout<<"Cannot open Csvfile";
        exit(1);
    }
    std::string str;
    // row
    for(int i=0; std::getline(_inputText,str); i++){
        std::vector<T> vec;
        std::string token;
        std::istringstream stream(str);

        if(i==0){
            //column
            for(int j=0; std::getline(stream,token,','); j++){
                if(j==0)continue;
                _columnNameVector.push_back(token);
            }
            continue;
        }
        //column
        for(int j=0; std::getline(stream,token,','); j++){
            if(j==0)_rowNameVector.push_back(token);
            else{
                try{
                    vec.push_back(boost::lexical_cast<T>(token));
                }catch(...){
                    std::cout<<"row:"<<i<<" column:"<<j<<" value:"<<token<<" error";
                    exit(0);
                }
            }
        }
        _extractedMatrix.push_back(vec);
    }

    _columnDimension = _columnNameVector.size();
    _rowDimension = _rowNameVector.size();
}//}}}

template<class T>
void CsvFileParser<T>::convertLog(){//{{{
    for(int i=0; i<_extractedMatrix.size(); i++){
        for(int j=0; j<_extractedMatrix[i].size(); j++){
            if(_extractedMatrix[i][j] > 10e-10){
                _extractedMatrix[i][j] = std::log(_extractedMatrix[i][j]);
            }else{
                _extractedMatrix[i][j] = -230.26;
            }
        }
    }
}//}}}
#endif
