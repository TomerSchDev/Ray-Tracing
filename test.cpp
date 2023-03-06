#include "json.hpp"
#include "iostream"
#include "fstream"
using json = nlohmann::json;
int main(int length, char **path){
    std::cout<<path[1]<<std::endl;
    std::ifstream dataFile(path[1]);
    if (dataFile)
        {
            auto jsonData = json::parse(dataFile);
            std::cout<<jsonData;
            dataFile.close();
        }
}