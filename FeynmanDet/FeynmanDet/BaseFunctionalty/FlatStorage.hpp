//
//  Stats-NZPaths.cpp
//  FeynmanDet
//
//  Created by Luis Paulo Santos on 24/02/2026.
//

#ifndef _FlatStorage
#define _FlatStorage

#include <vector>
#include <iostream>

class FlatStorage {
private:
    std::vector<unsigned int> data;
    std::vector<float> wR, wI;
    std::size_t element_size;   // size of each stored array

public:
    FlatStorage(std::size_t n)
        : element_size(n) {}

    void push_back(const unsigned long long* arr, float pR, float pI) {
        unsigned int arr_ui[128];
        
        for (int i=0; i< element_size; i++)
            arr_ui[i] = (unsigned int)arr[i];
        data.insert(data.end(), arr_ui, arr_ui + element_size);
        wR.push_back(pR);
        wI.push_back(pI);
    }

    unsigned int* operator[](std::size_t i) {
        return data.data() + i * element_size;
    }

    std::size_t size() const {
        return data.size() / element_size;
    }
    
    void file_save (FILE *f) {
        std::size_t const N = this->size();
        for (std::size_t i=0 ; i<N  ; i++) {
            fwrite((*this)[i], sizeof(unsigned int), element_size, f);
            float const pR = wR.at(i);
            float const pI = wI.at(i);
            fwrite(&pR, sizeof(float), 1, f);
            fwrite(&pI, sizeof(float), 1, f);
        }
    }
};

#endif


