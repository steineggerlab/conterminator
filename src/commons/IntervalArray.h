//
// Created by Martin Steinegger on 11/14/18.
//

#ifndef CONTERMINATOR_INTERVALARRAY_H
#define CONTERMINATOR_INTERVALARRAY_H

#include "MathUtil.h"
#include <utility>
#include <list>
#include <iostream>


class IntervalArray {

public:

    IntervalArray(){
        const int size = 1;
        array=(unsigned char *)calloc(size, sizeof(unsigned char));
        arraySizeInBytes = size;
        maxSizeInByte = size*8;
    }


    ~IntervalArray(){
        free(array);
        array = NULL;
    }

    void reset(){
        int ceilMax = MathUtil::ceilIntDivision(std::max(1,maxSizeInByte), 8);
        memset(array, 0, ceilMax * sizeof(unsigned char));
        maxSizeInByte = 0;
    }

    unsigned char getLowMask(unsigned int rest) {
        const unsigned char mask[8] = {0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE0, 0xC0,0x80};
        return mask[rest];
    }


    unsigned char getHighRest(unsigned int rest) {
        const unsigned char mask[8] = { 0xFF, 0x7F, 0x3F, 0x1F, 0x0F,0x07, 0x03, 0x01};
         return mask[rest];
    }


    void insert(int low, int high){
        if(low > high){
            std::swap(low, high);
        }
        //insert(array, low, high);
        maxSizeInByte = std::max(high, maxSizeInByte + 1);
        int ceilMax = MathUtil::ceilIntDivision(maxSizeInByte, 8);
        if(ceilMax >= arraySizeInBytes){
            int prevSize = arraySizeInBytes;
            arraySizeInBytes = std::max(arraySizeInBytes * 2, ceilMax+1);
            array = (unsigned char *) realloc(array,  arraySizeInBytes);
            memset(array+prevSize, 0,  arraySizeInBytes-prevSize);
        }
        bool lowFound = isSet(low);
        bool highFound =  isSet(high);;
        if((lowFound == true && highFound == true)){
            return;
        }
        unsigned int startPos=low/8;
        unsigned int startRest=low%8;
        unsigned int endPos=high/8;
        unsigned int endRest=high%8;
        for(size_t pos = startPos+1; pos < endPos; pos++ ){
            array[pos] = 0xFF;
        }
        unsigned char lowMask = getLowMask(startRest);
        unsigned char highMask = getHighRest(7-endRest);
        if(startPos == endPos){
            array[startPos] |= lowMask & highMask;
        }else{
            array[startPos] |= lowMask;
            array[endPos] |= highMask;
        }
    }


    bool doesOverlap(int low, int high)
    {

        if(low > high){
            std::swap(low, high);
        }
        bool lowFound;
        if(low >= maxSizeInByte){
            return false;
        }else{
            lowFound = isSet(low);
        }
        bool highFound;
        if(high >= maxSizeInByte){
            highFound = false;
        }else{
            highFound = isSet(high);;
        }
        return (lowFound || highFound);
    }

    bool isSet(int pos){
        unsigned int posIdx = pos/8;
        unsigned int posRest= pos%8;
        unsigned char value = array[posIdx];
        unsigned char check = (1U << posRest);
        return check & value;
    }


    void print(){
        bool started = false;
        for(int pos = 0; pos < maxSizeInByte; pos++){
            if(isSet(pos)  && started == false){
                started = true;
                std::cout << "[" << pos << ", ";
            }
            if(isSet(pos)==false && started== true){
                started = false;
                std::cout << pos - 1 << "]" << std::endl;
            }
        }
        if(started == true){
            std::cout << maxSizeInByte - 1 << "]" << std::endl;
        }
    }


private:
    unsigned char * array;
    int arraySizeInBytes;
    int maxSizeInByte;

};


#endif //CONTERMINATOR_INTERVALARRAY_H
