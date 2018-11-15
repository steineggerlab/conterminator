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

// Parts of the codes are taken from https://www.geeksforgeeks.org/interval-tree/

public:

    IntervalArray(){
        array=(unsigned char *)calloc(1024, sizeof(unsigned char));
        arraySize = 1024;
        maxSize = 1024*8;
    }


    ~IntervalArray(){
        free(array);
        array = NULL;
    }

    struct Node
    {
        Node() {}
        // interval data
        int low, high;
        // node data
        int max;
        Node *left, *right;
    };

    void reset(){
        memset(array, 0, (maxSize/8) * sizeof(unsigned char));
    }

    unsigned char getLowMask(unsigned int rest) {
        //const unsigned char mask[8] = {0x80, 0xC0, 0xE0, 0xF0, 0xF8, 0xFC, 0xFE, 0xFF};
        const unsigned char mask[8] = {0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE0, 0xC0,0x80};
        return mask[rest];
    }


    unsigned char getHighRest(unsigned int rest) {
        const unsigned char mask[8] = { 0xFF, 0x7F, 0x3F, 0x1F, 0x0F,0x07, 0x03, 0x01};

//        const unsigned char mask[8] = { 0x80, 0xC0, 0xE0, 0xF0, 0xF8, 0xFC, 0xFE, 0xFF };
         return mask[rest];
    }


    void insert(int low, int high){
        //insert(array, low, high);
        maxSize = std::max(high, maxSize);
        if((maxSize/8) >= arraySize){
            int ceil = MathUtil::ceilIntDivision(maxSize, 8);
            arraySize = std::max(arraySize *2, ceil);
            array = (unsigned char *) realloc(array,  arraySize);
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

        bool lowFound = isSet(low);
        bool highFound =  isSet(high);;
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
        for(int pos = 0; pos < maxSize; pos++){
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
            std::cout << maxSize - 1 << "]" << std::endl;
        }
    }


private:
    unsigned char * array;
    int arraySize;
    int maxSize;

};


#endif //CONTERMINATOR_INTERVALARRAY_H
