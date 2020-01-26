//
// Created by Martin Steinegger on 2020-01-25.
//

#ifndef CONTERMINATOR_KINGDOMEXPRESSION_H
#define CONTERMINATOR_KINGDOMEXPRESSION_H

#include "TaxonomyExpression.h"
#include "Debug.h"
#include <vector>

class KingdomExpression {


    std::vector<TaxonomyExpression*> taxTerms;

public:
    KingdomExpression(std::string & expression, NcbiTaxonomy &taxonomy){
        std::vector<std::string> expressionList = Util::split(expression, ",");
        for(size_t i = 0; i < expressionList.size(); i++){
            taxTerms.push_back(new TaxonomyExpression(expressionList[i],taxonomy));
        }
    }


    // this function returns the index of the term that fulfils the criteria
    // -1 means no term fulfils the criteria
    int isAncestorOf( unsigned int taxId){
        int index = -1;
        bool ancestor = false;
        for (size_t j = 0; j < taxTerms.size() && !ancestor; ++j) {
            ancestor |= taxTerms[j]->isAncestor(taxId);
            index = (ancestor) ? j : index;
        }
        return index;
    }

    ~KingdomExpression(){
        for (size_t i = 0; i < taxTerms.size(); i++){
            delete taxTerms[i];
        }
    }

    std::vector<TaxonomyExpression*> getTaxTerms(){
        return taxTerms;
    }
};


#endif //CONTERMINATOR_KINGDOMEXPRESSION_H
