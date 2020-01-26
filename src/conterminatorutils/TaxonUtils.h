//
// Created by Martin Steinegger on 2019-08-20.
//

#ifndef CONTERMINATOR_TAXONUTILS_H
#define CONTERMINATOR_TAXONUTILS_H
//
// Created by Martin Steinegger on 2019-07-17.
//

#include <string>
#include <algorithm>
#include "Util.h"
#include "DBWriter.h"
#include "NcbiTaxonomy.h"
#include "KingdomExpression.h"

class TaxonUtils{
public:
    struct RangEntry {
        unsigned int range;
        unsigned int kingdom;
        unsigned int species;
        unsigned int id;

        RangEntry() {};
        RangEntry(unsigned int range, unsigned int kingdom, unsigned int species, unsigned int id) :
                range(range), kingdom(kingdom), species(species), id(id) {}

        bool operator()(const RangEntry &lhs, const RangEntry &rhs) {
            if (lhs.range < rhs.range) return true;
            if (rhs.range < lhs.range) return false;
            if (lhs.kingdom < rhs.kingdom) return true;
            if (rhs.kingdom < lhs.kingdom) return false;
            if (lhs.species < rhs.species) return true;
            if (rhs.species < lhs.species) return false;
            if (lhs.id < rhs.id) return true;
            return false;
        }
    };

    struct TaxonInformation{
        TaxonInformation(unsigned int dbKey, int currTaxa, int ancestorTax, int start, int end, char * data) :
                dbKey(dbKey), currTaxa(currTaxa), termId(ancestorTax), start(start), end(end), data(data), overlaps(false), range(-1){}
        unsigned int dbKey;
        int currTaxa;
        int termId;
        int start;
        int end;
        char * data;
        bool overlaps;
        int range;
        static bool compareByRange(const TaxonInformation &first, const TaxonInformation &second) {
            if (first.range < second.range)
                return true;
            if (second.range < first.range)
                return false;
            if(first.termId < second.termId )
                return true;
            if(second.termId < first.termId )
                return false;
            if(first.currTaxa < second.currTaxa )
                return true;
            if(second.currTaxa < first.currTaxa )
                return false;
            return false;
        }

        static bool compareByDbKeyAndStart(const TaxonInformation &first, const TaxonInformation &second) {
            if (first.dbKey < second.dbKey)
                return true;
            if (second.dbKey < first.dbKey)
                return false;
            if(first.start < second.start )
                return true;
            if(second.start < first.start )
                return false;
            return false;
        }

        static bool compareByTaxAndStart(const TaxonInformation &first, const TaxonInformation &second) {
            if (first.termId < second.termId)
                return true;
            if (second.termId < first.termId)
                return false;
            if(first.start < second.start )
                return true;
            if(second.start < first.start )
                return false;
            if(first.currTaxa < second.currTaxa )
                return true;
            if(second.currTaxa < first.currTaxa )
                return false;
            return false;
        }
    };

    static bool compareToFirstInt(const std::pair<unsigned int, unsigned int>& lhs, const std::pair<unsigned int, unsigned int>& rhs){
        return (lhs.first <= rhs.first);
    }

    static unsigned int getTaxon(unsigned int id, std::vector<std::pair<unsigned int, unsigned int>> &mapping) {
        std::pair<unsigned int, unsigned int> val;
        val.first = id;
        std::vector<std::pair<unsigned int, unsigned int>>::iterator mappingIt = std::upper_bound(
                mapping.begin(), mapping.end(), val, compareToFirstInt);
        if (mappingIt->first != val.first) {
            return UINT_MAX;
        }
        return  mappingIt->second;
    }

    static std::vector<TaxonInformation> assignTaxonomy(std::vector<TaxonInformation>  &elements,
                                                        char *data, std::vector<std::pair<unsigned int, unsigned int>> & mapping,
                                                        NcbiTaxonomy & t, KingdomExpression & kingdomExpression,
                                                        std::vector<int> &blacklist,
                                                        size_t * taxaCounter, bool parseDbKey = false) {
        elements.clear();
        const char * entry[255];
        while (*data != '\0') {
            int termIndex;
            const size_t columns = Util::getWordsOfLine(data, entry, 255);
            if (columns == 0) {
                data = Util::skipLine(data);
                continue;
            }
            unsigned int id = Util::fast_atoi<unsigned int>(entry[0]);
            unsigned int taxon = getTaxon(id, mapping);;
            if (taxon == 0 || taxon == UINT_MAX) {
                goto next;
            }
            // remove blacklisted taxa
            for (size_t j = 0; j < blacklist.size(); ++j) {
                if (t.IsAncestor(blacklist[j], taxon)) {
                    goto next;
                }
            }
            termIndex = kingdomExpression.isAncestorOf(taxon);
            if(termIndex != -1) {
                taxaCounter[termIndex]++;
                int startPos = Util::fast_atoi<int>(entry[4]);
                int endPos   = Util::fast_atoi<int>(entry[5]);
                if(parseDbKey){
                    startPos = Util::fast_atoi<int>(entry[7]);
                    endPos   = Util::fast_atoi<int>(entry[8]);
                }
                elements.push_back(TaxonInformation(id, taxon, termIndex, std::min(startPos, endPos), std::max(startPos, endPos), data));
            }else{
                elements.push_back(TaxonInformation(id, taxon, 0, -1, -1, data));

            }
            next:
            data = Util::skipLine(data);
        }
        return elements;
    }
};


#endif //CONTERMINATOR_TAXONUTILS_H
