//
// Created by Martin Steinegger on 11/14/18.
//

#ifndef CONTERMINATOR_INTERVALTREE_H
#define CONTERMINATOR_INTERVALTREE_H


#include <utility>
#include <vector>
#include <iostream>


class IntervalTree {

// Parts of the codes are taken from https://www.geeksforgeeks.org/interval-tree/

public:

    IntervalTree(){
        root=NULL;
        nodesMaxelements = 1024;
        nodesBuffer = (Node*)malloc(nodesMaxelements* sizeof(Node));
        currNodeElement = 0;
    }


    ~IntervalTree(){
        root=NULL;
        free(nodesBuffer);
    }

    struct Node
    {
        Node(int low, int high) : low(low), high(high), left(NULL), right(NULL){}
        // interval data
        int low, high;
        // node data
        int max;
        Node *left, *right;
    };

    void reset(){
        root=NULL;
        currNodeElement=0;
    }

    Node * insert(int low, int high){
        root = insert(root, low, high);
        return root;
    }

    bool doesOverlap(int low, int high){
        return (overlapSearch(root, low, high) == NULL) ? false : true;
    }

    void print(){
        print(root);
    }


private:
    Node *root;
    Node * nodesBuffer;
    int nodesMaxelements;
    int currNodeElement;

    void print(Node *root)
    {
        if (root == NULL) return;

        print(root->left);

        std::cout << "[" << root->low << ", " << root->high << "]"
                  << " max = " << root->max << std::endl;

        print(root->right);
    }

    Node* createNode(int low, int high){
        if(currNodeElement >= nodesMaxelements){
            nodesMaxelements = nodesMaxelements * 2;
            realloc(nodesBuffer, nodesMaxelements * sizeof (Node));
        }

        Node * node = &nodesBuffer[currNodeElement];
        node->left=NULL;
        node->right=NULL;
        node->low  = low;
        node->high = high;
        node->max = high;

        currNodeElement++;
        return node;
    }

    Node * insert(Node *root, int low, int high) {
        // Tree is empty
        if (root == NULL)
            return createNode(low, high);

        // Get low value of interval at root
        int l = root->low;

        if (low < l){
            root->left = insert(root->left, low, high);
        }else{
            root->right = insert(root->right, low, high);
        }
        // Update the max value of this ancestor if needed
        if (root->max < high) {
            root->max = high;
        }

        return root;
    }


    bool checkOverlap(int i1Low, int i1High, int i2Low, int i2High)
    {
        return (i1Low <= i2High && i2Low<= i1High) ? true : false;
    }


    Node * overlapSearch(Node *root, int low, int high)
    {
        // Base Case, tree is empty
        if (root == NULL) {
            return NULL;
        }

        // If given interval overlaps with root
        if (checkOverlap(root->low, root->high, low, high))
            return root;

        // If left child of root is present and max of left child is
        // greater than or equal to given interval, then i may
        // overlap with an interval is left subtree
        if (root->left != NULL && root->left->max >= low)
            return overlapSearch(root->left, low, high);

        // Else interval can only overlap with right subtree
        return overlapSearch(root->right,  low, high);
    }

};


#endif //CONTERMINATOR_INTERVALTREE_H
